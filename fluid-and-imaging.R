## ----setup, echo=FALSE, message=FALSE, warning=FALSE------------------------------------------------
# devtools::install_url('https://cran.rstudio.com/src/contrib/Archive/calibFit/calibFit_2.1.0.tar.gz')
# remotes::install_github('atrihub/SRS')
# For ADNIMERGE, go to http://adni.loni.usc.edu/, https://adni.bitbucket.io/

library(Hmisc)
library(knitr)
library(tidyverse)
library(kableExtra)
library(gridExtra)
library(plotly)
library(calibFit)
library(SRS)
library(mixtools)
library(pROC)
library(party)
library(MoEClust)
library(nlme)
library(emmeans)
library(arsenal)
library(mvtnorm)
library(mice)

options(digits=2)

xaringanExtra::use_tile_view() # type 'o' to invoke tile view
theme_set(theme_bw())

# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
cbbPalette <-
    c("#0072B2", "#D55E00", "#E69F00",
      "#009E73", "#F0E442", "#999999",
      "#000000", "#56B4E9", "#CC79A7")
scale_colour_discrete <-
    function(...) scale_colour_manual(..., values = cbbPalette)
scale_fill_discrete <-
    function(...) scale_fill_manual(..., values = cbbPalette)
scale_colour_discrete <-
    function(...) scale_colour_manual(..., values = cbbPalette)
scale_fill_discrete <-
    function(...) scale_fill_manual(..., values = cbbPalette)

theme_table <- function(..., levs=2){
  theme_minimal(...) + 
    theme(
      panel.grid = element_blank(), 
      axis.text.x = element_blank(),
      axis.text.y = element_text(face='bold', color=cbbPalette[1:levs]),
      axis.title = element_blank())
}

center <- function(x) scale(x, scale = FALSE)






## ----echo=FALSE, fig.align='center', out.width='57%'------------------------------------------------
knitr::include_graphics("./images/atri.png")


## ----echo=FALSE, fig.align='center', out.width='47%'------------------------------------------------
knitr::include_graphics("./images/actc_logo.png")


## ----generate_batch_data----------------------------------------------------------------------------
# simulated data with batch effects
set.seed(20200225)

batch_data <- 
  tibble(
    batch = 1:10,
    Sigma = rgamma(n=10, shape=360/10, scale=10), # variance for each batch
    Mean = rnorm(n=10, mean=850, sd=200)) %>% # Mean for each batch
  group_by(batch) %>%
  nest() %>%
  mutate(Biomarker = map(data, ~ rnorm(n=50, .$Mean, .$Sigma))) %>%
  unnest(Biomarker) %>%
  unnest(data) %>%
  ungroup() %>%
  arrange(batch) %>%
  mutate(
    id = 1:(10*50),
    batch = as.factor(batch),
    Biomarker = ifelse(Biomarker<0, 0, Biomarker))


## ----batch_data_plot--------------------------------------------------------------------------------
ggplot(batch_data, aes(y=Biomarker, x=batch)) +
  geom_boxplot(outlier.shape=NA) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3, alpha=0.2)


## ----batch_data_summaries, results='asis', cache=FALSE----------------------------------------------
batch_data_sum <- batch_data %>% group_by(batch) %>%
  summarize(
    N=length(Biomarker),
    Mean=mean(Biomarker), 
    SD=sd(Biomarker),
    `SD/Mean = CV (%)`=SD/Mean*100)
batch_data_sum %>%
  kable(digits=2, format = 'html') %>%
  kable_styling(
    bootstrap_options=c('striped', 'condensed'),
    font_size=18, full_width=FALSE)


## ----echo=TRUE, results='markup'--------------------------------------------------------------------
anova(lm(Biomarker ~ batch, batch_data))


## ----batch_confounds--------------------------------------------------------------------------------
low_groups <- subset(batch_data_sum, Mean<median(batch_data_sum$Mean))$batch
batch_data <- batch_data %>%
  mutate(Group = ifelse(batch %in% low_groups, 'A', 'B'))
ggplot(batch_data, aes(y=Biomarker, x=batch)) +
  geom_boxplot(outlier.shape=NA) +
  geom_dotplot(aes(color=Group, fill=Group), 
    binaxis='y', stackdir='center', dotsize=0.3, alpha=0.5)


## ----batch_randomized-------------------------------------------------------------------------------
batch_data$Group <- sample(batch_data$Group, size=nrow(batch_data))
ggplot(batch_data, aes(y=Biomarker, x=batch)) +
  geom_boxplot(outlier.shape=NA) +
  geom_dotplot(aes(color=Group, fill=Group), 
    binaxis='y', stackdir='center', dotsize=0.3, alpha=0.5)


## ----randomization----------------------------------------------------------------------------------
data(srs_data)
# head(srs_data)

p.func.greedy.if.possible <- function(overallImbalance, treatmentCounts, maxCounts)
{
    cant.go <- treatmentCounts > maxCounts
    if(all(cant.go)) 
      stop("Randomization impossible. Probably need another treatment group.")

    number.of.treatments <- length(overallImbalance)
    k <- which(overallImbalance == min(overallImbalance))
    p.vec <- rep(0, number.of.treatments)
    p.vec[k] <- 1
    p.vec/sum(p.vec)
    p.vec[cant.go] <- 0

    if(all(p.vec == 0)){ # try less greedy
      number.of.treatments <- length(overallImbalance)
      p.star <- 2/3
      k <- which(overallImbalance == min(overallImbalance))
      if (length(k) > 1) {
          k <- sample(k, 1)
      }
      p.vec <- rep((1 - p.star)/(number.of.treatments - 1), number.of.treatments)
      p.vec[k] <- p.star
      p.vec      
      p.vec[cant.go] <- 0
    }
    
    p.vec
}

get.counts <- function(object)
{
  expt <- object@expt
  treatment.names <- expt@treatment.names
  factor.names <- expt@factor.names
  factor.level.names <- expt@factor.level.names
  treatment.names <- expt@treatment.names
  state.matrix <- object@stateTable
  tr.ratios <- object@tr.ratios
  
  tr.assignments <- object@tr.assignments
  tr.assignments$Treatment <- factor(tr.assignments$Treatment, 
    levels = treatment.names)
  tr.assignments$Counts <- factor(tr.assignments$Counts, 
    levels = factor.level.names[[which(factor.names == "Counts")]])
  tr.assignments <- with(tr.assignments, table(Counts, Treatment)) * 
    as.numeric(factor.level.names[[which(factor.names == "Counts")]])
  colSums(tr.assignments)
}  

expt <- ClinicalExperiment(number.of.factors = 3,
  factor.names = c('Counts', 'Group', 'Age'),
  number.of.factor.levels = c(2, 5, 2),
  factor.level.names = 
    list(c(4, 5), 1:5, c('young', 'old')),
  number.of.treatments = 13,
  treatment.names = as.character(1:13))

g.func <- function(imbalances)
{
    factor.weights <- c (1, 100, 1)
    imbalances %*% factor.weights
}

r.obj <- new("cPocockSimonRandomizer", expt, as.integer(20130827), 
  g.func=g.func, p.func = p.func.greedy.if.possible, max.counts = 30)

for(i in 1:nrow(srs_data)){
  r.obj <- randomize(r.obj, as.character(srs_data[i, "ID"]), 
     as.character(srs_data[i, expt@factor.names]))
}


## ----results='asis'---------------------------------------------------------------------------------
tr.assignments <- r.obj@tr.assignments %>%
  mutate(
    Treatment = factor(Treatment, levels = r.obj@expt@treatment.names),
    `Subject ID` = 1:nrow(r.obj@tr.assignments)
  ) %>%
  rename(
    Plate = Treatment,
    `Num. of samples` = Counts) %>%
  select(`Subject ID`, `Num. of samples`, Group, Age, Plate)

tr.assignments[1:10, ] %>%
  kable(digits=2, format = 'html') %>%
  kable_styling(
    bootstrap_options=c('striped', 'condensed'),
    font_size=18, full_width=FALSE)


## ---------------------------------------------------------------------------------------------------
tab <- with(tr.assignments, table(Plate, Age)) %>%
  as.data.frame() %>%
  pivot_wider(names_from='Age', values_from='Freq') %>%
  t() 
tab <- rbind(tab, `Num. samples`=as.numeric(get.counts(r.obj)))

tab %>%
  kable(digits=2, format = 'html') %>%
  kable_styling(
    bootstrap_options=c('striped', 'condensed'),
    font_size=18, full_width=FALSE)


## ----calibFit_fits, out.width='100%', fig.height=4, fig.width=4*(2)---------------------------------
data(HPLC)
data(ELISA)

linmodel <- lm(Response~Concentration, data=HPLC)
# The predicted response
HPLC$Fitted <- fitted(linmodel)

p1 <- ggplot(HPLC, aes(x=Concentration, y=Response)) +
  geom_point() +
  geom_line(aes(y=Fitted)) +
	xlab("Concentration (ng/ml)") +
	ylab("Response") +
	ggtitle("HPLC with ordinary least squares fit")

fplmodel <- with(ELISA,
  calib.fit(Concentration, Response, type="log.fpl")
)
# The predicted response
ELISA$Fitted <- fplmodel@fitted.values

p2 <- ggplot(ELISA, aes(x=log(Concentration), y=Response)) +
  geom_point() +
  geom_line(aes(y=Fitted)) +
	xlab("log(Concentration (ng/ml))") +
	ylab("Response") +
	ggtitle("ELISA with 4 parameter logistic fit")

grid.arrange(p1,p2,nrow=1)


## ----calibFit_residuals, out.width='100%', fig.height=4, fig.width=4*(2)----------------------------
p1 <- ggplot(HPLC, aes(x=Concentration, y=Response-Fitted)) +
  geom_point() +
  geom_hline(yintercept = 0) +
	xlab("Concentration (ng/ml)") +
	ylab("Residuals") +
	ggtitle("HPLC with ordinary least squares fit")

p2 <- ggplot(ELISA, aes(x=log(Concentration), y=Response-Fitted)) +
  geom_point() +
  geom_hline(yintercept = 0) +
	xlab("log(Concentration (ng/ml))") +
	ylab("Residuals") +
	ggtitle("ELISA with 4 parameter logistic fit")

grid.arrange(p1,p2,nrow=1)


## ----regression-movie, fig.show='animate', dev='jpeg'-----------------------------------------------
set.seed(20210524)
X <- rnorm(200, mean = 2, sd = 0.5)
e <- rnorm(200, mean = 0, sd = 0.5)
Y <-  X*2 + e 
lmfit <- lm(Y ~ X)
Fitted <- predict(lmfit)

Coef <- lmfit$coef + c(0, 1)

for(i in 1:5){
  Coef <- Coef + c(0, -0.21)
  Fitted <- as.matrix(cbind(1, X)) %*% Coef
  RSS <- sum((Fitted - X)^2)
  p <- qplot(x=X, y=Y) + 
    geom_point() + 
    geom_abline(intercept=Coef[1], slope=Coef[2]) +
    annotate("text", x = 1, y = 7, label = paste('RSS* =', round(RSS, digits = 0)))
  print(p)
}

RSS <- sum(lmfit$residuals^2)
p <- qplot(x=X, y=Y) + 
  geom_point() + 
  geom_abline(intercept=lmfit$coefficients[1], slope=lmfit$coefficients[2]) +
  annotate("text", x = 1, y = 7, label = paste('RSS* =', round(RSS, digits = 0)))
print(p)


## ----calib_fit--------------------------------------------------------------------------------------
cal.fpl <- with(ELISA, calib.fit(Concentration,Response,type="log.fpl"))
cal.lin.pom <- with(HPLC, calib.fit(Concentration,Response,type="lin.pom"))
cal.fpl.pom <- with(ELISA, calib.fit(Concentration,Response,type="log.fpl.pom"))

linpom.fit <- cal.lin.pom@fitted.values
fplpom.fit <- cal.fpl.pom@fitted.values

sig.lin <- cal.lin.pom@sigma
sig.fpl <- cal.fpl.pom@sigma

theta.lin <- cal.lin.pom@theta
theta.fpl <- cal.fpl.pom@theta

linpom.res <- cal.lin.pom@residuals*(1/((linpom.fit^theta.lin)*sig.lin))
fplpom.res <- cal.fpl.pom@residuals*(1/((fplpom.fit^theta.fpl)*sig.fpl))


## ----calibFit_pom_residuals, out.width='100%', fig.height=4, fig.width=4*(2)------------------------
p1 <- ggplot(HPLC, aes(x=linpom.fit, y=linpom.res)) +
  geom_point() +
  geom_hline(yintercept = 0) +
	xlab("Fitted Values (LS-POM)") +
	ylab("Standardized Residuals") +
	ggtitle("HPLC with least squares POM")

p2 <- ggplot(ELISA, aes(x=fplpom.fit, y=fplpom.res)) +
  geom_point() +
  geom_hline(yintercept = 0) +
	xlab("Fitted Values (FPL-POM)") +
	ylab("Standardized Residuals") +
	ggtitle("ELISA with 4 parameter logistic POM")

grid.arrange(p1,p2,nrow=1)


## ----calib_hplc_pom, fig.height=4.5, fig.width=4.5*2------------------------------------------------
par(mfrow=c(1,2))
ciu <- fitted(linmodel) + summary(linmodel)$sigma*qt(.975,linmodel$df)
cil <- fitted(linmodel) - summary(linmodel)$sigma*qt(.975,linmodel$df)

plot(HPLC[, c('Concentration', 'Response')], main = "HPLC data fit without POM",col="blue",pch=16)
lines(HPLC$Concentration,fitted(linmodel),col="lightblue")
lines(HPLC$Concentration,ciu,col="lightblue",lty=2)
lines(HPLC$Concentration,cil,col="grey",lty=2)

plot(cal.lin.pom,print=FALSE,main = "HPLC data fit with POM",xlab = "Concentration", ylab = "Response")


## ----calib_elisa_pom, fig.height=4.5, fig.width=4.5*2-----------------------------------------------
par(mfrow=c(1,2))
#par(mar = c(3.5,3.5,1.5,1.5))
plot(cal.fpl,print=FALSE,main = "ELISA fit without POM",xlab = "Concentration", ylab = "Response")

#par(mar = c(3.5,3.5,1.5,1.5))
plot(cal.fpl.pom,print=FALSE,main = "ELISA fit with POM",xlab = "Concentration", ylab = "Response")


## ----calibrated1, fig.height=4.5, fig.width=4.5, out.width='100%'-----------------------------------
calib.lin <- calib(cal.lin.pom, HPLC$Response)
plot(calib.lin, main="HPLC calribated with linear POM")


## ----calibrated2, fig.height=4.5, fig.width=4.5, out.width='100%'-----------------------------------
calib.fpl <- calib(cal.fpl.pom, ELISA$Response)
plot(calib.fpl, main="ELISA calribated with FPL POM")


## ----classification, fig.height=4, fig.width=6------------------------------------------------------
dd <- subset(ADNIMERGE::adnimerge, !is.na(ABETA)) %>%
  arrange(RID, EXAMDATE) %>%
  filter(!duplicated(RID)) %>%
  mutate(
    ABETA = as.numeric(gsub('<', '', gsub('>', '', ABETA))),
    TAU = as.numeric(gsub('<', '', gsub('>', '', TAU))),
    PTAU = as.numeric(gsub('<', '', gsub('>', '', PTAU))),
    DX2 = factor(DX, levels = c('CN', 'Dementia')))

ggplot(dd, aes(x=ABETA, y=TAU)) + 
  geom_point(aes(color=DX)) +
  scale_color_manual(values=c("#0072B2", "#E69F00", "#D55E00")) +
  theme(legend.position = c(0.80,0.75))


## ----classification_no_mci--------------------------------------------------------------------------
ggplot(subset(dd, DX!='MCI'), aes(x=ABETA, y=TAU)) + 
  geom_point(aes(color=DX2)) +
  scale_color_manual(values=c("#0072B2", "#D55E00"))


## ----ROC_abeta, fig.width=5, fig.height=5-----------------------------------------------------------
roc_abeta <- roc(DX2 ~ ABETA, subset(dd, DX!='MCI'))
# ggroc(roc_abeta)
plot(roc_abeta, print.thres="best", print.thres.best.method="youden")


## ----ROC_same, fig.width=5, fig.height=5------------------------------------------------------------
f1 <- function(x){
  unlist(lapply(x, function(y){
    if(y<0.5) return(0)
    if(y>=0.5) return(y) 
  }))
}

set.seed(20250423)
roc_same_auc <- tibble(x = seq(0,1, length = 10000)) %>%
  mutate(p1 = f1(x)) %>%
  rowwise() %>%
  mutate(Label1 = sample(c('A', 'B'), size = 1, prob = c(p1, 1-p1)))
roc_same_auc$Label2 <- rev(roc_same_auc$Label1)
roc_same_auc$Label2 <- ifelse(roc_same_auc$Label2 == 'B', 'A', 'B')

roc_1 <- roc(Label1 ~ x, roc_same_auc)
roc_2 <- roc(Label2 ~ x, roc_same_auc)
plot(roc_1, print.thres="best", print.thres.best.method="youden")
plot(roc_2, add=TRUE, col='orange', print.thres="best",
  print.thres.best.method="youden")


## ----ROC_abeta_tau, fig.width=5, fig.height=5-------------------------------------------------------
roc_tau <- roc(DX2 ~ TAU, subset(dd, DX!='MCI'))
roc_tauabeta <- roc(DX2 ~ I(TAU/ABETA), subset(dd, DX!='MCI'))

plot(roc_abeta, col='blue')
plot(roc_tau, add=TRUE, col='orange')
plot(roc_tauabeta, add=TRUE, col='red',
  print.thres="best", print.thres.best.method="youden")
legend("bottomright", 
  legend=c(expression(paste("A", beta)), "Tau", expression(paste("Tau/A", beta))),
       col=c("blue", "orange", "red"), lwd=2)



## ----echo=FALSE, eval=FALSE-------------------------------------------------------------------------
# roc.test(roc_abeta, roc_tau)
# roc.test(roc_abeta, roc_tauabeta)


## ----abeta_tau_scatter_youden, fig.height=2.5, fig.width=3*2----------------------------------------
ggplot(subset(dd, DX!='MCI'), aes(x=ABETA, y=TAU)) + 
  geom_point(aes(color=DX2)) +
  scale_color_manual(values=c("#0072B2", "#D55E00")) +
  geom_abline(intercept = 0, slope=0.394)


## ---------------------------------------------------------------------------------------------------
logistic.fit <- glm(DX2 ~ scale(ABETA) + scale(TAU), 
  data = subset(dd, DX!='MCI'), family=binomial)
summary(logistic.fit)$coef %>%
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var='Coefficient') %>%
  mutate(`Pr(>|z|)` = format.pval(round(`Pr(>|z|)`, digits = 3), eps = 0.001, digits=3)) %>%
  kable()


## ----logistic_pred_prob-----------------------------------------------------------------------------
ggplot(subset(dd, DX!='MCI'), aes(x=ABETA, y=TAU)) + 
  geom_point(aes(color=predict(logistic.fit, type='response'))) +
  scale_colour_gradient(low="#0072B2", high="#D55E00") +
  labs(color="Pred. prob. of Dementia") +
  geom_abline(intercept = 0, slope=0.394)


## ----ratio_gradient---------------------------------------------------------------------------------
ggplot(subset(dd, DX!='MCI'), aes(x=ABETA, y=TAU)) + 
  geom_point(aes(color=TAU/ABETA)) +
  scale_colour_gradient(low="#0072B2", high="#D55E00") +
  labs(color="Tau/Abeta ratio") +
  geom_line(data=data.frame(ABETA=80:1700, TAU=80:1700), 
            aes(color=TAU/ABETA), linetype='dashed') +
  geom_line(data=data.frame(ABETA=80:1700, TAU=1.5*(80:1700)), aes(color=TAU/ABETA)) +
  geom_line(data=data.frame(ABETA=80:1700, TAU=2*(80:1700)), aes(color=TAU/ABETA)) +
  geom_line(data=data.frame(ABETA=80:1700, TAU=1/2*(80:1700)), aes(color=TAU/ABETA)) +
  geom_line(data=data.frame(ABETA=80:1700, TAU=1/3*(80:1700)), aes(color=TAU/ABETA)) +
  geom_line(data=data.frame(ABETA=80:1700, TAU=1/4*(80:1700)), aes(color=TAU/ABETA)) +
  geom_line(data=data.frame(ABETA=80:1700, TAU=1/5*(80:1700)), aes(color=TAU/ABETA)) +
  coord_cartesian(xlim=c(200,1700), ylim=c(80,852))


## ----ratio_gradient_logistic------------------------------------------------------------------------
tau_fun <- function(abeta, p){
  # log(p/(1-p)) = y = a + b*abetaz + c*tauz
  # tauz = (y - a - b*abetaz)/c
  abetaz <- (abeta - attr(logistic.fit$model$`scale(ABETA)`, "scaled:center"))/
    attr(logistic.fit$model$`scale(ABETA)`, "scaled:scale")
  tauz <- (log(p/(1-p)) - logistic.fit$coef['(Intercept)']
    - logistic.fit$coef['scale(ABETA)']*abetaz)/
    logistic.fit$coef['scale(TAU)']
  tauz*attr(logistic.fit$model$`scale(TAU)`, "scaled:scale") +
    attr(logistic.fit$model$`scale(TAU)`, "scaled:center")
}
pd <- subset(dd, DX!='MCI') %>%
  mutate(prob = predict(logistic.fit, type='response'))
ggplot(pd, aes(x=ABETA, y=TAU, color=prob)) + 
  geom_point() +
  scale_colour_gradient(low="#0072B2", high="#D55E00") +
  labs(color="Pred. prob. of Dementia") +
  geom_line(data=data.frame(ABETA=80:1700, TAU=tau_fun(80:1700, p=0.01), prob=0.01)) +
  geom_line(data=data.frame(ABETA=80:1700, TAU=tau_fun(80:1700, p=0.1), prob=0.1)) +
  geom_line(data=data.frame(ABETA=80:1700, TAU=tau_fun(80:1700, p=0.25), prob=0.25)) +
  geom_line(data=data.frame(ABETA=80:1700, TAU=tau_fun(80:1700, p=0.5), prob=0.5)) +
  geom_line(data=data.frame(ABETA=80:1700, TAU=tau_fun(80:1700, p=0.75), prob=0.75)) +
  geom_line(data=data.frame(ABETA=80:1700, TAU=tau_fun(80:1700, p=0.9), prob=0.9)) +
  geom_line(data=data.frame(ABETA=80:1700, TAU=tau_fun(80:1700, p=0.99), prob=0.99)) +
  coord_cartesian(xlim=c(200,1700), ylim=c(80,852))


## ----ROC_logistic, fig.width=5, fig.height=5--------------------------------------------------------
roc_logistic.fit <- roc(DX2 ~ predict(logistic.fit, type='response'), subset(dd, DX!='MCI'))

plot(roc_abeta, col='blue')
plot(roc_tau, add=TRUE, col='orange')
plot(roc_tauabeta, add=TRUE, col='red')
plot(roc_logistic.fit, add=TRUE, col='purple')
legend("bottomright", 
  legend=c(expression(paste("A", beta)), "Tau", expression(paste("Tau/A", beta)), 'Logistic model'),
       col=c("blue", "orange", "red", "purple"), lwd=2)



## ----echo=FALSE, eval=FALSE-------------------------------------------------------------------------
# roc.test(roc_abeta, roc_logistic.fit)


## ---------------------------------------------------------------------------------------------------
logistic.fit2 <- glm(DX2 ~ scale(ABETA) + scale(TAU) + scale(I(AGE+Years.bl)) + as.factor(APOE4), family=binomial, data = subset(dd, DX!='MCI'))
summary(logistic.fit2)$coef %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var='Coefficient') %>%
  mutate(`Pr(>|z|)` = format.pval(round(`Pr(>|z|)`, digits = 3), eps = 0.001, digits=3)) %>%
  kable()


## ----eval=FALSE-------------------------------------------------------------------------------------
# roc_logistic.fit2 <- roc(DX2 ~ predict(logistic.fit2, type='response'), subset(dd, DX!='MCI'))
# roc.test(roc_abeta, roc_logistic.fit2)
# roc.test(roc_tauabeta, roc_logistic.fit2)
# auc(roc_logistic.fit2); ci(roc_logistic.fit2)


## ----tree1, fig.height=4.5, fig.width=5*(2)---------------------------------------------------------
tree.fit <- ctree(DX2 ~ ABETA + TAU, 
  data = subset(dd, DX!='MCI'), 
  controls = ctree_control(maxdepth=2))
plot(tree.fit)


## ----tree2------------------------------------------------------------------------------------------
ggplot(subset(dd, DX!='MCI'), aes(x=ABETA, y=TAU)) + 
  geom_point(aes(color=DX)) +
  scale_color_manual(values=c("#0072B2", "#D55E00")) +
  geom_vline(xintercept = 886) +
  #geom_hline(yintercept = 210) +
  #geom_hline(yintercept = 440) +
  geom_segment(aes(x = 0, y = 210, xend = 886, yend = 210)) +
  geom_segment(aes(x = 886, y = 440, xend = 1700, yend = 440))


## ----ROC_rf, fig.width=5, fig.height=5--------------------------------------------------------------
tree.fit2 <- ctree(DX2 ~ ABETA + TAU, data = subset(dd, DX!='MCI'))
tree.fit2.prob.dem <- predict(tree.fit2, type='prob') %>% 
  lapply(function(x) x[2]) %>% unlist()
rf.fit <- cforest(DX2 ~ ABETA + TAU, data = subset(dd, DX!='MCI'))
rf.fit.prob.dem <- predict(rf.fit, type='prob') %>% 
  lapply(function(x) x[2]) %>% unlist()
roc_tree.fit <- roc(DX2 ~ tree.fit2.prob.dem, subset(dd, DX!='MCI'))
roc_rf.fit <- roc(DX2 ~ rf.fit.prob.dem, subset(dd, DX!='MCI'))

plot(roc_abeta, col='blue')
plot(roc_tau, add=TRUE, col='orange')
plot(roc_tauabeta, add=TRUE, col='red')
plot(roc_logistic.fit, add=TRUE, col='purple')
plot(roc_tree.fit, add=TRUE, col='grey')
plot(roc_rf.fit, add=TRUE, col='black')
legend("bottomright", 
  legend=c(expression(paste("A", beta)), "Tau", expression(paste("Tau/A", beta)), 
      'Logistic model', 'Binary Tree', 'Random Forest'),
    col=c("blue", "orange", "red", "purple", "grey", "black"), lwd=2)


## ----echo=FALSE, eval=FALSE-------------------------------------------------------------------------
# roc.test(roc_abeta, roc_tree.fit)
# roc.test(roc_abeta, roc_rf.fit)


## ----density_Abeta, fig.height=2.25, fig.width=3*(2)------------------------------------------------
ggplot(subset(dd, DX!='MCI' & ABETA<1700), aes(x=ABETA)) + 
  geom_histogram(aes(y=..density..), alpha=0.5) +
  geom_density() +
  geom_rug(aes(color=DX))


## ----results='hide', cache=TRUE---------------------------------------------------------------------
#' Plot a Mixture Component
#' 
#' @param x Input data
#' @param mu Mean of component
#' @param sigma Standard deviation of component
#' @param lam Mixture weight of component
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

set.seed(1)
mixmdl <- normalmixEM(subset(dd, DX!='MCI' & ABETA<1700)$ABETA, k = 2)
## number of iterations= 121

set.seed(1)
mvmixmdl <- mvnormalmixEM(subset(dd, DX!='MCI' & ABETA<1700)[, c('ABETA', 'TAU')], k = 2)
## number of iterations= 91


## ----mixture_distribution_Abeta---------------------------------------------------------------------
data.frame(ABETA = mixmdl$x) %>%
  ggplot(aes(x=ABETA)) +
  geom_histogram(aes(y=..density..), alpha=0.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "#D55E00", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "#0072B2", lwd = 1.5) +
  ylab("Density")


## ---------------------------------------------------------------------------------------------------
post.df <- as.data.frame(cbind(x = mixmdl$x, mixmdl$posterior)) %>%
  arrange(x) %>% 
  rename(Abeta = x, `Prob. Abnormal` = comp.1, `Prob. Normal` = comp.2)
post.df %>% filter(Abeta > 1030 & Abeta < 1080) %>% kable()

# A reasonable cutoff might be:
# filter(post.df, `Prob. Abnormal` <= `Prob. Normal`)[1,]


## ----Bivariate_Density------------------------------------------------------------------------------
kd <- with(subset(dd, DX!='MCI' & ABETA<1700)[, c('ABETA', 'TAU')], MASS::kde2d(ABETA, TAU, n = 50))
fig <- plot_ly(x = kd$x, y = kd$y, z = kd$z) %>% add_surface() %>% layout(
  showlegend = FALSE,
  margin = list(l = 0, r = 0, b = 0, t = 0, pad = 0),
  scene = list(
      xaxis = list(title = "Abeta"),
      yaxis = list(title = "Tau"),
      zaxis = list(title = "Density")
    ))

htmlwidgets::saveWidget(fig, "bvdensity_csf_tau.html")


## ----bv_kernel_density------------------------------------------------------------------------------
kd <- with(subset(dd, DX!='MCI' & ABETA<1700)[, c('ABETA', 'TAU')], MASS::kde2d(ABETA, TAU, n = 50))
kdl <- expand.grid(i=1:50, j=1:50) %>%
  mutate(Abeta=NA, Tau=NA, density=NA)
for(r in 1:nrow(kdl)){
  kdl[r, 'Abeta'] <- kd$x[kdl$i[r]]
  kdl[r, 'Tau'] <- kd$y[kdl$j[r]]
  kdl[r, 'density'] <- kd$z[kdl$i[r], kdl$j[r]]
}

ggplot(kdl, aes(Abeta, Tau, z = density)) +
  geom_raster(aes(fill=density))+
  geom_contour(color='white')


## ----mvmix_post_prob, fig.width=5, fig.height=5, out.width='100%'-----------------------------------
post.df2 <- cbind(mvmixmdl$x, mvmixmdl$posterior) %>% 
  as.data.frame() %>% 
  rename(
    `Prob. Abnormal` = comp.1,
    `Prob. Normal` = comp.2
  )

ggplot(post.df2, aes(x=ABETA, y=TAU, color=`Prob. Abnormal`)) + 
  geom_point() +
  scale_colour_gradient(low="#0072B2", high="#D55E00") +
  geom_abline(intercept = 0, slope=0.394) +
  theme(legend.position=c(0.85, 0.82))


## ----mvmix_density, fig.width=5, fig.height=5, out.width='100%'-------------------------------------
plot(mvmixmdl, density = TRUE, alpha = c(0.01, 0.05, 0.10), 
  marginal = FALSE, whichplots=3, main2='', 
  xlab2='ABETA', ylab2='TAU',
  col2=c("#D55E00", "#0072B2"))


## ----echo=FALSE, fig.align='center', out.width='60%'------------------------------------------------
knitr::include_graphics("./images/moe.png")


## ----fig.height=5.9, fig.width=13-------------------------------------------------------------------
moedata <- dd %>%
  filter(DX!='MCI' & ABETA<1700) %>% 
  select(AGE, APOE4, ABETA, TAU, ADAS13) %>%
  mutate(APOE4 = as.factor(APOE4)) %>%
  na.omit()

moe23 <- MoE_clust(
  moedata %>% select(ABETA, TAU, ADAS13),
  G = 2:4,
  gating = ~ APOE4,
  expert = ~ AGE,
  network.data = moedata
)

plot(moe23, what="gpairs", jitter=TRUE)
# with(moedata, plot(moe23, what="gating", x.axis=APOE4, xlab="APOE4"))


## ---------------------------------------------------------------------------------------------------
dd2 <- subset(ADNIMERGE::adnimerge, VISCODE=='bl' & !is.na(Hippocampus) & DX=='CN') %>%
  mutate(
   APOEe4 = ifelse(APOE4>0, 1, 0),
   Sex = factor(PTGENDER, levels = c('Male', 'Female')))


## ----echo=TRUE--------------------------------------------------------------------------------------
lm_fit1_sex <- lm(I(Hippocampus/ICV*100) ~ Sex, data=dd2)


## ---------------------------------------------------------------------------------------------------
summary(lm_fit1_sex)$coef %>%
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var='Coefficient') %>%
  mutate(`Pr(>|t|)` = format.pval(round(`Pr(>|t|)`, digits = 3), eps = 0.001, digits=3)) %>%
  kable()


## ----echo=TRUE--------------------------------------------------------------------------------------
lm_fit2_sex <- lm(Hippocampus ~ ICV + Sex, data=dd2)


## ---------------------------------------------------------------------------------------------------
summary(lm_fit2_sex)$coef %>%
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var='Coefficient') %>%
  mutate(`Pr(>|t|)` = format.pval(round(`Pr(>|t|)`, digits = 3), eps = 0.001, digits=3)) %>%
  kable()


## ----hipp_sex_scatter, fig.height=4, fig.width=6----------------------------------------------------
ggplot(dd2, aes(y=Hippocampus, x=ICV, color=Sex)) +
  geom_point(alpha=0.5) +
  geom_smooth(method='lm') +
  xlab(expression(paste("ICV (", m*m^{3},")"))) +
  ylab(expression(paste("Hippocampus (", m*m^{3},")"))) +
  theme(legend.position = c(0.2, 0.8))


## ----hipp_sex_box, fig.height=4, fig.width=4*2------------------------------------------------------
p1 <- ggplot(dd2, aes(y=Hippocampus, x=Sex)) +
  geom_boxplot(outlier.shape=NA) +
  geom_dotplot(binaxis='y', stackdir='center', 
    dotsize=0.3, alpha=0.2, binwidth=100) +
  xlab('') + ylab(expression(paste("Hippocampus (", m*m^{3},")")))

p2 <- ggplot(dd2, aes(y=ICV, x=Sex)) +
  geom_boxplot(outlier.shape=NA) +
  geom_dotplot(binaxis='y', stackdir='center', 
    dotsize=0.3, alpha=0.2, binwidth=15000) +
  xlab('') + ylab(expression(paste("ICV (", m*m^{3},")")))

p3 <- ggplot(dd2, aes(y=(Hippocampus/ICV*100), x=Sex)) +
  geom_boxplot(outlier.shape=NA) +
  geom_dotplot(binaxis='y', stackdir='center', 
    dotsize=0.3, alpha=0.2, binwidth=0.0075) +
  xlab('') + ylab("Hippocampus (%ICV)")

p4 <- ggplot(dd2, 
  aes(y=Hippocampus-ICV*coef(lm_fit2_sex)['ICV'], x=Sex)) +
  geom_boxplot(outlier.shape=NA) +
  geom_dotplot(binaxis='y', stackdir='center', 
    dotsize=0.3, alpha=0.2, binwidth=100) +
  xlab('') + ylab(expression(paste("Hippocampus (", m*m^{3},", model-adjusted for ICV)")))

grid.arrange(p1,p2,p3,p4, nrow=1)


## ----hipp_adj_sex_scatter, fig.height=4, fig.width=6------------------------------------------------
ggplot(dd2, 
  aes(y=Hippocampus-ICV*coef(lm_fit2_sex)['ICV'], x=ICV, color=Sex)) +
  geom_point(alpha=0.5) +
  geom_smooth(method='lm') +
  xlab(expression(paste("ICV (", m*m^{3},")"))) +
  ylab(expression(paste("Hippocampus (", m*m^{3},", model-adjusted for ICV)"))) +
  theme(legend.position = c(0.2, 0.8))


## ---------------------------------------------------------------------------------------------------
ggplot() + 
  geom_segment(aes(x=1.009, y=0, xend=2.076, yend=100),
    arrow = arrow(length = unit(0.03, "npc"), ends='both')) +
  geom_abline(slope = 100/(2.076-1.009), 
    intercept=-1.009*100/(2.076-1.009), linetype='dashed') +
  geom_point(aes(x=1.009, y=0)) +
  geom_point(aes(x=2.076, y=100)) +
  xlab('PiB SUVr') + ylab('Centiloid') +
  xlim(1, 3) + ylim(-5, 200)
   






## ----echo=FALSE-------------------------------------------------------------------------------------
pibids <- unique(subset(ADNIMERGE::adnimerge, !is.na(PIB))$RID)
av45ids <- unique(subset(ADNIMERGE::adnimerge, !is.na(AV45))$RID)
set.seed(20210506)
holdout <- sample(intersect(pibids, av45ids), size = 10)

dd <- ADNIMERGE::adnimerge %>%
  filter(!RID %in% holdout) %>%
  arrange(RID, EXAMDATE) %>%
  select(RID, DX, PIB, AV45) %>%
  rename(PiB = PIB, Florbetapir = AV45) %>%
  mutate(PiB = as.numeric(PiB), Florbetapir = as.numeric(Florbetapir)) %>%
  group_by(RID) %>%
  fill(DX, .direction = "downup") %>%
  filter(!is.na(PiB) | !is.na(Florbetapir)) %>%
  filter(!duplicated(RID)) %>%
  pivot_longer(cols=c('PiB', 'Florbetapir'), names_to='Tracer', values_to='SUVR') %>%
  mutate(Tracer = factor(Tracer, levels = c('PiB', 'Florbetapir'))) %>%
  filter(!is.na(SUVR))
  
ggplot(dd, aes(x=SUVR)) +
  stat_ecdf(geom = "step") +
  facet_grid(.~Tracer, scales='free_x')


## ----echo=FALSE-------------------------------------------------------------------------------------
t1 <- with(dd, table(DX, Tracer))
t2 <- round(with(dd, prop.table(table(DX, Tracer), margin = 2))*100, 1)
tt <- t1
tt[,1] <- paste(t1[,1], paste0("(", t2[,1], "%)"))
tt[,2] <- paste(t1[,2], paste0("(", t2[,2], "%)"))
kable(tt)


## ---------------------------------------------------------------------------------------------------
invproptab <- 1/with(dd, prop.table(table(DX, Tracer), margin = 2))
kable(invproptab)


## ----echo=TRUE--------------------------------------------------------------------------------------
# Record the sampling adjustment weights in the data
dd <- dd %>% mutate(
  wt = case_when(
    DX == 'CN' & Tracer == 'PiB' ~ invproptab['CN', 'PiB'],
    DX == 'MCI' & Tracer == 'PiB' ~ invproptab['MCI', 'PiB'],
    DX == 'Dementia' & Tracer == 'PiB' ~ invproptab['Dementia', 'PiB'],
    DX == 'CN' & Tracer == 'Florbetapir' ~ invproptab['CN', 'Florbetapir'],
    DX == 'MCI' & Tracer == 'Florbetapir' ~ invproptab['MCI', 'Florbetapir'],
    DX == 'Dementia' & Tracer == 'Florbetapir' ~ invproptab['Dementia', 'Florbetapir']
  ))


## ----echo=FALSE-------------------------------------------------------------------------------------
ecdf.func <- function(x, weights=NULL, 
  type=c('i/n','(i-1)/(n-1)','i/(n+1)'), 
  normwt=TRUE, na.rm=TRUE,
  method = "linear", rule = 2, ...)
{
 type <- match.arg(type)
 ecdf.data <- Hmisc::wtd.Ecdf(x, weights=weights, type=type, 
  normwt=normwt, na.rm=na.rm)
 with(ecdf.data, approxfun(x, ecdf,
  na.rm=na.rm, rule = rule, ...))
}

inv.ecdf <- function(x, weights = NULL, probs = seq(0,1,by=0.01), 
  type = c("quantile", "(i-1)/(n-1)", "i/(n+1)", "i/n"), normwt = TRUE, 
  na.rm = TRUE, rule = 2, ...)
{
 type <- match.arg(type)
 inv.ecdf.data <- as.numeric(
    Hmisc::wtd.quantile(x, weights=weights, probs=probs, type=type,
     normwt=normwt, na.rm=na.rm))

 approxfun(probs, inv.ecdf.data, rule=rule, na.rm=na.rm, ...)
}


## ----echo=TRUE--------------------------------------------------------------------------------------
# Create adjusted ECDF functions (mapping SUVRs to Cumulative Probabilities)
# Hmisc::wtd.Ecdf returns a data.frame evaluating the ECDF at each observed value
PiB.ecdf <- with(subset(dd, Tracer == 'PiB'), ecdf.func(SUVR, weights=wt))
Fbp.ecdf <- with(subset(dd, Tracer == 'Florbetapir'), ecdf.func(SUVR, weights=wt))

# Create adjusted **inverse** ECDF functions 
# mapping Cumulative Probabilities (0 to 1 scale) to SUVRs
PiB.inv.ecdf <- with(subset(dd, Tracer == 'PiB'), 
 inv.ecdf(SUVR, weights=wt))
Fbp.inv.ecdf <- with(subset(dd, Tracer == 'Florbetapir'),  
 inv.ecdf(SUVR, weights=wt))


## ----echo=TRUE--------------------------------------------------------------------------------------
dd <- dd %>% mutate(
  `Adjusted cumulative probability` = case_when( # 
    Tracer == 'PiB' ~ PiB.ecdf(SUVR),
    Tracer == 'Florbetapir' ~ Fbp.ecdf(SUVR)),
  `Adjusted Z-score` = qnorm(`Adjusted cumulative probability`), # adjusted z-scores
  `Florbetapir to PiB adjusted SUVR` = case_when(
    Tracer == 'Florbetapir' ~ PiB.inv.ecdf(Fbp.ecdf(SUVR))),
  `PiB to Florbetapir adjusted SUVR` = case_when(
    Tracer == 'PiB' ~ Fbp.inv.ecdf(PiB.ecdf(SUVR))),
  CL = case_when(
    Tracer == 'PiB' ~ 100*(SUVR - 1.009)/1.067)) %>%
  arrange(Tracer, SUVR)


## ----weighted-ecdfs---------------------------------------------------------------------------------
ggplot(dd, aes(x=SUVR)) +
  stat_ecdf(geom = "step") +
  facet_grid(.~Tracer, scales='free_x') +
  geom_line(aes(x=SUVR, y=`Adjusted cumulative probability`), linetype='dashed') +
  ylab('Cumulative Probability')


## ---------------------------------------------------------------------------------------------------
dd %>% filter(Tracer == 'PiB') %>%
  ggplot(aes(x=CL, y=`Adjusted cumulative probability`)) +
  geom_line() +
  geom_point()


## ----pib-densities----------------------------------------------------------------------------------
dd %>% select(RID, DX, Tracer, SUVR, `Florbetapir to PiB adjusted SUVR`) %>%
  pivot_longer(c('SUVR', 'Florbetapir to PiB adjusted SUVR'), 
    names_to = 'Source', values_to = 'SUVR') %>%
  filter( (Tracer == 'PiB' & Source == 'SUVR') |
      (Tracer == 'Florbetapir' & Source == 'Florbetapir to PiB adjusted SUVR')) %>%
  mutate(Source = ifelse(Source == 'SUVR', 
    'Directly observed PiB SUVR', 'ECDF mapped Florbetapir to PiB SUVR')) %>%
ggplot(aes(x=SUVR, color=Source)) +
  geom_density(alpha=0.5)


## ----pib-densities-dx-------------------------------------------------------------------------------
dd %>% select(RID, DX, Tracer, SUVR, `Florbetapir to PiB adjusted SUVR`) %>%
  pivot_longer(c('SUVR', 'Florbetapir to PiB adjusted SUVR'), 
    names_to = 'Source', values_to = 'SUVR') %>%
  filter( (Tracer == 'PiB' & Source == 'SUVR') |
      (Tracer == 'Florbetapir' & Source == 'Florbetapir to PiB adjusted SUVR')) %>%
  mutate( Source = ifelse(Source == 'SUVR', 'Directly observed PiB SUVR', 'ECDF mapped Florbetapir to PiB SUVR')) %>%
  filter(!is.na(DX)) %>%
ggplot(aes(x=SUVR, color=Source)) +
  geom_density(alpha=0.5) +
  facet_grid(.~DX)


## ----z-score-densities------------------------------------------------------------------------------
ggplot() +
  geom_rug(data = subset(dd, Tracer == 'Florbetapir'), 
    aes(x=`Adjusted Z-score`, color = Tracer), alpha=0.5) +
  geom_rug(data = subset(dd, Tracer == 'PiB'), 
    aes(x=`Adjusted Z-score`, color = Tracer), alpha=0.5) +
  geom_density(data = subset(dd, Tracer == 'Florbetapir'), 
    aes(x=`Adjusted Z-score`, color = Tracer), alpha=0.5) +
  geom_density(data = subset(dd, Tracer == 'PiB'), aes(x=`Adjusted Z-score`, color = Tracer))




## ---------------------------------------------------------------------------------------------------
dd.validate <- full_join(
  ADNIMERGE::adnimerge %>% 
    arrange(RID, EXAMDATE) %>% 
    select(RID, EXAMDATE, PIB) %>%
    rename(PiB = PIB) %>%
    filter(RID %in% holdout) %>%
    filter(!is.na(PiB)) %>%
    filter(!duplicated(RID)),
  ADNIMERGE::adnimerge %>% 
    arrange(RID, EXAMDATE) %>% 
    select(RID, EXAMDATE, AV45) %>%
    rename(Florbetapir = AV45) %>%
    filter(RID %in% holdout) %>%
    filter(!is.na(Florbetapir)) %>%
    filter(!duplicated(RID)), by='RID', suffix = c('.PiB', '.Fbp')
  ) %>%
  mutate(
    `Years between scans` = as.numeric(EXAMDATE.Fbp - EXAMDATE.PiB)/365.25,
    `Royse et al linear map of PiB to Florbetapir` = 
      as.numeric(0.497 * PiB + 0.503),
    `Royse et al linear map of Florbetapir to PiB` = 
      as.numeric((Florbetapir - 0.503)/0.497),
    `Navitsky et al linear map of PiB to Florbetapir` = 
      as.numeric(0.536 * PiB + 0.502),
    `Navitsky et al linear map of Florbetapir to PiB` = 
      as.numeric((Florbetapir - 0.502)/0.536),
    `ECDF map of Florbetapir to PiB` = PiB.inv.ecdf(Fbp.ecdf(Florbetapir)),
    `ECDF map of PiB to Florbetapir` = Fbp.inv.ecdf(PiB.ecdf(PiB))
  )

fbb2pib_navitsky <- dd.validate %>%
  select(RID, `Years between scans`, PiB, 
    `Navitsky et al linear map of Florbetapir to PiB`,
    `ECDF map of Florbetapir to PiB`) %>%
  rename(
    `Navitsky et al linear map*` = `Navitsky et al linear map of Florbetapir to PiB`,
    `ECDF map` = `ECDF map of Florbetapir to PiB`) %>%
  pivot_longer(c('Navitsky et al linear map*', 'ECDF map'), 
    names_to = 'Method', values_to = 'Estimated PiB SUVR')

fbb2pib_navitsky_RMSE <- fbb2pib_navitsky %>% group_by(Method) %>%
  summarise(RMSE = sqrt(mean((PiB-`Estimated PiB SUVR`)^2))) %>%
  mutate(PiB = 2, `Estimated PiB SUVR` = 0.75,
    text = paste0('RMSE=', round(RMSE, digits=3)))

ggplot(fbb2pib_navitsky, aes(x=PiB, y=`Estimated PiB SUVR`)) +
  geom_point(aes(color=`Years between scans`)) +
  geom_abline(intercept = 0, slope=1) +
  facet_grid(.~Method) +
  xlab('Actual PiB SUVR') +
  geom_text(data=fbb2pib_navitsky_RMSE, aes(label=text)) +
  geom_hline(yintercept = 1.5, linetype='dashed') +
  geom_vline(xintercept = 1.5, linetype='dashed')


## ---------------------------------------------------------------------------------------------------
pib2fbb_navitsky <- dd.validate %>%
  select(RID, `Years between scans`, Florbetapir, 
    `Navitsky et al linear map of PiB to Florbetapir`,
    `ECDF map of PiB to Florbetapir`) %>%
  rename(
    `Navitsky et al linear map*` = `Navitsky et al linear map of PiB to Florbetapir`,
    `ECDF map` = `ECDF map of PiB to Florbetapir`) %>%
  pivot_longer(c('Navitsky et al linear map*', 'ECDF map'), 
    names_to = 'Method', values_to = 'Estimated Florbetapir SUVR')

pib2fbb_navitsky_RMSE <- pib2fbb_navitsky %>% group_by(Method) %>%
  summarise(RMSE = sqrt(mean((Florbetapir-`Estimated Florbetapir SUVR`)^2))) %>%
  mutate(Florbetapir = 1.5, `Estimated Florbetapir SUVR` = 0.75,
    text = paste0('RMSE=', round(RMSE, digits=3)))

ggplot(pib2fbb_navitsky, aes(x=Florbetapir, y=`Estimated Florbetapir SUVR`)) +
  geom_point(aes(color=`Years between scans`)) +
  geom_abline(intercept = 0, slope=1) +
  facet_grid(.~Method) +
  xlab('Actual Florbetapir SUVR') +
  geom_text(data=pib2fbb_navitsky_RMSE, aes(label=text)) +
  geom_hline(yintercept = 1.11, linetype='dashed') +
  geom_vline(xintercept = 1.11, linetype='dashed')


## ---------------------------------------------------------------------------------------------------
fbb2pib_royse <- dd.validate %>%
  select(RID, `Years between scans`, PiB, 
    `Royse et al linear map of Florbetapir to PiB`,
    `ECDF map of Florbetapir to PiB`) %>%
  rename(
    `Royse et al linear map*` = `Royse et al linear map of Florbetapir to PiB`,
    `ECDF map` = `ECDF map of Florbetapir to PiB`) %>%
  pivot_longer(c('Royse et al linear map*', 'ECDF map'), 
    names_to = 'Method', values_to = 'Estimated PiB SUVR')

fbb2pib_royse_RMSE <- fbb2pib_royse %>% group_by(Method) %>%
  summarise(RMSE = sqrt(mean((PiB-`Estimated PiB SUVR`)^2))) %>%
  mutate(PiB = 2, `Estimated PiB SUVR` = 0.75,
    text = paste0('RMSE=', round(RMSE, digits=3)))

ggplot(fbb2pib_royse, aes(x=PiB, y=`Estimated PiB SUVR`)) +
  geom_point(aes(color=`Years between scans`)) +
  geom_abline(intercept = 0, slope=1) +
  facet_grid(.~Method) +
  xlab('Actual PiB SUVR') +
  geom_text(data=fbb2pib_royse_RMSE, aes(label=text)) +
  geom_hline(yintercept = 1.5, linetype='dashed') +
  geom_vline(xintercept = 1.5, linetype='dashed')


## ---------------------------------------------------------------------------------------------------
pib2fbb_royse <- dd.validate %>%
  select(RID, `Years between scans`, Florbetapir, 
    `Royse et al linear map of PiB to Florbetapir`,
    `ECDF map of PiB to Florbetapir`) %>%
  rename(
    `Royse et al linear map*` = `Royse et al linear map of PiB to Florbetapir`,
    `ECDF map` = `ECDF map of PiB to Florbetapir`) %>%
  pivot_longer(c('Royse et al linear map*', 'ECDF map'), 
    names_to = 'Method', values_to = 'Estimated Florbetapir SUVR')

pib2fbb_royse_RMSE <- pib2fbb_royse %>% group_by(Method) %>%
  summarise(RMSE = sqrt(mean((Florbetapir-`Estimated Florbetapir SUVR`)^2))) %>%
  mutate(Florbetapir = 1.5, `Estimated Florbetapir SUVR` = 0.75,
    text = paste0('RMSE=', round(RMSE, digits=3)))

ggplot(pib2fbb_royse, aes(x=Florbetapir, y=`Estimated Florbetapir SUVR`)) +
  geom_point(aes(color=`Years between scans`)) +
  geom_abline(intercept = 0, slope=1) +
  facet_grid(.~Method) +
  xlab('Actual Florbetapir SUVR') +
  geom_text(data=pib2fbb_royse_RMSE, aes(label=text)) +
  geom_hline(yintercept = 1.11, linetype='dashed') +
  geom_vline(xintercept = 1.11, linetype='dashed')






## ----prep-data-for-combat, echo = FALSE-------------------------------------------------------------
dd.combat <- ADNIMERGE::adnimerge %>%
  arrange(RID, EXAMDATE) %>%
  mutate(Age = AGE + Years.bl) %>%
  select(RID, APOE4, Age, DX, MMSE, PIB, AV45, FBB) %>%
  group_by(RID) %>%
  fill(APOE4, DX, Age, MMSE, .direction = "downup") %>%
  mutate(across(PIB:FBB, as.numeric)) %>%
  pivot_longer(PIB:FBB, names_to = 'Tracer', values_to = 'SUVR') %>%
  filter(!is.na(SUVR) & !is.na(APOE4) & !is.na(DX) & !is.na(Age)) %>%
  mutate(
    log.SUVR = log(SUVR),
    Tracer = factor(Tracer, levels = c('PIB', 'AV45', 'FBB')),
    APOE4 = factor(APOE4, levels = 0:2))


## ----fit-model-for-combat---------------------------------------------------------------------------
combat.fit <- lme(log.SUVR ~ Tracer + APOE4 + DX + Age + MMSE, dd.combat,
  random = ~1|RID, weights = varIdent(form = ~ 1 | Tracer))


## ----echo = FALSE-----------------------------------------------------------------------------------
printCoefmat(summary(combat.fit)$tTable)


## ----echo = FALSE-----------------------------------------------------------------------------------
print(summary(combat.fit$modelStruct), digits = 4)


## ----echo = TRUE------------------------------------------------------------------------------------
dd.combat$fitted_log.SUVR <- predict(combat.fit) %>% as.numeric()

dd.combat$pred_log.SUVR_PiB <- predict(combat.fit, 
  newdata = dd.combat %>% mutate(Tracer = "PIB")) %>% as.numeric()


## ----echo = TRUE------------------------------------------------------------------------------------
(sdWeights <- coef(combat.fit$modelStruct$varStruct, unconstrained = FALSE))
dd.combat <- dd.combat %>% 
  mutate(
    resids = log.SUVR - fitted_log.SUVR,
    homog_resids = case_when(
      Tracer == 'PIB' ~ resids,
      Tracer == 'AV45' ~ resids * sdWeights['PIB'] ,
      Tracer == 'FBB' ~ resids * sdWeights['PIB'] / sdWeights['FBB']),
    ComBat_log.SUVR_PiB = homog_resids + pred_log.SUVR_PiB,
    fitted_SUVR = exp(fitted_log.SUVR),
    ComBat_SUVR_PiB = exp(ComBat_log.SUVR_PiB))


## ----echo = FALSE-----------------------------------------------------------------------------------
dd.combat %>% group_by(Tracer) %>%
  summarise(
    `SD of homog_resids` = sd(homog_resids),
    `SD of resids` = sd(resids))


## ---------------------------------------------------------------------------------------------------
ggplot(dd.combat, aes(x = SUVR, y =ComBat_SUVR_PiB, color = Tracer)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  xlab('Original SUVR') + ylab('ComBat SUVR') +
  scale_x_continuous(transform = 'log') +
  scale_y_continuous(transform = 'log')


## ----fig.height=5, fig.width=5*2.2------------------------------------------------------------------
p1 <- ggplot(dd.combat, aes(x = SUVR)) + 
  geom_density(aes(fill = Tracer, color = Tracer), alpha = 0.2) +
  scale_x_continuous(transform = 'log', limits = c(0.8,3)) +
  theme(legend.position = 'inside', legend.position.inside = c(0.9, 0.7)) +
  xlab('Original SUVR')
p2 <- ggplot(dd.combat, aes(x = ComBat_SUVR_PiB)) + 
  geom_density(aes(fill = Tracer, color = Tracer), alpha = 0.2) +
  scale_x_continuous(transform = 'log', limits = c(0.8,3)) +
  theme(legend.position = 'none') + 
  xlab('ComBat SUVR')
grid.arrange(p1,p2,nrow=2)

