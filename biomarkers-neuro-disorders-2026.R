## ----setup----------------------------------------------------------------------------------------------------------------
#| include: false

library(A4LEARN)
library(tidyverse)
library(arsenal)
library(kableExtra)
library(nlme)
library(emmeans)
library(splines)
library(clubSandwich)
library(broom)
library(mixtools)
library(patchwork)
# for sva::ComBat batch correction
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("sva")
library(sva)
library(reactable)
library(mmrm)
library(mice)

options(digits = 3)
theme_set(theme_bw())

cbbPalette <- c(
  "#E69F00", "#56B4E9", "#CC79A7",
  "#009E73", "#F0E442", "#999999",
  "#0072B2", "#D55E00", "#000000"
)

grpPalette <- c(LEARN = "#009E73", Placebo = "#0072B2", Solanezumab = "#D55E00")

scale_colour_discrete <- function(...) scale_colour_manual(..., values = cbbPalette)
scale_fill_discrete   <- function(...) scale_fill_manual(..., values = cbbPalette)

formatp <- function(x) case_when(
  x < 0.001 ~ "p<0.001",
  TRUE      ~ format.pval(x, digits = 3, eps = 0.001, nsmall = 3)
)

mean_ci <- function(x) {
  out <- Hmisc::smean.cl.normal(x)
  tibble(
    mean = unname(out[["Mean"]]),
    lower95 = unname(out[["Lower"]]),
    upper95 = unname(out[["Upper"]])
  )
}


## ----organize-a4learn-bl-data---------------------------------------------------------------------------------------------
#| echo: true

# Outcomes collected at Visit 1
V1OUTCOME <- A4LEARN::ADQS |>
  filter(VISITCD == "001", SUBSTUDY != "SF") |>
  select(BID, QSTESTCD, QSSTRESN) |>
  tidyr::pivot_wider(values_from = QSSTRESN, names_from = QSTESTCD)

# Outcomes collected at Visit 6
V6OUTCOME <- A4LEARN::ADQS |>
  filter(VISITCD == "006", SUBSTUDY != "SF") |>
  select(BID, QSTESTCD, QSSTRESN) |>
  tidyr::pivot_wider(values_from = QSSTRESN, names_from = QSTESTCD)

# ptau217 collected at Visit 1 or 6
V1.6PTAU217 <- A4LEARN::biomarker_pTau217 |>
    filter(TESTCD == "PTAU217", !is.na(ORRESRAW)) |>
    arrange(BID, VISCODE) |>
    filter(!duplicated(BID)) |>
    filter(VISCODE %in% c(1,6)) |>
    select(BID, PTAU217 = ORRESRAW)

# Collect baseline covariates, outcomes, ptau217
SUBJINFO <- A4LEARN::SUBJINFO |>
  filter(SUBSTUDY != "SF") |> 
  left_join(V6OUTCOME, by = "BID") |>
  left_join(V1OUTCOME |> 
    select(BID, CDRSB, CFITOTAL, CFISP, CFIPT, ADLPQPT, ADLPQSP),
    by = "BID") |> 
  left_join(V1.6PTAU217, by = "BID")


## ----organize-a4learn-long-data-------------------------------------------------------------------------------------------
#| echo: true

ADQS_PACC <- A4LEARN::ADQS |>
  filter(EPOCH != "SCREENING" | AVISIT == '006') |>
  filter(QSTESTCD == "PACC", !is.na(QSSTRESN)) |>
  rename(PACC = QSSTRESN, PACC_bl = QSBLRES, PACC_ch = QSCHANGE) |>
  select(BID, ASEQNCS, ASEQMMRM, SUBSTUDY, EPOCH, MITTFL, AVISIT, TX, ADURW, AGEYR, 
    AAPOEGNPRSNFLG, EDCCNTU, AMYLCENT, QSVERSION, PACC, PACC_bl, PACC_ch) |>
  mutate(
    Group = case_when(
      SUBSTUDY == "LEARN" ~ "LEARN",
      TRUE ~ TX) |> factor(levels = c("LEARN", "Placebo", "Solanezumab")),
    TX = factor(TX, levels = c("Placebo", "Solanezumab")),
    APOE4 = factor(AAPOEGNPRSNFLG),
    EDUC = EDCCNTU,
    VERSION = factor(QSVERSION),
    Week = case_when(
      as.numeric(AVISIT) < 200 ~ as.numeric(AVISIT)*4 - 24,
      TRUE ~ NA)) |> 
  filter(!is.na(Group)) |>
  left_join(V1.6PTAU217, by = "BID")

BASE_PACC <- ADQS_PACC |>
  arrange(BID, ADURW, ASEQNCS) |>
  group_by(BID) |>
  slice(1) |>
  ungroup()

ALL_BASE_PET <- SUBJINFO |>
  filter(!is.na(AMYLCENT))


## ----sim-batch-effect, results = 'hide', fig.show='hide'------------------------------------------------------------------
#| echo: true
set.seed(20200225)

batch_data <- 
  tibble(
    batch = 1:10, # simulate random batch-specific mean and variance
    Sigma = rgamma(n = 10, shape = 360 / 10, scale = 10),
    Mean  = rnorm(n = 10, mean = 850, sd = 200)) |>
  group_by(batch) |>
  nest() |> # creates "data" column with Mean and Sigma for each batch
  # simulate biomarker data per batch  
  mutate(Biomarker = purrr::map(data, ~rnorm(n = 50, .x$Mean, .x$Sigma))) |>
  unnest(Biomarker) |>
  unnest(data) |>
  ungroup() |>
  arrange(batch) |>
  mutate(batch = factor(batch))

ggplot(batch_data, aes(x = batch, y = Biomarker)) +
  geom_boxplot(outlier.shape = NA) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.3, alpha = 0.2)


## ----sim-batch-effect-fig-------------------------------------------------------------------------------------------------
#| echo: false
ggplot(batch_data, aes(x = batch, y = Biomarker)) +
  geom_boxplot(outlier.shape = NA) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.3, alpha = 0.2)


## ----anova-sim-batch-effect-----------------------------------------------------------------------------------------------
#| echo: true
anova(lm(Biomarker ~ batch, data = batch_data)) |>
  broom::tidy() |>
  mutate(p.value = formatp(p.value)) |>
  kable(digits = 3)


## ----gls-test-sim-batch-effect--------------------------------------------------------------------------------------------
#| echo: true
# Model 1 (Null): Assumes equal variance across all batches
equal_var <- gls(Biomarker ~ batch, data = batch_data)

# Model 2 (Alternative): Allows different variances for each batch
different_var <- gls(Biomarker ~ batch, data = batch_data,
  weights = varIdent(form = ~1 | batch))

# Perform the Likelihood Ratio Test
anova(equal_var, different_var) |>
  as_tibble(rownames = "model") |>
  mutate(`p-value` = formatp(`p-value`)) |>
  select(-call, -Model) |> 
  kable(digits = 3)


## ----combat---------------------------------------------------------------------------------------------------------------
#| echo: true
set.seed(20260603)
batch_data$Biomarker_ComBat <- sva::ComBat(
  dat = rbind(batch_data$Biomarker, jitter(batch_data$Biomarker, factor=1e-9)),
  batch = batch_data$batch,
  par.prior = TRUE,
  prior.plots = FALSE
)[1,]


## ----test-after-combat----------------------------------------------------------------------------------------------------
anova(lm(Biomarker_ComBat ~ batch, data = batch_data)) |>
  broom::tidy() |>
  mutate(p.value = formatp(p.value)) |>
  kable(digits = 3)

# Model 1 (Null): Assumes equal variance across all batches
fit_combat_equal_var <- gls(Biomarker_ComBat ~ batch, data = batch_data)

# Model 2 (Alternative): Allows different variances for each batch
fit_combat_diff_var <- gls(Biomarker_ComBat ~ batch, data = batch_data,
  weights = varIdent(form = ~1 | batch))

# Perform the Likelihood Ratio Test
anova(fit_combat_equal_var, fit_combat_diff_var) |>
  as_tibble(rownames = "model") |>
  mutate(`p-value` = formatp(`p-value`)) |>
  select(-call, -Model) |> 
  kable(digits = 3)


## ----combat-fig, echo=FALSE-----------------------------------------------------------------------------------------------
ggplot(batch_data, aes(x = batch, y = Biomarker_ComBat)) +
  geom_boxplot(outlier.shape = NA) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.3, alpha = 0.2) +
  ggtitle("ComBat adjusted biomarker by batch")


## ----a4learn-bl-char1-----------------------------------------------------------------------------------------------------
#| results: asis

A4labels <- c(
  AGEYR = "Age (years)",
  APOE4 = "APOEε4 carrier",
  EDUC = "Education (years)",
  PTAU217 = "Baseline pTau217 (U/ml)",
  AMYLCENT = "Baseline amyloid PET (CL)",
  PACC = "Baseline PACC"
)

arsenal::tableby(
  Group ~ AGEYR+ EDUC + PACC,
  data = BASE_PACC) |> 
summary(labelTranslations = A4labels, text = TRUE, digits = 1,
  numeric.stats = "meansd", cat.stats = "countpct")


## ----a4learn-bl-char2-----------------------------------------------------------------------------------------------------
#| results: asis

arsenal::tableby(
  Group ~ AMYLCENT + APOE4  + PTAU217,
  data = BASE_PACC) |> 
summary(labelTranslations = A4labels, text = TRUE, digits = 1,
  numeric.stats = "meansd", cat.stats = "countpct")


## ----pacc-box-violin-fig--------------------------------------------------------------------------------------------------
ggplot(BASE_PACC, aes(x = Group, y = PACC, fill = Group)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(fill = NA, alpha = 0.6, outlier.shape = NA) +
  ggbeeswarm::geom_beeswarm(alpha = 0.2, size = 1) +
  labs(x = "", y = "Baseline PACC") +
  theme(legend.position = "none") +
  scale_fill_manual(values = grpPalette)


## -------------------------------------------------------------------------------------------------------------------------
ggplot(BASE_PACC, aes(x = AMYLCENT, y = PACC, color = Group)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Baseline amyloid PET biomarker (CL)",
    y = "Baseline PACC") +
  scale_colour_manual(values = grpPalette)


## ----pacc-bl-model--------------------------------------------------------------------------------------------------------
fit_bl <- lm(PACC ~ scale(AMYLCENT) + scale(PTAU217) + scale(AGEYR) + APOE4 + 
  scale(EDUC), data = BASE_PACC)

broom::tidy(fit_bl) |>
  mutate(p.value = formatp(p.value)) |>
  kable(digits = 3)


## -------------------------------------------------------------------------------------------------------------------------
ggplot(ALL_BASE_PET, aes(x = SUVRCER)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "grey75", color = "white") +
  geom_density(color = cbbPalette[1], linewidth = 1.2) +
  labs(x = "All screening amyloid PET (SUVR)", y = "Density")


## ----ecdf-----------------------------------------------------------------------------------------------------------------
#| echo: true
ecdf_suvr <- ecdf(ALL_BASE_PET$SUVRCER)

ALL_BASE_PET <- ALL_BASE_PET |>
  mutate(
    SUVRCER_p = ecdf_suvr(SUVRCER), # 0-1 scale
    SUVRCER_z = qnorm(pmin(pmax(SUVRCER_p, 1e-6), 1 - 1e-6)),
    SUVRCER_pct = SUVRCER_p*100) # 0-100 scale


## ----pct-z-scores---------------------------------------------------------------------------------------------------------
#| echo: false
p1 <- ggplot(ALL_BASE_PET, aes(x = SUVRCER)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "grey75", color = "white") +
  geom_density(color = cbbPalette[1], linewidth = 1.2) +
  labs(x = "Amyloid PET (SUVR)", y = "Density")

p2 <- ggplot(ALL_BASE_PET |> arrange(SUVRCER), aes(x = SUVRCER, y = SUVRCER_pct*100)) +
  geom_line(color = cbbPalette[2], size = 2) +
  labs(x = "Amyloid PET (SUVR)", y = "Percentile (%)")

p3 <- ggplot(ALL_BASE_PET |> arrange(SUVRCER), aes(x = SUVRCER, y = SUVRCER_z)) +
  geom_line(color = cbbPalette[2], size = 2) +
  labs(x = "Amyloid PET (SUVR)", y = "Z-score")

p1 + p2 + p3


## ----pct-z-scores2--------------------------------------------------------------------------------------------------------
#| echo: false

p1 <- ggplot(ALL_BASE_PET, aes(x = SUVRCER_z)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "grey75", color = "white") +
  geom_density(color = cbbPalette[1], linewidth = 1.2) +
  labs(x = "ECDF-derived z-score for amyloid PET", y = "Density")

p2 <- ggplot(ALL_BASE_PET, aes(x = scale(SUVRCER))) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "grey75", color = "white") +
  geom_density(color = cbbPalette[1], linewidth = 1.2) +
  labs(x = "Usual z-score for amyloid PET", y = "Density")

p1 + p2


## ----pct-cl---------------------------------------------------------------------------------------------------------------
ggplot(ALL_BASE_PET, aes(x = AMYLCENT, y = SUVRCER_pct)) +
  geom_point() +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = 100, color = "red", linetype = "dashed", size = 1) +
  xlab("Amyloid PET (CL)") +
  ylab("Amyloid PET (Percentile)")


## ----mixture, include = FALSE---------------------------------------------------------------------------------------------
mix_centiloid <- mixtools::normalmixEM(ALL_BASE_PET$AMYLCENT, k = 2)

plot_mix_comp <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

mix_post <- mix_centiloid$posterior |>
  as_tibble() |>
  mutate(
    AMYLCENT = mix_centiloid$x,
    p1 = plot_mix_comp(AMYLCENT, 
      mu = mix_centiloid$mu[1], 
      sigma = mix_centiloid$sigma[1], 
      lam = mix_centiloid$lambda[1]),
    p2 = plot_mix_comp(AMYLCENT, 
      mu = mix_centiloid$mu[2], 
      sigma = mix_centiloid$sigma[2], 
      lam = mix_centiloid$lambda[2])) |>
  arrange(AMYLCENT)

mix_threshold <- with(mix_post |> filter(AMYLCENT>0, AMYLCENT<50), 
  approx(x = p1-p2, y = AMYLCENT, xout = 0))$y


## ----mixture-fig----------------------------------------------------------------------------------------------------------
ggplot(mix_post, aes(x=AMYLCENT)) +
  geom_histogram(aes(y = after_stat(density)), bins = 35, fill = "grey80") +
  geom_line(aes(y = p1), color = cbbPalette[2], linewidth = 1.2) +
  geom_line(aes(y = p2), color = cbbPalette[1], linewidth = 1.2) +
  geom_vline(xintercept = mix_threshold, color = "red", linetype = "dashed", size = 1) +
  annotate(
    "text",
    x = mix_threshold + 5, y = max(mix_post$p1) * 0.8,
    label = paste0("Mixture model threshold: ", round(mix_threshold, 1), " CL"),
    color = "red", size = 4, hjust = 0) +
  labs(x = "Amyloid PET (CL)", y = "Density")


## -------------------------------------------------------------------------------------------------------------------------
ggplot(ADQS_PACC, aes(x = ADURW, y = PACC, color = Group)) +
  geom_line(aes(group = BID), alpha = 0.12) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1.2) +
  labs(
    x = "Weeks from baseline",
    y = "PACC") +
  scale_color_manual(values = grpPalette)


## -------------------------------------------------------------------------------------------------------------------------
PACC_SUM <- ADQS_PACC |>
  filter(!is.na(Week), Week %in% seq(0, 456, 24)) |>
  group_by(Week, Group) |>
  summarise(
    n = n(),
    sd = sd(PACC),
    mean_ci(PACC),
    .groups = "drop")

PACC_SUM |>
  filter(Week %in% seq(0, 240, 48)) |>
  arrange(Week, Group) |>
  mutate(across(c(sd, mean, lower95, upper95), ~ round(., 2))) |>
  reactable()


## -------------------------------------------------------------------------------------------------------------------------
p1 <- ggplot(PACC_SUM, aes(x = Week, y = mean, color = Group, group = Group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = lower95, ymax = upper95), width = 0) +
  labs(
    x = "Weeks from baseline",
    y = "Mean PACC (95% CI)") +
  scale_color_manual(values = grpPalette)

t1 <- ggplot(PACC_SUM, aes(x = Week, y = Group, label = n, color = Group)) +
  geom_text(size = 3.5, show.legend = FALSE) +
  labs(x = "", y = "") +
  scale_color_manual(values = grpPalette) +
  theme_classic() +
  theme(
    # Strip away all the typical plot elements to make it look like a table
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank()) # Hide x-axis text (we rely on the main plot for week numbers)

p1 / t1 + 
  plot_layout(heights = c(4, 1))


## ----fit-lme--------------------------------------------------------------------------------------------------------------
fit_lme <- nlme::lme(
  PACC ~ splines::ns(ADURW, df = 2) * Group + AGEYR + APOE4 + EDUC + AMYLCENT,
  random = ~ ADURW | BID,
  data = ADQS_PACC,
  na.action = na.exclude,
  control = nlme::lmeControl(opt = "optim"))


## ----lme-coefs------------------------------------------------------------------------------------------------------------
coef_tab <- summary(fit_lme)$tTable |>
  as.data.frame() |>
  rownames_to_column("term") |>
  mutate(p.value = formatp(`p-value`)) |>
  select(term, Value, Std.Error, DF, `t-value`, p.value)

kable(coef_tab, digits = 3)


## ----lme-var-corr---------------------------------------------------------------------------------------------------------
VarCorr(fit_lme)


## ----lme-sub-plot---------------------------------------------------------------------------------------------------------
fit_lme_sub <- nlme::lme(
  PACC ~ splines::ns(ADURW, df = 2) * Group + AGEYR + APOE4 + EDUC + AMYLCENT,
  random = ~ ADURW | BID,
  data = ADQS_PACC |> filter(BID %in% unique(ADQS_PACC$BID)[1:12]),
  na.action = na.exclude,
  control = nlme::lmeControl(opt = "optim"))

plot(fit_lme_sub, fitted(.) ~ ADURW | BID, abline = 0, 
  xlab = 'Weeks from baseline')


## ----lme-means------------------------------------------------------------------------------------------------------------
visits <- seq(0, 456, 24)

pacc_emms <- ref_grid(fit_lme, 
  at = list(ADURW = visits),
  cov.reduce = AMYLCENT ~ Group) |>
  emmeans(specs = ~ Group | ADURW) |>
  as.data.frame()


## ----lme-fig--------------------------------------------------------------------------------------------------------------
p1 <- ggplot(pacc_emms, aes(x = ADURW, y = emmean, group = Group)) +
  geom_line(aes(color = Group), linewidth = 1.1) +
  geom_point(aes(color = Group), size = 1.8) +
  geom_ribbon(aes(fill = Group, ymin = lower.CL, ymax = upper.CL), alpha = 0.2) +
  labs(
    x = "Weeks from baseline",
    y = "Model-based mean PACC (95% CI)") +
  scale_color_manual(values = grpPalette) +
  scale_fill_manual(values = grpPalette) + 
  theme(legend.position = "none")

p2 <- ggplot(pacc_emms, aes(x = ADURW, y = emmean, group = Group)) +
  geom_line(data = ADQS_PACC, aes(group = BID, y = PACC,
    color = Group), alpha = 0.05) +
  geom_line(aes(color = Group), linewidth = 1.1) +
  geom_point(aes(color = Group), size = 1.8) +
  geom_ribbon(aes(fill = Group, ymin = lower.CL, ymax = upper.CL), alpha = 0.2) +
  labs(
    x = "Weeks from baseline",
    y = "Model-based mean PACC (95% CI)") +
  scale_color_manual(values = grpPalette) +
  scale_fill_manual(values = grpPalette) +
  guides(color = guide_legend(override.aes = list(alpha = 1)))

p1 + p2 + plot_layout(guides = "collect")


## ----lme-version----------------------------------------------------------------------------------------------------------
fit_lme_ver <- update(fit_lme, . ~ . + VERSION)

data.frame(
  Model = c("Without version effect", "With version effect"),
  AIC = c(AIC(fit_lme), AIC(fit_lme_ver)),
  BIC = c(BIC(fit_lme), BIC(fit_lme_ver))) |>
  arrange(AIC) |>
  kable(digits = 2)


## ----gls-spline-----------------------------------------------------------------------------------------------------------
ADQS_PACC_A4 <- ADQS_PACC |>
  filter(SUBSTUDY == "A4", !is.na(ASEQNCS), MITTFL == 1,
    EPOCH == "BLINDED TREATMENT" | AVISIT == "006") |>
  mutate(ASEQNCS_f = as.factor(ASEQNCS))

a4_spline_toeph <- gls(
    PACC ~ splines::ns(ADURW, df = 2) : TX + 
    AGEYR + APOE4 + EDUC + AMYLCENT + QSVERSION,
  data = ADQS_PACC_A4,
  weights = varIdent(form = ~ 1 | ASEQNCS),
  correlation = corARMA(form = ~ ASEQNCS | BID, p = 10))

ref_grid(a4_spline_toeph,
  at = list(ADURW = c(0,240), TX = levels(ADQS_PACC_A4$TX)),
  vcov. = clubSandwich::vcovCR(a4_spline_toeph, type = "CR2") |> as.matrix(), 
  data = ADQS_PACC_A4, 
  mode = "satterthwaite") |>
  emmeans(~ ADURW | TX) |>
  pairs(reverse = TRUE) |>
  as_tibble() |>
  select(-contrast) |>
  mutate(p.value = formatp(p.value)) |>
  kable(caption = "Mean PACC change from baseline at week 240 by treatment 
    group estimated from spline model.", digits = 2, booktabs = TRUE)

ref_grid(a4_spline_toeph, 
  at = list(ADURW = 240, TX = levels(ADQS_PACC_A4$TX)),
  vcov. = clubSandwich::vcovCR(a4_spline_toeph, type = "CR2") |> as.matrix(), 
  data = ADQS_PACC_A4, 
  mode = "satterthwaite") |>
  emmeans(specs = "TX", by = "ADURW") |>
  pairs(reverse = TRUE, adjust = "none") |>
  summary(infer = c(TRUE, TRUE)) |> 
  as_tibble() |>
  mutate(p.value = formatp(p.value)) |>
  kable(caption = "Mean PACC group change from baseline at week 240 by treatment 
    group estimated from spline model.", digits = 2)


## -------------------------------------------------------------------------------------------------------------------------
a4_spline_us <- mmrm(
  PACC ~ splines::ns(ADURW, df = 2) : TX + 
    AGEYR + APOE4 + EDUC + AMYLCENT + QSVERSION + us(ASEQNCS_f | BID),
  data = ADQS_PACC_A4,
  method = "Kenward-Roger")

ref_grid(a4_spline_us,
  at = list(ADURW = c(0,240), TX = levels(ADQS_PACC_A4$TX))) |>
  emmeans(~ ADURW | TX) |>
  pairs(reverse = TRUE) |>
  as_tibble() |>
  mutate(p.value = formatp(p.value)) |>
  select(-contrast) |>
  kable(caption = "Mean PACC change from baseline at week 240 by treatment 
    group estimated from spline model.", digits = 2)

ref_grid(a4_spline_us, 
  at = list(ADURW = 240, TX = levels(ADQS_PACC_A4$TX))) |>
  emmeans(specs = "TX", by = "ADURW") |>
  pairs(reverse = TRUE, adjust = "none") |> 
  summary(infer = c(TRUE, TRUE)) |> 
  as_tibble() |>
  mutate(p.value = formatp(p.value)) |>
  kable(caption = "Mean group difference in PACC change from baseline at week 240 estimated from spline model.", digits = 2)


## -------------------------------------------------------------------------------------------------------------------------
ADQS_PACC_A4_MMRM <- ADQS_PACC |>
  filter(SUBSTUDY == "A4", !is.na(ASEQMMRM), MITTFL == 1,
    EPOCH == "BLINDED TREATMENT" | AVISIT == "006") |>
  mutate(
    ASEQMMRM_f = as.factor(ASEQMMRM),
    AVISIT = as.factor(AVISIT))

a4_cat_us <- mmrm(
  PACC_ch ~ AVISIT * (TX + PACC_bl) +
    AGEYR + APOE4 + EDUC + AMYLCENT + QSVERSION + us(AVISIT | BID),
  data = ADQS_PACC_A4_MMRM,
  method = "Kenward-Roger")

ref_grid(a4_cat_us,
  at = list(AVISIT = '066', TX = levels(ADQS_PACC_A4_MMRM$TX))) |>
  emmeans(~ AVISIT | TX) |>
  summary(infer = c(TRUE, TRUE)) |> 
  as_tibble() |>
  mutate(p.value = formatp(p.value)) |>
  select(-AVISIT) |>
  kable(caption = "Mean PACC change from baseline at visit 66 by treatment 
    group estimated from categorical time MMRM.", digits = 2)

ref_grid(a4_cat_us, 
  at = list(AVISIT = '066', TX = levels(ADQS_PACC_A4_MMRM$TX))) |>
  emmeans(specs = "TX", by = "AVISIT") |>
  pairs(reverse = TRUE, adjust = "none") |> 
  summary(infer = c(TRUE, TRUE)) |> 
  as_tibble() |>
  mutate(p.value = formatp(p.value)) |>
  select(-AVISIT) |>
  kable(caption = "Mean group difference in PACC change from baseline at visit 66 estimated from categorical time MMRM.", digits = 2)


## -------------------------------------------------------------------------------------------------------------------------
set.seed(20240701)
ADQS_PACC_PSEUDO_MMRM <- ADQS_PACC_A4_MMRM[1:400,] |>
  mutate(
    TX = as.character(TX),
    TX = case_when(TX == "Solanezumab" ~ "Pseudo", TRUE ~ "Placebo"),
    # Add a linear PACC benefit for those randomize dto "Pseudo"
    PACC_ch = case_when(TX == "Pseudo" ~ PACC_ch + ADURW/240*4.5, TRUE ~ PACC))

pseudo_cat_us <- mmrm(
  PACC_ch ~ AVISIT * (TX + PACC_bl) +
    AGEYR + APOE4 + EDUC + AMYLCENT + QSVERSION + us(AVISIT | BID),
  data = ADQS_PACC_PSEUDO_MMRM,
  method = "Kenward-Roger")

ref_grid(pseudo_cat_us,
  at = list(AVISIT = '066', TX = unique(ADQS_PACC_PSEUDO_MMRM$TX))) |>
  emmeans(~ AVISIT | TX) |>
  summary(infer = c(TRUE, TRUE)) |> 
  as_tibble() |>
  mutate(p.value = formatp(p.value)) |>
  select(-AVISIT, -QSVERSION) |>
  kable(caption = "Mean PACC change from baseline at visit 66 by treatment 
    group estimated from categorical time MMRM.", digits = 2)

ref_grid(pseudo_cat_us, 
  at = list(AVISIT = '066', TX = unique(ADQS_PACC_PSEUDO_MMRM$TX))) |>
  emmeans(specs = "TX", by = "AVISIT") |>
  pairs(reverse = TRUE, adjust = "none") |> 
  summary(infer = c(TRUE, TRUE)) |> 
  as_tibble() |>
  mutate(p.value = formatp(p.value)) |>
  select(-AVISIT) |>
  kable(caption = "Mean group difference in PACC change from baseline at visit 66 estimated from categorical time MMRM.", digits = 2)


## ----results = 'hide', echo = FALSE---------------------------------------------------------------------------------------
# get default predictor matrix
trial_wide <- ADQS_PACC_PSEUDO_MMRM |>
  mutate(active = case_when(TX == "Pseudo" ~ 1, TRUE ~ 0)) |>
  select(BID, AVISIT, PACC_ch, active, PACC_bl, AGEYR, APOE4, 
    EDUC, AMYLCENT, QSVERSION, ADURW) |>
  pivot_wider(names_from = AVISIT, 
    values_from = c(PACC_ch, QSVERSION, ADURW)) |>
  assertr::verify(anyDuplicated(BID) == 0)

pm <- quickpred(trial_wide, mincor = 0.01)
pm[,"BID"] <- 0 # never use ID as predictor
# always use baseline predictors for PACC change:
pm[c("active", "AGEYR", "PACC_bl", "APOE4", "EDUC", "AMYLCENT"), 
  grepl("PACC_ch", colnames(pm))] <- 1
# always use version and time predictors for PACC change:
for(cc in unique(ADQS_PACC_PSEUDO_MMRM$AVISIT)){
  pm[paste0("PACC_ch_", cc), paste0("QSVERSION_", cc)] <- 1
  pm[paste0("PACC_ch_", cc), paste0("ADURW_", cc)] <- 1
}

meth <- make.method(trial_wide)


## ----results = 'hide', echo = TRUE----------------------------------------------------------------------------------------
trial_imp <- mice(trial_wide, m=5, predictorMatrix=pm, 
  seed = 20170714, maxit=100, method = meth, print = FALSE)


## ----echo = TRUE, size = 'scriptsize'-------------------------------------------------------------------------------------
trial_wide[,grepl("PACC_ch|BID", colnames(pm))] |> 
  filter(BID %in% c("B10350512", "B10656630", "B10892860"))
complete(trial_imp)[,grepl("PACC_ch|BID", colnames(pm))] |>
  filter(BID %in% c("B10350512", "B10656630", "B10892860"))


## ----echo = FALSE, size = 'scriptsize'------------------------------------------------------------------------------------
fits_mi <- with(data=trial_imp, lm(PACC_ch_066~active*scale(PACC_bl) +
  AGEYR + APOE4 + EDUC + AMYLCENT))
summary(fits_mi)


## ----echo = FALSE---------------------------------------------------------------------------------------------------------
summary(pool(fits_mi)) |>
  remove_rownames() |>
  column_to_rownames(var="term") |>
  mutate(p.value = formatp(p.value)) |>
  kable(digits = 3)


## -------------------------------------------------------------------------------------------------------------------------
adjust_dataset <- function(dat, delta) {
  treat_idx <- dat$active == 1
  # only modify imputed values in active arm
  idx <- missing_idx & treat_idx
  dat$PACC_ch_066[idx] <- dat$PACC_ch_066[idx] - delta
  dat
}

fit_models <- function(data_list) {
  lapply(data_list, function(dat) {
    lm(PACC_ch_066 ~ active * scale(PACC_bl) +
        AGEYR + APOE4 + EDUC + AMYLCENT,
      data = dat)})
}

tx_effect_mar <- pool(fits_mi)$pooled |> filter(term == "active") |> pull(estimate)
missing_idx <- is.na(trial_wide$PACC_ch_066)
completed_list <- complete(trial_imp, action = "all")

tipping_factors <- c(0, seq(0.5, 1, 0.1))
est_tipping <- lapply(tipping_factors, function(tf) {
  delta <- tf * tx_effect_mar
  adjusted_list <- lapply(completed_list, adjust_dataset, delta = delta)
  fits <- fit_models(adjusted_list)
  pooled <- summary(pool(as.mira(fits)))
  pooled |>
    dplyr::filter(term == "active") |>
    dplyr::mutate(
      tipping_factor = tf,
      delta = delta,
      .before = 1
    )
}) |>
  bind_rows()


## ----echo=FALSE-----------------------------------------------------------------------------------------------------------
bind_rows(est_tipping) |>
  mutate(p.value = formatp(p.value)) |>
  kable()


## -------------------------------------------------------------------------------------------------------------------------
sessionInfo()

