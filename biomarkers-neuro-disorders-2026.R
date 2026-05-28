## -----------------------------------------------------------------------------------------------------------
#| label: setup
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


## -----------------------------------------------------------------------------------------------------------
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

V1.6PTAU217 <- A4LEARN::biomarker_pTau217 |>
    filter(TESTCD == "PTAU217", !is.na(ORRESRAW)) |>
    arrange(BID, VISCODE) |>
    filter(!duplicated(BID)) |>
    filter(VISCODE %in% c(1,6)) |>
    select(BID, PTAU217 = ORRESRAW)

SUBJINFO <- A4LEARN::SUBJINFO |>
  filter(SUBSTUDY != "SF") |> 
  left_join(V6OUTCOME, by = "BID") |>
  left_join(V1OUTCOME |> 
    select(BID, CDRSB, CFITOTAL, CFISP, CFIPT, ADLPQPT, ADLPQSP),
    by = "BID") |> 
  left_join(V1.6PTAU217, by = "BID")


## -----------------------------------------------------------------------------------------------------------
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

ALL_BASE_PET <- A4LEARN::SUBJINFO |>
  filter(!is.na(AMYLCENT))


## -----------------------------------------------------------------------------------------------------------
#| echo: true
set.seed(20200225)

batch_data <- 
  tibble(
    batch = 1:10,
    Sigma = rgamma(n = 10, shape = 360 / 10, scale = 10),
    Mean  = rnorm(n = 10, mean = 850, sd = 200)
  ) |>
  group_by(batch) |>
  nest() |>
  mutate(Biomarker = purrr::map(data, ~rnorm(n = 50, .x$Mean, .x$Sigma))) |>
  unnest(Biomarker) |>
  unnest(data) |>
  ungroup() |>
  arrange(batch) |>
  mutate(
    id = row_number(),
    batch = factor(batch),
    Biomarker = pmax(Biomarker, 0)
  )

ggplot(batch_data, aes(x = batch, y = Biomarker)) +
  geom_boxplot(outlier.shape = NA) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.3, alpha = 0.2)


## -----------------------------------------------------------------------------------------------------------
#| echo: true
anova(lm(Biomarker ~ batch, data = batch_data)) |>
  broom::tidy() |>
  mutate(p.value = formatp(p.value)) |>
  kable(digits = 3)


## -----------------------------------------------------------------------------------------------------------
#| echo: true
# Model 1 (Null): Assumes equal variance across all batches
fit_equal_var <- gls(Biomarker ~ batch, data = batch_data)

# Model 2 (Alternative): Allows different variances for each batch
fit_diff_var <- gls(Biomarker ~ batch, data = batch_data,
  weights = varIdent(form = ~1 | batch))

# Perform the Likelihood Ratio Test
anova(fit_equal_var, fit_diff_var) |>
  as_tibble(rownames = "model") |>
  mutate(`p-value` = formatp(`p-value`)) |>
  select(-call, -Model) |> 
  kable(digits = 3)


## -----------------------------------------------------------------------------------------------------------
#| echo: true
set.seed(20260603)
batch_data$Biomarker_ComBat <- sva::ComBat(
  dat = rbind(batch_data$Biomarker, jitter(batch_data$Biomarker)),
  batch = batch_data$batch,
  par.prior = TRUE,
  prior.plots = FALSE
)[1,]


## -----------------------------------------------------------------------------------------------------------
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

ggplot(batch_data, aes(x = batch, y = Biomarker_ComBat)) +
  geom_boxplot(outlier.shape = NA) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.3, alpha = 0.2)


## -----------------------------------------------------------------------------------------------------------
#| results: asis

A4labels <- c(
  AGEYR = "Age (years)",
  APOE4 = "APOEε4 carrier",
  EDUC = "Education (years)",
  PTAU217 = "Baseline pTau217 (U/ml)",
  AMYLCENT = "Baseline amyloid PET (CL)",
  PACC = "Baseline PACC"
)

BASE_TAB <- arsenal::tableby(
  Group ~ AGEYR + APOE4 + EDUC + PTAU217 + AMYLCENT + PACC,
  data = BASE_PACC
)

summary(BASE_TAB, labelTranslations = A4labels, text = TRUE)


## -----------------------------------------------------------------------------------------------------------
ggplot(BASE_PACC, aes(x = Group, y = PACC, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  ggbeeswarm::geom_beeswarm(width = 0.15, alpha = 0.2, size = 1) +
  labs(x = "", y = "Baseline PACC") +
  theme(legend.position = "none") +
  scale_fill_manual(values = grpPalette)


## -----------------------------------------------------------------------------------------------------------
ggplot(BASE_PACC, aes(x = AMYLCENT, y = PACC, color = Group)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Baseline amyloid PET biomarker (CL)",
    y = "Baseline PACC") +
  scale_colour_manual(values = grpPalette)


## -----------------------------------------------------------------------------------------------------------
fit_bl <- lm(PACC ~ scale(AMYLCENT) + scale(PTAU217) + scale(AGEYR) + APOE4 + 
  scale(EDUC), data = BASE_PACC)

broom::tidy(fit_bl) |>
  mutate(p.value = formatp(p.value)) |>
  kable(digits = 3)


## -----------------------------------------------------------------------------------------------------------
ggplot(ALL_BASE_PET, aes(x = AMYLCENT)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "grey75", color = "white") +
  geom_density(color = cbbPalette[1], linewidth = 1.2) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = 100, color = "red", linetype = "dashed", size = 1) +
  labs(x = "All screening amyloid PET (CL)", y = "Density")


## -----------------------------------------------------------------------------------------------------------
#| echo: true
ecdf_suvr <- ecdf(ALL_BASE_PET$AMYLCENT)

ALL_BASE_PET <- ALL_BASE_PET |>
  mutate(
    AMYLCENT_pct = ecdf_suvr(AMYLCENT),
    AMYLCENT_z = qnorm(pmin(pmax(AMYLCENT_pct, 1e-6), 1 - 1e-6)))


## -----------------------------------------------------------------------------------------------------------
#| echo: false
p1 <- ggplot(ALL_BASE_PET, aes(x = AMYLCENT)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "grey75", color = "white") +
  geom_density(color = cbbPalette[1], linewidth = 1.2) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = 100, color = "red", linetype = "dashed", size = 1) +
  labs(x = "Amyloid PET (CL)", y = "Density")

p2 <- ggplot(ALL_BASE_PET |> arrange(AMYLCENT), aes(x = AMYLCENT, y = AMYLCENT_pct*100)) +
  geom_line(color = cbbPalette[2], size = 2) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = 100, color = "red", linetype = "dashed", size = 1) +
  labs(x = "Amyloid PET (CL)", y = "Percentile (%)")

p3 <- ggplot(ALL_BASE_PET |> arrange(AMYLCENT), aes(x = AMYLCENT, y = AMYLCENT_z)) +
  geom_line(color = cbbPalette[2], size = 2) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = 100, color = "red", linetype = "dashed", size = 1) +
  labs(x = "Amyloid PET (CL)", y = "Z-score")

p1 + p2 + p3


## -----------------------------------------------------------------------------------------------------------
#| echo: false

p1 <- ggplot(ALL_BASE_PET, aes(x = AMYLCENT_z)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "grey75", color = "white") +
  geom_density(color = cbbPalette[1], linewidth = 1.2) +
  labs(x = "ECDF-derived z-score for amyloid PET", y = "Density")

p2 <- ggplot(ALL_BASE_PET, aes(x = scale(AMYLCENT))) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "grey75", color = "white") +
  geom_density(color = cbbPalette[1], linewidth = 1.2) +
  labs(x = "Usual z-score for amyloid PET", y = "Density")

p1 + p2


## -----------------------------------------------------------------------------------------------------------
#| eval: true
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


## -----------------------------------------------------------------------------------------------------------
ggplot(ADQS_PACC, aes(x = ADURW, y = PACC, color = Group)) +
  geom_line(aes(group = BID), alpha = 0.12) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1.2) +
  labs(
    x = "Weeks from baseline",
    y = "PACC") +
  scale_color_manual(values = grpPalette)


## -----------------------------------------------------------------------------------------------------------
PACC_SUM <- ADQS_PACC |>
  filter(!is.na(Week), Week %in% seq(0, 456, 24)) |>
  group_by(Group, Week) |>
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


## -----------------------------------------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------------------------------------
fit_lme <- nlme::lme(
  PACC ~ ns(ADURW, df = 2) * Group + AGEYR + APOE4 + EDUC + AMYLCENT,
  random = ~ ADURW | BID,
  data = ADQS_PACC,
  na.action = na.exclude,
  control = nlme::lmeControl(opt = "optim"))


## -----------------------------------------------------------------------------------------------------------
coef_tab <- summary(fit_lme)$tTable |>
  as.data.frame() |>
  rownames_to_column("term") |>
  mutate(p.value = formatp(`p-value`)) |>
  select(term, Value, Std.Error, DF, `t-value`, p.value)

kable(coef_tab, digits = 3)


## -----------------------------------------------------------------------------------------------------------
visits <- seq(0, 456, 24)

pacc_emms <- ref_grid(fit_lme, 
  at = list(ADURW = visits),
  cov.reduce = AMYLCENT ~ Group) |>
  emmeans(specs = ~ Group | ADURW) |>
  as.data.frame()


## -----------------------------------------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------------------------------------
fit_lme_ver <- update(fit_lme, . ~ . + VERSION)

data.frame(
  Model = c("Without version effect", "With version effect"),
  AIC = c(AIC(fit_lme), AIC(fit_lme_ver)),
  BIC = c(BIC(fit_lme), BIC(fit_lme_ver))) |>
  arrange(AIC) |>
  kable(digits = 2)


## -----------------------------------------------------------------------------------------------------------
ADQS_PACC_A4 <- ADQS_PACC |>
  filter(SUBSTUDY == "A4", !is.na(ASEQNCS), MITTFL == 1,
    EPOCH == "BLINDED TREATMENT" | AVISIT == "006") |>
  mutate(ASEQNCS_f = as.factor(ASEQNCS))

a4_spline_gls <- gls(
    PACC ~ ns(ADURW, df = 2) : TX + 
    AGEYR + APOE4 + EDUC + AMYLCENT + QSVERSION,
  data = ADQS_PACC_A4,
  weights = varIdent(form = ~ 1 | ASEQNCS),
  correlation = corARMA(form = ~ ASEQNCS | BID, p = 10))

ref_grid(a4_spline_gls,
  at = list(ADURW = c(0,240), TX = levels(ADQS_PACC_A4$TX)),
  vcov. = clubSandwich::vcovCR(a4_spline_gls, type = "CR2") |> as.matrix(), 
  data = ADQS_PACC_A4, 
  mode = "satterthwaite") |>
  emmeans(~ ADURW | TX) |>
  pairs(reverse = TRUE) |>
  as_tibble() |>
  select(-contrast) |>
  mutate(p.value = formatp(p.value)) |>
  kable(caption = "Mean PACC change from baseline at week 240 by treatment 
    group estimated from spline model.", digits = 2, booktabs = TRUE)

ref_grid(a4_spline_gls, 
  at = list(ADURW = 240, TX = levels(ADQS_PACC_A4$TX)),
  vcov. = clubSandwich::vcovCR(a4_spline_gls, type = "CR2") |> as.matrix(), 
  data = ADQS_PACC_A4, 
  mode = "satterthwaite") |>
  emmeans(specs = "TX", by = "ADURW") |>
  pairs(reverse = TRUE, adjust = "none") |>
  summary(infer = c(TRUE, TRUE)) |> 
  as_tibble() |>
  mutate(p.value = formatp(p.value)) |>
  kable(caption = "Mean PACC group change from baseline at week 240 by treatment 
    group estimated from spline model.", digits = 2)


## -----------------------------------------------------------------------------------------------------------
a4_spline_mmrm <- mmrm(
  PACC ~ ns(ADURW, df = 2) : TX + 
    AGEYR + APOE4 + EDUC + AMYLCENT + QSVERSION + us(ASEQNCS_f | BID),
  data = ADQS_PACC_A4,
  method = "Kenward-Roger")

ref_grid(a4_spline_mmrm,
  at = list(ADURW = c(0,240), TX = levels(ADQS_PACC_A4$TX))) |>
  emmeans(~ ADURW | TX) |>
  pairs(reverse = TRUE) |>
  as_tibble() |>
  mutate(p.value = formatp(p.value)) |>
  select(-contrast) |>
  kable(caption = "Mean PACC change from baseline at week 240 by treatment 
    group estimated from spline model.", digits = 2)

ref_grid(a4_spline_mmrm, 
  at = list(ADURW = 240, TX = levels(ADQS_PACC_A4$TX))) |>
  emmeans(specs = "TX", by = "ADURW") |>
  pairs(reverse = TRUE, adjust = "none") |> 
  summary(infer = c(TRUE, TRUE)) |> 
  as_tibble() |>
  mutate(p.value = formatp(p.value)) |>
  kable(caption = "Mean group difference in PACC change from baseline at week 240 estimated from spline model.", digits = 2)


## -----------------------------------------------------------------------------------------------------------
ADQS_PACC_A4_MMRM <- ADQS_PACC |>
  filter(SUBSTUDY == "A4", !is.na(ASEQMMRM), MITTFL == 1,
    EPOCH == "BLINDED TREATMENT" | AVISIT == "006") |>
  mutate(
    ASEQMMRM_f = as.factor(ASEQMMRM),
    AVISIT = as.factor(AVISIT))

a4_cat_mmrm <- mmrm(
  PACC_ch ~ AVISIT * (TX + PACC_bl) +
    AGEYR + APOE4 + EDUC + AMYLCENT + QSVERSION + us(AVISIT | BID),
  data = ADQS_PACC_A4_MMRM,
  method = "Kenward-Roger")

ref_grid(a4_cat_mmrm,
  at = list(AVISIT = '066', TX = levels(ADQS_PACC_A4_MMRM$TX))) |>
  emmeans(~ AVISIT | TX) |>
  summary(infer = c(TRUE, TRUE)) |> 
  as_tibble() |>
  mutate(p.value = formatp(p.value)) |>
  select(-AVISIT) |>
  kable(caption = "Mean PACC change from baseline at visit 66 by treatment 
    group estimated from categorical time MMRM.", digits = 2)

ref_grid(a4_cat_mmrm, 
  at = list(AVISIT = '066', TX = levels(ADQS_PACC_A4_MMRM$TX))) |>
  emmeans(specs = "TX", by = "AVISIT") |>
  pairs(reverse = TRUE, adjust = "none") |> 
  summary(infer = c(TRUE, TRUE)) |> 
  as_tibble() |>
  mutate(p.value = formatp(p.value)) |>
  select(-AVISIT) |>
  kable(caption = "Mean group difference in PACC change from baseline at visit 66 estimated from categorical time MMRM.", digits = 2)


## -----------------------------------------------------------------------------------------------------------
sessionInfo()

