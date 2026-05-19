# =============================================================================
# experiment3_modeling.R
# Experiment 3 — Mixed-Effects Regression Models
#
# Purpose:
#   Fits two mixed-effects regression models on trial-level data to test
#   how headline morphing (credibility-boosted vs. engagement-boosted vs.
#   neutral) and instruction condition (accuracy vs. engagement) jointly
#   predict (1) whether a participant selects a headline and (2) how
#   credible they rate it. Outputs PNAS-style formatted tables via kableExtra.
#
# Input:  final_analysis_df.csv  — trial-level data produced by the Python
#                                   cleaning notebook
# Output: Rendered HTML tables (inline, via kbl/kable_classic)
# =============================================================================

# ── Libraries ─────────────────────────────────────────────────────────────────
# lme4 / lmerTest: mixed-effects models (lmerTest adds p-values to lmer)
# broom.mixed:     tidy() for extracting model coefficients as a data frame
# kableExtra:      PNAS-style formatted regression tables
# tidyverse:       data wrangling and piping
library(lme4)
library(stargazer)
library(dplyr)
library(kableExtra)
library(tidyverse)
library(tidyr)
library(broom.mixed)
library(ggplot2)
library(lmerTest)   # must load AFTER lme4; overrides lmer() to supply p-values


# ── Load data ─────────────────────────────────────────────────────────────────
# Trial-level dataframe: one row per participant × headline combination.
# Key columns used below:
#   Selected         — binary DV: did the participant add this headline to cart?
#   credibility_rating — continuous DV: participant's 1–7 credibility rating
#   MorphCredibility — dummy: 1 if headline was credibility-boosted, else 0
#   MorphEngagement  — dummy: 1 if headline was engagement-boosted, else 0
#   Instruction      — dummy: 0 = accuracy instruction, 1 = engagement instruction
#   PROLIFIC_PID     — participant identifier (used as random effect)
#   GOID             — headline identifier (used as random effect)
df <- read.csv("/Users/Rhea/Documents/responsibility_lab/online_news/final_analysis_df.csv",
               stringsAsFactors = FALSE)

# Rename ID columns to intuitive grouping labels for lme4 random-effect syntax
df$participant <- df$PROLIFIC_PID
df$headline    <- df$GOID

# ── Sanity checks ─────────────────────────────────────────────────────────────
# Confirm expected sample size, trial count, and baseline selection rate
cat("N participants:", length(unique(df$participant)), "\n")
cat("N headlines:   ", length(unique(df$headline)), "\n")
cat("N trials:      ", nrow(df), "\n")
cat("Selection rate:", round(mean(df$Selected), 3), "\n\n")


# =============================================================================
# MODEL 1 — Predicting Headline Selection (Binary DV)
# -----------------------------------------------------------------------------
# A mixed-effects logistic regression (glmer, binomial family).
# Fixed effects: morph type, instruction condition, and their interactions.
# Random effects: by-participant and by-headline intercepts, capturing
#   individual differences in selection tendency and headline-level variance.
# Optimizer: bobyqa with increased iterations to ensure convergence.
# =============================================================================
model <- glmer(
  Selected ~ MorphCredibility + MorphEngagement + Instruction
  + MorphCredibility:Instruction + MorphEngagement:Instruction
  + (1 | participant) + (1 | headline),
  data    = df,
  family  = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

# ── Extract fixed-effect statistics manually ──────────────────────────────────
# (glmer does not produce p-values via summary(); computed here from z-scores)
coefs  <- fixef(model)               # point estimates (log-odds)
ses    <- sqrt(diag(vcov(model)))    # standard errors
zvals  <- coefs / ses                # z-statistics
pvals  <- 2 * pnorm(abs(zvals), lower.tail = FALSE)  # two-tailed p-values

# ── Extract random-effect variance components ─────────────────────────────────
# VarCorr gives variance (vcov) and SD (sdcor) for each grouping factor
vc       <- as.data.frame(VarCorr(model))
var_part <- round(vc[vc$grp == "participant", "vcov"],  3)
sd_part  <- round(vc[vc$grp == "participant", "sdcor"], 3)
var_head <- round(vc[vc$grp == "headline",    "vcov"],  3)
sd_head  <- round(vc[vc$grp == "headline",    "sdcor"], 3)

# ── Build formatted output table ──────────────────────────────────────────────
# tidy() extracts fixed effects as a clean data frame.
# Predictor terms are relabeled for readability, estimates rounded,
# and significance stars appended to formatted p-values.
model_tbl <- tidy(model, effects = "fixed", conf.int = FALSE) %>%
  mutate(
    # Recode raw term names to human-readable predictor labels
    term = recode(term,
                  "(Intercept)"                  = "Intercept (neutral headline, accuracy instructions)",
                  "MorphCredibility"             = "Credibility-boosted",
                  "MorphEngagement"              = "Engagement-boosted",
                  "Instruction"                  = "Instruction",
                  "MorphCredibility:Instruction" = "Credibility-boosted × Instruction",
                  "MorphEngagement:Instruction"  = "Engagement-boosted × Instruction"
    ),
    estimate  = round(estimate,   3),
    std.error = round(std.error,  3),
    statistic = round(statistic,  3),
    # Compute significance stars BEFORE overwriting p.value
    stars   = symnum(p.value,
                     cutpoints = c(0, .001, .01, .05, 1),
                     symbols   = c("***", "**", "*", "")),
    # Format p-values: "< .001" for very small values, 3 decimal places otherwise
    p.value = paste0(
      case_when(
        p.value < .001 ~ "< .001",
        TRUE           ~ formatC(p.value, digits = 3, format = "f")
      ),
      stars
    )
  ) %>%
  select(
    `Predictor` = term,
    `β`         = estimate,
    `SE`        = std.error,
    `z`         = statistic,   # z-statistic for logistic model
    `p`         = p.value
  )

# ── Render PNAS-style HTML table ──────────────────────────────────────────────
model_tbl %>%
  kbl(
    align    = c("l", "r", "r", "r", "r"),
    caption  = "Mixed-effects logistic regression predicting headline selection.",
    booktabs = TRUE
  ) %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman") %>%
  row_spec(0, bold = TRUE) %>%
  row_spec(nrow(model_tbl), extra_css = "border-bottom: 1px solid black;") %>%
  column_spec(1, width = "12em") %>%
  column_spec(2:5, width = "5em") %>%
  add_footnote(
    c("* p < .05  ** p < .01  *** p < .001",
      "Includes random effects for participants and headlines.",
      "Credibility-boosted = 1 and engagement-boosted = 0 for credibility variant headlines.",
      "Credibility-boosted = 0 and engagement-boosted = 1 for engagement variant headlines.",
      "Instruction = 1 for the engagement instruction condition, and 0 for the accuracy instruction condition."
      ),
    notation = "none"
  )


# =============================================================================
# MODEL 2 — Predicting Perceived Credibility (Continuous DV)
# -----------------------------------------------------------------------------
# A mixed-effects linear regression (lmer, REML estimation).
# Same fixed-effect structure as Model 1, but the outcome is a 1–7
# credibility rating rather than a binary selection decision.
# lmerTest overrides lmer() to produce Satterthwaite p-values in summary().
# =============================================================================
model_linear <- lmer(
  credibility_rating ~ MorphCredibility + MorphEngagement + Instruction
  + MorphCredibility:Instruction + MorphEngagement:Instruction
  + (1 | participant) + (1 | headline),
  data    = df,
  REML    = TRUE   # REML is preferred for unbiased variance component estimation
)
summary(model_linear)

# ── Build formatted output table (same pipeline as Model 1) ───────────────────
model_tbl <- tidy(model_linear, effects = "fixed", conf.int = FALSE) %>%
  mutate(
    term = recode(term,
                  "(Intercept)"                        = "Intercept (neutral headline, accuracy instructions)",
                  "MorphCredibility"                   = "Credibility-boosted",
                  "MorphEngagement"                    = "Engagement-boosted",
                  "Instruction"                        = "Instruction",
                  "MorphCredibility:Instruction"       = "Credibility-boosted × Instruction",
                  "MorphEngagement:Instruction"        = "Engagement-boosted × Instruction"
    ),
    estimate  = round(estimate,   3),
    std.error = round(std.error,  3),
    statistic = round(statistic,  3),
    stars = as.character(symnum(p.value,   # cast to character for safe paste0()
                   cutpoints = c(0, .001, .01, .05, 1),
                   symbols   = c("***", "**", "*", ""))),
    p.value = paste0(
      case_when(
        p.value < .001 ~ "< .001",
        TRUE           ~ formatC(p.value, digits = 3, format = "f")
      ),
      stars
    )
  ) %>%
  select(
    `Predictor` = term,
    `β`         = estimate,
    `SE`        = std.error,
    `t`         = statistic,   # t-statistic for linear model (vs. z above)
    `p`         = p.value
  )

# ── Render PNAS-style HTML table ──────────────────────────────────────────────
model_tbl %>%
  kbl(
    align    = c("l", "r", "r", "r", "r"),
    caption  = "Mixed-effects linear regression predicting credibility ratings.",
    booktabs = TRUE
  ) %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman") %>%
  row_spec(0, bold = TRUE) %>%
  row_spec(nrow(model_tbl), extra_css = "border-bottom: 1px solid black;") %>%
  column_spec(1, width = "12em") %>%
  column_spec(2:5, width = "5em") %>%
  add_footnote(
    c("* p < .05  ** p < .01  *** p < .001",
      "Includes random effects for participants and headlines.",
      "Credibility-boosted = 1 and engagement-boosted = 0 for credibility variant headlines.",
      "Credibility-boosted = 0 and engagement-boosted = 1 for engagement variant headlines.",
      "Instruction = 1 for the engagement instruction condition, and 0 for the accuracy instruction condition."
    ),
    notation = "none"
  )
# Table caption: illustrates how credibility ratings vary by headline variant
# type, instruction condition, and their interactions.
