library(lme4)
library(stargazer)
library(dplyr)
library(kableExtra)
library(tidyverse)
library(tidyr)
library(broom.mixed)
library(ggplot2)
library(lmerTest)

# ── Load data ──────────────────────────────────────────────────────────────
df <- read.csv("/Users/Rhea/Documents/responsibility_lab/online_news/final_analysis_df.csv", stringsAsFactors = FALSE)

# Random effect grouping columns
df$participant <- df$PROLIFIC_PID
df$headline    <- df$GOID

# ── Sanity check ───────────────────────────────────────────────────────────
cat("N participants:", length(unique(df$participant)), "\n")
cat("N headlines:   ", length(unique(df$headline)), "\n")
cat("N trials:      ", nrow(df), "\n")
cat("Selection rate:", round(mean(df$Selected), 3), "\n\n")


# ── Fit model ──────────────────────────────────────────────────────────────
model <- glmer(
  Selected ~ MorphCredibility + MorphEngagement + Instruction
  + MorphCredibility:Instruction + MorphEngagement:Instruction
  + (1 | participant) + (1 | headline),
  data    = df,
  family  = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

# ── Extract fixed effects, SEs, z-scores, p-values ────────────────────────
coefs  <- fixef(model)
ses    <- sqrt(diag(vcov(model)))
zvals  <- coefs / ses
pvals  <- 2 * pnorm(abs(zvals), lower.tail = FALSE)

# ── Extract random effect variances/SDs ───────────────────────────────────
vc       <- as.data.frame(VarCorr(model))
var_part <- round(vc[vc$grp == "participant", "vcov"],  3)
sd_part  <- round(vc[vc$grp == "participant", "sdcor"], 3)
var_head <- round(vc[vc$grp == "headline",    "vcov"],  3)
sd_head  <- round(vc[vc$grp == "headline",    "sdcor"], 3)


# Extract tidy model output
model_tbl <- tidy(model, effects = "fixed", conf.int = FALSE) %>%
  mutate(
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
    # Compute stars from raw p.value FIRST, then format
    stars   = symnum(p.value,
                     cutpoints = c(0, .001, .01, .05, 1),
                     symbols   = c("***", "**", "*", "")),
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
    `z`         = statistic,
    `p`         = p.value
  )

# Render PNAS-style table
model_tbl %>%
  kbl(
    align   = c("l", "r", "r", "r", "r"),
    caption = "Mixed-effects logistic regression predicting headline selection.",
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
# caption: Table illustrates headline selection likelihood by headline variant type, 
# instruction condition, and the interactions of these independent variables.


# ---- MODEL #2: PREDICTING PERCEIVED CREDIBILITY -----------------------
model_linear <- lmer(
  credibility_rating ~ MorphCredibility + MorphEngagement + Instruction
  + MorphCredibility:Instruction + MorphEngagement:Instruction
  + (1 | participant) + (1 | headline),
  data    = df,
  REML    = TRUE
)
summary(model_linear)

# Extract tidy model output
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
    stars = symnum(p.value,
                   cutpoints = c(0, .001, .01, .05, 1),
                   symbols   = c("***", "**", "*", "")),
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
    `t`         = statistic,
    `p`         = p.value
  )

# Render PNAS-style table
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

