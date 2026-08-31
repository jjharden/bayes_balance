# ==============================================================================
# Title: A Falsification Framework for Investigating Observed Covariate Balance
# Author: Jeffrey J. Harden
# Affiliation: Department of Political Science, University of Notre Dame
# Email: jharden@nd.edu
# Journal: Observational Studies
# Script: bayes_balance.R
# Description: Methodological foundation, core functions, and full diagnostic 
#              workflow (bootstrapping, summaries, tables, and ridge plots) 
#              using cobalt::bal.tab() and the LaLonde (1986) dataset.
# ==============================================================================

# ------------------------------------------------------------------------------
# INTRODUCTION & TUTORIAL FOR NEW USERS
# ------------------------------------------------------------------------------
# Welcome to the Bayesian Bootstrap Balance Falsification Framework!
#
# WHY USE THIS FRAMEWORK?
# Standard balance assessment relies on point estimates (e.g., mean differences)
# and null hypothesis significance testing (NHST). This common "folk definition"
# can hide substantial imbalance in higher distributional moments (variance, 
# overlap) and over-relies on p-values.
#
# My framework uses the Bayesian Bootstrap (Rubin 1981) to generate posterior
# distributions for balance statistics across multiple moments. With posteriors, 
# you can move past NHST and make intuitive probabilistic statements, such as:
#   1. "What is the probability that an observed covariate's balance statistic
#       falls inside an Equivalence Range (ER)?"
#   2. "What is the probability that matching improved balance by at least 5 
#       percentage points?"
#
# HOW THIS SCRIPT WORKS WITH 'COBALT':
# Rather than writing custom metrics for every dataset, this framework uses 
# cobalt::bal.tab() inside the bootstrap function. 'cobalt' is a useful package 
# for balance diagnostics in R. Passing bootstrap weights directly into bal.tab()
# allows you to easily extract Standardized Mean Differences (SMD), Variance 
# Ratios (VR), Kolmogorov-Smirnov (KS) statistics, and Overlap (OVL) statistics.
#
# WORKFLOW OVERVIEW:
# 1. Prepare data & run your adjustment method (e.g., MatchIt::matchit()).
# 2. Define a single-pass cobalt bootstrap extraction function.
# 3. Run bayesboot() to generate posterior distributions across all covariates.
# 4. Compute target metrics: Pr(In ER) and Pr(delta improvement).
# 5. Output summary tables and publication-ready ggridges/ggplot2 visualizations.
# ------------------------------------------------------------------------------

# Load required packages
library(MatchIt)
library(cobalt)
library(bayesboot)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(ggpp)

# ==============================================================================
# 1. CORE HELPER FUNCTIONS
# ==============================================================================

#' Probability of Balance Improvement
#' 
#' Calculates the posterior probability that covariate adjustment improved balance
#' by at least a specified threshold (delta).
#' 
#' @param ps_un Vector of posterior draws for unadjusted data
#' @param ps_adj Vector of posterior draws for adjusted data
#' @param dist Minimum meaningful improvement threshold (default = 0.05 for OVL)
#' @return Scalar probability between 0 and 1 rounded to 3 decimal places
pr_improve <- function(ps_un, ps_adj, dist = .05) {
  round(mean(ps_un - ps_adj > dist), 3)
}

# Define arguments for bal.tab to extract multiple distributional moments
bt_args <- list(un = TRUE, abs = FALSE,
                stats = c("mean.diffs", "variance.ratios",
                          "ks.statistics", "ovl.coefficients"))

#' Bootstrap Extraction Function via cobalt
#' 
#' Evaluates bal.tab across bootstrap draws and extracts balance statistics 
#' for all covariates into a single wide data frame.
#' 
#' @param dat Dataset passed to bayesboot
#' @param w Vector of observation weights from bayesboot
#' @param m.out The matchit object to evaluate
#' @param args List of arguments for bal.tab
bt_boot <- function(dat, w, m.out, args){
  bt <- bal.tab(m.out, un = args$un, abs = args$abs,
                stats = args$stats, s.weights = w)
  
  result <- as.data.frame(cbind(t(bt$Balance[ , "Diff.Un"]), t(bt$Balance[ , "Diff.Adj"]),
                                t(bt$Balance[ , "V.Ratio.Un"]), t(bt$Balance[ , "V.Ratio.Adj"]),
                                t(bt$Balance[ , "KS.Un"]), t(bt$Balance[ , "KS.Adj"]),
                                t(bt$Balance[ , "OVL.Un"]), t(bt$Balance[ , "OVL.Adj"])))
  
  colnames(result) <- c(paste(rownames(bt$Balance), "Diff.Un", sep = "_"),
                        paste(rownames(bt$Balance), "Diff.Adj", sep = "_"),
                        paste(rownames(bt$Balance), "V.Ratio.Un", sep = "_"),
                        paste(rownames(bt$Balance), "V.Ratio.Adj", sep = "_"),
                        paste(rownames(bt$Balance), "KS.Un", sep = "_"),
                        paste(rownames(bt$Balance), "KS.Adj", sep = "_"),
                        paste(rownames(bt$Balance), "OVL.Un", sep = "_"),
                        paste(rownames(bt$Balance), "OVL.Adj", sep = "_"))
  
  result <- result %>% mutate_all(~ replace(., is.na(.), 0))
  return(result)
}

# ==============================================================================
# 2. DATA PREPARATION & MATCHING (LALONDE [1986] EXAMPLE)
# ==============================================================================

# Load LaLonde (1986) data from cobalt
data("lalonde", package = "cobalt")

# Define covariates matching the dataset schema
covariates <- c("age", "educ", "race", "married", "nodegree", "re74", "re75")
fml <- as.formula(paste("treat ~", paste(covariates, collapse = " + ")))

# Run 1:1 nearest neighbor matching
m.out <- matchit(fml, data = lalonde, method = "nearest")

# ==============================================================================
# 3. BAYESIAN BOOTSTRAP BALANCE ANALYSIS
# ==============================================================================

set.seed(13355)
iterations <- 5000

cat("Running Bayesian Bootstrap across all covariates...\n")

# Generate posterior distributions for unadjusted and adjusted balance statistics 
# across all covariates simultaneously.
boot.out <- bayesboot(data = lalonde, statistic = bt_boot,
                      R = iterations, .progress = "text", use.weights = TRUE,
                      m.out = m.out, args = bt_args)

# Convert bayesboot output to data.frame for downstream data manipulation
boot.out <- as.data.frame(boot.out)

# ==============================================================================
# 4. SUMMARY TABLE GENERATION
# ==============================================================================

# Extract balance row names generated by cobalt, excluding the distance measure
bt_temp <- bal.tab(m.out, un = TRUE, stats = "ovl.coefficients")
bal_covariates <- setdiff(rownames(bt_temp$Balance), "distance")

# Construct posterior summary table (tab) containing point estimates and 95% CIs
summary_rows <- list()
for (cov in bal_covariates) {
  unadj_draws <- boot.out[[paste0(cov, "_OVL.Un")]]
  adj_draws <- boot.out[[paste0(cov, "_OVL.Adj")]]
  
  summary_rows[[cov]] <- data.frame(
    ov_pm_un = mean(unadj_draws),
    ov_lo_un = quantile(unadj_draws, 0.025),
    ov_hi_un = quantile(unadj_draws, 0.975),
    ov_pm_adj = mean(adj_draws),
    ov_lo_adj = quantile(adj_draws, 0.025),
    ov_hi_adj = quantile(adj_draws, 0.975)
  )
}

tab <- bind_rows(summary_rows)
rownames(tab) <- bal_covariates

# Print diagnostic summary
summary_df <- data.frame(
  Covariate = bal_covariates,
  Unadjusted_OVL = sprintf("%.3f [%.3f, %.3f]", tab$ov_pm_un, tab$ov_lo_un, tab$ov_hi_un),
  Adjusted_OVL = sprintf("%.3f [%.3f, %.3f]", tab$ov_pm_adj, tab$ov_lo_adj, tab$ov_hi_adj),
  Pr_5pp_Improvement = sapply(bal_covariates, function(cov) {
    pr_improve(boot.out[[paste0(cov, "_OVL.Un")]], boot.out[[paste0(cov, "_OVL.Adj")]])
  })
)

cat("\n======================================================================\n")
cat("          BAYESIAN BALANCE DIAGNOSTIC TABLE (LaLonde Data)           \n")
cat("======================================================================\n")
print(summary_df, row.names = FALSE)

# ==============================================================================
# 5. VISUALIZATION (Figure 5 Replication)
# ==============================================================================

# Prepare posterior draw long format data frames (stripping suffixes)
d_ovl1 <- boot.out %>%
  pivot_longer(cols = ends_with("OVL.Un"), names_to = "variable", values_to = "posterior") %>%
  mutate(var = sub("_OVL\\.Un$", "", variable), adjust = "Unadjusted") %>%
  select(var, posterior, adjust) %>%
  filter(var != "distance")

d_ovl2 <- boot.out %>%
  pivot_longer(cols = ends_with("OVL.Adj"), names_to = "variable", values_to = "posterior") %>%
  mutate(var = sub("_OVL\\.Adj$", "", variable), adjust = "Adjusted") %>%
  select(var, posterior, adjust) %>%
  filter(var != "distance")

d_ovl <- bind_rows(d_ovl1, d_ovl2)

# Order the covariates
ordered_vars <- c("age", "educ", "married", "nodegree", 
                  "race_black", "race_hispan", "race_white", 
                  "re74", "re75")

d_ovl$var <- factor(d_ovl$var, levels = ordered_vars)

# Format tab data frame
tab_ovl <- data.frame(
  var = rep(rownames(tab), times = 2),
  pm = c(tab$ov_pm_un, tab$ov_pm_adj),
  lo = c(tab$ov_lo_un, tab$ov_lo_adj),
  hi = c(tab$ov_hi_un, tab$ov_hi_adj),
  adjust = rep(c("Unadjusted", "Adjusted"), each = nrow(tab))
)

# Factor levels
tab_ovl$var <- factor(tab_ovl$var, levels = ordered_vars)

# Calculate numeric position with vertical nudge (0.25) and discrete dodging (+/- 0.1)
tab_ovl$y_num <- as.numeric(tab_ovl$var)
tab_ovl$y_pos <- tab_ovl$y_num + 0.25 + ifelse(tab_ovl$adjust == "Adjusted", -0.1, 0.1)

# Compute probability of improvement labels
imp_labs <- list()
txt <- "Pr(5pp improvement) = "
imp_labs[[1]] <- paste0(txt, pr_improve(boot.out$age_OVL.Un, boot.out$age_OVL.Adj))
imp_labs[[2]] <- paste0(txt, pr_improve(boot.out$educ_OVL.Un, boot.out$educ_OVL.Adj))
imp_labs[[3]] <- paste0(txt, pr_improve(boot.out$married_OVL.Un, boot.out$married_OVL.Adj))
imp_labs[[4]] <- paste0(txt, pr_improve(boot.out$nodegree_OVL.Un, boot.out$nodegree_OVL.Adj))
imp_labs[[5]] <- paste0(txt, pr_improve(boot.out$race_black_OVL.Un, boot.out$race_black_OVL.Adj))
imp_labs[[6]] <- paste0(txt, pr_improve(boot.out$race_hispan_OVL.Un, boot.out$race_hispan_OVL.Adj))
imp_labs[[7]] <- paste0(txt, pr_improve(boot.out$race_white_OVL.Un, boot.out$race_white_OVL.Adj))
imp_labs[[8]] <- paste0(txt, pr_improve(boot.out$re74_OVL.Un, boot.out$re74_OVL.Adj))
imp_labs[[9]] <- paste0(txt, pr_improve(boot.out$re75_OVL.Un, boot.out$re75_OVL.Adj))

# Set visualization theme
theme_set(theme_bw(base_size = 20))

# Ridge plot
ggplot(d_ovl, aes(x = posterior, y = var, fill = adjust)) +
  geom_density_ridges(alpha = .4, scale = 1) +
  labs(
    x = "OVL statistic posterior distributions",
    y = "Covariates"
  ) +
  scale_y_discrete(labels = c("Age", "Education\n (years)", "Married", "No high\n school",
                              "Black", "Hispanic", "White", "Earnings\n (1974)", "Earnings\n (1975)"),
                   expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 1, .1), expand = c(0, 0)) +
  scale_fill_cyclical(
    values = c("#0000ff", "#ff0000", "#8080ff", "#ff8080"),
    name = NULL, guide = "legend"
  ) +
  geom_segment(data = tab_ovl, 
               aes(x = lo, xend = hi, y = y_pos, yend = y_pos, color = adjust),
               inherit.aes = FALSE) +
  geom_point(data = tab_ovl, 
             aes(x = pm, y = y_pos, color = adjust),
             inherit.aes = FALSE) +
  scale_color_manual(name = NULL, values = c("#0000ff", "#ff0000")) + 
  annotate("text", x = .45, y = 1.5, label = imp_labs[[1]]) +
  annotate("text", x = .35, y = 2.5, label = imp_labs[[2]]) +
  annotate("text", x = .525, y = 3.5, label = imp_labs[[3]]) +
  annotate("text", x = .325, y = 4.5, label = imp_labs[[4]]) +
  annotate("text", x = .15, y = 5.5, label = imp_labs[[5]]) +
  annotate("text", x = .35, y = 6.5, label = imp_labs[[6]]) +
  annotate("text", x = .38, y = 7.5, label = imp_labs[[7]]) +
  annotate("text", x = .575, y = 8.5, label = imp_labs[[8]]) +
  annotate("text", x = .4, y = 9.5, label = imp_labs[[9]]) +
  coord_cartesian(clip = "off") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  guides(color = "none")