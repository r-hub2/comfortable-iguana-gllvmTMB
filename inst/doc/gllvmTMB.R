## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.path = "",
  fig.width = 6, fig.height = 4,
  fig.align = "center",
  message = FALSE, warning = FALSE
)

## ----first-model-sketch, eval = FALSE-----------------------------------------
# fit <- gllvmTMB(
#   traits(length, mass, wing, tarsus, bill) ~ 1 +
#     latent(1 | individual, d = 2),
#   data = df_wide,
#   unit = "individual",
#   family = gaussian()
# )

## ----setup--------------------------------------------------------------------
library(gllvmTMB)

## ----example-data-------------------------------------------------------------
example_path <- system.file(
  "extdata", "examples", "morphometrics-example.rds",
  package = "gllvmTMB"
)
if (!nzchar(example_path)) {
  stop("Install gllvmTMB to use its morphometrics example data.")
}
morph <- readRDS(example_path)
df <- morph$data_long
df_wide <- morph$data_wide
truth <- morph$truth
trait_names <- truth$trait_names
n_traits <- length(trait_names)
Lambda_true <- truth$Lambda

## ----example-data-preview-----------------------------------------------------
morph$story$question
head(df_wide, 3)
head(df, 6)

## ----fit-call-----------------------------------------------------------------
fit <- gllvmTMB(
  traits(length, mass, wing, tarsus, bill) ~ 1 +
    latent(1 | individual, d = 2),
  data = df_wide,
  unit = "individual",
  family = gaussian()
)

# A finite likelihood confirms that the model fitted; diagnostics come next.
as.numeric(logLik(fit))

## ----fit-long-----------------------------------------------------------------
fit_long <- gllvmTMB(
  value ~ 0 + trait + latent(0 + trait | individual, d = 2),
  data = df,
  trait = "trait",
  unit = "individual",
  family = gaussian()
)

# Same fit:
all.equal(as.numeric(logLik(fit)), as.numeric(logLik(fit_long)))

## ----fit-health---------------------------------------------------------------
health <- check_gllvmTMB(fit)
health[, c("component", "status", "message", "action")]

## ----residual-table-----------------------------------------------------------
rq <- residuals(
  fit,
  type = "simulation_rank",
  nsim = 199,
  seed = 1,
  condition_on_RE = FALSE
)
head(rq[, c("trait", "family", "observed", "residual", "status")], 8)
stopifnot(diff(range(rq$residual, na.rm = TRUE)) > 1)

## ----residual-qq, fig.width = 6.2, fig.height = 4.2, out.extra='class="wide-scientific-figure"', fig.cap = "Marginal simulation-rank residual Q-Q check for the fitted Gaussian morphometrics model. New random effects are generated for each fitted-model draw. The plot is a diagnostic display, not interval calibration or a posterior predictive check."----
predictive_check(
  fit,
  type = "rq_qq",
  residual_type = "simulation_rank",
  nsim = 199,
  seed = 1,
  condition_on_RE = FALSE
)

## ----sigma-table--------------------------------------------------------------
sigma_rows <- extract_Sigma_table(fit, level = "unit")
sigma_rows

## ----sigma-truth--------------------------------------------------------------
round(truth$Sigma, 2)

## ----ord, fig.cap = "Individual scores on the two latent trait axes."---------
ord <- extract_ordination(fit, level = "unit")
plot(ord$scores[, 1], ord$scores[, 2], pch = 19, col = "steelblue",
     xlab = "LV 1", ylab = "LV 2",
     main = "Individuals on latent trait axes")
abline(h = 0, v = 0, lty = 2, col = "grey70")

## ----communality--------------------------------------------------------------
extract_communality(fit, "unit")

## ----cor----------------------------------------------------------------------
corr_rows <- extract_correlations(fit, tier = "unit")
corr_rows

## ----cor-plot, fig.width = 7, fig.height = 3.8, fig.cap = "Pairwise point estimates of trait correlations from the first fitted model. No interval bounds are requested in this first-fit view."----
plot_correlations(corr_rows, sort = "magnitude")

## ----cor-matrix, fig.width = 7.4, fig.height = 5.1, out.extra='class="wide-scientific-figure"', fig.cap = "Point estimates of trait correlations from the first fitted model. The matrix is descriptive; interval routes are target-specific and opt-in."----
plot_correlations(
  corr_rows,
  style = "heatmap",
  matrix_layout = "by_level",
  label_type = "estimate",
  title = "First-model trait correlations"
)

## ----total-variance-profile---------------------------------------------------
profile_trait <- which.max(diag(truth$Sigma))
total_variance_ci <- profile_ci_total_variance(
  fit,
  tier = "unit",
  trait_idx = profile_trait
)
total_variance_ci[, c(
  "trait", "estimate", "lower", "upper", "interval_status"
)]

## ----loadings-----------------------------------------------------------------
Lambda_hat <- extract_loadings(fit, level = "unit")
rownames(Lambda_hat) <- paste0("trait_", seq_len(n_traits))
colnames(Lambda_hat) <- paste0("LV", seq_len(ncol(Lambda_hat)))
round(Lambda_hat, 3)

## ----cov----------------------------------------------------------------------
list(true_LL  = round(Lambda_true %*% t(Lambda_true), 2),
     fitted_LL = round(Lambda_hat  %*% t(Lambda_hat),  2))

