# =============================================================================
# workflow_all_in_one.R
# SINGLE self-contained script: simulate -> impute -> qgcomp, WQS, BKMR -> compare.
# No source() of other files, no relative paths. Open it and Run Source (Ctrl+Shift+S),
# or in an R console:  source("workflow_all_in_one.R")
#
# Install once:
#   install.packages(c("mice","MASS","dplyr","ggplot2","qgcomp","gWQS","bkmr"))
# Each method is wrapped so a missing package does NOT stop the rest.
#
# The three methods are the canonical exposure-mixture trio:
#   WQS  (weighted quantile sum, single-direction index)   -> gWQS
#   qgcomp (quantile g-computation, signed index)           -> qgcomp
#   BKMR (flexible kernel response surface)                 -> bkmr
# all combined with multiple imputation (MICE) and pooled correctly.
# =============================================================================

suppressPackageStartupMessages({
  library(MASS); library(mice); library(dplyr); library(qgcomp); library(bkmr)
})
set.seed(2026)

# ---- CONFIG -----------------------------------------------------------------
DEMO  <- TRUE                 # TRUE = quick first run; FALSE = real (m=20, iter high)
m_imp <- if (DEMO) 5    else 20
iter  <- if (DEMO) 2000 else 5000
n <- 250; p <- 5
exp_nms   <- paste0("X", 1:p)
cov_names <- c("Age", "Sex", "Smoking", "BMI")

# ---- Pooling helpers --------------------------------------------------------
barnard_rubin_df <- function(m, b, t, dfcom = Inf) {
  lambda <- (1 + 1/m) * b / t; lambda[lambda < 1e-4] <- 1e-4
  dfold <- (m - 1) / lambda^2
  dfobs <- (dfcom + 1)/(dfcom + 3) * dfcom * (1 - lambda)
  ifelse(is.infinite(dfcom), dfold, dfold * dfobs / (dfold + dfobs))
}
pool_scalar <- function(est, v, dfcom = Inf, level = 0.95) {
  m <- length(est); qbar <- mean(est); ubar <- mean(v); b <- var(est)
  t <- ubar + (1 + 1/m) * b; df <- barnard_rubin_df(m, b, t, dfcom)
  riv <- (1 + 1/m) * b / ubar; fmi <- (riv + 2/(df + 3)) / (riv + 1)
  tcrit <- qt(1 - (1 - level)/2, pmax(df, 1e-3))
  data.frame(estimate = qbar, se = sqrt(t),
             lower = qbar - tcrit*sqrt(t), upper = qbar + tcrit*sqrt(t),
             df = df, fmi = fmi)
}
make_common_breaks <- function(imp, exp_nms, q = 4) {
  s <- complete(imp, "long")
  setNames(lapply(exp_nms, function(v) {
    br <- unique(quantile(s[[v]], probs = seq(0,1,by=1/q), na.rm = TRUE))
    br[1] <- -Inf; br[length(br)] <- Inf; br }), exp_nms)
}

# ---- 0. Synthetic data + MICE ----------------------------------------------
Sigma <- 0.5 ^ abs(outer(1:p, 1:p, "-"))
Z <- MASS::mvrnorm(n, rep(0, p), Sigma); colnames(Z) <- exp_nms
Age <- rnorm(n,50,8); Sex <- rbinom(n,1,.5); Smoking <- rbinom(n,1,.3); BMI <- rnorm(n,27,4)
h <- 0.5*sin(Z[,1]) + 0.4*Z[,2] - 0.3*Z[,4]          # X1 non-linear; X3,X5 null
y <- as.numeric(h + 0.02*(Age-50) + .3*Sex + .4*Smoking + rnorm(n))
dat <- data.frame(y, Age, Sex, Smoking, BMI, Z)
for (v in c(exp_nms,"Age","BMI")) dat[sample(n, round(.05*n)), v] <- NA
imp   <- mice(dat, m = m_imp, method = "pmm", seed = 1234, printFlag = FALSE)
dfcom <- n - (1 + p + length(cov_names))

# =============================================================================
# 1. qgcomp  (signed index; no directional constraint)
# =============================================================================
cat("\n==================== qgcomp ====================\n")
brk  <- make_common_breaks(imp, exp_nms, 4)
form <- reformulate(c(exp_nms, cov_names), response = "y")
qg_fits <- lapply(seq_len(m_imp), function(i)
  qgcomp.glm.noboot(form, expnms = exp_nms, data = complete(imp, i),
                    family = gaussian(), q = NULL, breaks = brk))
qg_psi <- pool_scalar(sapply(qg_fits, function(f) f$psi),
                      sapply(qg_fits, function(f) f$var.psi),
                      dfcom = qg_fits[[1]]$fit$df.residual)
qg_wm <- t(sapply(qg_fits, function(f) {
  v <- c(f$pos.weights, -f$neg.weights)[exp_nms]; replace(v, is.na(v), 0) }))
qg_w <- data.frame(exposure = exp_nms, weight = round(colMeans(qg_wm), 3))
cat(sprintf("[df check] dfcom=%.0f  df=%.1f  (m-1=%d)  fmi=%.2f\n",
            qg_fits[[1]]$fit$df.residual, qg_psi$df, m_imp - 1, qg_psi$fmi))
print(round(qg_psi, 4)); print(qg_w)

# =============================================================================
# 2. WQS  (weighted quantile sum; gWQS). Single-direction index (positive here)
#    -> the directional-homogeneity assumption is the contrast with qgcomp/BKMR.
#    Pooled with Rubin on the WQS index coefficient. Optional: skips if missing.
# =============================================================================
cat("\n==================== WQS (gWQS) ====================\n")
wqs_b <- NULL; wqs_w <- NULL
if (requireNamespace("gWQS", quietly = TRUE)) {
  tryCatch({
    wqs_fits <- lapply(seq_len(m_imp), function(i)
      gWQS::gwqs(reformulate(c("wqs", cov_names), response = "y"),
                 mix_name = exp_nms, data = complete(imp, i), q = 4,
                 validation = 0, b = 100, b1_pos = TRUE,
                 family = "gaussian", seed = 1000 + i))
    bse <- t(sapply(wqs_fits, function(f) {
      co <- summary(f)$coefficients; co["wqs", c("Estimate", "Std. Error")] }))
    wqs_b <- pool_scalar(bse[, 1], bse[, 2]^2, dfcom = dfcom)
    wmat <- sapply(wqs_fits, function(f) {
      fw <- f$final_weights; fw$mean_weight[match(exp_nms, fw$mix_name)] })
    wqs_w <- data.frame(exposure = exp_nms, weight = round(rowMeans(wmat), 3))
    cat(sprintf("[df check] df=%.1f  fmi=%.2f\n", wqs_b$df, wqs_b$fmi))
    print(round(wqs_b, 4)); print(wqs_w)
  }, error = function(e) cat("[WQS skipped]", conditionMessage(e), "\n"))
} else cat("[WQS skipped] package 'gWQS' not installed.\n")

# =============================================================================
# 3. BKMR  (flexible surface; PIPs averaged; overall effect approx AND exact)
# =============================================================================
cat("\n==================== BKMR ====================\n")
s <- complete(imp, "long")
z_center <- sapply(s[, exp_nms], mean); z_scale <- sapply(s[, exp_nms], sd)
z_grid <- setNames(lapply(exp_nms, function(v) {
  zs <- (s[[v]] - z_center[v]) / z_scale[v]
  seq(quantile(zs,.05), quantile(zs,.95), length.out = 50) }), exp_nms)
Zstd <- scale(s[, exp_nms], center = z_center, scale = z_scale)
qpoint <- function(q) matrix(apply(Zstd, 2, quantile, probs = q), nrow = 1,
                             dimnames = list(NULL, exp_nms))
bk_fits <- vector("list", m_imp)
for (i in seq_len(m_imp)) {
  cat("  fitting imputation", i, "of", m_imp, "\n")
  d  <- complete(imp, i)
  Zi <- sweep(sweep(as.matrix(d[, exp_nms]), 2, z_center, "-"), 2, z_scale, "/")
  colnames(Zi) <- exp_nms
  Xi <- model.matrix(reformulate(cov_names), d)[, -1]
  set.seed(1000 + i)
  bk_fits[[i]] <- kmbayes(y = d$y, Z = Zi, X = Xi, iter = iter,
                          family = "gaussian", varsel = TRUE, verbose = FALSE)
}
pips <- data.frame(variable = ExtractPIPs(bk_fits[[1]])$variable,
                   PIP = Reduce("+", lapply(bk_fits, function(f) ExtractPIPs(f)$PIP)) / m_imp)
qs_seq <- seq(0.25, 0.75, by = 0.05)
risk_all <- bind_rows(lapply(seq_along(bk_fits), function(i) {
  r <- OverallRiskSummaries(bk_fits[[i]], qs = qs_seq, q.fixed = 0.5, method = "approx")
  r$imp <- i; r }))
bk_approx <- bind_rows(lapply(qs_seq, function(q) {
  s2 <- risk_all[risk_all$quantile == q, ]
  cbind(quantile = q, pool_scalar(s2$est, s2$sd^2, dfcom = dfcom)[, c("estimate","lower","upper")]) }))
bk_exact <- bind_rows(lapply(qs_seq, function(q) {
  Znew <- rbind(qpoint(q), qpoint(0.50))
  Xnew <- matrix(0, nrow = nrow(Znew), ncol = length(cov_names))
  draws <- unlist(lapply(bk_fits, function(f) {
    sp <- SamplePred(f, Znew = Znew, Xnew = Xnew, sel = NULL); as.numeric(sp[,1]-sp[,2]) }))
  data.frame(quantile = q, estimate = mean(draws),
             lower = quantile(draws,.025), upper = quantile(draws,.975)) }))
cat("\n-- PIPs --\n"); print(pips)
cat("\n-- overall approx --\n"); print(round(bk_approx, 4))
cat("\n-- overall exact  --\n"); print(round(bk_exact, 4))

# =============================================================================
# 4. COMPARE  (direction + drivers; estimands related, not identical)
# =============================================================================
cat("\n==================== COMPARISON ====================\n")
bk75 <- bk_approx[which.min(abs(bk_approx$quantile - .75)), ]
rows <- list("qgcomp (psi, +1 quartile)"  = c(qg_psi$estimate, qg_psi$lower, qg_psi$upper),
             "BKMR approx (75th vs 50th)"  = c(bk75$estimate, bk75$lower, bk75$upper))
if (!is.null(wqs_b))
  rows[["WQS (index coef)"]] <- c(wqs_b$estimate, wqs_b$lower, wqs_b$upper)
overall_tab <- data.frame(method = names(rows),
  estimate = round(sapply(rows, `[`, 1), 4),
  lower = round(sapply(rows, `[`, 2), 4),
  upper = round(sapply(rows, `[`, 3), 4),
  excludes_0 = sapply(rows, function(v) (v[2] > 0) | (v[3] < 0)), row.names = NULL)

imp_tab <- data.frame(exposure = exp_nms,
                      qgcomp_weight = qg_w$weight[match(exp_nms, qg_w$exposure)],
                      bkmr_PIP = pips$PIP[match(exp_nms, pips$variable)])
if (!is.null(wqs_w)) imp_tab$wqs_weight <- wqs_w$weight[match(exp_nms, wqs_w$exposure)]

cat("\n-- overall mixture effect --\n");  print(overall_tab, row.names = FALSE)
cat("\n-- per-exposure importance --\n"); print(imp_tab, row.names = FALSE)
cat("\nReading: compare DIRECTION and which exposures each method flags",
    "(true drivers = X1, X2, X4). WQS assumes one direction; qgcomp allows",
    "both; BKMR adds the non-linear shape.\n")

dir.create("results", showWarnings = FALSE)
write.csv(overall_tab, "results/comparison_overall.csv",    row.names = FALSE)
write.csv(imp_tab,     "results/comparison_importance.csv", row.names = FALSE)
cat("\nDone. Tables written to results/.\n")
