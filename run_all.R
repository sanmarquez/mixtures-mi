# =============================================================================
# run_all.R  —  run the whole workflow on the same synthetic data, then compare.
# Run from the repository root:  source("run_all.R")
# Each step is isolated: if one method's package is missing, the rest continue.
# Order matters: 04 (compare) reads the CSVs written by 01-03, so it runs last.
# =============================================================================

run_step <- function(label, path) tryCatch({
  message("\n==================== ", label, " ====================")
  source(path, local = new.env())
}, error = function(e) message("[", label, " failed] ", conditionMessage(e)))

run_step("1  qgcomp",                       "R/01_qgcomp_mi.R")
run_step("2  WQS (gWQS)",                   "R/02_wqs_mi.R")
run_step("3  BKMR",                         "R/03_bkmr_mi.R")
run_step("4  causalbkmr cross-check (opt)", "R/05_bkmr_causalbkmr.R")
run_step("5  COMPARE methods",              "R/04_compare_methods.R")
# Optional: BWQS (needs a Stan-compatible BWQS build) -> R/06_bwqs_optional.R

message("\nDone. Tables in results/, figures in figures/.")
