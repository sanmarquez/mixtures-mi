# Exposure-mixture models on multiply-imputed data

Reproducible example of the canonical environmental-mixture trio — **Weighted
Quantile Sum (WQS)**, **quantile g-computation (qgcomp)** and **Bayesian Kernel
Machine Regression (BKMR)** — each combined with **multiple imputation (MICE)**
and **pooled correctly**, on one shared synthetic dataset.

> R · CRAN-only for the three main methods · runs out of the box · one synthetic
> dataset shared across all methods (fixed seeds → fully reproducible).

---

## The problem

Studies of chemical mixtures (metals, EDCs, air pollutants…) face two hard
problems *at the same time*:

1. **Correlated exposures with a joint effect.** Standard one-exposure-at-a-time
   regression cannot answer "what does the *mixture* do, and which components
   drive it?" A family of mixture methods exists (WQS, qgcomp, BKMR), each with
   different assumptions — and they can disagree.
2. **Missing data.** Exposures and covariates are routinely incomplete. Multiple
   imputation (MICE) handles this, but *combining* a mixture model across imputed
   datasets is where analyses quietly go wrong — especially for the Bayesian and
   index methods, where a naïve pooling gives wrong confidence intervals.

This repo shows how to do both correctly, and how to **read the three methods
together** instead of trusting one in isolation.

## Objectives

- Fit WQS, qgcomp and BKMR on the **same** multiply-imputed data.
- **Pool each correctly** across imputations (the part most tutorials skip):
  Rubin's rules with the small-sample degrees-of-freedom correction for the
  index/scalar quantities, posterior-draw pooling for the Bayesian surface,
  averaging for inclusion probabilities.
- **Compare and triangulate**: agreement across methods — and *where* they
  disagree — is what makes a mixture finding credible.

## Methods at a glance

| Method | Estimand | Key assumption | Pooling across imputations |
|---|---|---|---|
| **WQS** (`gWQS`) | one weighted-index coefficient | single effect direction | Rubin on the index coefficient (with `dfcom`) |
| **qgcomp** | signed mixture scalar `psi` | linear, but direction free | Rubin on `psi` (with `dfcom`) |
| **BKMR** (`bkmr`) | flexible response surface + PIPs | none (non-linear, interactions) | PIPs averaged; effect by Rubin (approx) or pooled posterior draws (exact) |

## Repository structure

```
mixtures-mi/
├── R/
│   ├── 00_simulate_data.R    # shared synthetic cohort + MICE (same seed = same data)
│   ├── utils_pooling.R       # Rubin pooling (dfcom), common breaks, common grid
│   ├── 01_qgcomp_mi.R        # quantile g-computation + MI
│   ├── 02_wqs_mi.R           # Weighted Quantile Sum (gWQS) + MI
│   ├── 03_bkmr_mi.R          # BKMR + MI  (overall effect: approx AND exact)
│   ├── 04_compare_methods.R  # compares the three, prints a reading
│   ├── 05_bkmr_causalbkmr.R  # OPTIONAL: cross-check BKMR pooling vs causalbkmr
│   └── 06_bwqs_optional.R    # OPTIONAL: Bayesian WQS (needs a Stan-compatible build)
├── workflow_all_in_one.R     # everything in ONE script (easiest to run)
├── run_all.R                 # same, but by sourcing the modular R/ files
├── LICENSE · .gitignore · README.md
```

## Requirements

```r
install.packages(c("mice", "MASS", "dplyr", "ggplot2", "qgcomp", "gWQS", "bkmr"))
# optional extras:
# install.packages("rstan");  remotes::install_github("ElenaColicino/bwqs")   # 06
# remotes::install_github("zc2326/causalbkmr")                                # 05
```

R ≥ 4.1. The three main methods use only CRAN packages.

## How to run

Easiest — open `workflow_all_in_one.R` and Run Source (it is self-contained,
no relative paths):

```r
source("workflow_all_in_one.R")
```

Or the modular version from the repo root:

```r
source("run_all.R")
```

`DEMO <- TRUE` (top of the script) runs quickly with `m = 5` imputations; set it
to `FALSE` for the real settings (`m = 20`, more MCMC iterations).

## The synthetic data (so the results are checkable)

`n = 250`, five correlated exposures `X1…X5`. In the true model: **X1** has a
non-linear effect, **X2** a positive and **X4** a negative linear effect; **X3**
and **X5** have none. About 5% of values are set missing. So each method can be
checked against a known truth.

## Example output (reproducible demo, `m = 5`, seed 2026)

**Overall mixture effect** — all three agree on a clear positive joint effect:

| method | estimate | 95% CI | excludes 0 |
|---|---|---|---|
| qgcomp (`psi`, +1 quartile) | 0.465 | 0.284 – 0.645 | ✅ |
| BKMR (75th vs 50th) | 0.434 | 0.227 – 0.641 | ✅ |
| WQS (index coefficient) | 0.635 | 0.503 – 0.767 | ✅ |

**Per-exposure importance** — and here the methods *differ informatively*:

| exposure | qgcomp weight | BKMR PIP | WQS weight |
|---|---|---|---|
| X1 | +0.58 | 1.00 | 0.68 |
| X2 | +0.36 | 0.95 | 0.32 |
| **X4** | **−0.72** | **0.94** | **0.00** |
| X3 | −0.28 | 0.43 | 0.00 |
| X5 | +0.06 | 0.35 | 0.01 |

**The reading.** All three recover X1 and X2. The interesting row is **X4**, a
true *negative* driver: qgcomp (−0.72) and BKMR (PIP 0.94) both detect it, but
**WQS misses it entirely (weight 0)** because its single-direction assumption
cannot represent a component acting the opposite way. BKMR's curve also shows the
joint response is not perfectly linear — the shape the index methods assume away.
This is exactly why the three are run together: agreement where they share
assumptions, and a clear, interpretable divergence where they don't.

## How the pooling is done (the decisions worth reviewing)

- **Common quantile scale / grid.** Breaks (WQS, qgcomp) and the BKMR grid are
  fixed once on the stacked imputations, so a "quartile" or a grid point means
  the same thing in every imputed dataset before results are combined.
- **Degrees of freedom.** Rubin pooling uses the Barnard–Rubin small-sample
  correction and is given the complete-data df (`dfcom`); without it the pooled
  df collapse and the CIs come out far too wide. Each script prints a `[df check]`.
- **PIPs are averaged**, not Rubin-pooled (they are probabilities, not estimates
  with a sampling variance).
- **BKMR exact vs approx.** The overall effect is pooled both by Rubin on the
  per-fit summaries (*approx*) and by pooling the full posterior draws (*exact*);
  in the example they agree, which is the robustness check.
- **Number of imputations.** `m = 20` for real use (rule of thumb `m ≥ 100·FMI`);
  `m = 5` in the demo. Too-small `m` adds Monte-Carlo noise, not bias.

## Optional cross-checks

- `05_bkmr_causalbkmr.R` reproduces the BKMR pooling with the official
  `causalbkmr` package (off-CRAN; its `exact` route can fail on some versions —
  take the exact result from `03`).
- `06_bwqs_optional.R` adds Bayesian WQS. The `BWQS` package ships a pre-2.32
  Stan model that fails to compile on `rstan`/`StanHeaders ≥ 2.32`, so it is
  optional; `02_wqs_mi.R` (frequentist WQS) is the robust default.

## References

- Carrico C, Gennings C, Wheeler DC, Factor-Litvak P. *Characterization of WQS
  regression for highly correlated data.* J Agric Biol Environ Stat. 2015;20:100–120.
- Keil AP, et al. *A quantile-based g-computation approach to the effects of
  exposure mixtures.* Environ Health Perspect. 2020;128(4):047004.
- Bobb JF, et al. *Bayesian kernel machine regression…* Biostatistics. 2015;16(3):493–508.
- Bobb JF, et al. *Statistical software for … BKMR.* Environ Health. 2018;17:67.
- van Buuren S, Groothuis-Oudshoorn K. *mice…* J Stat Softw. 2011;45(3).
- Barnard J, Rubin DB. *Small-sample degrees of freedom with multiple
  imputation.* Biometrika. 1999;86(4):948–955.

## License

MIT — see `LICENSE` (add your name/institution).
