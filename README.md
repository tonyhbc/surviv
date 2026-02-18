

# **{PACKAGE}**: Instrumental variable methods for survival analysis

Causal inference for **time-to-event outcomes** under **unmeasured confounding**, with a focus on IV-based Cox model estimators and longitudinal/sequential trial emulation.

---

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("{OWNER}/{REPO}")
```

Load the package:

```r
library({PACKAGE})  # replace with the actual package name
```

> **Tip:** Replace `{OWNER}`, `{REPO}`, and `{PACKAGE}` throughout this README.

---

## Why this package?

In many observational survival studies, treatment is confounded by factors that are not fully measured. **Instrumental variables (IVs)** can enable causal estimation when a valid instrument exists (e.g., provider preference, site-level variation, policy/eligibility thresholds), provided standard IV assumptions hold.

This package collects several **IV estimators for Cox-type survival models**, spanning classic single-shot IV estimators and modern longitudinal/sequential designs.

---

## Methods implemented

Below are the **four major IV approaches** provided in this package.

### 1) Cox IV via orthogonality / method-of-moments (OMOM)

An IV estimator for the Cox proportional hazards model based on **orthogonality conditions / method-of-moments** ideas (e.g., MacKenzie and related work).

- **Target:** log hazard ratio (Cox PH scale)
- **Use when:** you have a valid IV and want a Cox-model effect without explicitly modeling frailty

**Main function:** `coxiv_omom()` *(replace with your actual function name if different)*

---

### 2) Two-stage residual inclusion with frailty (TSRI-F), constant treatment effect

A **two-stage residual inclusion (2SRI)** approach for Cox models that incorporates an **individual frailty term** to represent unmeasured confounding.

- **Target:** constant log hazard ratio
- **Use when:** unmeasured confounding induces subject-level heterogeneity well captured by frailty

**Main function:** `coxiv_tsrif(..., tvareff = FALSE)` *(or equivalent in your API)*

---

### 3) TSRI-F with a time-varying (piecewise) treatment effect

An extension of TSRI-F allowing a **piecewise-constant** (change-point) treatment effect, e.g., early vs late effect.

- **Target:** two log hazard ratios before/after a prespecified change time
- **Use when:** treatment effects plausibly differ over time (early/late effect)

**Main function:** `coxiv_tsrif(..., tvareff = TRUE, tchange = ...)` *(or equivalent in your API)*

---

### 4) Sequential 2SRI Cox (Seq-2SRI) for longitudinal observational data

A **sequential trial emulation** approach for treatment initiation over follow-up, combined with a **2SRI-style control-function** adjustment. This leverages repeated eligibility/landmarking to estimate a causal effect using stacked “emulated trials”.

- **Target:** causal log hazard ratio under a target-trial style estimand
- **Use when:** treatment can be initiated during follow-up and you want a sequential/landmark-style analysis

**Main workflow (typical):** `seqem()` → `coxiv_seq()` *(replace if your function names differ)*

---

## Methods at a glance

| Method                 | Handles unmeasured confounding via      |        Time-varying treatment effect | Longitudinal / sequential design | Main entry point                    |
| ---------------------- | --------------------------------------- | -----------------------------------: | -------------------------------: | ----------------------------------- |
| OMOM Cox IV            | orthogonality / moments                 | optional (depends on implementation) |                                ✗ | `coxiv_omom()`                      |
| TSRI-F Cox (constant)  | residual inclusion + frailty            |                                    ✗ |                                ✗ | `coxiv_tsrif(..., tvareff = FALSE)` |
| TSRI-F Cox (piecewise) | residual inclusion + frailty            |                                    ✓ |                                ✗ | `coxiv_tsrif(..., tvareff = TRUE)`  |
| Seq-2SRI Cox           | sequential emulation + control-function | optional (depends on implementation) |                                ✓ | `seqem()` / `coxiv_seq()`           |

---

## Quick start (minimal skeleton)

> These are intentionally minimal “shapes” of calls—see package documentation/vignettes for fully runnable examples.

### TSRI-F Cox

```r
# fit <- coxiv_tsrif(
#   formula     = Surv(time, event) ~ A + X1 + X2,
#   trtformula  = A ~ Z + X1 + X2,
#   data        = dat,
#   tvareff     = FALSE
# )
```

### Sequential 2SRI Cox (trial emulation)

```r
# em <- seqem(
#   data  = long_data,
#   id    = "id",
#   tvtrt = "A"
#   # ... additional longitudinal settings ...
# )
#
# fit <- coxiv_seq(
#   formula    = Surv(start, stop, event) ~ A + X1 + X2,
#   trtformula = A ~ Z + X1 + X2,
#   data       = em,
#   id         = "id",
#   tvtrt      = "A",
#   iv         = "Z"
# )
```

---

## Included example data

The package includes a **small, analysis-ready example dataset** intended for demonstrations and unit tests.

- **Type:** longitudinal cohort / registry-style survival data (start–stop or visit-based format)
- **Purpose:** demonstrate IV-based causal estimation of treatment effects on time-to-event outcomes
- **Key elements:** a treatment/exposure of interest, baseline measured covariates, and at least one candidate instrument (an exogenous source of treatment variation)

> **Edit the study description below to match your dataset (2–4 sentences):**  
> “This dataset represents a de-identified cohort of patients followed for a clinical event over time. Treatment initiation may occur during follow-up, and an instrument provides quasi-random variation in treatment assignment. The dataset is included to illustrate the workflow and reproduce package examples; it is not intended for clinical inference.”

To see what datasets are available:

```r
data(package = "{PACKAGE}")
```

Load the example dataset (replace `example_data` with your dataset name):

```r
# data(example_data, package = "{PACKAGE}")
```

---

## Assumptions (brief)

All IV methods rely on standard IV assumptions:

- **Relevance:** the instrument predicts treatment
- **Exclusion restriction:** the instrument affects the outcome only through treatment
- **Independence:** the instrument is independent of unmeasured confounders (possibly conditional on measured covariates)

Each estimator also has method-specific modeling assumptions (e.g., Cox PH structure, frailty specification, or sequential emulation/weighting assumptions).

---

## Citation

If you use this package in academic work, please cite:

- `{Your manuscript / preprint / software paper citation here}`

---

## Contributing

Issues and pull requests are welcome.  
For bug reports, please include a minimal reproducible example and `sessionInfo()`.

---

## License

{Specify license here (e.g., MIT, GPL-3). If you already have a LICENSE file, this can be brief.}
