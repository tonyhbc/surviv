

# `surviv`: Instrumental variable methods for survival analysis

Instrumental variable methods for **time-to-event outcomes** under **unmeasured confounding**, with a focus on IV-based Cox model estimators under various real-world data setting.

---

## Installation

Install the development version from GitHub and load the package:

```r
install.packages("remotes")
remotes::install_github("tonyhbc/surviv")
library({surviv})
```

---

## Analytic goal at a glance

In many survival follow-up studies, treatment is confounded by factors that are not fully measured. **Instrumental variables (IVs)** can enable causal estimation and mitigate bias due to unmeasured confounding when a valid instrument exists (e.g., randomization, provider preference, site-level variation, policy/eligibility enactment), provided standard IV assumptions hold.

This package collects several **IV estimators for Cox-type survival models**, spanning classic baseline IV estimators and time-varying designs for flexible real-world analytics.

---

## Methods implemented

Below are the **four major IV approaches** provided in this package.

### 1) Cox IV via orthogonality method-of-moments (Cox OMOM)

An IV estimator for the Cox proportional hazards model based on self-orthogonality of IV relative to unmeasured confounder to construct MOM-stype estimating equations. [MacKenzie et al. 2014](https://link.springer.com/article/10.1007/s10742-014-0117-x).

- **Estimate:** marginal log hazard ratio
- **Use when:** one has a baseline IV, a numeric/binary baseline treatment, and want a marginal hazard ratio estimand from Cox model.
  
**Main function:** `coxiv_omom()`.

---

### 2) Two-stage residual inclusion with frailty (TSRI-F)

A **two-stage residual inclusion (2SRI)** approach for Cox models that incorporates an **individual frailty term** to absorb unobserved heterogeneity and estimation noise.

- **Estimate:** conditional log hazard ratio
- **Use when:** one has a baseline IV, a numeric/binary baseline treatment, and assume a __constant treatment effect__ over time.

**Main function:** `coxiv_tsrif(..., tvareff = FALSE)`.
  
---

### 3) TSRI-F with a change-point effect
An extension of TSRI-F allowing a **piecewise-constant** (change-point) treatment effect, i.e., early vs late effect.

- **Estimate:** early and late conditional log hazard ratio (before and after a prespecified change point).
- **Use when:** a baseline IV, a numeric/binary baseline treatment, and believes treatment effect is distinct between early and late follow-up period.
- 
**Main function:** `coxiv_tsrif(..., tvareff = TRUE, tchange = t_c)`.

---

### 4) Sequential 2SRI Cox (Seq-2SRI) for longitudinal follow-up data

A **sequential target trial emulation** (STE) approach for treatment initiation over follow-up, combined with a **2SRI control-function** adjustment. This leverages repeated eligibility/landmarking to estimate a treatment effect using stacked “emulated trials” for a time-dependent binary treatment that might initiate any time during follow-up.

- **Estimate:** conditional log hazard ratio (possible to distinguish by treatment start time)
- **Use when:** one has a time-dependent binary treatment that is monotonic (once initiate, always sustain), confounders collected at baseline, a baseline valid IV, and possibly a treatment effect that changes depend on initiation time.
  
**Main workflow (typical):** `seqem()` → `coxiv_seq()` *(replace if your function names differ)*

---

## Methods at a glance

| Method                    | Handles unmeasured confounding via       |        Time-varying treatment effect |              Longitudinal design |     Function call                   |
| ----------------------    | ---------------------------------------  | -----------------------------------: | -------------------------------: | ----------------------------------- |
| Cox OMOM                  | orthogonality / moments                  |                                    ✗ |                                ✗ | `coxiv_omom()`                      |
| TSRI-F Cox (constant)     | residual inclusion + frailty             |                                    ✗ |                                ✗ | `coxiv_tsrif(..., tvareff = FALSE)` |
| TSRI-F Cox (change-point) | residual inclusion + frailty             |                                    ✓ |                                ✗ | `coxiv_tsrif(..., tvareff = TRUE)`  |
| Seq-2SRI Cox              | sequential emulation + residual inclusion| optional (depends on implementation) |                                ✓ | `seqem()` / `coxiv_seq()`           |

---

## Example data

The package includes two example datasets: `VitD` and `vascular` for demonstration purpose.

**VitD**: VitD intake and mortality study.
- *Type:* a classical survival data with baseline covariates, IV, and right-censored survival time. Mirrored from the [`ivtools`](https://cran.r-project.org/web/packages/ivtools/) package.
- *Purpose:* for demonstration of estimators `coxiv_omom()` and `coxiv_tsrif()`.
- *Key elements:* the treatment is vitamin D intake and the IV is mutation status in Flaggrin gene.

**vascular**: Vascular reintervention surgery during follow-up after EVAR procedure 
- *Type:* longitudinal cohort / registry-style survival data (start–stop format) with a time-dependent treatment and baseline IV.
- *Purpose:* for demonstration of estimators `seqem()`, `seqcox()`, and `coxiv_seq()`.
- *Key elements:* the treatment is the time-dependent vascular reintervention surgery that could take place any time during follow-up, including time zero, and assumed to remain treated once initiated. The IV is the center-level
    prevalence of reintervention surgery. The study cohort is a group of abdominal aortic aneurysm patients who just went through a first-line EVAR procedure (time zero).
  
To access the data:

```r
data(VitD, package = "surviv")
data(vascular, package = "surviv")
```

---

## Citation

If you use this package in academic work, please cite:

- ``

---

## Contributing

Issues and pull requests are welcome.  
For bug reports, please include a minimal reproducible example and `sessionInfo()`.

---

## License

