#' Sequential trial emulation for time-dependent treatment data
#'
#' @description
#' Construct a stacked "sequence of trials" dataset from longitudinal
#' survival data in counting-process form, following the sequential
#' trials emulation framework of Gran et al. and Keogh et al.
#'
#' The function:
#' \enumerate{
#'   \item Identifies trial start times: the earliest cohort start time plus
#'         all times at which any individual initiates treatment (0 -> 1).
#'   \item Augments each individual's trajectory by splitting intervals at
#'         these trial start times so that each at-risk subject has an
#'         explicit row starting at each trial time.
#'   \item For each trial, includes individuals who are under observation at
#'         the trial start and have not yet started treatment, and constructs
#'         trial-specific baseline versions of treatment and covariates
#'         (suffix `.base`), new follow-up times `start.new`, `stop.new`, and
#'         row index `rownum` within `(id, trial)`.
#'   \item Implements artificial censoring by restricting to rows where the
#'         current treatment equals the baseline treatment in that trial.
#' }
#'
#' @param data A `data.frame` in counting-process form with potentially
#'   multiple rows per individual.
#' @param start Name of the start time variable.
#' @param stop Name of the stop time variable.
#' @param event Name of the event indicator (0/1).
#' @param tvtrt Name of the time-varying treatment variable (binary 0/1).
#' @param id Name of the individual identifier.
#' @param covs Optional character vector of covariate names in `data` for which
#'   trial-specific baseline versions should be created. For each `v` in
#'   `covs`, `seqem()` creates `v.base` representing the value of `v` at
#'   the start of each emulated trial for each individual. A baseline version
#'   of `tvtrt` is always created as `tvtrt.base`.
#' @param coarsen Optional character to coarsen the time scale.
#'
#' @return An object of class `"seqem"` with components:
#' \itemize{
#'   \item `data_seqorig` – stacked emulated trials **before** artificial
#'         censoring, with columns:
#'         `id`, `trial`, `trial_time`, `rownum`, `visit`,
#'         original `start`, `stop`, `event`, `tvtrt`, covariates,
#'         `start.new`, `stop.new`, `tvtrt.base`, and `{cov}.base` for
#'         covariates in `covs`.
#'   \item `data_seqem` – stacked emulated trials **after** artificial
#'         censoring, obtained by restricting to rows with
#'         `tvtrt == tvtrt.base`. This is the dataset typically used for
#'         sequential Cox / MSM fitting with IPACW.
#'   \item `summary` – a character string describing the STE object (used by
#'         `print.seqem()`).
#'   \item `trial_summary` – a per-trial summary `data.frame` with one row
#'         per emulated trial and columns:
#'         `trial` (trial index), `t0` (trial start time),
#'         `n` (unique individuals), `rows` (person-time rows),
#'         `treated` (unique baseline-treated individuals),
#'         and `events` (terminal events) in that trial.
#' }
#'
#' Attributes on `data_seqorig` and `data_seqem`:
#' \itemize{
#'   \item `"trial_times"` – vector of trial start times.
#'   \item `"id_var"`, `"start_var"`, `"stop_var"`, `"event_var"`,
#'         `"tvtrt_var"` – names of the corresponding variables.
#' }
#'
#' @export
seqem <- function(data,
                  start,
                  stop,
                  event,
                  tvtrt,
                  id,
                  covs = NULL,
                  coarsen = c("none", "floor", "ceiling"),
                  cbin_width = NULL) {

  coarsen <- match.arg(coarsen)

  ## ---- Basic checks ------------------------------------------------------

  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.")
  }
  for (arg in list(start = start, stop = stop, event = event,
                   tvtrt = tvtrt, id = id)) {
    nm  <- names(arg)
    val <- arg[[1L]]
    if (!is.character(val) || length(val) != 1L) {
      stop(sprintf("Argument '%s' must be a single character string.", nm))
    }
    if (!(val %in% names(data))) {
      stop(sprintf("Variable '%s' not found in 'data'.", val))
    }
  }

  if (!is.null(covs)) {
    if (!is.character(covs)) {
      stop("`covs` must be NULL or a character vector of variable names.")
    }
    missing_covs <- setdiff(covs, names(data))
    if (length(missing_covs) > 0L) {
      stop("The following covariates in `covs` are not in `data`: ",
           paste(missing_covs, collapse = ", "))
    }
  }

  if (coarsen != "none") {
    if (is.null(cbin_width) || !is.numeric(cbin_width) || length(cbin_width) != 1L || cbin_width <= 0) {
      stop("When `coarsen != \"none\"`, `cbin_width` must be a single positive numeric.")
    }
  }

  df <- data

  # Event as integer 0/1
  df[[event]] <- as.integer(df[[event]])

  # Sort by id and start time
  df <- df[order(df[[id]], df[[start]]), ]

  # At most one event per id
  ev_by_id <- tapply(df[[event]], df[[id]], function(x) sum(x, na.rm = TRUE))
  if (any(ev_by_id > 1L, na.rm = TRUE)) {
    stop("More than one event (", event, "==1) detected for some individuals. ",
         "The sequential trial emulation assumes a single terminal event per id.")
  }

  # tvtrt to numeric 0/1 if logical/factor
  if (is.logical(df[[tvtrt]])) {
    df[[tvtrt]] <- as.integer(df[[tvtrt]])
  } else if (is.factor(df[[tvtrt]])) {
    lev <- levels(df[[tvtrt]])
    if (length(lev) != 2L) {
      stop("Time-varying treatment '", tvtrt, "' is a factor with more than two levels. ",
           "Currently only binary treatment (0/1) is supported.")
    }
    df[[tvtrt]] <- as.integer(df[[tvtrt]] == lev[2L])
  }

  # Positive-length intervals
  if (any(df[[stop]] <= df[[start]])) {
    stop("Found rows with 'stop' <= 'start'. Ensure intervals are of positive length.")
  }

  ## ---- 0. Optional discretization / time coarsening (NEW) ---------------

  if (coarsen != "none") {

    origin <- min(df[[start]], na.rm = TRUE)

    snap_floor <- function(t) origin + floor((t - origin) / cbin_width) * cbin_width
    snap_ceil  <- function(t) origin + ceiling((t - origin) / cbin_width) * cbin_width

    snap_start <- if (coarsen == "floor") snap_floor else snap_ceil

    .coarsen_one_id <- function(df_i) {
      df_i <- df_i[order(df_i[[start]]), , drop = FALSE]

      entry <- min(df_i[[start]], na.rm = TRUE)
      exit  <- max(df_i[[stop]],  na.rm = TRUE)

      # True event time (if any): the stop time of the event row
      has_event <- any(df_i[[event]] == 1L, na.rm = TRUE)
      t_event   <- if (has_event) min(df_i[[stop]][df_i[[event]] == 1L], na.rm = TRUE) else Inf

      # Snap ALL row start times (these correspond to *any* TD covariate changes)
      s_raw <- df_i[[start]]
      s_new <- snap_start(s_raw)

      # Never create time before entry; keep the first start at the true entry
      s_new[1L] <- entry
      if (coarsen == "floor") {
        s_new <- pmax(s_new, entry)
      }

      # If multiple original rows map to the same snapped start, keep the LAST one
      # (i.e., latest observed values within that bin before/at the snapped boundary).
      grp <- split(seq_len(nrow(df_i)), s_new)
      keep_idx <- vapply(grp, max, integer(1))
      keep_idx <- keep_idx[order(as.numeric(names(keep_idx)))]

      df_c  <- df_i[keep_idx, , drop = FALSE]
      s_c   <- s_new[keep_idx]

      # Rebuild stop times from the next snapped start; last stop is min(exit, t_event)
      stop_c <- c(s_c[-1L], min(exit, t_event))

      # Assign new times, wipe events then set event once if applicable
      df_c[[start]] <- s_c
      df_c[[stop]]  <- stop_c
      df_c[[event]] <- 0L

      # Drop non-positive intervals (can happen if last snapped start >= exit/t_event)
      ok <- df_c[[stop]] > df_c[[start]]
      df_c <- df_c[ok, , drop = FALSE]
      if (nrow(df_c) == 0L) return(NULL)

      # Place the event in the interval containing t_event (and truncate after it)
      if (is.finite(t_event)) {
        j <- which(df_c[[start]] < t_event & df_c[[stop]] >= t_event)[1L]
        if (!is.na(j)) {
          df_c[[stop]][j]  <- t_event
          df_c[[event]][j] <- 1L
          df_c <- df_c[seq_len(j), , drop = FALSE]
        } else {
          # If no interval contains it (edge cases), ignore; original data should avoid this.
        }
      }

      df_c
    }

    ids <- unique(df[[id]])
    out_list <- vector("list", length(ids))
    for (ii in seq_along(ids)) {
      this_id <- ids[ii]
      df_i <- df[df[[id]] == this_id, , drop = FALSE]
      out_list[[ii]] <- .coarsen_one_id(df_i)
    }
    df2 <- do.call(rbind, Filter(Negate(is.null), out_list))
    rownames(df2) <- NULL
    df2 <- df2[order(df2[[id]], df2[[start]]), , drop = FALSE]

    # Re-check positive-length after coarsening
    if (any(df2[[stop]] <= df2[[start]])) {
      stop("Coarsening produced rows with 'stop' <= 'start'. ",
           "This should not happen; please report this example.")
    }

    df <- df2
  }

  ## ---- 1. Trial start times (baseline + treatment initiations) ----------

  prev_trt <- ave(df[[tvtrt]], df[[id]],
                  FUN = function(z) c(0, head(z, -1L)))
  switch_rows    <- (df[[tvtrt]] == 1L) & (prev_trt == 0L)
  trt_init_times <- df[[start]][switch_rows]

  min_start   <- min(df[[start]], na.rm = TRUE)
  trial_times <- sort(unique(c(min_start, trt_init_times)))
  K           <- length(trial_times)
  if (K < 1L) stop("No valid trial start times found.")

  ## ---- 2. Augment each individual's trajectory at trial_times ----------

  .split_subject <- function(df_id, start_name, stop_name, event_name, cuts_subject) {
    if (length(cuts_subject) == 0L) return(df_id)

    out     <- vector("list", length = 0L)
    idx_out <- 1L

    for (i in seq_len(nrow(df_id))) {
      row <- df_id[i, , drop = FALSE]
      s0  <- row[[start_name]]
      s1  <- row[[stop_name]]
      ev  <- row[[event_name]]

      cuts_row <- cuts_subject[cuts_subject > s0 & cuts_subject < s1]

      if (length(cuts_row) == 0L) {
        out[[idx_out]] <- row
        idx_out <- idx_out + 1L
      } else {
        pts   <- c(s0, cuts_row, s1)
        n_seg <- length(pts) - 1L
        for (j in seq_len(n_seg)) {
          new_row               <- row
          new_row[[start_name]] <- pts[j]
          new_row[[stop_name]]  <- pts[j + 1L]
          if (j < n_seg) {
            new_row[[event_name]] <- 0L
          } else {
            new_row[[event_name]] <- ev
          }
          out[[idx_out]] <- new_row
          idx_out <- idx_out + 1L
        }
      }
    }

    do.call(rbind, out)
  }

  ids         <- unique(df[[id]])
  df_aug_list <- vector("list", length(ids))
  names(df_aug_list) <- as.character(ids)

  for (ii in seq_along(ids)) {
    this_id <- ids[ii]
    df_i    <- df[df[[id]] == this_id, , drop = FALSE]
    df_i    <- df_i[order(df_i[[start]]), ]

    s_min <- min(df_i[[start]])
    s_max <- max(df_i[[stop]])

    cuts_i <- trial_times[trial_times > s_min & trial_times < s_max]

    df_aug_list[[ii]] <- .split_subject(df_i,
                                        start_name = start,
                                        stop_name  = stop,
                                        event_name = event,
                                        cuts_subject = cuts_i)
  }

  df_aug <- do.call(rbind, df_aug_list)
  rownames(df_aug) <- NULL

  df_aug <- df_aug[order(df_aug[[id]], df_aug[[start]]), ]
  df_aug$visit <- ave(rep(1L, nrow(df_aug)), df_aug[[id]], FUN = cumsum)

  ## ---- 3. Precompute entry, exit, first treatment times -----------------

  id_levels <- sort(unique(df_aug[[id]]))

  entry_time <- tapply(df_aug[[start]], df_aug[[id]], min)
  exit_time  <- tapply(df_aug[[stop]],  df_aug[[id]], max)

  first_trt_time <- tapply(df_aug[[start]][df_aug[[tvtrt]] == 1L],
                           df_aug[[id]][df_aug[[tvtrt]] == 1L],
                           min)

  entry_time_all <- setNames(entry_time[match(id_levels, names(entry_time))],
                             id_levels)
  exit_time_all  <- setNames(exit_time[match(id_levels, names(exit_time))],
                             id_levels)
  first_trt_all  <- setNames(rep(Inf, length(id_levels)), id_levels)
  if (length(first_trt_time)) {
    first_trt_all[names(first_trt_time)] <- first_trt_time
  }

  ## ---- 4. Construct sequence of emulated trials (uncensored) -----------

  dat_seq_list <- vector("list", K)
  covs_for_baseline <- unique(c(covs, tvtrt))

  for (k in seq_along(trial_times)) {
    t0          <- trial_times[k]
    trial_index <- k - 1L  # 0-based as in Gran/Keogh

    elig_ids <- id_levels[
      entry_time_all <= t0 &
        exit_time_all  > t0 &
        first_trt_all  >= t0
    ]
    if (!length(elig_ids)) next

    idx_k   <- df_aug[[id]] %in% elig_ids & df_aug[[start]] >= t0
    trial_df <- df_aug[idx_k, , drop = FALSE]
    if (!nrow(trial_df)) next

    trial_df$trial      <- trial_index
    trial_df$trial_time <- t0

    trial_df <- trial_df[order(trial_df[[id]], trial_df[[start]]), ]

    trial_df$rownum <- ave(seq_len(nrow(trial_df)),
                           interaction(trial_df[[id]], trial_df$trial, drop = TRUE),
                           FUN = seq_along)

    # Baseline versions of tvtrt and covariates at trial start
    first_rows <- !duplicated(trial_df[[id]])
    for (v in covs_for_baseline) {
      if (!v %in% names(trial_df)) {
        stop("Baseline covariate '", v, "' not found in trial data. ",
             "Check that it exists in the original data.")
      }
      base_vals <- trial_df[[v]][first_rows]
      names(base_vals) <- as.character(trial_df[[id]][first_rows])
      v_base <- base_vals[as.character(trial_df[[id]])]

      trial_df[[paste0(v, ".base")]] <- v_base
    }

    trial_df$start.new <- trial_df[[start]] - t0
    trial_df$stop.new  <- trial_df[[stop]]  - t0

    dat_seq_list[[k]] <- trial_df
  }

  dat_seq_list <- Filter(Negate(is.null), dat_seq_list)
  if (length(dat_seq_list) == 0L) {
    stop("No eligible individuals for any emulated trial after applying ",
         "entry/exit and 'not previously treated' criteria.")
  }

  dat_seq_orig <- do.call(rbind, dat_seq_list)
  rownames(dat_seq_orig) <- NULL

  ## ---- 5. Artificial censoring (MANDATORY) -----------------------------

  tv_base_name <- paste0(tvtrt, ".base")
  if (!tv_base_name %in% names(dat_seq_orig)) {
    stop("Baseline treatment variable '", tv_base_name,
         "' not found in stacked trials.")
  }

  keep    <- dat_seq_orig[[tvtrt]] == dat_seq_orig[[tv_base_name]]
  dat_seq <- dat_seq_orig[keep, , drop = FALSE]

  ## ---- 6. Attach attributes & build summary ----------------------------

  attr(dat_seq_orig, "trial_times") <- trial_times
  attr(dat_seq_orig, "id_var")      <- id
  attr(dat_seq_orig, "start_var")   <- start
  attr(dat_seq_orig, "stop_var")    <- stop
  attr(dat_seq_orig, "event_var")   <- event
  attr(dat_seq_orig, "tvtrt_var")   <- tvtrt
  attr(dat_seq_orig, "coarsen")     <- coarsen
  attr(dat_seq_orig, "cbin_width")  <- cbin_width

  attr(dat_seq, "trial_times") <- trial_times
  attr(dat_seq, "id_var")      <- id
  attr(dat_seq, "start_var")   <- start
  attr(dat_seq, "stop_var")    <- stop
  attr(dat_seq, "event_var")   <- event
  attr(dat_seq, "tvtrt_var")   <- tvtrt
  attr(dat_seq, "coarsen")     <- coarsen
  attr(dat_seq, "cbin_width")  <- cbin_width

  ## ---- 6b. Trial-level summary dataset ---------------------------------

  trial_var <- "trial"
  trial_summary <- NULL

  if (trial_var %in% names(dat_seq)) {
    trial_times_attr <- attr(dat_seq, "trial_times")
    id_var_attr      <- attr(dat_seq, "id_var")
    event_var_attr   <- attr(dat_seq, "event_var")
    tvtrt_var_attr   <- attr(dat_seq, "tvtrt_var")

    trials_present <- sort(unique(dat_seq[[trial_var]]))

    if (!is.null(id_var_attr) &&
        !is.null(event_var_attr) &&
        id_var_attr %in% names(dat_seq) &&
        event_var_attr %in% names(dat_seq)) {

      ids_per_trial <- tapply(dat_seq[[id_var_attr]], dat_seq[[trial_var]],
                              function(z) length(unique(z)))
      rows_per_trial <- tapply(seq_len(nrow(dat_seq)), dat_seq[[trial_var]], length)
      events_per_trial <- tapply(dat_seq[[event_var_attr]], dat_seq[[trial_var]],
                                 function(z) sum(z, na.rm = TRUE))

      # treated: number of unique ids with baseline treatment = 1 in each trial
      tv_base_name_attr <- if (!is.null(tvtrt_var_attr)) paste0(tvtrt_var_attr, ".base") else NULL
      treated_raw <- NULL
      if (!is.null(tv_base_name_attr) && tv_base_name_attr %in% names(dat_seq)) {
        treated_raw <- tapply(
          dat_seq[[id_var_attr]][dat_seq[[tv_base_name_attr]] == 1],
          dat_seq[[trial_var]][dat_seq[[tv_base_name_attr]] == 1],
          function(z) length(unique(z))
        )
      }

      align_vec <- function(v) v[match(trials_present, as.numeric(names(v)))]

      ids_aligned    <- align_vec(ids_per_trial)
      rows_aligned   <- align_vec(rows_per_trial)
      events_aligned <- align_vec(events_per_trial)
      if (!is.null(treated_raw)) {
        treated_aligned <- align_vec(treated_raw)
        treated_aligned[is.na(treated_aligned)] <- 0L
      } else {
        treated_aligned <- rep(NA_integer_, length(trials_present))
      }

      t0s <- if (!is.null(trial_times_attr)) {
        trial_times_attr[pmin(trials_present + 1L, length(trial_times_attr))]
      } else {
        rep(NA_real_, length(trials_present))
      }

      trial_summary <- data.frame(
        trial   = trials_present,
        t0      = t0s,
        n       = as.numeric(ids_aligned),
        rows    = as.numeric(rows_aligned),
        treated = as.numeric(treated_aligned),
        events  = as.numeric(events_aligned)
      )
      rownames(trial_summary) <- NULL
    }
  }

  summary_text <- .seqem_summary_text(dat_seq)

  out <- list(
    data_seqorig   = dat_seq_orig,
    data_seqem     = dat_seq,
    summary        = summary_text,
    trial_summary  = trial_summary
  )
  class(out) <- "seqem"
  out
}


# Internal helper: build multi-line summary string for a seqem object
# Internal helper: build multi-line summary string for a seqem object
.seqem_summary_text <- function(dat) {
  trial_times <- attr(dat, "trial_times")
  id_var      <- attr(dat, "id_var")
  start_var   <- attr(dat, "start_var")
  stop_var    <- attr(dat, "stop_var")
  event_var   <- attr(dat, "event_var")
  tvtrt_var   <- attr(dat, "tvtrt_var")
  coarsen     <- attr(dat, "coarsen")
  cbin_width  <- attr(dat, "cbin_width")

  n_rows <- nrow(dat)
  n_ids  <- if (!is.null(id_var) && id_var %in% names(dat)) {
    length(unique(dat[[id_var]]))
  } else NA_integer_

  if (!("trial" %in% names(dat))) {
    return("Object lacks 'trial' variable; not a valid STE dataset.\n")
  }

  trials_present   <- sort(unique(dat[["trial"]]))
  n_trials_present <- length(trials_present)
  n_trials_total   <- if (!is.null(trial_times)) length(trial_times) else n_trials_present

  buf <- character()

  buf <- c(buf,
           "Sequential trial emulation data (seqem)\n",
           "--------------------------------------------------\n")
  buf <- c(buf, sprintf("Rows in data_seqem : %d\n", n_rows))
  if (!is.na(n_ids)) {
    buf <- c(buf, sprintf("Ids in data_seqem  : %d\n", n_ids))
  }
  line_trials <- sprintf("Trials (present)   : %d", n_trials_present)
  if (!is.null(trial_times)) {
    line_trials <- paste0(
      line_trials,
      sprintf("  (%d trial start time%s in 'trial_times')",
              n_trials_total,
              if (n_trials_total == 1L) "" else "s")
    )
  }
  buf <- c(buf, line_trials, "\n")

  if (!is.null(coarsen)) {
    if (coarsen == "none") {
      buf <- c(buf, "Time coarsening      : none\n")
    } else {
      buf <- c(buf, sprintf("Time coarsening      : %s (bin width = %s)\n",
                            coarsen, as.character(cbin_width)))
    }
  }

  if (!is.null(id_var)) {
    buf <- c(buf, "Core variables:\n")
    buf <- c(buf, sprintf("  id    : %s\n", id_var))
    buf <- c(buf, sprintf("  start : %s\n", start_var))
    buf <- c(buf, sprintf("  stop  : %s\n", stop_var))
    buf <- c(buf, sprintf("  event : %s\n", event_var))
    buf <- c(buf, sprintf("  tvtrt : %s\n", tvtrt_var))
    buf <- c(buf, "  trial : integer trial index (0-based)\n")

    tv_base_name <- if (!is.null(tvtrt_var)) paste0(tvtrt_var, ".base") else NULL
    if (!is.null(tv_base_name) && tv_base_name %in% names(dat)) {
      buf <- c(buf,
               sprintf("  Baseline treatment variable: %s\n", tv_base_name))
    }
    if (all(c("start.new", "stop.new") %in% names(dat))) {
      buf <- c(buf, "  start.new, stop.new present.\n")
    }
  }

  # Per-trial summary (top 10)
  if (n_trials_present > 0L &&
      !is.null(id_var) && !is.null(event_var) &&
      id_var %in% names(dat) && event_var %in% names(dat)) {

    buf <- c(buf, "\nPer-trial summary (present trials):\n")

    ids_per_trial <- tapply(dat[[id_var]], dat[["trial"]],
                            function(z) length(unique(z)))
    rows_per_trial <- tapply(seq_len(nrow(dat)), dat[["trial"]], length)
    events_per_trial <- tapply(dat[[event_var]], dat[["trial"]],
                               function(z) sum(z, na.rm = TRUE))

    # treated: unique ids with baseline treatment = 1 in each trial
    tv_base_name <- if (!is.null(tvtrt_var)) paste0(tvtrt_var, ".base") else NULL
    treated_raw <- NULL
    if (!is.null(tv_base_name) && tv_base_name %in% names(dat)) {
      treated_raw <- tapply(
        dat[[id_var]][dat[[tv_base_name]] == 1],
        dat[["trial"]][dat[[tv_base_name]] == 1],
        function(z) length(unique(z))
      )
    }

    align_vec <- function(v) v[match(trials_present, as.numeric(names(v)))]

    ids_per_trial    <- align_vec(ids_per_trial)
    rows_per_trial   <- align_vec(rows_per_trial)
    events_per_trial <- align_vec(events_per_trial)

    if (!is.null(treated_raw)) {
      treated_per_trial <- align_vec(treated_raw)
      treated_per_trial[is.na(treated_per_trial)] <- 0L
    } else {
      treated_per_trial <- rep(NA_integer_, n_trials_present)
    }

    t0s <- if (!is.null(trial_times)) {
      trial_times[pmin(trials_present + 1L, length(trial_times))]
    } else {
      rep(NA_real_, n_trials_present)
    }

    max_show <- min(n_trials_present, 10L)
    for (j in seq_len(max_show)) {
      tr     <- trials_present[j]
      t0_str <- if (is.na(t0s[j])) "NA" else format(t0s[j], digits = 6, trim = TRUE)
      buf <- c(buf,
               sprintf("  trial %d (t0 = %s): rows = %d, ids = %d, treated = %d, events = %d\n",
                       tr, t0_str,
                       rows_per_trial[j], ids_per_trial[j],
                       treated_per_trial[j], events_per_trial[j]))
    }
    if (n_trials_present > 10L) {
      extra <- n_trials_present - 10L
      buf <- c(buf,
               sprintf("  ... (%d additional trial%s not shown)\n",
                       extra, if (extra == 1L) "" else "s"))
    }
  }

  paste0(buf, collapse = "")
}

