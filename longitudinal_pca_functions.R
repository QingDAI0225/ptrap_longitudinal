# Longitudinal CLR-PCA analysis functions

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(mgcv)
})

# -----------------------------
# small utilities
# -----------------------------
rmse_vec <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
mae_vec  <- function(obs, pred) mean(abs(obs - pred), na.rm = TRUE)
r2_vec   <- function(obs, pred) {
  ss_res <- sum((obs - pred)^2, na.rm = TRUE)
  ss_tot <- sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
  if (!is.finite(ss_tot) || isTRUE(all.equal(ss_tot, 0))) return(NA_real_)
  1 - ss_res / ss_tot
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

pct_overlap <- function(x, y) {
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  denom <- length(union(x, y))
  if (denom == 0) return(NA_real_)
  length(intersect(x, y)) / denom
}

assert_has_columns <- function(df, required_cols, df_name = deparse(substitute(df))) {
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(df_name, " is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  invisible(TRUE)
}

resolve_feature_id_col <- function(clr_tbl, feature_id_col = "id") {
  if (!is.null(feature_id_col) && feature_id_col %in% colnames(clr_tbl)) {
    return(feature_id_col)
  }
  NULL
}

resolve_analysis_axes <- function(axis_df = NULL, axes = NULL, n_pc = NULL, max_axes = 2) {
  if (!is.null(axes)) return(axes)
  
  if (!is.null(axis_df)) {
    pc_cols <- grep("^PC\\d+$", colnames(axis_df), value = TRUE)
    pc_cols <- pc_cols[order(as.integer(sub("^PC", "", pc_cols)))]
    return(head(pc_cols, max_axes))
  }
  
  if (!is.null(n_pc)) {
    return(paste0("PC", seq_len(min(n_pc, max_axes))))
  }
  
  c("PC1", "PC2")
}

compute_k_for_gam <- function(axis_df, group_col = NULL, k_time = 6, min_k = 3) {
  assert_has_columns(axis_df, "day_num", "axis_df")
  
  n_day <- dplyr::n_distinct(axis_df$day_num[is.finite(axis_df$day_num)])
  if (!is.finite(n_day) || n_day < 2) {
    stop("Need at least 2 distinct day_num values to fit a temporal GAM.")
  }
  
  k_use <- min(k_time, n_day - 1)
  
  if (!is.null(group_col) && group_col %in% colnames(axis_df)) {
    n_min_group_day <- axis_df %>%
      dplyr::filter(is.finite(day_num)) %>%
      dplyr::group_by(.data[[group_col]]) %>%
      dplyr::summarise(n_day = dplyr::n_distinct(day_num), .groups = "drop") %>%
      dplyr::pull(n_day) %>%
      min(na.rm = TRUE)
    
    k_use <- min(k_use, n_min_group_day - 1)
  }
  
  if (!is.finite(k_use) || k_use < min_k) {
    stop("Insufficient unique time points to fit GAM with requested grouping structure.")
  }
  
  k_use
}

# -----------------------------
# PCA helpers
# -----------------------------
get_clr_sample_matrix <- function(clr_tbl, id_cols, sample_names = NULL,
                                  feature_id_col = "id") {
  feature_id_col <- resolve_feature_id_col(clr_tbl, feature_id_col)
  sample_cols <- setdiff(colnames(clr_tbl), id_cols)
  
  if (!is.null(sample_names)) {
    missing_samples <- setdiff(sample_names, sample_cols)
    if (length(missing_samples) > 0) {
      stop("Missing sample columns: ", paste(missing_samples, collapse = ", "))
    }
    sample_cols <- sample_names
  }
  
  X <- as.matrix(clr_tbl[, sample_cols, drop = FALSE])
  storage.mode(X) <- "double"
  
  if (!is.null(feature_id_col)) {
    rownames(X) <- clr_tbl[[feature_id_col]]
  }
  
  t(X)
}

fit_pca_from_clr <- function(clr_tbl, id_cols, sample_names = NULL, n_pc = 5,
                             feature_id_col = "id") {
  feature_id_col <- resolve_feature_id_col(clr_tbl, feature_id_col)
  X <- get_clr_sample_matrix(
    clr_tbl,
    id_cols = id_cols,
    sample_names = sample_names,
    feature_id_col = feature_id_col
  )
  
  pca <- prcomp(X, center = TRUE, scale. = FALSE)
  keep_pc <- seq_len(min(n_pc, ncol(pca$x)))
  
  scores <- as.data.frame(pca$x[, keep_pc, drop = FALSE]) %>%
    rownames_to_column("Sample.Name")
  colnames(scores)[-1] <- paste0("PC", seq_len(ncol(scores) - 1))
  
  loadings <- as.data.frame(pca$rotation[, keep_pc, drop = FALSE]) %>%
    rownames_to_column("feature")
  colnames(loadings)[-1] <- paste0("PC", seq_len(ncol(loadings) - 1))
  
  if (!is.null(feature_id_col)) {
    feat_map <- clr_tbl %>% select(any_of(c(feature_id_col, "Taxon"))) %>% distinct()
    loadings <- loadings %>% left_join(feat_map, by = c("feature" = feature_id_col))
  }
  
  list(pca = pca, scores = scores, loadings = loadings)
}

align_loading_sign <- function(ref_loadings, new_loadings, axes = c("PC1", "PC2")) {
  out <- new_loadings
  sign_tbl <- tibble(axis = axes, sign_flip = 1)
  
  for (ax in axes) {
    ref_ax <- ref_loadings %>% select(feature, all_of(ax))
    new_ax <- out %>% select(feature, all_of(ax))
    tmp <- inner_join(ref_ax, new_ax, by = "feature", suffix = c("_ref", "_new"))
    r <- safe_cor(tmp[[paste0(ax, "_ref")]], tmp[[paste0(ax, "_new")]])
    if (is.finite(r) && r < 0) {
      out[[ax]] <- -out[[ax]]
      sign_tbl$sign_flip[sign_tbl$axis == ax] <- -1
    }
  }
  list(loadings = out, sign_tbl = sign_tbl)
}

fit_pca_from_matrix <- function(X, n_pc = 5, feat_map = NULL) {
  stopifnot(is.matrix(X))
  storage.mode(X) <- "double"
  
  pca <- prcomp(X, center = TRUE, scale. = FALSE)
  keep_pc <- seq_len(min(n_pc, ncol(pca$x)))
  
  scores <- as.data.frame(pca$x[, keep_pc, drop = FALSE]) %>%
    tibble::rownames_to_column("Sample.Name")
  colnames(scores)[-1] <- paste0("PC", seq_len(ncol(scores) - 1))
  
  loadings <- as.data.frame(pca$rotation[, keep_pc, drop = FALSE]) %>%
    tibble::rownames_to_column("feature")
  colnames(loadings)[-1] <- paste0("PC", seq_len(ncol(loadings) - 1))
  
  if (!is.null(feat_map)) {
    loadings <- loadings %>% left_join(feat_map, by = "feature")
  }
  
  list(pca = pca, scores = scores, loadings = loadings)
}

get_top_features <- function(loadings_df, axis = "PC1", top_k = 20) {
  loadings_df %>%
    mutate(abs_loading = abs(.data[[axis]])) %>%
    arrange(desc(abs_loading)) %>%
    slice_head(n = top_k) %>%
    pull(feature)
}

build_axis_df <- function(clr_tbl, metadata_batch, id_cols, n_pc = 5,
                          X_all = NULL, feat_map = NULL) {
  assert_has_columns(metadata_batch, c("Sample.Name", "ptrap", "date"), "metadata_batch")
  
  if (is.null(X_all)) {
    fit <- fit_pca_from_clr(
      clr_tbl = clr_tbl,
      id_cols = id_cols,
      sample_names = metadata_batch$Sample.Name,
      n_pc = n_pc,
      feature_id_col = if ("id" %in% id_cols) "id" else NULL
    )
  } else {
    X_use <- X_all[metadata_batch$Sample.Name, , drop = FALSE]
    fit <- fit_pca_from_matrix(
      X = X_use,
      n_pc = n_pc,
      feat_map = feat_map
    )
  }
  
  axis_df <- metadata_batch %>%
    left_join(fit$scores, by = "Sample.Name") %>%
    arrange(ptrap, date) %>%
    mutate(
      ptrap = as.integer(ptrap),
      ptrap_f = factor(ptrap),
      day_num = as.numeric(date - min(date, na.rm = TRUE))
    )
  
  list(axis_df = axis_df, fit = fit)
}

prepare_pca_cache <- function(clr_tbl, id_cols, sample_names, feature_id_col = "id") {
  feature_id_col <- resolve_feature_id_col(clr_tbl, feature_id_col)
  
  X_all <- get_clr_sample_matrix(
    clr_tbl = clr_tbl,
    id_cols = id_cols,
    sample_names = sample_names,
    feature_id_col = feature_id_col
  )
  
  feat_map <- NULL
  if (!is.null(feature_id_col)) {
    feat_map <- clr_tbl %>%
      select(any_of(c(feature_id_col, "Taxon"))) %>%
      distinct() %>%
      rename(feature = all_of(feature_id_col))
  }
  
  list(
    X_all = X_all,
    feat_map = feat_map,
    feature_id_col = feature_id_col
  )
}

add_treatment_labels <- function(axis_df,
                                 treated_ptraps = c(1, 2),
                                 control_ptrap = 3) {
  axis_df %>%
    mutate(
      treatment = case_when(
        ptrap %in% treated_ptraps ~ "Treated",
        ptrap == control_ptrap    ~ "Control",
        TRUE ~ "Other"
      ),
      treatment = factor(treatment, levels = c("Control", "Treated", "Other"))
    )
}

# -----------------------------
# trajectory summaries and GAMs
# -----------------------------
summarize_axis_shift <- function(axis_df, axes = c("PC1", "PC2")) {
  axis_df %>%
    arrange(ptrap, date) %>%
    group_by(ptrap, treatment) %>%
    summarise(
      across(
        all_of(axes),
        list(initial = first,
             final = last,
             delta = ~ dplyr::last(.) - dplyr::first(.)),
        .names = "{.col}_{.fn}"
      ),
      n_timepoints = n(),
      .groups = "drop"
    )
}

fit_axis_gams <- function(axis_df, axes = c("PC1", "PC2"),
                          mode = c("treatment", "ptrap"),
                          k_time = 6) {
  mode <- match.arg(mode)
  out <- vector("list", length(axes))
  names(out) <- axes
  
  if (mode == "treatment") {
    k_use <- compute_k_for_gam(axis_df, group_col = "treatment", k_time = k_time)
    formula_builder <- function(ax) {
      as.formula(
        paste0(
          ax,
          " ~ treatment + s(day_num, by = treatment, k = ", k_use, ") + s(ptrap_f, bs = 're')"
        )
      )
    }
  } else {
    k_use <- compute_k_for_gam(axis_df, group_col = "ptrap_f", k_time = k_time)
    formula_builder <- function(ax) {
      as.formula(paste0(ax, " ~ ptrap_f + s(day_num, by = ptrap_f, k = ", k_use, ")"))
    }
  }
  
  for (ax in axes) {
    out[[ax]] <- mgcv::gam(formula_builder(ax), data = axis_df, method = "REML")
  }
  out
}

fit_treatment_gams <- function(axis_df, axes = c("PC1", "PC2"), k_time = 6) {
  fit_axis_gams(axis_df, axes = axes, mode = "treatment", k_time = k_time)
}

fit_ptrap_gams <- function(axis_df, axes = c("PC1", "PC2"), k_time = 6) {
  fit_axis_gams(axis_df, axes = axes, mode = "ptrap", k_time = k_time)
}

augment_gam_preds <- function(axis_df, gam_list, suffix = "gam") {
  out <- axis_df
  for (nm in names(gam_list)) {
    out[[paste0(nm, "_", suffix, "_pred")]] <- predict(gam_list[[nm]], newdata = out)
  }
  out
}

summarize_gams <- function(gam_list) {
  purrr::imap_dfr(gam_list, function(fit, nm) {
    sm <- summary(fit)
    s_tbl <- as.data.frame(sm$s.table)
    
    if (nrow(s_tbl) == 0) {
      smooth_terms <- NA_character_
      smooth_p_values <- NA_character_
    } else {
      s_tbl$term <- rownames(sm$s.table)
      smooth_terms <- paste(s_tbl$term, collapse = "; ")
      smooth_p_values <- paste(signif(s_tbl$`p-value`, 3), collapse = "; ")
    }
    
    tibble(
      axis = nm,
      deviance_explained = sm$dev.expl,
      adj_r_sq = sm$r.sq,
      smooth_terms = smooth_terms,
      smooth_p_values = smooth_p_values
    )
  })
}

plot_axis_trajectory <- function(axis_df, axis = "PC1",
                                 pred_col = paste0(axis, "_ptrap_gam_pred")) {
  assert_has_columns(axis_df, c("date", axis, pred_col, "ptrap", "treatment"), "axis_df")
  
  ggplot(axis_df, aes(x = date, y = .data[[axis]], color = factor(ptrap), shape = treatment)) +
    geom_point(size = 2) +
    geom_line(aes(y = .data[[pred_col]], group = ptrap), linewidth = 0.8) +
    labs(x = "Date", y = axis, color = "P-trap", shape = "Group",
         title = paste0(axis, " trajectory by p-trap")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_treatment_smooth <- function(axis_df, gam_list, axis = "PC1", n_grid = 200,
                                  dataset_label = NULL) {
  assert_has_columns(axis_df, c("date", "day_num", axis, "ptrap_f", "treatment"), "axis_df")
  
  stopifnot(axis %in% names(gam_list))
  fit <- gam_list[[axis]]
  
  date0 <- min(axis_df$date, na.rm = TRUE)
  
  newdat <- tidyr::expand_grid(
    day_num = seq(min(axis_df$day_num, na.rm = TRUE),
                  max(axis_df$day_num, na.rm = TRUE),
                  length.out = n_grid),
    treatment = factor(c("Control", "Treated"),
                       levels = levels(axis_df$treatment))
  ) %>%
    dplyr::mutate(
      ptrap_f = factor(levels(axis_df$ptrap_f)[1], levels = levels(axis_df$ptrap_f)),
      date = date0 + day_num
    )
  
  pred <- predict(
    fit,
    newdata = newdat,
    se.fit = TRUE,
    exclude = "s(ptrap_f)"
  )
  
  newdat <- newdat %>%
    dplyr::mutate(
      fit = as.numeric(pred$fit),
      se = as.numeric(pred$se.fit),
      lower = fit - 1.96 * se,
      upper = fit + 1.96 * se
    )
  
  plot_title <- if (is.null(dataset_label)) {
    paste0(axis, " treatment-level smooth trajectories")
  } else {
    paste0(dataset_label, ": ", axis, " treatment-level smooth trajectories")
  }
  
  ggplot(axis_df, aes(x = date, y = .data[[axis]], color = treatment)) +
    geom_point(size = 1.8, alpha = 0.75) +
    geom_ribbon(
      data = newdat,
      aes(x = date, ymin = lower, ymax = upper, fill = treatment),
      inherit.aes = FALSE,
      alpha = 0.18,
      color = NA
    ) +
    geom_line(
      data = newdat,
      aes(x = date, y = fit, color = treatment),
      inherit.aes = FALSE,
      linewidth = 1.1
    ) +
    labs(
      x = "Date",
      y = axis,
      color = "Group",
      fill = "Group",
      title = plot_title
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# -----------------------------
# treated-to-treated concordance
# -----------------------------
calc_treated_concordance <- function(axis_df,
                                     treated_ptraps = c(1, 2),
                                     axes = c("PC1", "PC2"),
                                     use_pred = TRUE) {
  stopifnot(length(treated_ptraps) == 2)
  
  collapse_one_per_day <- function(df) {
    pred_cols <- grep("_ptrap_gam_pred$", colnames(df), value = TRUE)
    
    df %>%
      group_by(date, ptrap) %>%
      summarise(
        across(all_of(c(axes, pred_cols)), ~ mean(.x, na.rm = TRUE)),
        .groups = "drop"
      )
  }
  
  d1 <- axis_df %>% filter(ptrap == treated_ptraps[1]) %>% collapse_one_per_day()
  d2 <- axis_df %>% filter(ptrap == treated_ptraps[2]) %>% collapse_one_per_day()
  
  merged <- inner_join(d1, d2, by = "date", suffix = c("_1", "_2"))
  
  if (nrow(merged) < 3) {
    return(tibble(axis = axes, obs_cor = NA_real_, pred_cor = NA_real_, n_overlap_dates = nrow(merged)))
  }
  
  purrr::map_dfr(axes, function(ax) {
    obs_cor <- safe_cor(merged[[paste0(ax, "_1")]], merged[[paste0(ax, "_2")]])
    pred_cor <- NA_real_
    
    pred_cols <- c(paste0(ax, "_ptrap_gam_pred_1"), paste0(ax, "_ptrap_gam_pred_2"))
    if (use_pred && all(pred_cols %in% colnames(merged))) {
      pred_cor <- safe_cor(merged[[pred_cols[1]]], merged[[pred_cols[2]]])
    }
    
    tibble(axis = ax, obs_cor = obs_cor, pred_cor = pred_cor, n_overlap_dates = nrow(merged))
  })
}

# -----------------------------
# prediction on PC scores
# -----------------------------
make_lag_pairs <- function(axis_df, axes = c("PC1", "PC2")) {
  stopifnot(all(axes %in% colnames(axis_df)))
  
  axis_df %>%
    arrange(ptrap, date) %>%
    group_by(ptrap) %>%
    mutate(
      next_date = lead(date),
      dt_days = as.numeric(next_date - date),
      across(all_of(axes), lead, .names = "{.col}_next"),
      pair_idx = row_number()
    ) %>%
    ungroup() %>%
    filter(!is.na(next_date))
}

run_next_step_prediction <- function(axis_df,
                                     axes = c("PC1", "PC2"),
                                     train_frac = 0.7,
                                     min_train_pairs = 4) {
  pair_df <- make_lag_pairs(axis_df, axes = axes)
  
  if (nrow(pair_df) == 0) {
    return(list(
      pred_df = pair_df,
      fits = list(),
      metrics_by_ptrap = tibble(),
      metrics_by_treatment = tibble()
    ))
  }
  
  pred_df <- pair_df %>%
    group_by(ptrap) %>%
    mutate(
      n_pairs = n(),
      train_n_raw = pmax(min_train_pairs, floor(train_frac * n_pairs)),
      train_n = if_else(n_pairs >= 2, pmin(train_n_raw, n_pairs - 1L), 0L),
      is_train = pair_idx <= train_n
    ) %>%
    ungroup()
  
  for (ax in axes) {
    pred_df[[paste0(ax, "_pred_naive")]] <- pred_df[[ax]]
  }
  
  train_df <- pred_df %>% filter(is_train)
  test_df  <- pred_df %>% filter(!is_train)
  
  if (nrow(train_df) < length(axes) + 3 || nrow(test_df) == 0) {
    for (ax in axes) {
      pred_df[[paste0(ax, "_pred_dynamic")]] <- NA_real_
    }
    return(list(
      pred_df = pred_df,
      fits = list(),
      metrics_by_ptrap = tibble(),
      metrics_by_treatment = tibble()
    ))
  }
  
  fits <- list()
  rhs <- paste(c(axes, "dt_days", "treatment", "ptrap_f"), collapse = " + ")
  for (ax in axes) {
    form <- as.formula(paste0(ax, "_next ~ ", rhs))
    fits[[ax]] <- lm(form, data = train_df)
    pred_df[[paste0(ax, "_pred_dynamic")]] <- predict(fits[[ax]], newdata = pred_df)
  }
  
  metric_row <- function(df_sub) {
    out <- list(n_test_pairs = nrow(df_sub))
    for (ax in axes) {
      obs <- df_sub[[paste0(ax, "_next")]]
      p_n <- df_sub[[paste0(ax, "_pred_naive")]]
      p_d <- df_sub[[paste0(ax, "_pred_dynamic")]]
      out[[paste0(ax, "_RMSE_naive")]]   <- rmse_vec(obs, p_n)
      out[[paste0(ax, "_RMSE_dynamic")]] <- rmse_vec(obs, p_d)
      out[[paste0(ax, "_MAE_naive")]]    <- mae_vec(obs, p_n)
      out[[paste0(ax, "_MAE_dynamic")]]  <- mae_vec(obs, p_d)
      out[[paste0(ax, "_R2_naive")]]     <- r2_vec(obs, p_n)
      out[[paste0(ax, "_R2_dynamic")]]   <- r2_vec(obs, p_d)
    }
    tibble::as_tibble(out)
  }
  
  metrics_by_ptrap <- pred_df %>%
    filter(!is_train) %>%
    group_by(ptrap, treatment) %>%
    group_modify(~ metric_row(.x)) %>%
    ungroup()
  
  metrics_by_treatment <- pred_df %>%
    filter(!is_train) %>%
    group_by(treatment) %>%
    group_modify(~ metric_row(.x)) %>%
    ungroup()
  
  list(
    pred_df = pred_df,
    fits = fits,
    metrics_by_ptrap = metrics_by_ptrap,
    metrics_by_treatment = metrics_by_treatment
  )
}

# -----------------------------
# stability: leave-one-week-out
# -----------------------------
run_week_stability <- function(clr_tbl, metadata_batch, id_cols,
                               ref_fit, top_k = 20,
                               axes = c("PC1", "PC2"),
                               X_all = NULL, feat_map = NULL) {
  assert_has_columns(metadata_batch, c("Sample.Name", "week"), "metadata_batch")
  
  weeks <- metadata_batch %>%
    filter(!is.na(week)) %>%
    distinct(week) %>%
    arrange(week) %>%
    pull(week)
  
  ref_loadings <- ref_fit$loadings %>% select(feature, any_of(axes), any_of("Taxon"))
  selection_list <- list()
  fold_stats <- list()
  
  for (wk in weeks) {
    keep_meta <- metadata_batch %>% filter(week != wk)
    if (nrow(keep_meta) < 4) next
    
    if (is.null(X_all)) {
      fit_i <- fit_pca_from_clr(
        clr_tbl,
        id_cols = id_cols,
        sample_names = keep_meta$Sample.Name,
        n_pc = max(2, length(axes)),
        feature_id_col = if ("id" %in% id_cols) "id" else NULL
      )
    } else {
      X_sub <- X_all[keep_meta$Sample.Name, , drop = FALSE]
      fit_i <- fit_pca_from_matrix(
        X = X_sub,
        n_pc = max(2, length(axes)),
        feat_map = feat_map
      )
    }
    
    aligned <- align_loading_sign(ref_loadings, fit_i$loadings, axes = axes)
    ld_i <- aligned$loadings
    
    fold_stats[[as.character(wk)]] <- purrr::map_dfr(axes, function(ax) {
      ref_ax <- ref_loadings %>% select(feature, all_of(ax))
      new_ax <- ld_i %>% select(feature, all_of(ax))
      tmp <- inner_join(ref_ax, new_ax, by = "feature", suffix = c("_ref", "_new"))
      
      top_ref <- get_top_features(ref_loadings, axis = ax, top_k = top_k)
      top_new <- get_top_features(ld_i, axis = ax, top_k = top_k)
      
      tibble(
        omitted_week = wk,
        axis = ax,
        loading_cor = safe_cor(tmp[[paste0(ax, "_ref")]], tmp[[paste0(ax, "_new")]]),
        top_overlap = pct_overlap(top_ref, top_new),
        n_keep_samples = nrow(keep_meta)
      )
    })
    
    selection_list[[as.character(wk)]] <- purrr::map_dfr(axes, function(ax) {
      ld_i_local <- ld_i
      if (!"Taxon" %in% colnames(ld_i_local)) {
        ld_i_local$Taxon <- NA_character_
      }
      
      ld_i_local %>%
        mutate(abs_loading = abs(.data[[ax]])) %>%
        arrange(desc(abs_loading)) %>%
        slice_head(n = top_k) %>%
        transmute(
          omitted_week = wk,
          axis = ax,
          feature,
          Taxon,
          loading = .data[[ax]]
        )
    })
  }
  
  fold_stats <- bind_rows(fold_stats)
  selection_tbl <- bind_rows(selection_list)
  
  if (nrow(fold_stats) == 0 || nrow(selection_tbl) == 0) {
    return(list(
      fold_stats = fold_stats,
      selection_tbl = selection_tbl,
      selection_freq = tibble()
    ))
  }
  
  n_folds <- dplyr::n_distinct(fold_stats$omitted_week)
  selection_freq <- selection_tbl %>%
    group_by(axis, feature, Taxon) %>%
    summarise(
      n_selected = n(),
      selection_freq = n() / n_folds,
      mean_loading = mean(loading, na.rm = TRUE),
      sd_loading = sd(loading, na.rm = TRUE),
      sign_consistency = mean(sign(loading) == sign(mean(loading, na.rm = TRUE)), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(axis, desc(selection_freq), desc(abs(mean_loading)))
  
  list(
    fold_stats = fold_stats,
    selection_tbl = selection_tbl,
    selection_freq = selection_freq
  )
}

run_overlap_concordance <- function(asv_mat, pairs_list) {
  assert_has_columns(pairs_list, c("pair_id", "Samp25", "Samp35"), "pairs_list")

  res <- map_dfr(seq_len(nrow(pairs_list)), function(i) {
    samp25 <- pairs_list$Samp25[i]
    samp35 <- pairs_list$Samp35[i]
    pid <- pairs_list$pair_id[i]

    if (!(samp25 %in% colnames(asv_mat) && samp35 %in% colnames(asv_mat))) {
      return(tibble(
        pair_id = pid,
        spearman_rho = NA_real_,
        pearson_r = NA_real_,
        available = FALSE
      ))
    }

    x <- asv_mat[, samp25]
    y <- asv_mat[, samp35]
    tibble(
      pair_id = pid,
      spearman_rho = safe_cor(x, y, method = "spearman"),
      pearson_r = safe_cor(x, y, method = "pearson"),
      available = TRUE
    )
  })

  summary_tbl <- res %>%
    summarise(
      n_pairs = n(),
      n_available = sum(available, na.rm = TRUE),
      mean_spearman = mean(spearman_rho, na.rm = TRUE),
      median_spearman = median(spearman_rho, na.rm = TRUE),
      mean_pearson = mean(pearson_r, na.rm = TRUE),
      median_pearson = median(pearson_r, na.rm = TRUE)
    )

  list(pairwise = res, summary = summary_tbl)
}

# -----------------------------
# top-level runner for one batch
# -----------------------------
run_longitudinal_pca_analysis <- function(clr_tbl,
                           metadata_batch,
                           id_cols,
                           batch_label = "batch",
                           treated_ptraps = c(1, 2),
                           control_ptrap = 3,
                           n_pc = 5,
                           axes = NULL,
                           top_k = 20,
                           train_frac = 0.7,
                           min_train_pairs = 4,
                           k_time = 6) {
  
  pca_cache <- prepare_pca_cache(
    clr_tbl = clr_tbl,
    id_cols = id_cols,
    sample_names = metadata_batch$Sample.Name,
    feature_id_col = if ("id" %in% id_cols) "id" else NULL
  )
  
  built <- build_axis_df(
    clr_tbl = clr_tbl,
    metadata_batch = metadata_batch,
    id_cols = id_cols,
    n_pc = n_pc,
    X_all = pca_cache$X_all,
    feat_map = pca_cache$feat_map
  )
  
  axis_df <- built$axis_df %>%
    add_treatment_labels(treated_ptraps = treated_ptraps, control_ptrap = control_ptrap)
  
  analysis_axes <- resolve_analysis_axes(axis_df = axis_df, axes = axes, n_pc = n_pc, max_axes = 2)
  
  treat_gams <- fit_treatment_gams(axis_df, axes = analysis_axes, k_time = k_time)
  ptrap_gams <- fit_ptrap_gams(axis_df, axes = analysis_axes, k_time = k_time)
  
  axis_df <- axis_df %>%
    augment_gam_preds(treat_gams, suffix = "treat_gam") %>%
    augment_gam_preds(ptrap_gams, suffix = "ptrap_gam")
  
  pred_obj <- run_next_step_prediction(
    axis_df,
    axes = analysis_axes,
    train_frac = train_frac,
    min_train_pairs = min_train_pairs
  )
  
  stab_obj <- run_week_stability(
    clr_tbl = clr_tbl,
    metadata_batch = metadata_batch,
    id_cols = id_cols,
    ref_fit = built$fit,
    top_k = top_k,
    axes = analysis_axes,
    X_all = pca_cache$X_all,
    feat_map = pca_cache$feat_map
  )
  
  treated_conc <- calc_treated_concordance(
    axis_df,
    treated_ptraps = treated_ptraps,
    axes = analysis_axes,
    use_pred = TRUE
  )
  
  list(
    batch_label = batch_label,
    analysis_axes = analysis_axes,
    axis_df = axis_df,
    pca_fit = built$fit,
    axis_shift = summarize_axis_shift(axis_df, axes = analysis_axes),
    treatment_gams = treat_gams,
    ptrap_gams = ptrap_gams,
    treatment_gam_summary = summarize_gams(treat_gams),
    ptrap_gam_summary = summarize_gams(ptrap_gams),
    prediction = pred_obj,
    stability = stab_obj,
    treated_concordance = treated_conc
  )
}
