# ================================
# Radar Plot Functions
# ================================

#' Create Radar Plot with Statistical Tests
#'
#' Creates radar (spider) charts with optional statistical testing and p-values.
#' Supports two modes: groups mode (with statistics) and items mode (without statistics).
#'
#' @param data Data frame containing the data to plot
#' @param value_vars Character vector of column names to plot on radar axes
#' @param group_var Character. Grouping variable name (only for mode="groups")
#' @param mode Character. One of "auto", "groups", or "items". Auto-detects based on group_var
#' @param items_layout Character. For items mode: "grid" (one radar per item),
#'   "single" (one chosen item), or "overlay_all" (all items on one radar)
#' @param item_name Character. For items_layout="single", specify which item to plot
#' @param show_average Logical. Show average line in items mode
#' @param engine Character. For groups mode: "ggplot" or "fmsb"
#' @param p_method Character. Statistical test: "auto", "t.test", "anova", or "none"
#' @param max_scale Numeric. Maximum value for radar axis scale
#' @param colors Character vector of colors for groups/items
#' @param title Character. Plot title
#' @param text_radius Numeric. For fmsb groups mode: p-value text radius
#' @param mfrow Numeric vector. Items grid layout, e.g., c(3,4)
#' @param facet_ncol Integer. Facet columns (unused, placeholder)
#' @param show_ci Logical. Show confidence intervals for groups mode (default: TRUE)
#' @param ci_level Numeric. Confidence level for intervals (default: 0.95)
#' @param ci_alpha Numeric. Transparency for CI ribbons in ggplot (default: 0.15)
#'
#' @return Invisibly returns a list containing plot object (if ggplot), p-values,
#'   means, confidence intervals, engine type, mode, test type, and max_scale
#'
#' @export
#' @examples
#' \dontrun{
#' # Groups mode with t-test and confidence intervals
#' data <- data.frame(
#'   group = rep(c("Control", "Treatment"), each = 20),
#'   var1 = rnorm(40, 5, 2),
#'   var2 = rnorm(40, 7, 2),
#'   var3 = rnorm(40, 6, 2)
#' )
#' 
#' # With 95% confidence intervals (default)
#' result <- radar_with_p(
#'   data = data,
#'   value_vars = c("var1", "var2", "var3"),
#'   group_var = "group",
#'   mode = "groups",
#'   engine = "ggplot",
#'   p_method = "auto",
#'   show_ci = TRUE
#' )
#' 
#' # With custom 99% confidence intervals
#' result <- radar_with_p(
#'   data = data,
#'   value_vars = c("var1", "var2", "var3"),
#'   group_var = "group",
#'   mode = "groups",
#'   engine = "fmsb",
#'   show_ci = TRUE,
#'   ci_level = 0.99
#' )
#' 
#' # Without confidence intervals
#' result <- radar_with_p(
#'   data = data,
#'   value_vars = c("var1", "var2", "var3"),
#'   group_var = "group",
#'   mode = "groups",
#'   show_ci = FALSE
#' )
#' 
#' # Items mode - grid layout
#' items_data <- data.frame(
#'   item1 = c(5, 7, 6, 8, 5),
#'   item2 = c(7, 6, 8, 7, 6),
#'   item3 = c(6, 8, 7, 9, 7)
#' )
#' rownames(items_data) <- paste("Item", 1:5)
#' 
#' radar_with_p(
#'   data = items_data,
#'   value_vars = colnames(items_data),
#'   mode = "items",
#'   items_layout = "grid",
#'   show_average = TRUE
#' )
#' }
radar_with_p <- function(
    data,
    value_vars,
    group_var = NULL,
    mode = c("auto", "groups", "items"),
    items_layout = c("grid", "single", "overlay_all"),
    item_name = NULL,
    show_average = TRUE,
    engine = c("ggplot", "fmsb"),
    p_method = c("auto", "t.test", "anova", "none"),
    max_scale = NULL,
    colors = NULL,
    title = NULL,
    text_radius = 1.35,
    mfrow = NULL,
    facet_ncol = 4,
    show_ci = TRUE,
    ci_level = 0.95,
    ci_alpha = 0.15
){
  
  # Check required packages
  required_pkgs <- c("fmsb")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required. Install with: install.packages('", pkg, "')")
    }
  }
  
  # Match arguments
  mode <- match.arg(mode)
  engine <- match.arg(engine)
  p_method <- match.arg(p_method)
  items_layout <- match.arg(items_layout)
  
  # Helper function
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  # Validate inputs
  stopifnot(all(value_vars %in% names(data)))
  
  # -------------------------
  # Decide mode automatically
  # -------------------------
  if (mode == "auto") {
    mode <- if (!is.null(group_var) && group_var %in% names(data)) "groups" else "items"
  }
  
  # -------------------------
  # ITEMS MODE (no p-values)
  # -------------------------
  if (mode == "items") {
    mat <- data[, value_vars, drop = FALSE]
    if (!is.null(rownames(mat)) && any(duplicated(rownames(mat)))) {
      warning("Duplicated rownames detected; consider making them unique for clearer legends.")
    }
    
    # Bounds + average
    colMax <- function(x) apply(x, 2, max, na.rm = TRUE)
    colMin <- function(x) apply(x, 2, min, na.rm = TRUE)
    max_row <- colMax(mat)
    min_row <- colMin(mat)
    avg_row <- colMeans(mat, na.rm = TRUE)
    
    # Axis scale
    if (is.null(max_scale)) {
      mx <- suppressWarnings(max(as.matrix(mat), na.rm = TRUE))
      if (!is.finite(mx)) mx <- 1
      max_scale <- max(1, ceiling(mx / 5) * 5)
    }
    grid_breaks <- pretty(c(0, max_scale), n = 5)
    seg_count <- length(grid_breaks) - 1
    
    # Colors
    n_items <- nrow(mat)
    if (items_layout == "overlay_all") {
      if (is.null(colors) || length(colors) < n_items) {
        colors <- grDevices::hcl.colors(n_items, palette = "Dark 3")
      }
      line_cols <- colors
      fill_cols <- rep(NA, n_items)
      avg_cols <- c("grey30")
      avg_fills <- c(grDevices::adjustcolor("grey60", 0.25))
    } else {
      # grid/single: default one item color
      if (is.null(colors)) colors <- "steelblue4"
    }
    
    # Draw
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar), add = TRUE)
    
    if (items_layout == "grid") {
      # Sensible default grid
      if (is.null(mfrow)) {
        n <- n_items
        ncol <- if (n <= 12) 4 else 4
        nrow <- ceiling(n / ncol)
        mfrow <- c(nrow, ncol)
      }
      par(mfrow = mfrow, mar = rep(0.8, 4))
      
      for (i in seq_len(n_items)) {
        poly <- if (isTRUE(show_average)) {
          rbind(max_row, min_row, avg_row, mat[i, , drop = FALSE])
        } else {
          rbind(max_row, min_row, mat[i, , drop = FALSE])
        }
        
        fmsb::radarchart(
          poly,
          seg = seg_count,
          caxislabels = grid_breaks,
          axistype = 1,
          cglcol = "grey85", cglty = 1, cglwd = 0.7,
          axislabcol = "grey30", vlcex = 0.85,
          pcol = if (isTRUE(show_average)) c("grey40", colors) else colors,
          pfcol = if (isTRUE(show_average)) c(grDevices::adjustcolor("grey60", 0.25), NA) else NA,
          plwd = if (isTRUE(show_average)) c(2, 2) else 2,
          pty = 32,
          title = rownames(mat)[i] %||% paste("Item", i)
        )
        
        if (isTRUE(show_average)) {
          legend("bottomright", bty = "n",
                 legend = c("Average", rownames(mat)[i] %||% paste("Item", i)),
                 lwd = c(2, 2), col = c("grey40", colors), cex = 0.8)
        }
      }
      
    } else if (items_layout == "single") {
      idx <- if (!is.null(item_name)) match(item_name, rownames(mat)) else 1L
      if (is.na(idx)) stop(sprintf("item_name '%s' not found in rownames(data).", item_name))
      
      par(mfrow = c(1, 1), mar = c(1.5, 1.5, 2, 1.5))
      
      poly <- if (isTRUE(show_average)) {
        rbind(max_row, min_row, avg_row, mat[idx, , drop = FALSE])
      } else {
        rbind(max_row, min_row, mat[idx, , drop = FALSE])
      }
      
      fmsb::radarchart(
        poly,
        seg = seg_count,
        caxislabels = grid_breaks,
        axistype = 1,
        cglcol = "grey80", cglty = 1, cglwd = 0.8,
        axislabcol = "grey30", vlcex = 0.95,
        pcol = if (isTRUE(show_average)) c("grey40", colors) else colors,
        pfcol = if (isTRUE(show_average)) c(grDevices::adjustcolor("grey60", 0.25), NA) else NA,
        plwd = if (isTRUE(show_average)) c(2, 2) else 2,
        pty = 32,
        title = sprintf("%s%s",
                        rownames(mat)[idx] %||% paste("Item", idx),
                        if (isTRUE(show_average)) " vs Average" else "")
      )
      
      if (isTRUE(show_average)) {
        legend("topright", bty = "n",
               legend = c("Average", rownames(mat)[idx] %||% paste("Item", idx)),
               lwd = c(2, 2), col = c("grey40", colors), cex = 0.9)
      }
      
    } else { # overlay_all
      par(mfrow = c(1, 1), mar = c(1.5, 1.5, 2, 1.5))
      
      poly <- if (isTRUE(show_average)) {
        rbind(rep(max_scale, ncol(mat)), rep(0, ncol(mat)), avg_row, mat)
      } else {
        rbind(rep(max_scale, ncol(mat)), rep(0, ncol(mat)), mat)
      }
      
      pcol <- if (isTRUE(show_average)) c("grey30", line_cols) else line_cols
      pfcol <- if (isTRUE(show_average)) c(grDevices::adjustcolor("grey60", 0.25), fill_cols) else fill_cols
      plwd <- if (isTRUE(show_average)) c(3, rep(2, n_items)) else rep(2, n_items)
      
      fmsb::radarchart(
        poly,
        seg = seg_count,
        caxislabels = grid_breaks,
        axistype = 1,
        cglcol = "grey80", cglty = 1, cglwd = 0.8,
        axislabcol = "grey30", vlcex = 0.9,
        pcol = pcol,
        pfcol = pfcol,
        plwd = plwd,
        title = title %||% if (isTRUE(show_average)) "All items vs average" else "All items"
      )
      
      legend("topright", bty = "n",
             legend = c(if (isTRUE(show_average)) "Average", rownames(mat)),
             lwd = c(if (isTRUE(show_average)) 3, rep(2, n_items)),
             col = c(if (isTRUE(show_average)) "grey30", line_cols),
             cex = 0.9)
    }
    
    return(invisible(list(
      engine = "fmsb",
      mode = "items",
      items_layout = items_layout,
      show_average = show_average,
      max_scale = max_scale
    )))
  }
  
  # -------------------------
  # GROUPS MODE (stats + p-values)
  # -------------------------
  
  # Check additional packages for groups mode
  if (engine == "ggplot") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package 'ggplot2' is required for engine='ggplot'. Install with: install.packages('ggplot2')")
    }
    if (!requireNamespace("tidyr", quietly = TRUE)) {
      stop("Package 'tidyr' is required. Install with: install.packages('tidyr')")
    }
  }
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Install with: install.packages('dplyr')")
  }
  
  if (p_method != "none") {
    if (!requireNamespace("rstatix", quietly = TRUE)) {
      stop("Package 'rstatix' is required for statistical tests. Install with: install.packages('rstatix')")
    }
  }
  
  stopifnot(!is.null(group_var), group_var %in% names(data))
  
  # Prepare data
  df <- data %>%
    dplyr::select(dplyr::all_of(c(group_var, value_vars))) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(value_vars), as.numeric)) %>%
    tidyr::drop_na(dplyr::all_of(c(group_var, value_vars)))
  
  grp_fac <- factor(df[[group_var]])
  n_groups <- nlevels(grp_fac)
  
  # Determine test mode
  counts <- table(grp_fac)
  has_reps <- all(counts >= 2)
  test_mode <- switch(p_method,
                      "t.test" = if (n_groups == 2 && has_reps) "t.test" else "none",
                      "anova" = if (n_groups >= 2 && has_reps) "anova" else "none",
                      "auto" = if (n_groups == 2 && has_reps) "t.test" else if (n_groups >= 3 && has_reps) "anova" else "none",
                      "none" = "none")
  
  # Format p-values
  fmt_p <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 0.001) "<0.001***"
    else if (p < 0.01) sprintf("%.3f**", p)
    else if (p < 0.05) sprintf("%.3f*", p)
    else "ns"
  }
  
  # Calculate p-values per variable
  p_tbl <- purrr::map_dfr(value_vars, function(v) {
    p <- NA_real_
    if (test_mode == "t.test" && n_groups == 2) {
      p <- tryCatch(t.test(df[[v]] ~ grp_fac)$p.value, error = function(e) NA_real_)
    } else if (test_mode == "anova") {
      aov_res <- tryCatch(
        rstatix::anova_test(data = df, formula = as.formula(paste(v, "~", group_var))),
        error = function(e) NULL
      )
      if (!is.null(aov_res)) {
        p <- tryCatch(aov_res$p[aov_res$Effect %in% c(group_var, "group")][1], error = function(e) NA_real_)
        if (is.na(p)) p <- tryCatch(aov_res$p[1], error = function(e) NA_real_)
      }
    }
    dplyr::tibble(variable = v, p = p, p_label = fmt_p(p))
  }) %>%
    dplyr::mutate(axis_with_p = if (test_mode == "none") variable else paste0(variable, "\n(", p_label, ")"))
  
  # Group means
  means_wide <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_var))) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(value_vars), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  
  # Calculate confidence intervals
  ci_data <- NULL
  if (show_ci) {
    z_score <- qnorm((1 + ci_level) / 2)
    
    ci_data <- df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_var))) %>%
      dplyr::summarise(
        dplyr::across(
          dplyr::all_of(value_vars),
          list(
            mean = ~mean(.x, na.rm = TRUE),
            se = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))),
            n = ~sum(!is.na(.x))
          ),
          .names = "{.col}_{.fn}"
        ),
        .groups = "drop"
      )
    
    # Reshape CI data for easier joining
    ci_long <- ci_data %>%
      tidyr::pivot_longer(
        cols = -dplyr::all_of(group_var),
        names_to = c("variable", ".value"),
        names_pattern = "(.+)_(mean|se|n)"
      ) %>%
      dplyr::mutate(
        ci_lower = mean - z_score * se,
        ci_upper = mean + z_score * se,
        ci_lower = pmax(ci_lower, 0)  # Don't go below 0
      ) %>%
      dplyr::select(dplyr::all_of(group_var), variable, ci_lower, ci_upper, n)
  }
  
  # Axis scale (consider CI upper bounds if showing CIs)
  if (is.null(max_scale)) {
    if (show_ci && !is.null(ci_data)) {
      mx <- ci_long %>%
        dplyr::pull(ci_upper) %>%
        max(na.rm = TRUE)
    } else {
      mx <- means_wide %>%
        dplyr::select(dplyr::all_of(value_vars)) %>%
        as.matrix() %>%
        max(na.rm = TRUE)
    }
    if (!is.finite(mx)) mx <- 1
    max_scale <- max(1, ceiling(mx / 5) * 5)
  }
  
  # Colors
  if (is.null(colors)) {
    cols <- grDevices::rainbow(nlevels(grp_fac), s = 0.65, v = 0.9)
    if (engine == "ggplot") names(cols) <- levels(grp_fac)
    colors <- cols
  }
  
  if (engine == "ggplot") {
    plot_long <- means_wide %>%
      tidyr::pivot_longer(cols = dplyr::all_of(value_vars), names_to = "variable", values_to = "value") %>%
      dplyr::left_join(dplyr::select(p_tbl, variable, axis_with_p), by = "variable") %>%
      dplyr::mutate(
        axis_with_p = factor(axis_with_p, levels = unique(p_tbl$axis_with_p)),
        !!group_var := factor(.data[[group_var]], levels = levels(grp_fac))
      )
    
    # Add CI data if available
    if (show_ci && !is.null(ci_long)) {
      plot_long <- plot_long %>%
        dplyr::left_join(ci_long, by = c(group_var, "variable"))
    }
    
    g <- ggplot2::ggplot(
      plot_long,
      ggplot2::aes(x = axis_with_p, y = value, group = .data[[group_var]],
                   color = .data[[group_var]], fill = .data[[group_var]])
    )
    
    # Add CI ribbons if available
    if (show_ci && !is.null(ci_long)) {
      g <- g + 
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
          alpha = ci_alpha,
          color = NA
        )
    }
    
    g <- g +
      ggplot2::geom_polygon(alpha = 0.18, linewidth = 1) +
      ggplot2::geom_point(size = 2.4, show.legend = FALSE) +
      ggplot2::coord_polar(start = 0) +
      ggplot2::scale_y_continuous(
        limits = c(0, max_scale),
        breaks = seq(0, max_scale, by = max(5, round(max_scale / 6)))
      ) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::labs(
        title = title %||% paste0("Radar", 
                                  if (test_mode != "none") " with p-values" else "",
                                  if (show_ci) paste0(" (", ci_level*100, "% CI)") else ""),
        x = NULL, y = NULL, color = "Group", fill = "Group",
        caption = if (test_mode != "none") "* p<0.05  ** p<0.01  *** p<0.001" else NULL
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        legend.position = "right",
        panel.grid.major = ggplot2::element_line(color = "grey85"),
        axis.text.x = ggplot2::element_text(size = 11, lineheight = 0.9),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16)
      )
    
    print(g)
    return(invisible(list(
      plot = g,
      p_values = p_tbl,
      means = means_wide,
      ci = if (show_ci) ci_long else NULL,
      ci_level = if (show_ci) ci_level else NULL,
      engine = "ggplot",
      mode = "groups",
      test = test_mode,
      max_scale = max_scale
    )))
  }
  
  # fmsb engine (groups)
  mat <- means_wide
  rownames(mat) <- mat[[group_var]]
  mat <- mat[, value_vars, drop = FALSE]
  
  max_min <- rbind(
    rep(max_scale, length(value_vars)),
    rep(0, length(value_vars))
  )
  colnames(max_min) <- colnames(mat)
  plot_mat <- rbind(max_min, mat)
  
  base_cols <- if (!is.null(names(colors))) unname(colors) else colors
  line_cols <- base_cols
  fill_cols <- grDevices::adjustcolor(base_cols, alpha.f = 0.35)
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mar = c(1, 1, 2, 1))
  
  fmsb::radarchart(
    plot_mat,
    axistype = 1,
    pcol = line_cols,
    pfcol = fill_cols,
    plwd = 2,
    cglcol = "grey60", cglty = 1, cglwd = 0.8,
    axislabcol = "grey30",
    vlcex = 1.05
  )
  
  # Add confidence interval bands if requested
  if (show_ci && !is.null(ci_long)) {
    k <- length(value_vars)
    angles <- seq(from = pi / 2, to = -3 * pi / 2, length.out = k + 1)[-(k + 1)]
    
    # Draw CI bands for each group
    for (grp_idx in seq_len(nrow(mat))) {
      grp_name <- rownames(mat)[grp_idx]
      ci_grp <- ci_long[ci_long[[group_var]] == grp_name, ]
      
      if (nrow(ci_grp) > 0) {
        # Match order to value_vars
        ci_grp <- ci_grp[match(value_vars, ci_grp$variable), ]
        
        # Convert to polar coordinates for lower bound
        x_lower <- cos(angles) * ci_grp$ci_lower / max_scale
        y_lower <- sin(angles) * ci_grp$ci_lower / max_scale
        
        # Convert to polar coordinates for upper bound
        x_upper <- cos(angles) * ci_grp$ci_upper / max_scale
        y_upper <- sin(angles) * ci_grp$ci_upper / max_scale
        
        # Draw CI as shaded polygon between lower and upper bounds
        x_poly <- c(x_lower, rev(x_upper))
        y_poly <- c(y_lower, rev(y_upper))
        
        polygon(x_poly, y_poly, 
                col = grDevices::adjustcolor(base_cols[grp_idx], alpha.f = 0.15),
                border = NA)
      }
    }
    
    # Redraw the main lines on top of CI bands
    for (grp_idx in seq_len(nrow(mat))) {
      vals <- as.numeric(mat[grp_idx, ])
      x_vals <- cos(angles) * vals / max_scale
      y_vals <- sin(angles) * vals / max_scale
      
      # Close the polygon
      x_vals <- c(x_vals, x_vals[1])
      y_vals <- c(y_vals, y_vals[1])
      
      lines(x_vals, y_vals, col = line_cols[grp_idx], lwd = 2)
    }
  }
  
  # P-values around the chart
  if (test_mode != "none") {
    k <- length(value_vars)
    angles <- seq(from = pi / 2, to = -3 * pi / 2, length.out = k + 1)[-(k + 1)]
    for (i in seq_len(k)) {
      text(cos(angles[i]) * text_radius,
           sin(angles[i]) * text_radius,
           labels = p_tbl$p_label[i],
           cex = 0.8, col = "darkred", xpd = NA)
    }
  }
  
  legend("topright",
         legend = rownames(plot_mat)[-(1:2)],
         bty = "n", pch = 20, col = fill_cols,
         text.col = "black", cex = 1.0, pt.cex = 2)
  
  title(main = title %||% paste0("Radar", 
                                 if (test_mode != "none") " with p-values" else "",
                                 if (show_ci) paste0(" (", ci_level*100, "% CI)") else ""))
  
  invisible(list(
    plot = NULL,
    p_values = p_tbl,
    means = means_wide,
    ci = if (show_ci) ci_long else NULL,
    ci_level = if (show_ci) ci_level else NULL,
    engine = "fmsb",
    mode = "groups",
    test = test_mode,
    max_scale = max_scale
  ))
}
