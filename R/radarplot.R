# ================================
# Packages
# ================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rstatix)
  library(fmsb)
  library(purrr)
  library(grDevices)  # adjustcolor, palettes
})

`%||%` <- function(a, b) if (!is.null(a)) a else b  # small helper

# ================================
# One function to rule them all
#   - mode = "groups"  : multi-row per group (stats + p-values; engines: ggplot/fmsb)
#   - mode = "items"   : rows are items (no stats), layouts:
#         items_layout = "grid"        (one radar per item)
#                        "single"      (one chosen item)
#                        "overlay_all" (all items on one radar)
# ================================
radar_with_p <- function(
    data,
    value_vars,
    group_var = NULL,                # only for mode="groups"
    mode = c("auto", "groups", "items"),
    items_layout = c("grid", "single", "overlay_all"),
    item_name = NULL,                # for items_layout="single"
    show_average = TRUE,             # for items layouts
    engine = c("ggplot", "fmsb"),    # groups-mode engine; items uses fmsb
    p_method = c("auto", "t.test", "anova", "none"),
    max_scale = NULL,
    colors = NULL,
    title = NULL,
    text_radius = 1.35,              # fmsb groups-mode: p-value text radius
    mfrow = NULL,                    # items grid layout, e.g., c(3,4)
    facet_ncol = 4                   # (unused here; placeholder if you add ggplot facets)
){
  mode     <- match.arg(mode)
  engine   <- match.arg(engine)
  p_method <- match.arg(p_method)
  items_layout <- match.arg(items_layout)

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
    colMax <- function(x) apply(x, 2, max)
    colMin <- function(x) apply(x, 2, min)
    max_row <- colMax(mat)
    min_row <- colMin(mat)
    avg_row <- colMeans(mat)

    # axis scale
    if (is.null(max_scale)) {
      mx <- suppressWarnings(max(as.matrix(mat), na.rm = TRUE))
      if (!is.finite(mx)) mx <- 1
      max_scale <- max(1, ceiling(mx / 5) * 5)
    }
    grid_breaks <- pretty(c(0, max_scale), n = 5)
    seg_count   <- length(grid_breaks) - 1

    # Colors
    n_items <- nrow(mat)
    if (items_layout == "overlay_all") {
      if (is.null(colors) || length(colors) < n_items) {
        colors <- grDevices::hcl.colors(n_items, palette = "Dark 3")
      }
      line_cols <- colors
      fill_cols <- rep(NA, n_items)
      avg_cols  <- c("grey30")
      avg_fills <- c(grDevices::adjustcolor("grey60", 0.25))
    } else {
      # grid/single: default one item color
      if (is.null(colors)) colors <- "steelblue4"
    }

    # Draw
    opar <- par(no.readonly = TRUE); on.exit(par(opar), add = TRUE)

    if (items_layout == "grid") {
      # sensible default grid like 3x4 if <=12 items, else 4 cols
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
          cglcol = "grey85", cglty = 1, cglwd = 0.7,  # polygon grid
          axislabcol = "grey30", vlcex = 0.85,
          pcol  = if (isTRUE(show_average)) c("grey40", colors) else colors,
          pfcol = if (isTRUE(show_average)) c(grDevices::adjustcolor("grey60", 0.25), NA) else NA,
          plwd  = if (isTRUE(show_average)) c(2, 2) else 2,
          pty   = 32,
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
      par(mfrow = c(1,1), mar = c(1.5, 1.5, 2, 1.5))
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
        pcol  = if (isTRUE(show_average)) c("grey40", colors) else colors,
        pfcol = if (isTRUE(show_average)) c(grDevices::adjustcolor("grey60", 0.25), NA) else NA,
        plwd  = if (isTRUE(show_average)) c(2, 2) else 2,
        pty   = 32,
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
      par(mfrow = c(1,1), mar = c(1.5, 1.5, 2, 1.5))
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
        pcol  = pcol,
        pfcol = pfcol,
        plwd  = plwd,
        title = title %||% if (isTRUE(show_average)) "All items vs average" else "All items"
      )
      legend("topright", bty = "n",
             legend = c(if (isTRUE(show_average)) "Average", rownames(mat)),
             lwd = c(if (isTRUE(show_average)) 3, rep(2, n_items)),
             col = c(if (isTRUE(show_average)) "grey30", line_cols),
             cex = 0.9)
    }

    return(invisible(list(
      engine = "fmsb", mode = "items", items_layout = items_layout,
      show_average = show_average, max_scale = max_scale
    )))
  }

  # -------------------------
  # GROUPS MODE (stats + p-values)
  # -------------------------
  stopifnot(!is.null(group_var), group_var %in% names(data))

  df <- data %>%
    select(all_of(c(group_var, value_vars))) %>%
    mutate(across(all_of(value_vars), as.numeric)) %>%
    drop_na(all_of(c(group_var, value_vars)))

  grp_fac  <- factor(df[[group_var]])
  n_groups <- nlevels(grp_fac)

  # Determine test mode (fallback to "none" if not enough replicates)
  counts <- table(grp_fac)
  has_reps <- all(counts >= 2)
  test_mode <- switch(p_method,
                      "t.test" = if (n_groups == 2 && has_reps) "t.test" else "none",
                      "anova"  = if (n_groups >= 2 && has_reps) "anova"  else "none",
                      "auto"   = if (n_groups == 2 && has_reps) "t.test" else if (n_groups >= 3 && has_reps) "anova" else "none",
                      "none"   = "none")

  fmt_p <- function(p){
    if (is.na(p)) return("NA")
    if (p < 0.001) "<0.001***"
    else if (p < 0.01) sprintf("%.3f**", p)
    else if (p < 0.05) sprintf("%.3f*",  p)
    else "ns"
  }

  # p-values per variable
  p_tbl <- map_dfr(value_vars, function(v){
    p <- NA_real_
    if (test_mode == "t.test" && n_groups == 2) {
      p <- tryCatch(t.test(df[[v]] ~ grp_fac)$p.value, error = function(e) NA_real_)
    } else if (test_mode == "anova") {
      aov_res <- tryCatch(rstatix::anova_test(data = df, formula = as.formula(paste(v, "~", group_var))), error = function(e) NULL)
      if (!is.null(aov_res)) {
        p <- tryCatch(aov_res$p[aov_res$Effect %in% c(group_var, "group")][1], error = function(e) NA_real_)
        if (is.na(p)) p <- tryCatch(aov_res$p[1], error = function(e) NA_real_)
      }
    }
    tibble(variable = v, p = p, p_label = fmt_p(p))
  }) %>%
    mutate(axis_with_p = if (test_mode == "none") variable else paste0(variable, "\n(", p_label, ")"))

  # group means (plotted)
  means_wide <- df %>%
    group_by(.data[[group_var]]) %>%
    summarise(across(all_of(value_vars), ~mean(.x, na.rm = TRUE)), .groups = "drop")

  # axis scale
  if (is.null(max_scale)) {
    mx <- means_wide %>% select(all_of(value_vars)) %>% as.matrix() %>% max(na.rm = TRUE)
    if (!is.finite(mx)) mx <- 1
    max_scale <- max(1, ceiling(mx / 5) * 5)
  }

  # colors
  if (is.null(colors)) {
    cols <- grDevices::rainbow(nlevels(grp_fac), s = 0.65, v = 0.9)
    if (engine == "ggplot") names(cols) <- levels(grp_fac)
    colors <- cols
  }

  if (engine == "ggplot") {
    plot_long <- means_wide %>%
      pivot_longer(cols = all_of(value_vars), names_to = "variable", values_to = "value") %>%
      left_join(select(p_tbl, variable, axis_with_p), by = "variable") %>%
      mutate(
        axis_with_p = factor(axis_with_p, levels = unique(p_tbl$axis_with_p)),
        !!group_var := factor(.data[[group_var]], levels = levels(grp_fac))
      )

    g <- ggplot(
      plot_long,
      aes(x = axis_with_p, y = value, group = .data[[group_var]],
          color = .data[[group_var]], fill = .data[[group_var]])
    ) +
      geom_polygon(alpha = 0.18, linewidth = 1) +
      geom_point(size = 2.4, show.legend = FALSE) +
      coord_polar(start = 0) +
      scale_y_continuous(
        limits = c(0, max_scale),
        breaks = seq(0, max_scale, by = max(5, round(max_scale/6)))
      ) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      labs(
        title = title %||% paste0("Radar", if (test_mode != "none") " with p-values" else ""),
        x = NULL, y = NULL, color = "Group", fill = "Group",
        caption = if (test_mode != "none") "* p<0.05  ** p<0.01  *** p<0.001" else NULL
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "right",
        panel.grid.major = element_line(color = "grey85"),
        axis.text.x = element_text(size = 11, lineheight = .9),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
      )

    print(g)
    return(invisible(list(plot = g, p_values = p_tbl, means = means_wide,
                          engine = "ggplot", mode = "groups",
                          test = test_mode, max_scale = max_scale)))
  }

  # fmsb engine (groups)
  mat <- means_wide
  rownames(mat) <- mat[[group_var]]
  mat <- mat[, value_vars, drop = FALSE]

  max_min <- rbind(rep(max_scale, length(value_vars)),
                   rep(0,         length(value_vars)))
  colnames(max_min) <- colnames(mat)
  plot_mat <- rbind(max_min, mat)

  base_cols <- if (!is.null(names(colors))) unname(colors) else colors
  line_cols <- base_cols
  fill_cols <- grDevices::adjustcolor(base_cols, alpha.f = 0.35)

  old_par <- par(no.readonly = TRUE); on.exit(par(old_par), add = TRUE)
  par(mar = c(1,1,2,1))

  fmsb::radarchart(
    plot_mat,
    axistype = 1,
    pcol  = line_cols,
    pfcol = fill_cols,
    plwd  = 2,
    cglcol = "grey60", cglty = 1, cglwd = 0.8,   # polygon grid
    axislabcol = "grey30",
    vlcex = 1.05
  )

  # p-values around the chart (if applicable)
  if (test_mode != "none") {
    k <- length(value_vars)
    angles <- seq(from = pi/2, to = -3*pi/2, length.out = k + 1)[-(k + 1)]
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

  title(main = title %||% paste0("Radar", if (test_mode != "none") " with p-values" else ""))

  invisible(list(plot = NULL, p_values = p_tbl, means = means_wide,
                 engine = "fmsb", mode = "groups",
                 test = test_mode, max_scale = max_scale))
}

# ================================
# Example A — Two groups → fmsb (t-tests)
# ================================
set.seed(123)
exA <- data.frame(
  group = rep(c("Control", "Treatment"), each = 50),
  Clarity    = c(rnorm(50, 6, 1.5),  rnorm(50, 7.5, 1.5)),
  Speed      = c(rnorm(50, 7, 1.5),  rnorm(50, 5,   1.5)),
  Accuracy   = c(rnorm(50, 8, 1.5),  rnorm(50, 8.5, 1.5)),
  Innovation = c(rnorm(50, 5, 1.5),  rnorm(50, 7,   1.5)),
  Strategy   = c(rnorm(50, 6.5,1.5), rnorm(50, 6.8, 1.5))
)

resA <- radar_with_p(
  data = exA,
  value_vars = names(exA)[-1],
  group_var = "group",
  mode = "groups",
  engine = "fmsb",
  p_method = "auto",   # auto chooses t.test for 2 groups (if replicates exist)
  title = "Control vs Treatment (fmsb)"
)
# resA$p_values holds the p-value table

# ================================
# Example B — ≥3 groups → ggplot (ANOVA)
# ================================
set.seed(42)
npg <- 20
groups <- c("Healthy", "Mild", "Moderate", "Severe")
gene_data <- data.frame(
  patient_id = paste0("P", 1:(npg * 4)),
  severity   = factor(rep(groups, each = npg), levels = groups, ordered = TRUE)
) %>%
  mutate(
    sev_num = as.numeric(severity),
    Gene_A = rnorm(n(), 10 + (sev_num * 1.5), 2),
    Gene_B = rnorm(n(),  5 + (sev_num * 2.0), 2.2),
    Gene_C = rnorm(n(), 25 - (sev_num * 2.0), 2.8),
    Gene_D = rnorm(n(), 15 + c(5, 0, 1, 6)[sev_num], 2),
    Gene_E = rnorm(n(), 12, 2.6),
    Gene_F = rnorm(n(), 20, 3)
  ) %>% select(-sev_num)

sev_cols <- c(Healthy="#4CAF50", Mild="#FFC107", Moderate="#FF9800", Severe="#F44336")

resB <- radar_with_p(
  data = gene_data,
  value_vars = c("Gene_A","Gene_B","Gene_C","Gene_D","Gene_E","Gene_F"),
  group_var = "severity",
  mode = "groups",
  engine = "ggplot",
  p_method = "auto",             # auto chooses ANOVA for ≥3 groups (if replicates exist)
  max_scale = 30,                # optional
  colors = sev_cols,
  title  = "Gene Expression by Severity (ggplot)"
)
# print(resB$plot)  # if needed

# ================================
# Countries matrix (items mode demos C–F)
# ================================
countries <- read.csv(text="
country,healthcare,education,innovation,environment,safety,infrastructure
Switzerland,19,18,20,19,18,19
Germany,17,16,18,16,15,17
Czechia,12,14,13,12,,13
Poland,11,13,10,9,12,
Spain,15,16,,16,15,14
Italy,14,15,12,13,13,12
Netherlands,18,17,19,,16,16
Sweden,19,18,,19,18,
France,16,16,16,15,14,15
Austria,17,16,17,18,17,17
Portugal,14,15,12,16,16,13
", stringsAsFactors = FALSE)
countries[is.na(countries)] <- 0
rownames(countries) <- countries$country
countries$country <- NULL

# ================================
# Example C — One radar per country (grid, no average)
# ================================
radar_with_p(
  data = countries,
  value_vars = colnames(countries),
  mode = "items",
  items_layout = "grid",
  show_average = FALSE,
  engine = "fmsb",   # ignored in items mode; fmsb used
  title = "Countries — one radar each (no average)"
)

# ================================
# Example D — One radar per country vs average (grid with average overlay)
# ================================
radar_with_p(
  data = countries,
  value_vars = colnames(countries),
  mode = "items",
  items_layout = "grid",
  show_average = TRUE,
  engine = "fmsb",
  title = "Countries — each vs average"
)

# ================================
# Example E — Single country (with average overlay)
# ================================
radar_with_p(
  data = countries,
  value_vars = colnames(countries),
  mode = "items",
  items_layout = "single",
  item_name = "Switzerland",  # change to any rowname
  show_average = TRUE,
  engine = "fmsb",
  title = "Switzerland vs average"
)

# ================================
# Example F — All countries on one radar (optional average overlay)
# ================================
radar_with_p(
  data = countries,
  value_vars = colnames(countries),
  mode = "items",
  items_layout = "overlay_all",
  show_average = TRUE,
  engine = "fmsb",
  title = "All countries vs average"
)
