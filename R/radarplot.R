  * Automatically decides **t-test** (if there are exactly 2 groups) or **one-way ANOVA** (≥3 groups), or you can force either.
* Can draw with either **base `fmsb::radarchart`** or a clean **`ggplot2` polar** radar.
* Adds nicely formatted p-values to the axis labels (ggplot) or around the chart (fmsb).
* Returns the plot (for ggplot) **and** a p-value table invisibly.

```r
# ---- Packages you need ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix)
library(fmsb)

# ---- One function to rule them all ----
radar_with_p <- function(
    data,
    group_var,
    value_vars,
    engine = c("ggplot", "fmsb"),
    p_method = c("auto", "t.test", "anova"),
    max_scale = NULL,         # if NULL, auto-compute from data
    colors = NULL,            # named vector ok for ggplot; for fmsb a vector is fine
    title = NULL,
    text_radius = 1.35        # only used for fmsb: where p-values are written
){
  engine   <- match.arg(engine)
  p_method <- match.arg(p_method)

  # basic checks
  stopifnot(group_var %in% names(data))
  stopifnot(all(value_vars %in% names(data)))
  df <- data %>%
    select(all_of(c(group_var, value_vars))) %>%
    mutate(across(all_of(value_vars), as.numeric)) %>%
    drop_na(all_of(c(group_var, value_vars)))

  # group info
  grp <- df[[group_var]]
  grp_fac <- as.factor(grp)
  n_groups <- nlevels(grp_fac)

  # choose test
  test_mode <- switch(p_method,
                      "t.test" = "t.test",
                      "anova"  = "anova",
                      "auto"   = if (n_groups == 2) "t.test" else "anova"
  )

  # p-value formatter
  fmt_p <- function(p){
    if (is.na(p)) return("NA")
    if (p < 0.001) "<0.001***"
    else if (p < 0.01) sprintf("%.3f**", p)
    else if (p < 0.05) sprintf("%.3f*",  p)
    else "ns"
  }

  # compute p-values per variable
  p_tbl <- purrr::map_dfr(value_vars, function(v){
    if (test_mode == "t.test" && n_groups == 2) {
      p <- tryCatch(
        t.test(df[[v]] ~ grp_fac)$p.value,
        error = function(e) NA_real_
      )
    } else {
      # one-way ANOVA via rstatix
      aov_res <- rstatix::anova_test(data = df, formula = as.formula(paste(v, "~", group_var)))
      # rstatix keeps a row with term equal to group_var; pick that p
      p <- tryCatch(aov_res$p[aov_res$Effect %in% c(group_var, "group")][1], error = function(e) NA_real_)
      if (is.na(p)) {
        # fallback if structure differs
        p <- tryCatch(aov_res$p[1], error = function(e) NA_real_)
      }
    }
    tibble(variable = v, p = p, p_label = fmt_p(p))
  }) %>%
    mutate(axis_with_p = paste0(variable, "\n(", p_label, ")"))

  # group means
  means_wide <- df %>%
    group_by(.data[[group_var]]) %>%
    summarise(across(all_of(value_vars), ~mean(.x, na.rm = TRUE)), .groups = "drop")

  # auto max scale if needed (nice round up to /5)
  if (is.null(max_scale)) {
    mx <- means_wide %>% select(all_of(value_vars)) %>% as.matrix() %>% max(na.rm = TRUE)
    if (!is.finite(mx)) mx <- 1
    max_scale <- max(1, ceiling(mx / 5) * 5)
  }

  # default colors if none given
  if (is.null(colors)) {
    # base palette without extra deps
    cols <- grDevices::rainbow(n_groups, s = 0.6, v = 0.9)
    if (engine == "ggplot") {
      names(cols) <- levels(grp_fac)
      colors <- cols
    } else {
      colors <- cols
    }
  }

  # --- Draw ---
  if (engine == "ggplot") {
    plot_long <- means_wide %>%
      pivot_longer(cols = all_of(value_vars), names_to = "variable", values_to = "value") %>%
      left_join(select(p_tbl, variable, axis_with_p), by = "variable") %>%
      mutate(
        axis_with_p = factor(axis_with_p, levels = p_tbl$axis_with_p),
        !!group_var := factor(.data[[group_var]], levels = levels(grp_fac))
      )

    g <- ggplot(
      plot_long,
      aes(x = axis_with_p, y = value,
          group = .data[[group_var]],
          color = .data[[group_var]],
          fill  = .data[[group_var]])
    ) +
      geom_polygon(alpha = 0.18, linewidth = 1, show.legend = TRUE) +
      geom_point(size = 2.5, show.legend = FALSE) +
      coord_polar(start = 0) +
      scale_y_continuous(limits = c(0, max_scale), breaks = seq(0, max_scale, by = max(5, round(max_scale/6)))) +
      labs(
        title = title %||% "Radar with p-values",
        x = NULL, y = NULL, fill = "Group", color = "Group",
        caption = "* p<0.05  ** p<0.01  *** p<0.001"
      ) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "right",
        panel.grid.major = element_line(color = "grey85"),
        axis.text.x = element_text(size = 11, lineheight = .9),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
      )

    print(g)
    invisible(list(plot = g, p_values = p_tbl, means = means_wide, engine = "ggplot",
                   test = test_mode, max_scale = max_scale))

  } else { # engine == "fmsb"
    # build matrix for fmsb: first two rows are max/min
    mat <- means_wide
    rownames(mat) <- mat[[group_var]]
    mat <- mat[, value_vars, drop = FALSE]
    max_min <- rbind(rep(max_scale, length(value_vars)),
                     rep(0,         length(value_vars)))
    colnames(max_min) <- colnames(mat)
    plot_mat <- rbind(max_min, mat)

    # semi-transparent fill for groups
    base_cols <- if (!is.null(names(colors))) unname(colors) else colors
    line_cols <- base_cols
    fill_cols <- grDevices::adjustcolor(base_cols, alpha.f = 0.35)

    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    par(mar = c(1,1,2,1))

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

    # p-values text around the chart
    k <- length(value_vars)
    angles <- seq(from = pi/2, to = -3*pi/2, length.out = k + 1)[-(k + 1)]
    for (i in seq_len(k)) {
      text(cos(angles[i]) * text_radius,
           sin(angles[i]) * text_radius,
           labels = p_tbl$p_label[i],
           cex = 0.8, col = "darkred", xpd = NA)
    }

    legend("topright",
           legend = rownames(plot_mat)[-(1:2)],
           bty = "n", pch = 20, col = fill_cols,
           text.col = "black", cex = 1.0, pt.cex = 2)

    title(main = title %||% "Radar with p-values")

    invisible(list(plot = NULL, p_values = p_tbl, means = means_wide, engine = "fmsb",
                   test = test_mode, max_scale = max_scale))
  }
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
```

## Examples

### 1) Your first dataset → `fmsb` + automatic t-tests (2 groups)

```r
set.seed(123)
data <- data.frame(
  group = rep(c("Control", "Treatment"), each = 50),
  Clarity = c(rnorm(50, mean = 6, sd = 1.5), rnorm(50, mean = 7.5, sd = 1.5)),
  Speed = c(rnorm(50, mean = 7, sd = 1.5), rnorm(50, mean = 5, sd = 1.5)),
  Accuracy = c(rnorm(50, mean = 8, sd = 1.5), rnorm(50, mean = 8.5, sd = 1.5)),
  Innovation = c(rnorm(50, mean = 5, sd = 1.5), rnorm(50, mean = 7, sd = 1.5)),
  Strategy = c(rnorm(50, mean = 6.5, sd = 1.5), rnorm(50, mean = 6.8, sd = 1.5))
)

res1 <- radar_with_p(
  data = data,
  group_var = "group",
  value_vars = colnames(data)[-1],
  engine = "fmsb",            # draw with fmsb
  p_method = "auto",           # auto = t.test for 2 groups
  title = "Control vs Treatment"
)
# res1$p_values contains the p-value table
```

### 2) Gene-severity example → `ggplot` + automatic ANOVA (≥3 groups)

```r
set.seed(42)
n_per_group <- 20
groups <- c("Healthy", "Mild", "Moderate", "Severe")
gene_data <- data.frame(
  patient_id = paste0("Patient_", 1:(n_per_group * 4)),
  severity = factor(rep(groups, each = n_per_group), levels = groups, ordered = TRUE)
) %>%
  mutate(
    sev_numeric = as.numeric(severity),
    Gene_A = rnorm(n(), mean = 10 + (sev_numeric * 1.5), sd = 2),
    Gene_B = rnorm(n(), mean = 5 + (sev_numeric * 2), sd = 2.5),
    Gene_C = rnorm(n(), mean = 25 - (sev_numeric * 2), sd = 3),
    Gene_D = rnorm(n(), mean = 15 + c(5, 0, 1, 6)[sev_numeric], sd = 2),
    Gene_E = rnorm(n(), mean = 12, sd = 3),
    Gene_F = rnorm(n(), mean = 20, sd = 3)
  ) %>% select(-sev_numeric)

severity_palette <- c("Healthy"="#4CAF50","Mild"="#FFC107","Moderate"="#FF9800","Severe"="#F44336")

res2 <- radar_with_p(
  data = gene_data,
  group_var = "severity",
  value_vars = c("Gene_A","Gene_B","Gene_C","Gene_D","Gene_E","Gene_F"),
  engine = "ggplot",           # draw with ggplot
  p_method = "auto",           # auto = ANOVA for >=3 groups
  max_scale = 30,              # optional; otherwise auto
  colors = severity_palette,
  title = "Gene Expression Profile by Disease Severity"
)
# print(res2$plot) if you didn't already see it; res2$p_values has the stats
```

That’s it—one function, both styles, automatic stats. If you want it to **always** use t-tests or **always** ANOVA, just set `p_method = "t.test"` or `"anova"` explicitly.
