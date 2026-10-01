# =============================================================================
# Fig7.R
#
# Draws panels B, C, D and E of Figure 7 from four small tables. Nothing else is
# needed: put Fig7B.csv, Fig7C.csv, Fig7D.csv and Fig7E.csv beside this file and
# run it.
#
#   Fig7B   the measured change in translation efficiency of each variant against
#           the change measured by luciferase assay
#   Fig7C   the features of the model, their coefficients and 95% intervals, and
#           how often stability selection kept each one
#   Fig7D   the sequences held out from fitting: observed against predicted TE
#   Fig7E   predicted TE against wild-type luciferase activity
#
# The drawing code is the code that made the published figures, so a panel drawn
# here is the same panel. Each section is marked with the panel it draws, and
# every figure is written as a .pdf and a .png into FIG_DIR.
# =============================================================================

rm(list = ls())
suppressPackageStartupMessages({ library(tidyverse); library(ggrepel); library(patchwork) })

# Where the csv files are. Set this if they are not beside this script; left
# empty, the places below are tried in turn.
DATA_DIR <- ""

if (!nzchar(DATA_DIR)) {
  cands <- c(".", "figure_data", file.path("..", "figure_data"),
             "public_5UTR_TE_model/figure_data",
             file.path("polysome_redo", "wt_lasso", "public_5UTR_TE_model", "figure_data"),
             file.path("/mnt/cloudBackup/GoogleDrive_thejovialjoy/Projects/rna_utranstab",
                       "polysome_redo", "wt_lasso", "public_5UTR_TE_model", "figure_data"))
  hit <- cands[file.exists(file.path(cands, "Fig7B.csv"))]
  if (!length(hit))
    stop("Fig7B.csv not found. Set DATA_DIR at the top of this script.\n  looked in:\n    ",
         paste(normalizePath(cands, mustWork = FALSE), collapse = "\n    "))
  DATA_DIR <- hit[1]
}
DATA_DIR <- normalizePath(DATA_DIR)

FIG_DIR <- file.path(DATA_DIR, "figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

# the shape of the coefficient panel, as the published figure used it
COEF_ROW_MM  <- 3.3      # height of one feature row, mm
COEF_BAR_MM  <- 62       # width of the coefficient panel, mm
COEF_FREQ_MM <- 20       # width of the selection-probability panel, mm
COEF_GAP_MM  <- 9        # space between the two, mm
COEF_BAND    <- 0.52     # thickness of an interval, as a share of the row
COEF_DOT     <- 0.40     # size of an estimate, as a share of the row
COEF_LEGEND  <- FALSE
CLASS_PT     <- 10       # size of the class names beside the ribbon
NEIGHBOUR_R  <- 0.02     # density radius in the observed-against-predicted figure

# the figure font, chosen once and used by every figure below
FONT <- "sans"
if (requireNamespace("systemfonts", quietly = TRUE)) {
  fams <- unique(systemfonts::system_fonts()$family)
  for (ff in c("Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"))
    if (ff %in% fams) { FONT <- ff; break }
}

cat("[Data]    <-", DATA_DIR, "\n")
cat("[Figures] ->", normalizePath(FIG_DIR), "| font:", FONT, "\n\n")

MM <- 1/25.4                                   # millimetres to inches

# the figure font, chosen once and used by every figure below
FONT <- "sans"
if (requireNamespace("systemfonts", quietly = TRUE)) {
  fams <- unique(systemfonts::system_fonts()$family)
  for (ff in c("Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"))
    if (ff %in% fams) { FONT <- ff; break }
}

BASE <- 8                                      # base font size, pt
AXIS_TITLE_PT <- BASE + 1.5                    # axis titles, which carry long labels
# and read small at the base size
cat("\n[Figures] ->", FIG_DIR, "| font:", FONT, "\n")

# palette: low saturation, harmonised, legible in greyscale and for the common
# forms of colour blindness (blue against terracotta rather than red/green)
C_UP   <- "#2F5D8A"   # steel blue    - increases TE
C_DN   <- "#C2593F"   # terracotta    - decreases TE
C_WT   <- "#1F7A6D"   # deep teal     - wild-type luciferase
C_LIB  <- "#2F5D8A"   # library constructs
C_FL   <- "#8A5A9E"   # full-length constructs
INK    <- "#1F1F1F"; AXIS <- "#333333"; SOFT <- "#5A5A5A"
LINE   <- "#1F1F1F"; BAND <- "#C9D1DC"; IDL <- "#A0A0A0"; QUAD <- "#EEF2F6"

theme_pub <- function(base = BASE) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(axis.line  = element_line(linewidth = .3, colour = AXIS),
          axis.ticks = element_line(linewidth = .3, colour = AXIS),
          axis.ticks.length = unit(1.1, "mm"),
          axis.text  = element_text(size = base - .5, colour = AXIS),
          axis.title = element_text(size = AXIS_TITLE_PT, colour = INK),
          axis.title.x = element_text(margin = margin(t = 4)),
          axis.title.y = element_text(margin = margin(r = 4)),
          plot.title = element_text(size = base - .5, colour = SOFT, hjust = 0,
                                    margin = margin(b = 5)),
          plot.title.position = "plot",
          legend.background = element_blank(), legend.key = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = base - .5, colour = INK),
          legend.key.size = unit(3, "mm"),
          plot.margin = margin(5, 6, 4, 4))
}
TXT <- (BASE - .5)/.pt                          # geom text size matching axis text

# P values written as a journal would: three decimals down to 0.001, then
# mantissa x 10^exponent
p_expr <- function(p) {
  if (!is.finite(p)) return(quote(italic(P) == "NA"))
  if (p <= .Machine$double.xmin) return(quote(italic(P) < 10^-300))
  if (p >= 1e-3) return(bquote(italic(P) == .(sprintf("%.3f", p))))
  e <- floor(log10(p)); m <- p/10^e
  if (round(m, 1) >= 10) { m <- m/10; e <- e + 1 }
  bquote(italic(P) == .(sprintf("%.1f", m)) %*% 10^.(e))
}
title_spearman <- function(rho, p)
  bquote("Spearman" ~ rho == .(sprintf("%.3f", rho)) * "," ~~ .(p_expr(p)))
TE_lab <- function(prefix) bquote(.(prefix) ~ log[2] * "(Poly / (Free + 40S))")

# feature labels: k-mers named by their length in the RNA alphabet (3-mer CGG),
# like the RBP motifs beside them


save_fig <- function(p, name, w_mm, h_mm) {
  base <- file.path(FIG_DIR, name); w <- w_mm*MM; h <- h_mm*MM
  ggsave(paste0(base, ".pdf"), p, width = w, height = h, device = cairo_pdf)
  ggsave(paste0(base, ".png"), p, width = w, height = h, dpi = 600, bg = "white")
  svg_ok <- FALSE
  if (requireNamespace("svglite", quietly = TRUE)) {
    svg_ok <- tryCatch({
      ggsave(paste0(base, ".svg"), p, width = w, height = h, device = svglite::svglite)
      TRUE
    }, error = function(e) {
      if (grDevices::dev.cur() > 1) try(grDevices::dev.off(), silent = TRUE)
      message("  svglite failed (", conditionMessage(e), "); using grDevices::svg")
      FALSE
    })
  }
  if (!svg_ok)
    ggsave(paste0(base, ".svg"), p, width = w, height = h, device = grDevices::svg)
  cat(sprintf("  %-36s %3.0f x %3.0f mm  (.pdf .png .svg)\n", name, w_mm, h_mm))
}
## ---- coefficients, grouped by feature class ------------------------------------
# Features are grouped by what they describe, so the figure says what KIND of
# sequence information the model uses. Each class is marked by a colour ribbon
# with its name on the left and a faint tint behind its rows; class colours are
# earth tones so they never compete with the blue/terracotta that encodes the
# direction of the effect.
# One coefficient (the CG k-mer) is several times larger than the rest and would
# compress every other bar into a sliver, so the axis is cut just past the
# second-largest effect and the cut bar carries a break mark.
feat_class <- function(id) dplyr::case_when(
  grepl("^kmer_", id)                    ~ "k-mer",
  grepl("^RBP_",  id)                    ~ "RBP motif",
  grepl("^ARE_",  id)                    ~ "ARE",
  grepl("^miRNA_|^miR_", id)             ~ "miRNA site",
  grepl("^uAUG_", id)                    ~ "Upstream AUG",
  grepl("^RG4|stemLoop|struc|psudoK|pseudoK|KISS|freeEnergy", id) ~ "Structure",
  grepl("way$", id)                      ~ "Conservation",
  TRUE                                   ~ "Composition")
CLASS_LEVELS <- c("Composition", "k-mer", "RBP motif", "ARE", "miRNA site",
                  "Upstream AUG", "Structure", "Conservation")
CLASS_COL  <- c(Composition = "#9A7F55", `k-mer` = "#5E8571",
                `RBP motif` = "#6E6A96", ARE = "#9A6653",
                `miRNA site` = "#4E7B9C", `Upstream AUG` = "#B08347",
                Structure = "#7D8B6A", Conservation = "#8A6E8F")
CLASS_TINT <- c(Composition = "#F6F1E8", `k-mer` = "#EEF4F0",
                `RBP motif` = "#F1F0F6", ARE = "#F7EFEC",
                `miRNA site` = "#EDF2F6", `Upstream AUG` = "#F8F2E9",
                Structure = "#F1F3EE", Conservation = "#F4F0F5")

prep_coef <- function(tbl, n_top) {
  tbl %>% dplyr::arrange(dplyr::desc(abs(coef))) %>%
    { if (is.finite(n_top)) dplyr::slice_head(., n = n_top) else . } %>%
    dplyr::mutate(label = pretty_feat(label, feature_id),
                  label = ifelse(duplicated(label) | duplicated(label, fromLast = TRUE),
                                 paste0(label, " [", feature_id, "]"), label),
                  class = factor(feat_class(feature_id), levels = CLASS_LEVELS),
                  direction = factor(direction, levels = c("increases", "decreases"))) %>%
    dplyr::filter(!is.na(class)) %>%
    dplyr::arrange(class, coef) %>%                   # within a class, by effect
    dplyr::mutate(label = factor(label, levels = unique(label)))
}

# Widths are set in millimetres from the content itself -- the longest feature
# name and the longest class name -- so no label can be clipped however many rows
# the figure has. Text width is estimated from the font size: about 0.55 em per
# character for Arial, with a margin.
char_mm <- function(pt, bold = FALSE) pt*0.3528*(if (bold) .62 else .56)

# The shape of the coefficient figure. A tall row height with narrow panels
# gives a figure close to square, which sits beside the other panels of a
# composite figure; the wide, short default reads better on its own.
# Two shapes are written for each coefficient figure: the near-square one, which
# sits beside another panel in a composite figure, and the long one, which reads
# better on its own. Each is a row height and three widths, in mm.
COEF_SHAPES <- list(
  square = list(row = 4.6, bar = 46, freq = 13, gap = 4, suffix = ""),
  wide   = list(row = 3.3, bar = 62, freq = 20, gap = 9, suffix = "_wide"))

COEF_ROW_MM  <- 4.6      # height of one feature row, mm
COEF_BAR_MM  <- 46       # width of the coefficient panel, mm
COEF_FREQ_MM <- 13       # width of the selection-probability panel, mm
COEF_GAP_MM  <- 4        # space between the two, mm
COEF_BAND    <- 0.52     # thickness of an interval, as a share of the row
COEF_DOT     <- 0.40     # size of an estimate, as a share of the row
COEF_LEGEND  <- FALSE    # the colours are explained by the axis, so the legend
# under the figure is off by default

make_coef_fig <- function(tbl, with_freq, row_mm = COEF_ROW_MM, txt = BASE - .5,
                          bar_mm = COEF_BAR_MM, freq_mm = COEF_FREQ_MM,
                          gap_mm = COEF_GAP_MM) {
  lo <- min(0, min(tbl$ci_low)); hi <- max(0, max(tbl$ci_high))
  xlim <- c(lo - .04*(hi - lo), hi + .06*(hi - lo))
  
  present <- levels(droplevels(tbl$class))
  facet <- facet_grid(class ~ ., scales = "free_y", space = "free_y")
  no_strip <- theme(strip.background = element_blank(), strip.text = element_blank(),
                    panel.spacing.y = unit(1.8, "mm"))
  tints <- lapply(present, function(k)
    geom_rect(data = data.frame(class = factor(k, levels = CLASS_LEVELS)),
              aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = CLASS_TINT[[k]]))
  
  ## ---- the class ribbon, sized to its longest name ----
  # The class names are written the same way for every class: along the ribbon
  # when every one of them fits, beside it otherwise. Deciding this once rather
  # than class by class keeps a short class from being turned while its
  # neighbour lies flat.
  rib_lab <- tbl %>% dplyr::group_by(class) %>%
    dplyr::summarise(n = dplyr::n(),
                     mid = levels(droplevels(label))[ceiling(dplyr::n()/2)],
                     .groups = "drop") %>%
    dplyr::mutate(name = as.character(class),
                  len_mm = nchar(name)*char_mm(CLASS_PT, bold = TRUE),
                  nudge = ifelse(n %% 2 == 0, .5, 0))
  rib_lab$vertical <- all(rib_lab$n*row_mm >= rib_lab$len_mm + 2)
  TILE <- 1.6                                             # ribbon width, mm
  horiz_mm <- if (!rib_lab$vertical[1]) max(rib_lab$len_mm) else 0
  rib_mm <- max(8, horiz_mm + TILE + 3.5)
  ribbon <- ggplot(tbl, aes(y = label)) + facet +
    lapply(present, function(k)
      geom_tile(data = tbl %>% dplyr::filter(class == k),
                aes(x = rib_mm - TILE/2, y = label),
                width = TILE, height = 1, fill = CLASS_COL[[k]])) +
    lapply(seq_len(nrow(rib_lab)), function(i) {
      r <- rib_lab[i, ]
      geom_text(data = r, aes(y = mid), inherit.aes = FALSE,
                x = if (r$vertical) rib_mm - TILE - 2.4 else rib_mm - TILE - 1.4,
                label = r$name, angle = if (r$vertical) 90 else 0,
                hjust = if (r$vertical) .5 else 1, vjust = .5, nudge_y = r$nudge,
                size = CLASS_PT/.pt, family = FONT, fontface = "bold",
                colour = CLASS_COL[[r$name]])
    }) +
    scale_x_continuous(limits = c(0, rib_mm), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void(base_family = FONT) + no_strip +
    theme(plot.margin = margin(0, 0, 0, 0))
  
  ## ---- the coefficients: estimate and 95% interval ----
  # The interval is drawn as a band that fades outwards rather than as a whisker:
  # each interval is cut into NSLICE pieces whose opacity falls with the distance
  # from the estimate, so the eye is drawn to the middle of the interval.
  NSLICE <- 48
  slices <- purrr::map_dfr(seq_len(nrow(tbl)), function(i) {
    r <- tbl[i, ]
    e <- seq(r$ci_low, r$ci_high, length.out = NSLICE + 1)
    mid <- (utils::head(e, -1) + e[-1])/2
    half <- max((r$ci_high - r$ci_low)/2, 1e-12)
    tibble::tibble(label = r$label, class = r$class, direction = r$direction,
                   x0 = utils::head(e, -1), x1 = e[-1],
                   a = (1 - abs(mid - r$coef)/half)^1.6)
  })
  bar <- ggplot(tbl, aes(coef, label)) + facet + tints +
    geom_vline(xintercept = 0, colour = AXIS, linewidth = .3) +
    geom_segment(data = slices, aes(x = x0, xend = x1, y = label, yend = label,
                                    colour = direction, alpha = a),
                 inherit.aes = FALSE, linewidth = row_mm*COEF_BAND, lineend = "butt") +
    geom_point(aes(fill = direction), shape = 21, size = row_mm*COEF_DOT,
               colour = "white", stroke = .45) +
    scale_colour_manual(values = c(increases = C_UP, decreases = C_DN), guide = "none") +
    scale_alpha_continuous(range = c(.06, .55), guide = "none") +
    scale_fill_manual(values = c(increases = C_UP, decreases = C_DN),
                      labels = c(increases = "Increases TE", decreases = "Decreases TE")) +
    scale_x_continuous(breaks = scales::breaks_pretty(4), expand = c(0, 0)) +
    coord_cartesian(xlim = xlim, clip = "off") +
    labs(x = "Coefficient (standardised)  ·  95% interval", y = NULL) +
    theme_pub() + no_strip +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_text(colour = INK, size = txt, margin = margin(r = 2)),
          axis.title.x = element_text(margin = margin(t = 1.5)),
          legend.position = if (COEF_LEGEND) "bottom" else "none",
          legend.justification = "center", legend.margin = margin(2, 0, 0, 0))
  
  ## ---- absolute widths, so the names always have room ----
  lab_mm <- max(nchar(as.character(tbl$label)))*char_mm(txt) + 3
  panels <- list(ribbon, bar); w <- c(rib_mm, bar_mm)
  if (with_freq) {
    # from zero, so the length of each bar is the selection probability itself
    # the value is written on the bar: inside it when the bar is long enough to
    # hold the text, just outside it when it is not
    freq_lab <- tbl %>%
      dplyr::mutate(txt = sprintf("%.0f%%", 100*sel_freq),
                    inside = sel_freq > .28,
                    x = ifelse(inside, sel_freq - .02, sel_freq + .02),
                    hj = ifelse(inside, 1, 0),
                    col = ifelse(inside, "white", SOFT))
    freq <- ggplot(tbl, aes(sel_freq, label)) + facet + tints +
      geom_col(width = .68, fill = "#8E9AA8") +
      geom_text(data = freq_lab, aes(x = x, y = label, label = txt, hjust = hj,
                                     colour = col),
                inherit.aes = FALSE, size = (txt - 1.2)/.pt, family = FONT,
                fontface = "bold", show.legend = FALSE) +
      scale_colour_identity() +
      scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1),
                         labels = scales::percent_format(accuracy = 1),
                         expand = expansion(mult = c(0, .02))) +
      labs(x = "Selection\nprobability", y = NULL) +
      theme_pub() + no_strip +
      theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
            axis.text.y = element_blank())
    # an explicit gap, so the tick labels of the two x-axes cannot meet
    panels <- c(panels, list(plot_spacer(), freq)); w <- c(w, gap_mm, freq_mm)
  }
  p <- wrap_plots(panels, nrow = 1, widths = unit(w, "mm")) &
    theme(plot.margin = margin(4, 3, 3, 2))
  list(plot = p,
       width_mm  = sum(w) + lab_mm + (length(w) - 1)*3 + 14,
       height_mm = nrow(tbl)*row_mm + (length(present) - 1)*1.8 +
         if (COEF_LEGEND) 30 else 22)
}

# the x-axis label and the builder the mutant/WT panel uses
TE_X      <- bquote(Delta * log[2] * "(Poly / (Free + 40S))")

mk_mtwt <- function(df, xvar, xlab, stats, subtitle = NULL,
                    colour_by = NULL, one_colour = C_LIB, labs_ = NULL) {
  p <- ggplot(df, aes(.data[[xvar]], luciferase_log2_mt_over_wt)) +
    # the two quadrants where the MPRA and the luciferase agree on the direction
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = QUAD) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = QUAD) +
    geom_hline(yintercept = 0, colour = "#B8BEC6", linewidth = .3) +
    geom_vline(xintercept = 0, colour = "#B8BEC6", linewidth = .3) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = LINE,
                fill = BAND, alpha = .45, linewidth = .5)
  p <- p + if (is.null(colour_by))
    geom_point(shape = 21, size = 2.6, fill = one_colour, colour = "white", stroke = .45)
  else
    geom_point(aes(fill = .data[[colour_by]]), shape = 21, size = 2.6,
               colour = "white", stroke = .45, alpha = .72)
  p <- p +
    ggrepel::geom_text_repel(aes(label = Gene), size = TXT, family = FONT,
                             fontface = "italic", colour = INK, box.padding = .45,
                             point.padding = .3, min.segment.length = Inf, seed = 41,
                             max.overlaps = Inf) +
    labs(title = title_spearman(stats$spearman, stats$spearman_p), subtitle = subtitle,
         x = xlab, y = expression(log[2] ~ "translation efficiency (mutant / WT)")) +
    theme_pub() + theme(plot.subtitle = element_text(size = BASE - 1, colour = SOFT))
  if (is.null(colour_by)) return(p)
  p <- p + scale_fill_manual(values = c(Library = C_LIB, `Full-length` = C_FL),
                             labels = if (is.null(labs_)) ggplot2::waiver() else labs_)
  if (is.null(labs_)) p + guides(fill = "none")
  else p + theme(legend.position = "top", legend.justification = "right",
                 legend.direction = "vertical", legend.title = element_blank(),
                 legend.background = element_blank(), legend.key = element_blank(),
                 legend.margin = margin(0, 0, 1, 0), legend.spacing.y = unit(0, "mm"),
                 legend.key.height = unit(3.0, "mm"), legend.key.width = unit(3.0, "mm"),
                 legend.text = element_text(size = BASE - 1, colour = INK),
                 plot.title = element_text(margin = margin(b = 0))) +
    guides(fill = guide_legend(ncol = 1, override.aes = list(size = 2.4, alpha = 1)))
}

## =============================================================================
## Panel C  --  the features of the model
## =============================================================================
cf_all <- readr::read_csv(file.path(DATA_DIR, "Fig7C.csv"), show_col_types = FALSE) %>%
  dplyr::transmute(feature_id, label,
                   coef = coefficient, ci_low, ci_high,
                   sel_freq = selection_probability,
                   class = factor(feat_class(feature_id), levels = CLASS_LEVELS),
                   direction = factor(direction, levels = c("increases", "decreases"))) %>%
  dplyr::filter(!is.na(class)) %>%
  dplyr::arrange(class, coef) %>%             # within a class, by effect
  dplyr::mutate(label = factor(label, levels = unique(label)))

for (with_freq in c(TRUE, FALSE)) {
  g <- make_coef_fig(cf_all, with_freq)
  save_fig(g$plot, paste0("Fig7C", if (with_freq) "_with_probability" else ""),
           g$width_mm, g$height_mm)
}
cat("  (the sequence logos are drawn by model_5UTR_TE.R, from the ATtRACT matrices)\n")


## =============================================================================
## Panel B  --  measured change in TE against the luciferase ratio
## =============================================================================
b <- readr::read_csv(file.path(DATA_DIR, "Fig7B.csv"), show_col_types = FALSE) %>%
  dplyr::rename(luciferase_log2_mt_over_wt = luciferase_log2_mt_over_WT)
n_by  <- b %>% dplyr::count(construct)
lab_b <- setNames(sprintf("%s (n = %d)", n_by$construct, n_by$n), n_by$construct)
st_b  <- list(spearman   = cor(b$delta_TE, b$luciferase_log2_mt_over_wt, method = "spearman"),
              spearman_p = suppressWarnings(cor.test(b$delta_TE, b$luciferase_log2_mt_over_wt,
                                                     method = "spearman"))$p.value)
save_fig(mk_mtwt(b, "delta_TE", TE_X, st_b, colour_by = "construct", labs_ = lab_b),
         "Fig7B", 86, 88)


## =============================================================================
## Panel D  --  the held-out sequences, coloured by how crowded they are
## =============================================================================
fig_best <- readr::read_csv(file.path(DATA_DIR, "Fig7D.csv"), show_col_types = FALSE) %>%
  dplyr::rename(predicted = predicted_TE, observed = observed_TE)
rho2 <- cor(fig_best$predicted, fig_best$observed, method = "spearman")
p2_spearman <- suppressWarnings(cor.test(fig_best$predicted, fig_best$observed,
                                         method = "spearman"))$p.value
n_fig2 <- nrow(fig_best)

leg_theme <- theme(legend.position = "right",
                   legend.title = element_text(size = BASE - .5, colour = INK,
                                               margin = margin(b = 6), lineheight = 1.1),
                   legend.text = element_text(size = BASE - 1, colour = SOFT),
                   legend.margin = margin(0, 0, 0, 3))
bar <- guide_colourbar(barwidth = unit(2.4, "mm"), barheight = unit(20, "mm"),
                       ticks.colour = "white", title.position = "top")
# sample size above the word Density
leg_title <- bquote(atop(italic(n) == .(format(n_fig2, big.mark = ",")), "Density"))

sx <- (fig_best$observed  - min(fig_best$observed))/diff(range(fig_best$observed))
sy <- (fig_best$predicted - min(fig_best$predicted))/diff(range(fig_best$predicted))
d2 <- outer(sx, sx, "-")^2 + outer(sy, sy, "-")^2
fb <- fig_best %>%
  dplyr::mutate(neighbours = rowSums(d2 <= NEIGHBOUR_R^2) - 1L) %>%   # not itself
  dplyr::arrange(neighbours)                 # densest drawn last, on top
p2 <- ggplot(fb, aes(predicted, observed)) +
  geom_point(aes(colour = neighbours), size = 1.1, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "black",
              fill = "grey70", alpha = .45, linewidth = .5) +
  scale_colour_viridis_c(option = "magma", end = .97, name = leg_title,
                         breaks = scales::breaks_pretty(3), guide = bar) +
  labs(title = title_spearman(rho2, p2_spearman),
       x = TE_lab("Predicted"), y = TE_lab("Observed")) +
  theme_pub() + leg_theme
save_fig(p2, "Fig7D", 98, 80)   # wider by the colour bar


## =============================================================================
## Panel E  --  predicted TE against wild-type luciferase
## =============================================================================
lw <- readr::read_csv(file.path(DATA_DIR, "Fig7E.csv"), show_col_types = FALSE) %>%
  dplyr::rename(luc_log2 = luciferase_log2, rep_lo = replicate_lo, rep_hi = replicate_hi)
rho3 <- cor(lw$predicted_TE, lw$luc_log2, method = "spearman")
p3_spearman <- suppressWarnings(cor.test(lw$predicted_TE, lw$luc_log2,
                                         method = "spearman"))$p.value
p3 <- ggplot(lw, aes(predicted_TE, luc_log2)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = LINE,
              fill = BAND, alpha = .45, linewidth = .5) +
  geom_errorbar(aes(ymin = rep_lo, ymax = rep_hi), width = 0, linewidth = .35,
                colour = C_WT) +
  geom_point(shape = 21, size = 2.6, fill = C_WT, colour = "white", stroke = .45) +
  ggrepel::geom_text_repel(aes(label = Gene), size = TXT, family = FONT,
                           fontface = "italic", colour = INK, box.padding = .45,
                           point.padding = .3, min.segment.length = Inf, seed = 41) +
  labs(title = title_spearman(rho3, p3_spearman), x = TE_lab("Predicted"),
       y = expression(log[2] ~ "luciferase activity")) +
  theme_pub()
save_fig(p3, "Fig7E", 84, 80)

cat("\n=== Done ===\n")