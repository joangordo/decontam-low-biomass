# Summary figure for the README: what each cleaning choice removes.
#
# Reads the cached simulation results and writes a light and a dark version
# to figures/ (GitHub shows the one that matches the reader's theme).
# Run from the project root:  Rscript analysis/readme_figure.R

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

sim <- readRDS("analysis/results/simulation_results.rds")$results

# Very-low-biomass studies (sample DNA about equal to contaminant DNA) with
# realistic, leaky controls: 4 numbers of controls x 20 studies = 80 per point.
d <- sim |>
  filter(leak > 0, biomass_ratio == 1,
         method %in% c("frequency", "prevalence", "in negative control",
                       "abundance >1%")) |>
  group_by(conc_noise, method, threshold) |>
  summarise(contam = mean(contam_reads_removed),
            taxa = mean(true_taxa_lost),
            reads = mean(true_reads_lost), .groups = "drop") |>
  mutate(
    safe = reads <= 1 & taxa <= 2,          # the safety rule used in Part 2
    family = case_when(method == "frequency" ~ "frequency",
                       method == "prevalence" ~ "prevalence",
                       TRUE ~ "filter"),
    family = factor(family, c("frequency", "prevalence", "filter")),
    panel = factor(ifelse(conc_noise < 0.3, "Precise DNA quantification",
                          "Imprecise DNA quantification"),
                   c("Precise DNA quantification", "Imprecise DNA quantification")),
    taxa_plot = pmax(taxa, 0.3)
  )

# Label positions are set by hand so nothing overlaps (x is on a log scale).
P <- levels(d$panel)
lab <- tibble::tribble(
  ~panel, ~x,   ~y,    ~hjust, ~label,                                ~family,      ~bold,
  P[1],   0.41, 0.700, 0.5,    "0.10",                                "frequency",  FALSE,
  P[1],   1.44, 0.925, 0.5,    "0.25",                                "frequency",  FALSE,
  P[1],   23.8, 0.910, 0.5,    "0.50",                                "frequency",  FALSE,
  P[1],   2.2,  0.815, 0,      "frequency method",                    "frequency",  TRUE,
  P[1],   0.85, 0.106, 1,      "0.10",                                "prevalence", FALSE,
  P[1],   4.8,  0.290, 0,      "0.25",                                "prevalence", FALSE,
  P[1],   14,   0.765, 0,      "0.50",                                "prevalence", FALSE,
  P[1],   1.3,  0.040, 0,      "prevalence method",                   "prevalence", TRUE,
  P[2],   0.61, 0.020, 0.5,    "0.10",                                "frequency",  FALSE,
  P[2],   2.29, 0.125, 0.5,    "0.25",                                "frequency",  FALSE,
  P[2],   23.8, 0.375, 0.5,    "0.50",                                "frequency",  FALSE,
  P[2],   5,    0.195, 0,      "frequency method",                    "frequency",  TRUE,
  P[2],   0.76, 0.165, 0.5,    "0.10",                                "prevalence", FALSE,
  P[2],   3.6,  0.350, 1,      "0.25",                                "prevalence", FALSE,
  P[2],   10.2, 0.800, 1,      "0.50",                                "prevalence", FALSE,
  P[2],   4.5,  0.600, 1,      "prevalence method",                   "prevalence", TRUE
)
filters <- tibble::tribble(
  ~x, ~y,    ~hjust, ~label,
  48, 1.040, 0.5,    "remove every taxon\nseen in a blank",
  45, 0.073, 1,      "keep taxa\nabove 1%"
)
filters <- tidyr::crossing(panel = P, filters)
lab$panel <- factor(lab$panel, P); filters$panel <- factor(filters$panel, P)

themes <- list(
  light = list(bg = "#ffffff", ink = "#1f2328", ink2 = "#59636e", grid = "#e1e0d9",
               pal = c(frequency = "#2a78d6", prevalence = "#eb6834", filter = "#1baf7a")),
  dark = list(bg = "#0d1117", ink = "#f0f6fc", ink2 = "#9198a1", grid = "#2c2c2a",
              pal = c(frequency = "#3987e5", prevalence = "#d95926", filter = "#199e70"))
)

dir.create("figures", showWarnings = FALSE)
for (mode in names(themes)) {
  th <- themes[[mode]]
  p <- ggplot(d, aes(taxa_plot, contam / 100, colour = family)) +
    annotate("rect", xmin = 0.25, xmax = 2, ymin = -Inf, ymax = Inf,
             fill = th$ink2, alpha = 0.08) +
    annotate("text", x = 0.28, y = 1.04, hjust = 0, size = 3, colour = th$ink2,
             label = "loses at most 2 real taxa") +
    geom_path(data = filter(d, family != "filter"),
              aes(group = interaction(family, panel)), linewidth = 0.7) +
    geom_point(aes(shape = safe), size = 3, stroke = 1.1, fill = th$bg) +
    geom_text(data = filter(lab, !bold), aes(x, y, label = label, hjust = hjust),
              inherit.aes = FALSE, size = 3, colour = th$ink2) +
    geom_text(data = filter(lab, bold), aes(x, y, label = label, hjust = hjust, colour = family),
              inherit.aes = FALSE, size = 3.4, fontface = "bold") +
    geom_text(data = filters, aes(x, y, label = label, hjust = hjust), inherit.aes = FALSE,
              size = 3, colour = th$ink2, lineheight = 0.9) +
    facet_wrap(~panel) +
    scale_x_log10(breaks = c(0.3, 1, 3, 10, 30), labels = c("0.3", "1", "3", "10", "30"),
                  limits = c(0.25, 90)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1.08),
                       breaks = seq(0, 1, 0.25)) +
    scale_colour_manual(values = th$pal, guide = "none") +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                       labels = c(`TRUE` = "meets the safety rule: loses at most 1% of real reads and 2 real taxa",
                                  `FALSE` = "removes real biology"),
                       breaks = c("TRUE", "FALSE"), name = NULL) +
    labs(
      title = "Cleaning very-low-biomass 16S data: what you remove vs what you lose",
      subtitle = "Simulated studies with realistic (leaky) negative controls; each point is the mean of 80 studies",
      x = "Real taxa removed per study (log scale; further right is worse)",
      y = "Contaminant reads removed (higher is better)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.background = element_rect(fill = th$bg, colour = NA),
      panel.background = element_rect(fill = th$bg, colour = NA),
      panel.grid.major = element_line(colour = th$grid, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      text = element_text(colour = th$ink),
      axis.text = element_text(colour = th$ink2),
      axis.title = element_text(colour = th$ink2, size = 10),
      strip.text = element_text(colour = th$ink, face = "bold", size = 11, hjust = 0),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(colour = th$ink2, size = 10),
      legend.position = "bottom",
      legend.text = element_text(colour = th$ink2, size = 9),
      legend.key = element_rect(fill = th$bg, colour = NA),
      plot.margin = margin(12, 16, 8, 12)
    ) +
    guides(shape = guide_legend(override.aes = list(colour = th$ink2)))
  ggsave(sprintf("figures/readme-tradeoff-%s.png", mode), p, width = 10, height = 5.2,
         dpi = 160, bg = th$bg)
}
message("Saved figures/readme-tradeoff-light.png and figures/readme-tradeoff-dark.png")
