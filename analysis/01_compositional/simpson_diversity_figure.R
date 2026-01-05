# ============================================================================
# Simpson Diversity Figure with Significance Brackets
# Publication-quality figure for HMSC EV miRNA compositional analysis
# ============================================================================
#
# USAGE:
#   source(here("analysis", "01_compositional", "simpson_diversity_figure.R"))
#
# INPUT:
#   results/01_compositional/diversity_metrics_by_sample.csv
#
# OUTPUT:
#   results/01_compositional/simpson_diversity_significance.png
#   results/01_compositional/simpson_diversity_significance.pdf
#
# ============================================================================

library(here)
library(dplyr)
library(ggplot2)

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------

input_file <- here("results", "01_compositional", "diversity_metrics_by_sample.csv")
output_dir <- here("results", "01_compositional")

# ----------------------------------------------------------------------------
# Load Data
# ----------------------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, 
       "\nRun CompositionalAnalysis.Rmd first to generate diversity metrics.")
}

diversity_df <- read.csv(input_file)

# Set treatment factor levels in biological order
# Relabel CF_Asp_Neg to CF for publication
diversity_df$Treatment <- factor(diversity_df$Treatment, 
                                  levels = c("UNSTIM", "HC", "CF_Asp_Neg"),
                                  labels = c("UNSTIM", "HC", "CF"))

cat("Loaded diversity data:", nrow(diversity_df), "samples\n")
cat("Treatment groups:", paste(levels(diversity_df$Treatment), collapse = ", "), "\n\n")

# ----------------------------------------------------------------------------
# Calculate Tukey HSD Results for Annotations
# ----------------------------------------------------------------------------

# Run ANOVA and Tukey HSD on Simpson diversity
simpson_aov <- aov(Simpson ~ Treatment, data = diversity_df)
simpson_tukey <- TukeyHSD(simpson_aov)

# Extract results into dataframe
tukey_df <- as.data.frame(simpson_tukey$Treatment)
tukey_df$comparison <- rownames(tukey_df)

# Parse comparison names into group1 and group2
tukey_results <- tukey_df %>%
  mutate(
    # Split "HC-UNSTIM" into group1="HC", group2="UNSTIM"
    group1 = sapply(strsplit(comparison, "-"), `[`, 1),
    group2 = sapply(strsplit(comparison, "-"), `[`, 2),
    p_value = `p adj`
  ) %>%
  dplyr::select(group1, group2, p_value)

# Add significance labels
tukey_results$label <- case_when(
  tukey_results$p_value < 0.001 ~ "***",
  tukey_results$p_value < 0.01  ~ "**",
  tukey_results$p_value < 0.05  ~ "*",
  TRUE ~ "ns"
)

# Report ANOVA result
anova_pval <- summary(simpson_aov)[[1]][["Pr(>F)"]][1]
cat("Simpson diversity ANOVA p-value:", format(anova_pval, digits = 4), "\n")
cat("Tukey HSD results:\n")
print(tukey_results)

# Only show significant comparisons
show_ns <- FALSE
if(!show_ns) {
  tukey_results <- tukey_results %>% filter(p_value < 0.05)
}
# ----------------------------------------------------------------------------
# Calculate bracket positions
# ----------------------------------------------------------------------------

# Get y-axis maximum for positioning brackets
y_max <- max(diversity_df$Simpson)
y_range <- diff(range(diversity_df$Simpson))

# Map treatment names to x positions
treatment_positions <- c("UNSTIM" = 1, "HC" = 2, "CF" = 3)

# Calculate bracket heights (stagger multiple brackets)
tukey_results <- tukey_results %>%
  mutate(
    x1 = treatment_positions[group1],
    x2 = treatment_positions[group2],
    x_mid = (x1 + x2) / 2,
    # Stagger brackets based on span width
    y_bracket = y_max + y_range * (0.08 + 0.12 * (row_number() - 1))
  )

# ----------------------------------------------------------------------------
# Create Publication-Quality Figure
# ----------------------------------------------------------------------------

# Color palette (colorblind-friendly)
treatment_colors <- c("UNSTIM" = "#66C2A5",   # Teal
                      "HC" = "#FC8D62",        # Coral
                      "CF" = "#8DA0CB")        # Periwinkle

p_simpson <- ggplot(diversity_df, aes(x = Treatment, y = Simpson, fill = Treatment)) +
  
  # Boxplot
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6, 
               color = "gray30", size = 0.5) +
  
 # Individual points
  geom_jitter(width = 0.15, alpha = 0.7, size = 2.5, shape = 21, 
              color = "gray30", stroke = 0.5) +
  
  # Mean diamond
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5, 
               fill = "white", color = "black", stroke = 1) +
  
  # Colors
  scale_fill_manual(values = treatment_colors) +
  
  # Axis labels
  labs(
    x = NULL,
    y = "Simpson Diversity Index (1-D)"
  ) +
  
  # Theme
  theme_classic(base_size = 12) +
  theme(
    # Remove legend (treatments labeled on x-axis)
    legend.position = "none",
    
    # Axis formatting
    axis.text.x = element_text(size = 11, color = "black", face = "bold"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    
    # Plot margins
    plot.margin = margin(t = 30, r = 15, b = 10, l = 10)
  ) +
  
  # Y-axis scale (add room for brackets)
  scale_y_continuous(
    limits = c(min(diversity_df$Simpson) - 0.02, 
               y_max + y_range * (0.08 + 0.12 * nrow(tukey_results) + 0.05)),
    expand = c(0, 0)
  )

# Add significance brackets
if(nrow(tukey_results) > 0) {
  
  bracket_height <- y_range * 0.02  # Height of vertical segments
  
  # Build dataframes for all bracket segments (avoids ggplot2 lazy evaluation issue)
  # Left vertical segments
  left_segs <- data.frame(
    x = tukey_results$x1,
    xend = tukey_results$x1,
    y = tukey_results$y_bracket,
    yend = tukey_results$y_bracket + bracket_height
  )
  
  # Horizontal segments
  horiz_segs <- data.frame(
    x = tukey_results$x1,
    xend = tukey_results$x2,
    y = tukey_results$y_bracket + bracket_height,
    yend = tukey_results$y_bracket + bracket_height
  )
  
  # Right vertical segments
  right_segs <- data.frame(
    x = tukey_results$x2,
    xend = tukey_results$x2,
    y = tukey_results$y_bracket + bracket_height,
    yend = tukey_results$y_bracket
  )
  
  # Text labels
  text_labels <- data.frame(
    x = tukey_results$x_mid,
    y = tukey_results$y_bracket + bracket_height + y_range * 0.02,
    label = tukey_results$label
  )
  
  # Add all bracket components using data parameter
  p_simpson <- p_simpson +
    geom_segment(data = left_segs, aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE, color = "black", size = 0.5) +
    geom_segment(data = horiz_segs, aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE, color = "black", size = 0.5) +
    geom_segment(data = right_segs, aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE, color = "black", size = 0.5) +
    geom_text(data = text_labels, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = 5, fontface = "bold")
}

# Display the plot
print(p_simpson)

# ----------------------------------------------------------------------------
# Save Figures
# ----------------------------------------------------------------------------

# PNG for Google Docs (high resolution)
ggsave(file.path(output_dir, "simpson_diversity_significance.png"), p_simpson, 
       width = 5, height = 5, dpi = 300, bg = "white")

# PDF for publication
ggsave(file.path(output_dir, "simpson_diversity_significance.pdf"), p_simpson, 
       width = 5, height = 5, device = cairo_pdf)

cat("\nFigures saved:\n")
cat("  -", file.path(output_dir, "simpson_diversity_significance.png"), "(300 dpi)\n")
cat("  -", file.path(output_dir, "simpson_diversity_significance.pdf"), "(vector)\n")

# ----------------------------------------------------------------------------
# Summary Statistics for Figure Legend
# ----------------------------------------------------------------------------

cat("\n=== Summary Statistics for Figure Legend ===\n\n")

summary_stats <- diversity_df %>%
  group_by(Treatment) %>%
  summarise(
    n = n(),
    Mean = mean(Simpson),
    SD = sd(Simpson),
    Median = median(Simpson),
    .groups = 'drop'
  )

print(summary_stats, digits = 3)

cat("\nSuggested figure legend text:\n")
cat("--------------------------------------------------------------------------------\n")
cat("Simpson diversity index (1-D) of miRNA cargo in HMSC-derived extracellular\n")
cat("vesicles by treatment condition. Boxplots show median (line), interquartile\n")
cat("range (box), and 1.5×IQR (whiskers). White diamonds indicate group means.\n")
cat("Individual samples shown as points. **p < 0.01 (Tukey HSD post-hoc test\n")
cat("following one-way ANOVA, p = 0.0008).\n")
cat("--------------------------------------------------------------------------------\n")
