print(paste(date(), "Starting script"))
source("/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/06b0_sources.R")

save_plot <- function(out_path, in_plot, in_scale = 1) {

  ggsave(
    out_path,
    in_plot,
    bg = "white",
    scale = in_scale,
    width = 178,
    height = 178 / 2,
    units = "mm"
  )

}

break_plot <- function(
  in_plot,
  break_point = 6000,
  title_ = NULL,
  tag_ = NULL,
  x_ = NULL,
  y_ = NULL
) {

  lower_plot <- in_plot +
    coord_cartesian(
      ylim = c(0, break_point),
      expand = FALSE
    ) +
    labs(
      title = NULL,
      x = x_,
      y = y_
    )

  # Second figure - setting y-axis to a higher range
  upper_plot <- in_plot +
    coord_cartesian(
      ylim = c(break_point, NA),
      expand = FALSE
    ) +
    guides(x = "none") +  # Remove x-axis
    labs(
      title = title_,
      tag = tag_,
      x = NULL,
      y = NULL
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_blank(),
      legend.position = "none",  # Remove legend
      panel.background = element_blank(),
      axis.line = element_blank()
    )

  # Combine the two plots using patchwork
  upper_plot / lower_plot
}
svg_to_patchwork <- function(cairo_svg_path) {
  # Parse the cairo svg into a Picture object
  pic <- grImport2::readPicture(cairo_svg_path)

  # Convert Picture -> grob
  g <- grImport2::pictureGrob(pic)
  patchwork::wrap_elements(full = g, clip = FALSE)
}


##########################################################

current_dir <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/"

print(paste(date(), "Processing: 06b1_f1_acc.R"))
source(file.path(current_dir, "06b1_f1_acc.R"))
gc()

print(paste(date(), "Processing: 06b2_GATK_vs_JAX.R"))
source(file.path(current_dir, "06b2_GATK_vs_JAX.R"))
gc()

print(paste(date(), "Processing: 06b3_visualize_overlap.R"))
source(file.path(current_dir, "06b3_visualize_overlap.R"))
gc()

print(paste(date(), "Processing: 06b4_visualize_snps_reads.R"))
source(file.path(current_dir, "06b4_visualize_snps_reads.R"))
gc()

print(paste(date(), "Processing: 06b5_visualize_models.R"))
source(file.path(current_dir, "06b5_visualize_models.R"))
gc()

print(paste(date(), "Processing: 06b6_resources_plots.R"))
source(file.path(current_dir, "06b6_resources_plots.R"))
gc()

print(paste(date(), "Processing: 06b7_get_pca.R"))
source(file.path(current_dir, "06b7_get_pca.R"))
gc()

##########################################################

# Figure 1

# local / manual
fig1a <- svg_to_patchwork(fig1a_path) +
  labs(tag = "A")

# 06b7_get_pca.R
fig1b <- wgs_pca_plot +
  labs(tag = "B")

# Combine with fig1b in patchwork
fig1 <- wrap_plots(
  fig1a,
  fig1b,
  widths = c(3, 1)
)

#########################################################

# Figure 2
# 06b2_GATK_vs_JAX.R
fig2a <- wrap_plots(
  rna_unfiltered_venn +
    labs(tag = "A"),
  rna_filtered_venn,
  wgs_unfiltered_venn,
  wgs_filtered_venn,
  ncol = 1,
  nrow = 4,
  widths = c(1, 1),
  heights = c(1, 1)
)
# 06b2_GATK_vs_JAX.R
fig2b <- wrap_plots(
  gatk_rna_boxplot +
    labs(tag = "B"),
  jax_rna_boxplot,
  gatk_wgs_boxplot,
  jax_wgs_boxplot,
  ncol = 2,
  widths = c(1, 1),
  heights = c(1, 1)
)
# 06b5_visualize_models.R
fig2c <- break_plot(
  model_bar_plots$models,
  break_point = 300000,
  title_ = "Models",
  y_ = "Total number of model SNPs"
)
fig2c[[1]] <- fig2c[[1]] + labs(tag = "C")

# 06b3_visualize_overlap
fig2d <- wrap_plots(
  rna_boxplot_a +
    labs(tag = "D"),
  rna_boxplot_a_p,
  wgs_boxplot_a,
  wgs_boxplot_a_p,
  ncol = 2,
  widths = c(1, 1),
  heights = c(1, 1)
)

fig2 <- wrap_plots(
  free(fig2a), fig2b, fig2c, fig2d,
  nrow = 1,
  widths = c(1, 1, 1, 2)
)

#########################################################

# Figure 3

# 06b1_f1_acc.R
out_f1_plot[[1]] <- out_f1_plot[[1]] + labs(tag = "A")
fig3a <- out_f1_plot

# 06b6_resources_plots.R
resources_out_plot[[1]] <- resources_out_plot[[1]] + labs(tag = "B")
fig3b <- resources_out_plot

fig3 <- wrap_plots(
  fig3a,
  fig3b
)

#########################################################

# Figure S1

# 06b4_visualize_snps_reads.R
out_by_reads[[1]][[1]] <- out_by_reads[[1]][[1]] + labs(tag = "A")
figS1a <- out_by_reads

# 06b4_visualize_snps_reads.R
rna_wgs_scatterplot[[1]] <- rna_wgs_scatterplot[[1]] + labs(tag = "B")
figS1b <- rna_wgs_scatterplot

figS1 <- wrap_plots(
  figS1a,
  figS1b,
  widths = c(2, 1)
)

#########################################################

# Figure S2

# 
figS2a <- paper_euler

# 
figS2b <- 

figS2 <- wrap_plots(
  figS2a,
  figS2b,
  widths = c(2, 1)
)

#########################################################

# Save

# Figure 1
save_plot(
  file.path(out_dir, "fig1.png"),
  fig1,
  in_scale = 3
)
save_plot(
  file.path(out_dir, "fig1.svg"),
  fig1,
  in_scale = 3
)

# Figure 2
save_plot(
  file.path(out_dir, "fig2.png"),
  fig2,
  in_scale = 2.2
)
save_plot(
  file.path(out_dir, "fig2.svg"),
  fig2,
  in_scale = 2.2
)

# Figure 3
save_plot(
  file.path(out_dir, "fig3.png"),
  fig3,
  in_scale = 2.2
)
save_plot(
  file.path(out_dir, "fig3.svg"),
  fig3,
  in_scale = 2.2
)

# Figure S1
save_plot(
  file.path(out_dir, "figS1.png"),
  figS1,
  in_scale = 4
)
save_plot(
  file.path(out_dir, "figS1.svg"),
  figS1,
  in_scale = 4
)

# Figure S2
save_plot(
  file.path(out_dir, "figS2a.png"),
  figS2a,
  in_scale = 4
)
save_plot(
  file.path(out_dir, "figS2.svg"),
  figS2,
  in_scale = 4
)
