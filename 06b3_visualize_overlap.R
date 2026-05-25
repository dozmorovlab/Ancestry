source("/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/06b0_sources.R")

rna_tools_paths <- list.files(rna_tools_dir, full.names = TRUE)
wgs_tools_paths <- list.files(wgs_tools_dir, full.names = TRUE)
rna_papers_paths <- list.files(rna_papers_dir, full.names = TRUE)
wgs_papers_paths <- list.files(wgs_papers_dir, full.names = TRUE)
model_tools_paths <- list.files(model_tools_dir, full.names = TRUE)
model_papers_paths <- list.files(model_papers_dir, full.names = TRUE)

process_samples <- function(in_dt, sample_name = "Sample") {

  sample_map <- c(
    "VCU-CC-121" = "VCU-OC-121",
    "VCU-CC-122" = "VCU-CO-122",
    "VCU-CC-126" = "VCU-CO-126",
    "VCU-CC-146" = "VCU-CO-146",
    "VCU-CC-147" = "VCU-CO-147",
    "VCU-CO-097" = "VCU-LC-097",
    "VCU-LC-099" = "VCU-CO-099"
  )

  exclude_list <- c(
    "VCU-OC-113",
    "VCU-PC-124",
    "VCU-PC-127",
    "VCU-CO-078"
  )

  in_dt <- as.data.table(in_dt)

  # Sample adjustment
  in_dt[[sample_name]] <- ifelse(
    in_dt[[sample_name]] %in% names(sample_map),
    sample_map[in_dt[[sample_name]]],
    in_dt[[sample_name]]
  )
  in_dt <- in_dt[! get(sample_name) %in% exclude_list]
  in_dt
}

paths <- c(
  rna_tools_paths,
  wgs_tools_paths,
  rna_papers_paths,
  wgs_papers_paths
)
model_paths <- c(
  model_tools_paths,
  model_papers_paths
)
tools <- gsub(".*\\.", "", gsub("\\.aims.bed?", "", paths))
samples <- gsub("\\..*", "", gsub("_.*", "", basename(paths)))
models <- gsub("\\.bed", "", basename(model_paths))
bed_dt <- data.table(
  sample = c(
    samples,
    rep_len("models", length(model_paths))
  ),
  tool = c(
    tools,
    models
  ),
  seq_type = c(
    rep_len("RNA", length(rna_tools_paths)),
    rep_len("WGS", length(wgs_tools_paths)),
    rep_len("RNA", length(rna_papers_paths)),
    rep_len("WGS", length(wgs_papers_paths)),
    rep_len("models", length(model_paths))
  ),
  aim_type = c(
    rep_len("tools", length(rna_tools_paths)),
    rep_len("tools", length(wgs_tools_paths)),
    rep_len("papers", length(rna_papers_paths)),
    rep_len("papers", length(wgs_papers_paths)),
    rep_len("tools", length(model_tools_paths)),
    rep_len("papers", length(model_papers_paths))
  ),
  path = c(paths, model_paths)
)

bed_dt <- process_samples(bed_dt, sample_name = "sample")

# bed_dt[, n_snps := vapply(path, function(p) {
#   nrow(tryCatch(fread(p), error = function(e) data.table()))
# }, integer(1))]

plan(multisession, workers = max(1, parallel::detectCores() - 1))

bed_dt[, n_snps := future_vapply(path, function(p) {
  tryCatch(
    countLines(p),
    error = function(e) 0L
  )
}, integer(1))]

#setnames(out_dt, c("tool"), c("AIMs Source"))

rna_tools <- bed_dt[seq_type == "RNA" & aim_type == "tools"]
wgs_tools <- bed_dt[seq_type == "WGS" & aim_type == "tools"]
rna_papers <- bed_dt[seq_type == "RNA" & aim_type == "papers"]
wgs_papers <- bed_dt[seq_type == "WGS" & aim_type == "papers"]
model_tools <- bed_dt[seq_type == "models" & aim_type == "tools"]
model_papers <- bed_dt[seq_type == "models" & aim_type == "papers"]

rna_tools[, n_snps_percent := n_snps / setNames(model_tools$n_snps, model_tools$tool)[tool]]
wgs_tools[, n_snps_percent := n_snps / setNames(model_tools$n_snps, model_tools$tool)[tool]]
rna_papers[, n_snps_percent := n_snps / setNames(model_papers$n_snps, model_papers$tool)[tool]]
wgs_papers[, n_snps_percent := n_snps / setNames(model_papers$n_snps, model_papers$tool)[tool]]

# Rename to official
rna_tools[, tool := official_names[tool]]
wgs_tools[, tool := official_names[tool]]
rna_papers[, tool := official_names[tool]]
wgs_papers[, tool := official_names[tool]]

rna_tools[, color := main_colors[tool]]
wgs_tools[, color := main_colors[tool]]
rna_papers[, color := main_colors[tool]]
wgs_papers[, color := main_colors[tool]]

rna_boxplot_a <- ggplot(
  rna_tools,
  aes(
    x = reorder(tool, n_snps, FUN = median, decreasing = TRUE),
    y = n_snps,
    fill = color
  )
) +
  geom_boxplot() +
  labs(
    title = "RNA-Seq",
    subtitle = NULL,
    x = NULL,
    y = "Number of sample SNPs\noverlaping model SNPs"
  ) +
  theme_minimal() +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_identity() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

rna_boxplot_b <- ggplot(
  rna_papers,
  aes(
    x = reorder(tool, n_snps, FUN = median, decreasing = TRUE),
    y = n_snps,
    fill = color
  )
) +
  geom_boxplot() +
  labs(
    title = "RNA-Seq",
    subtitle = NULL,
    x = NULL,
    y = "Number of sample SNPs\noverlaping model SNPs"
  ) +
  theme_minimal() +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_identity() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

wgs_boxplot_a <- ggplot(
  wgs_tools,
  aes(
    x = reorder(tool, n_snps, FUN = median, decreasing = TRUE),
    y = n_snps,
    fill = color
  )
) +
  geom_boxplot() +
  labs(
    title = "WGS",
    subtitle = NULL,
    x = NULL,
    y = "Number of sample SNPs overlaping model SNPs"
  ) +
  theme_minimal() +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_identity() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

wgs_boxplot_b <- ggplot(
  wgs_papers,
  aes(
    x = reorder(tool, n_snps, FUN = median, decreasing = TRUE),
    y = n_snps,
    fill = color
  )
) +
  geom_boxplot() +
  labs(
    title = "WGS",
    subtitle = NULL,
    x = NULL,
    y = "Number of sample SNPs overlaping model SNPs"
  ) +
  theme_minimal() +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_identity() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

# Percents

rna_boxplot_a_p <- ggplot(
  rna_tools,
  aes(
    x = reorder(tool, n_snps_percent, FUN = median, decreasing = TRUE),
    y = n_snps_percent,
    fill = color
  )
) +
  geom_boxplot() +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "Percent coverage of model SNPs"
  ) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_identity() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

rna_boxplot_b_p <- ggplot(
  rna_papers,
  aes(
    x = reorder(tool, n_snps_percent, FUN = median, decreasing = TRUE),
    y = n_snps_percent,
    fill = color
  )
) +
  geom_boxplot() +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "Percent coverage of model SNPs"
  ) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_identity() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

wgs_boxplot_a_p <- ggplot(
  wgs_tools,
  aes(
    x = reorder(tool, n_snps_percent, FUN = median, decreasing = TRUE),
    y = n_snps_percent,
    fill = color
  )
) +
  geom_boxplot() +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "Percent coverage of model SNPs"
  ) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_identity() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

wgs_boxplot_b_p <- ggplot(
  wgs_papers,
  aes(
    x = reorder(tool, n_snps_percent, FUN = median, decreasing = TRUE),
    y = n_snps_percent,
    fill = color
  )
) +
  geom_boxplot() +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "Percent coverage of model SNPs"
  ) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_identity() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

################################################

tools_plot <- rna_boxplot_a + wgs_boxplot_a
papers_plot <- rna_boxplot_b + wgs_boxplot_b

tools_plot_p <- rna_boxplot_a_p + wgs_boxplot_a_p
papers_plot_p <- rna_boxplot_b_p + wgs_boxplot_b_p

# ggsave(
#   file.path(out_dir, "tools_overlap.png"),
#   tools_plot,
#   width = 6,
#   height = 3,
#   dpi = 600,
#   bg = "white",
#   scale = 2
# )
# ggsave(
#   file.path(out_dir, "papers_overlap.png"),
#   papers_plot,
#   width = 6,
#   height = 3,
#   dpi = 600,
#   bg = "white",
#   scale = 2
# )

# ggsave(
#   file.path(out_dir, "tools_overlap_p.png"),
#   tools_plot_p,
#   width = 6,
#   height = 3,
#   dpi = 600,
#   bg = "white",
#   scale = 2
# )
# ggsave(
#   file.path(out_dir, "papers_overlap_p.png"),
#   papers_plot_p,
#   width = 6,
#   height = 3,
#   dpi = 600,
#   bg = "white",
#   scale = 2
# )

out_dt <- bed_dt[, .(
  `N Samples` = .N,
  Min         = as.numeric(min(n_snps, na.rm = TRUE)),
  Q1          = as.numeric(quantile(n_snps, 0.25, na.rm = TRUE)),
  Median      = as.numeric(median(n_snps, na.rm = TRUE)),
  Mean        = as.numeric(mean(n_snps, na.rm = TRUE)),
  Q3          = as.numeric(quantile(n_snps, 0.75, na.rm = TRUE)),
  Max         = as.numeric(max(n_snps, na.rm = TRUE))
), by = .(tool, seq_type)]

write.csv(
  out_dt,
  file.path(out_dir, "overlap_table.csv")
)

##########################################################################


