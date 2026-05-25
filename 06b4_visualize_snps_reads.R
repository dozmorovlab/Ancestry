source("/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/06b0_sources.R")

rna_reads <- fread(rna_reads_dir)
wgs_reads <- fread(wgs_reads_dir)

rna_reads$Sample <- gsub("(_|.merged).*", "", basename(rna_reads$file))
wgs_reads$Sample <- gsub("(_|.merged).*", "", basename(wgs_reads$file))

rna_reads <- process_samples(rna_reads)
wgs_reads <- process_samples(wgs_reads)

rna_dt <- process_gatk_vs_jax(in_RNA_csv)
wgs_dt <- process_gatk_vs_jax(in_WGS_csv)

rna_snps <- rna_dt[
  Type %in% c("Shared", "GATK Unique"),
  .(Total_SNPs = sum(SNPs)),
  by = Sample
]
wgs_snps <- wgs_dt[
  Type %in% c("Shared", "GATK Unique"),
  .(Total_SNPs = sum(SNPs)),
  by = Sample
]

rna_dt <- merge(rna_snps, rna_reads, by = "Sample", all = TRUE)
wgs_dt <- merge(wgs_snps, wgs_reads, by = "Sample", all = TRUE)

rna_scatterplot <- ggscatter(
  rna_dt,
  x = "total_reads",
  y = "Total_SNPs",
  add = "reg.line",
  conf.int = TRUE,
  cor.coef = TRUE,
  cor.method = "pearson",
  title = "RNA-Seq",
  subtitle = "Pearson"
) +
  labs(
    x = "Reads",
    y = "SNPs"
  ) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma)


wgs_scatterplot <- ggscatter(
  wgs_dt,
  x = "total_reads",
  y = "Total_SNPs",
  add = "reg.line",
  conf.int = TRUE,
  cor.coef = TRUE,
  cor.method = "pearson",
  title = "WGS",
  subtitle = "Pearson"
) +
  labs(
    x = "Reads",
    y = "SNPs"
  ) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma)


rna_wgs_scatterplot <- rna_scatterplot / wgs_scatterplot

# ggsave(
#   file.path(out_dir, "pearson_scatter.png"),
#   rna_wgs_scatterplot,
#   width = 3,
#   height = 6,
#   dpi = 600,
#   scale = 2
# )
#######################################################
rna_snps_plot <- ggplot(
  rna_snps,
  aes(
    y = Total_SNPs,
    x = reorder(Sample, Total_SNPs)
  )
) +
  geom_col() +
  theme_minimal() +
  labs(
    title = "RNA-Seq",
    subtitle = "SNPs",
    x = NULL,
    y = NULL
  ) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma)

wgs_snps_plot <- ggplot(
  wgs_snps,
  aes(
    y = Total_SNPs,
    x = reorder(Sample, Total_SNPs)
  )
) +
  geom_col() +
  theme_minimal() +
  labs(
    title = "WGS",
    subtitle = "SNPs",
    x = NULL,
    y = NULL
  ) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma)

snps_sample_plot <- rna_snps_plot | wgs_snps_plot

rna_snps$Type <- "RNA"
wgs_snps$Type <- "WGS"

snps_df <- rbindlist(list(rna_snps, wgs_snps))

snps_box_plot <- ggplot(
  snps_df,
  aes(
    y = Total_SNPs,
    x = factor(Type, levels = c("WGS", "RNA"))
  )
) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    x = NULL,
    y = NULL
  ) +
  scale_y_continuous(labels = scales::comma)

# ggsave(
#   file.path(out_dir, "snps.png"),
#   snps_sample_plot,
#   width = 3,
#   height = 3,
#   dpi = 600,
#   scale = 3
# )

#######################################################################

rna_snps_plot_by_reads <- ggplot(
  rna_dt,
  aes(
    y = Total_SNPs,
    x = reorder(Sample, total_reads)
  )
) +
  geom_col() +
  theme_minimal() +
  labs(
    title = "RNA-Seq",
    subtitle = "SNPs",
    x = NULL,
    y = NULL
  ) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma)

wgs_snps_plot_by_reads <- ggplot(
  wgs_dt,
  aes(
    y = Total_SNPs,
    x = reorder(Sample, total_reads)
  )
) +
  geom_col() +
  theme_minimal() +
  labs(
    title = "WGS",
    subtitle = "SNPs",
    x = NULL,
    y = NULL
  ) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma)

#########
rna_reads_plot <- ggplot(
  rna_dt,
  aes(
    y = total_reads,
    x = reorder(Sample, total_reads)
  )
) +
  geom_col() +
  theme_minimal() +
  labs(
    title = "RNA-Seq",
    subtitle = "Reads",
    x = NULL,
    y = NULL
  ) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma)

wgs_reads_plot <- ggplot(
  wgs_dt,
  aes(
    y = total_reads,
    x = reorder(Sample, total_reads)
  )
) +
  geom_col() +
  theme_minimal() +
  labs(
    title = "WGS",
    subtitle = "Reads",
    x = NULL,
    y = NULL
  ) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma)
#########
rna_reads_plot_by_snps <- ggplot(
  rna_dt,
  aes(
    y = total_reads,
    x = reorder(Sample, Total_SNPs)
  )
) +
  geom_col() +
  theme_minimal() +
  labs(
    title = "RNA-Seq",
    subtitle = "Reads",
    x = NULL,
    y = NULL
  ) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma)

wgs_reads_plot_by_snps <- ggplot(
  wgs_dt,
  aes(
    y = total_reads,
    x = reorder(Sample, Total_SNPs)
  )
) +
  geom_col() +
  theme_minimal() +
  labs(
    title = "WGS",
    subtitle = "Reads",
    x = NULL,
    y = NULL
  ) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma)


rna_out_by_snps <- rna_snps_plot + rna_reads_plot_by_snps
wgs_out_by_snps <- wgs_snps_plot + wgs_reads_plot_by_snps
rna_out_by_reads <- rna_reads_plot + rna_snps_plot_by_reads
wgs_out_by_reads <- wgs_reads_plot + wgs_snps_plot_by_reads

out_by_snps <- rna_out_by_snps / wgs_out_by_snps
out_by_reads <- rna_out_by_reads / wgs_out_by_reads

# ggsave(
#   file.path(out_dir, "snps_by_reads.png"),
#   out_by_reads,
#   width = 3,
#   height = 3,
#   dpi = 600,
#   scale = 5
# )

# ggsave(
#   file.path(out_dir, "reads_by_snps.png"),
#   out_by_snps,
#   width = 3,
#   height = 3,
#   dpi = 600,
#   scale = 5
# )
