source("/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/06b0_sources.R")

# /lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/resources/get_resources.sh
usage_dt <- fread(
  usage_path,
  colClasses = list(character = c("Start", "End"))
)

end_dates <- c(
  # Admixture_wgs_old = "2025-08-07T14:51:36",
  # Admixture_rna_old = "2025-08-07T08:32:17",

  Admixture_wgs = "2025-09-15T11:25:48",
  Admixture_rna = "2025-09-15T11:25:42",

  AEon_rna = "2025-07-31T10:21:47",
  AEon_wgs = "2025-07-29T06:44:22",

  EthSEQ_rna = "2025-07-30T10:13:56",
  EthSEQ_wgs = "2025-07-29T06:44:36",

  gnomAD_rna = "2025-07-30T14:22:06",
  gnomAD_wgs = "2025-07-29T06:44:42",

  `JAX_rna` = "2025-07-30T14:21:52",
  `JAX_wgs` = "2025-07-19T11:42:12",

  RAIDS_rna = "2025-08-05T06:31:05",
  RAIDS_wgs_1 = "2025-08-08T09:38:50"
)

# Test
# usage_dt[
#   (seq_len(.N) %in% grep("2025-07", End)) &
#     (seq_len(.N) %in% grep("EthSEQ", JobName)) &
#     (seq_len(.N) %in% grep("UNLIMITED", Timelimit, invert = TRUE))
# ][, .(Start, End, MaxVMSize, MaxDiskWrite, MaxDiskRead)]

select_dt <- usage_dt[
  End %in% end_dates &
    (seq_len(.N) %in% grep("UNLIMITED", Timelimit, invert = TRUE))
]
select_dt[, Name := names(end_dates)]
select_dt[, Tool := gsub("_.*", "", Name)]
select_dt[, Type := sub(".*(rna|wgs).*", "\\1", Name)]

plot_dt <- select_dt[, .(Name, Tool, Type, MaxRSS, CPUTimeRAW)]
plot_dt[, MaxRSS_GB := as.integer(gsub("K", "", MaxRSS)) / 1024 / 1024]

# plot_dt <- merge(
#   plot_dt[seq_len(.N) %in% grep("rna", Name)],
#   plot_dt[seq_len(.N) %in% grep("wgs", Name)],
#   by = "Tool",
#   suffixes = c(" RNA", " WGS")
# )

rna_wgs_names <- c(
  rna = "RNA-Seq",
  wgs = "WGS"
)

plot_dt[, Tool := official_names[Tool]]
plot_dt[, color := main_colors[Tool]]
plot_dt[, Type := rna_wgs_names[Type]]

mem_plot <- ggplot(
  plot_dt,
  aes(
    x = interaction(Type, Tool, sep = ": "),
    y = `MaxRSS_GB`,
    group = Type,
    fill = color
  )
) +
  geom_col(position = "dodge", color = "black") +
  theme_minimal() +
  labs(
    x = NULL,
    y = "GB",
    title = "Max Memory Use"
  ) +
  scale_fill_identity() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  scale_y_continuous(labels = scales::comma)

cpu_plot <- ggplot(
  plot_dt,
  aes(
    x = interaction(Type, Tool, sep = ": "),
    y = `CPUTimeRAW` / 60 / 60,
    group = Type,
    fill = color
  )
) +
  geom_col(position = "dodge", color = "black") +
  theme_minimal() +
  labs(
    x = NULL,
    y = "hours",
    title = "CPU Time"
  ) +
  scale_fill_identity() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  scale_y_continuous(labels = scales::comma)

resources_out_plot <- mem_plot + cpu_plot

# ggsave(
#   file.path(out_dir, "resources.png"),
#   resources_out_plot,
#   bg = "white"
# )

# ggsave(
#   file.path(out_dir, "memory_use.png"),
#   mem_plot,
#   bg = "white"
# )
# TotalCPU
# MaxRSS
# MaxVMSize
# MaxDiskWrite
# ConsumedEnergy / ConsumedEnergyRaw
