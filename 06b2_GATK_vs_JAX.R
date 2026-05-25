source("/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/06b0_sources.R")

rna_dt <- process_gatk_vs_jax(in_RNA_csv)
wgs_dt <- process_gatk_vs_jax(in_WGS_csv)
rna_unfiltered_dt <- process_gatk_vs_jax(in_RNA_unfiltered_csv)
wgs_unfiltered_dt <- process_gatk_vs_jax(in_WGS_unfiltered_csv)

########################################################
# RNA
boxplot_plot_rna <- ggplot(rna_dt, aes(x = Type, y = SNPs)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "RNA",
    y = "SNP count per sample",
    x = NULL
  )

# Calculate group means and standard errors
summary_dt_rna <- aggregate(SNPs ~ Type, data = rna_dt, function(x) {
  m <- mean(x, na.rm = TRUE)
  se <- sd(x, na.rm = TRUE) / sqrt(length(x))
  c(mean = m, se = se)
})

# 'SNPs' column is a matrix, separate it into two columns:
summary_dt_rna <- cbind(
  summary_dt_rna[, "Type", drop = FALSE],
  as.data.frame(summary_dt_rna$SNPs)
)


# Plot
err_plot_rna <- ggplot(summary_dt_rna, aes(x = Type, y = mean)) +
  geom_col() +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.2
  ) +
  theme_minimal() +
  labs(
    title = "RNA",
    y = "Mean SNPs +- standard error",
    x = NULL
  )

rna_plot <- err_plot_rna

# WGS
boxplot_plot_wgs <- ggplot(wgs_dt, aes(x = Type, y = SNPs)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "WGS",
    y = "SNP count per sample",
    x = NULL
  )

# Calculate group means and standard errors
summary_dt_wgs <- aggregate(SNPs ~ Type, data = wgs_dt, function(x) {
  m <- mean(x, na.rm = TRUE)
  se <- sd(x, na.rm = TRUE) / sqrt(length(x))
  c(mean = m, se = se)
})

summary_dt_wgs <- cbind(
  summary_dt_wgs[, "Type", drop = FALSE],
  as.data.frame(summary_dt_wgs$SNPs)
)

# Plot
err_plot_wgs <- ggplot(summary_dt_wgs, aes(x = Type, y = mean)) +
  geom_col() +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.2
  ) +
  theme_minimal() +
  labs(
    title = "WGS",
    y = "Mean SNPs +- standard error",
    x = NULL
  )

wgs_plot <- err_plot_wgs

##################################################################

# Summary

GATK_wgs_snps <- summary_dt_wgs[summary_dt_wgs$Type == "GATK Unique", "mean"] +
  summary_dt_wgs[summary_dt_wgs$Type == "Shared", "mean"]
JAX_wgs_snps <- summary_dt_wgs[summary_dt_wgs$Type == "JAX Unique", "mean"] +
  summary_dt_wgs[summary_dt_wgs$Type == "Shared", "mean"]
GATK_rna_snps <- summary_dt_rna[summary_dt_rna$Type == "GATK Unique", "mean"] +
  summary_dt_rna[summary_dt_rna$Type == "Shared", "mean"]
JAX_rna_snps <- summary_dt_rna[summary_dt_rna$Type == "JAX Unique", "mean"] +
  summary_dt_rna[summary_dt_rna$Type == "Shared", "mean"]

print(paste0("WGS vs RNA: ", GATK_wgs_snps / GATK_rna_snps, " times fewer SNPs in RNA than WGS"))
print(paste0("WGS GATK vs JAX: ", GATK_wgs_snps / JAX_wgs_snps, " times fewer SNPs in JAX than GATK"))
print(paste0("RNA GATK vs JAX: ", GATK_rna_snps / JAX_rna_snps, " times fewer SNPs in JAX than GATK"))

#####################################################################
# eulerr
fit_rna <- euler(c(
  GATK = round(
    summary_dt_rna[summary_dt_rna$Type == "GATK Unique", "mean"],
    1
  ),
  JAX = round(
    summary_dt_rna[summary_dt_rna$Type == "JAX Unique", "mean"],
    1
  ),
  "GATK&JAX" = round(
    summary_dt_rna[summary_dt_rna$Type == "Shared", "mean"],
    1
  )
))

fit_wgs <- euler(c(
  GATK = round(
    summary_dt_wgs[summary_dt_wgs$Type == "GATK Unique", "mean"],
    1
  ),
  JAX = round(
    summary_dt_wgs[summary_dt_wgs$Type == "JAX Unique", "mean"],
    1
  ),
  "GATK&JAX" = round(
    summary_dt_wgs[summary_dt_wgs$Type == "Shared", "mean"],
    1
  )
))

rna_e <- as.ggplot(
  ~eulerr:::plot.eulergram(
    plot(
      fit_rna,
      quantities = TRUE,
    )
  ),
  scale = 1
)

wgs_e <- as.ggplot(
  ~eulerr:::plot.eulergram(
    plot(
      fit_wgs,
      quantities = TRUE,
    )
  ),
  scale = 1
)

eulerr_plot <- rna_e + wgs_e

ggsave(
  file.path(out_dir, "GATK_vs_JAX.eulerr.png"),
  eulerr_plot,
  scale = 3,
  height = 4,
  width = 4
)

##################################################################

# Unfiltered

# RNA
boxplot_plot_rna_unfiltered <- ggplot(rna_unfiltered_dt, aes(x = Type, y = SNPs)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "rna_unfiltered",
    y = "SNP count per sample",
    x = NULL
  )

# Calculate group means and standard errors
summary_dt_rna_unfiltered <- aggregate(SNPs ~ Type, data = rna_unfiltered_dt, function(x) {
  m <- mean(x, na.rm = TRUE)
  se <- sd(x, na.rm = TRUE) / sqrt(length(x))
  c(mean = m, se = se)
})

# 'SNPs' column is a matrix, separate it into two columns:
summary_dt_rna_unfiltered <- cbind(
  summary_dt_rna_unfiltered[, "Type", drop = FALSE],
  as.data.frame(summary_dt_rna_unfiltered$SNPs)
)


# Plot
err_plot_rna_unfiltered <- ggplot(summary_dt_rna_unfiltered, aes(x = Type, y = mean)) +
  geom_col() +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.2
  ) +
  theme_minimal() +
  labs(
    title = "rna_unfiltered",
    y = "Mean SNPs +- standard error",
    x = NULL
  )

rna_unfiltered_plot <- err_plot_rna_unfiltered

# WGS
boxplot_plot_wgs_unfiltered <- ggplot(wgs_unfiltered_dt, aes(x = Type, y = SNPs)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "wgs_unfiltered",
    y = "SNP count per sample",
    x = NULL
  )

# Calculate group means and standard errors
summary_dt_wgs_unfiltered <- aggregate(SNPs ~ Type, data = wgs_unfiltered_dt, function(x) {
  m <- mean(x, na.rm = TRUE)
  se <- sd(x, na.rm = TRUE) / sqrt(length(x))
  c(mean = m, se = se)
})

summary_dt_wgs_unfiltered <- cbind(
  summary_dt_wgs_unfiltered[, "Type", drop = FALSE],
  as.data.frame(summary_dt_wgs_unfiltered$SNPs)
)

# Plot
err_plot_wgs_unfiltered <- ggplot(summary_dt_wgs_unfiltered, aes(x = Type, y = mean)) +
  geom_col() +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.2
  ) +
  theme_minimal() +
  labs(
    title = "wgs_unfiltered",
    y = "Mean SNPs +- standard error",
    x = NULL
  )

wgs_unfiltered_plot <- err_plot_wgs_unfiltered

# eulerr
fit_unfiltered_rna <- euler(c(
  GATK = round(
    summary_dt_rna_unfiltered[summary_dt_rna_unfiltered$Type == "GATK Unique", "mean"],
    1
  ),
  JAX = round(
    summary_dt_rna_unfiltered[summary_dt_rna_unfiltered$Type == "JAX Unique", "mean"],
    1
  ),
  "GATK&JAX" = round(
    summary_dt_rna_unfiltered[summary_dt_rna_unfiltered$Type == "Shared", "mean"],
    1
  )
))

fit_unfiltered_wgs <- euler(c(
  GATK = round(
    summary_dt_wgs_unfiltered[summary_dt_wgs_unfiltered$Type == "GATK Unique", "mean"],
    1
  ),
  JAX = round(
    summary_dt_wgs_unfiltered[summary_dt_wgs_unfiltered$Type == "JAX Unique", "mean"],
    1
  ),
  "GATK&JAX" = round(
    summary_dt_wgs_unfiltered[summary_dt_wgs_unfiltered$Type == "Shared", "mean"],
    1
  )
))

rna_unfiltered_e <- as.ggplot(
  ~eulerr:::plot.eulergram(
    plot(
      fit_unfiltered_rna,
      quantities = TRUE,
    )
  ),
  scale = 1
)

wgs_unfiltered_e <- as.ggplot(
  ~eulerr:::plot.eulergram(
    plot(
      fit_unfiltered_wgs,
      quantities = TRUE,
    )
  ),
  scale = 1
)

eulerr_all_plot <- (rna_e + wgs_e) / (rna_unfiltered_e + wgs_unfiltered_e)

ggsave(
  file.path(out_dir, "GATK_vs_JAX_unfiltered.eulerr.png"),
  eulerr_all_plot,
  scale = 3,
  height = 4,
  width = 4
)

###############

# Your numbers (can be decimals)
fmt_val <- function(df, type, in_title = "") {
  row <- df[df$Type == type, ]
  paste0(
    in_title, "\n",
    format(round(row$mean, 1), big.mark = ",", scientific = FALSE),
    "\n± ",
    format(round(row$se, 1), big.mark = ",", scientific = FALSE)
  )
}

make_venn <- function(
  a, b, ab,
  color_a = "#1a8ab8", color_b = "white", color_ab = NULL,
  in_title = NULL, text_size = 4
) {
  circles <- data.frame(
    x = c(1, 2),       # circle centers
    y = c(1, 1),
    r = c(1, 1),
    set = c("A", "B")
  )

  # Basic plot
  ggplot() +
    geom_circle(
      data = circles,
      aes(x0 = x, y0 = y, r = r, fill = set),
      alpha = 0.3,
      color = "black"
    ) +
    annotate("text", x = 0.5, y = 1, label = a, size = text_size) +
    annotate("text", x = 2.5, y = 1, label = b, size = text_size) +
    annotate("text", x = 1.5, y = 1, label = ab, size = text_size) +
    coord_fixed() +
    theme_void() +
    scale_fill_manual(values = c(color_a, color_ab, color_b)) +
    theme(legend.position = "none") +
    labs(title = in_title)
}

rna_unfiltered_venn  <- make_venn(
  a = fmt_val(summary_dt_rna_unfiltered, "GATK Unique", "GATK"),
  b = fmt_val(summary_dt_rna_unfiltered, "JAX Unique", "JAX"),
  ab = fmt_val(summary_dt_rna_unfiltered, "Shared"),
  in_title = "RNA-Seq\nUnfiltered",
  text_size = 3,
  color_a = main_colors[["RNA GATK Unfiltered"]],
  color_b = main_colors[["RNA JAX Unfiltered"]]
)
wgs_unfiltered_venn  <- make_venn(
  a = fmt_val(summary_dt_wgs_unfiltered, "GATK Unique", "GATK"),
  b = fmt_val(summary_dt_wgs_unfiltered, "JAX Unique", "JAX"),
  ab = fmt_val(summary_dt_wgs_unfiltered, "Shared"),
  in_title = "WGS\nUnfiltered",
  text_size = 3,
  color_a = main_colors[["WGS GATK Unfiltered"]],
  color_b = main_colors[["WGS JAX Unfiltered"]]
)
rna_filtered_venn  <- make_venn(
  a = fmt_val(summary_dt_rna, "GATK Unique", "GATK"),
  b = fmt_val(summary_dt_rna, "JAX Unique", "JAX"),
  ab = fmt_val(summary_dt_rna, "Shared"),
  in_title = "\nFiltered",
  text_size = 3,
  color_a = main_colors[["RNA GATK Filtered"]],
  color_b = main_colors[["RNA JAX Filtered"]]
)
wgs_filtered_venn  <- make_venn(
  a = fmt_val(summary_dt_wgs, "GATK Unique", "GATK"),
  b = fmt_val(summary_dt_wgs, "JAX Unique", "JAX"),
  ab = fmt_val(summary_dt_wgs, "Shared"),
  in_title = "\nFiltered",
  text_size = 3,
  color_a = main_colors[["WGS GATK Filtered"]],
  color_b = main_colors[["WGS JAX Filtered"]]
)

#########################################################################

summary_dt_rna$Overlap <- summary_dt_rna$Type
summary_dt_rna$Type <- "RNA"
summary_dt_rna_unfiltered$Overlap <- summary_dt_rna_unfiltered$Type
summary_dt_rna_unfiltered$Type <- "RNA unfiltered"
summary_dt_wgs$Overlap <- summary_dt_wgs$Type
summary_dt_wgs$Type <- "WGS"
summary_dt_wgs_unfiltered$Overlap <- summary_dt_wgs_unfiltered$Type
summary_dt_wgs_unfiltered$Type <- "WGS unfiltered"


summary_all <- rbindlist(list(
  summary_dt_rna, summary_dt_rna_unfiltered,
  summary_dt_wgs, summary_dt_wgs_unfiltered
))

write.table(
  summary_all,
  file.path(out_dir, "GATK_vs_JAX.summary.csv"),
  sep = ",",
  row.names = FALSE
)

##########################################################

rna_dt <- process_gatk_vs_jax(in_RNA_csv)
wgs_dt <- process_gatk_vs_jax(in_WGS_csv)

rna_unfiltered_dt <- process_gatk_vs_jax(in_RNA_unfiltered_csv)
wgs_unfiltered_dt <- process_gatk_vs_jax(in_WGS_unfiltered_csv)

rna_snps <- rna_dt[
  Type %in% c("Shared", "GATK Unique"),
  .(Total_SNPs = sum(SNPs)),
  by = Sample
]
rna_snps$Type <- "RNA GATK Filtered"
wgs_snps <- wgs_dt[
  Type %in% c("Shared", "GATK Unique"),
  .(Total_SNPs = sum(SNPs)),
  by = Sample
]
wgs_snps$Type <- "WGS GATK Filtered"
rna_unfiltered_snps <- rna_unfiltered_dt[
  Type %in% c("Shared", "GATK Unique"),
  .(Total_SNPs = sum(SNPs)),
  by = Sample
]
rna_unfiltered_snps$Type <- "RNA GATK Unfiltered"
wgs_unfiltered_snps <- wgs_unfiltered_dt[
  Type %in% c("Shared", "GATK Unique"),
  .(Total_SNPs = sum(SNPs)),
  by = Sample
]
wgs_unfiltered_snps$Type <- "WGS GATK Unfiltered"


rna_snps_jax <- rna_dt[
  Type %in% c("Shared", "JAX Unique"),
  .(Total_SNPs = sum(SNPs)),
  by = Sample
]
rna_snps_jax$Type <- "RNA JAX Filtered"
wgs_snps_jax <- wgs_dt[
  Type %in% c("Shared", "JAX Unique"),
  .(Total_SNPs = sum(SNPs)),
  by = Sample
]
wgs_snps_jax$Type <- "WGS JAX Filtered"
rna_unfiltered_snps_jax <- rna_unfiltered_dt[
  Type %in% c("Shared", "JAX Unique"),
  .(Total_SNPs = sum(SNPs)),
  by = Sample
]
rna_unfiltered_snps_jax$Type <- "RNA JAX Unfiltered"
wgs_unfiltered_snps_jax <- wgs_unfiltered_dt[
  Type %in% c("Shared", "JAX Unique"),
  .(Total_SNPs = sum(SNPs)),
  by = Sample
]
wgs_unfiltered_snps_jax$Type <- "WGS JAX Unfiltered"


########################################################################
snps_p_dt <- rbindlist(list(
  rna_snps, rna_unfiltered_snps, wgs_snps, wgs_unfiltered_snps,
  rna_snps_jax, rna_unfiltered_snps_jax, wgs_snps_jax, wgs_unfiltered_snps_jax
))
snps_p_dt[, Type := factor(Type, levels = c(
  "RNA GATK Unfiltered", "RNA GATK Filtered", "RNA JAX Unfiltered", "RNA JAX Filtered",
  "WGS GATK Unfiltered", "WGS GATK Filtered", "WGS JAX Unfiltered", "WGS JAX Filtered"
))]

snps_p_dt[, color := main_colors[Type]]

x_dict <- c(
  "RNA GATK Unfiltered" = "Unfiltered",
  "RNA GATK Filtered" = "Filtered",
  "RNA JAX Unfiltered" = "Unfiltered",
  "RNA JAX Filtered" = "Filtered",
  "WGS GATK Unfiltered" = "Unfiltered",
  "WGS GATK Filtered" = "Filtered",
  "WGS JAX Unfiltered" = "Unfiltered",
  "WGS JAX Filtered" = "Filtered"
)
x_order <- c("Unfiltered", "Filtered")

gatk_rna_boxplot <- ggplot(
  data = snps_p_dt[Type %in% c("RNA GATK Unfiltered", "RNA GATK Filtered")],
  mapping = aes(
    x = factor(x_dict[Type], levels = x_order),
    y = Total_SNPs,
    fill = color
  )
) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "RNA-Seq\nGATK",
    x = NULL,
    y = "Number of sample SNPs"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  scale_fill_identity() +
  scale_y_continuous(labels = scales::comma)

gatk_wgs_boxplot <- ggplot(
  data = snps_p_dt[Type %in% c("WGS GATK Unfiltered", "WGS GATK Filtered")],
  mapping = aes(
    x = factor(x_dict[Type], levels = x_order),
    y = Total_SNPs,
    fill = color
  )
) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "WGS\nGATK",
    x = NULL,
    y = "Number of sample SNPs"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  scale_fill_identity() +
  scale_y_continuous(labels = scales::comma)

jax_rna_boxplot <- ggplot(
  data = snps_p_dt[Type %in% c("RNA JAX Unfiltered", "RNA JAX Filtered")],
  mapping = aes(
    x = factor(x_dict[Type], levels = x_order),
    y = Total_SNPs,
    fill = color
  )
) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "\nJAX",
    x = NULL,
    y = "Number of sample SNPs"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  scale_fill_identity() +
  scale_y_continuous(labels = scales::comma)

jax_wgs_boxplot <- ggplot(
  data = snps_p_dt[Type %in% c("WGS JAX Unfiltered", "WGS JAX Filtered")],
  mapping = aes(
    x = factor(x_dict[Type], levels = x_order),
    y = Total_SNPs,
    fill = color
  )
) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "\nJAX",
    x = NULL,
    y = "Number of sample SNPs"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  scale_fill_identity() +
  scale_y_continuous(labels = scales::comma)

all_boxplot_sep <-  wrap_plots(list(
  gatk_rna_boxplot, gatk_wgs_boxplot, jax_rna_boxplot, jax_wgs_boxplot
))

all_boxplot_lin <- ggplot(
  data = snps_p_dt,
  mapping = aes(
    x = factor(x_dict[Type], levels = x_order),
    y = Total_SNPs,
    fill = color
  )
) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "SNPs",
    x = NULL,
    y = "Number of sample SNPs"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  scale_fill_identity() +
  scale_y_continuous(labels = scales::comma)

all_boxplot_log <- ggplot(
  data = snps_p_dt,
  mapping = aes(
    x = factor(x_dict[Type], levels = x_order),
    y = Total_SNPs,
    fill = color
  )
) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "SNPs",
    x = NULL,
    y = "Number of sample SNPs"
  ) +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  scale_fill_identity() +
  scale_y_continuous(labels = scales::comma)

# ggsave(
#   file.path(out_dir, "GATK_vs_JAX_snps_sep.png"),
#   all_boxplot_sep,
#   bg = "white",
#   scale = 2,
#   height = 4,
#   width = 4
# )
# ggsave(
#   file.path(out_dir, "GATK_vs_JAX_snps_lin.png"),
#   all_boxplot_lin,
#   bg = "white",
#   scale = 2,
#   height = 4,
#   width = 4
# )
# ggsave(
#   file.path(out_dir, "GATK_vs_JAX_snps_log.png"),
#   all_boxplot_log,
#   bg = "white",
#   scale = 2,
#   height = 4,
#   width = 4
# )
