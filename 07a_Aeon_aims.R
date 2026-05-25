options(scipen = 999)

out_dir <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/aims_tools/models"
dir.create(out_dir, recursive = TRUE)

# URL
url <- "https://raw.githubusercontent.com/GenomicRisk/aeon/main/aeon_ancestry/refs/g1k_allele_freqs.txt"

g1k_freqs <- read.delim(url, sep = "\t", header = TRUE)

write.table(
  g1k_freqs[, c("CHROM", "START", "STOP", "VAR_ID")],
  file = file.path(out_dir, "Aeon.bed"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
