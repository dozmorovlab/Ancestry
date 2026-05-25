library(GenomicRanges)

options(scipen = 999)

snpwt_path <- "/lustre/home/harrell_lab/JAX/ref_data/ancestry_panel_v2.snpwt"

out_dir <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/aims_tools/models"
dir.create(out_dir, recursive = TRUE)

snpwt_df <- read.delim(
  snpwt_path,
  sep = " ",
  skip = 5,
  header = FALSE
)

snpid <- snpwt_df[, 1]
chrom <- sapply(strsplit(snpid, ":"), `[`, 1)
pos <- as.numeric(sapply(strsplit(snpid, ":"), `[`, 2))

gr <- GRanges(
  seqnames = chrom,
  ranges = IRanges(start = pos, end = pos),
  snp.id = snpid
)
seqlevelsStyle(gr) <- "UCSC"

out_df <- as.data.frame(gr)[c("seqnames", "start", "end", "snp.id")]
out_df$start <- out_df$start - 1

write.table(
  out_df,
  file = file.path(out_dir, "JAX_SNPWeights.bed"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
