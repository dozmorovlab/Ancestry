library(SNPRelate)
library(GenomicRanges)

options(scipen = 999)

out_dir <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/aims_tools/models"
dir.create(out_dir, recursive = TRUE)

gds <- snpgdsOpen(file.path(
  "/lustre/home/harrell_lab/bulkRNASeq/Ancestry/References/Population/matGeno1000g.gds"
))

chrom <- read.gdsn(index.gdsn(gds, "snp.chromosome"))
pos <- read.gdsn(index.gdsn(gds, "snp.position"))
snpid <- read.gdsn(index.gdsn(gds, "snp.id"))
snpgdsClose(gds)


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
  file = file.path(out_dir, "RAIDS.bed"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
