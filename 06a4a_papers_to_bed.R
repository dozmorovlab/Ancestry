library(data.table)
library(openxlsx)
library(httr)
library(ggplot2)
library(patchwork)
library(biomaRt)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(ggupset)

out_bed_dir <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/aims_tools/papers/bed"
dir.create(out_bed_dir)

# Kosoy R, Nassir R, Tian C, White PA, Butler LM, Silva G, et al. Ancestry informative marker sets for determining continental origin and admixture proportions in common populations in America.
kosoy_path <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/aims_tools/papers/humu_20822_sm_tables2.txt"

# Nassir R, Kosoy R, Tian C, White PA, Butler LM, Silva G, et al. An ancestry informative marker set for determining continental origin: validation and extension using human genome diversity panels.
nassir_path <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/aims_tools/papers/12863_2008_697_MOESM1_ESM.xls"

# Torres JB, Stone AC, Kittles R. An anthropological genetic perspective on Creolization in the Anglophone Caribbean.
torres_path <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/aims_tools/papers/ajpa22261-sup-0001-suppinfo1.xlsx"

kosoy_dt <- fread(kosoy_path, skip = 1)
nassir_dt <- as.data.table(readxl::read_xls(nassir_path, skip = 1))
torres_dt <- as.data.table(readxl::read_xlsx(torres_path, skip = 3))

rs_dt_long <- data.table(
  name = c(
    rep_len("Kosoy", length(kosoy_dt$`NCBI SNP Reference`)),
    rep_len("Nassir", length(nassir_dt$NCBI)),
    rep_len("Torres", length(torres_dt$`SNP rs`))
  ),
  rs = c(
    kosoy_dt$`NCBI SNP Reference`,
    nassir_dt$NCBI,
    torres_dt$`SNP rs`
  )
)

rs_dt_wide <- dcast(
  unique(rs_dt_long[, .(rs, name)]),  # remove potential duplicates
  rs ~ name,
  fun.aggregate = function(x) as.logical(length(x)),
  value.var = "name"
)

#######################################################

rs_list <- unique(rs_dt_wide$rs)

snp_mart <- useEnsembl(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp", mirror = "useast")
snp_info <- getBM(
  attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end', 'allele', 'ensembl_gene_stable_id'),
  filters = 'snp_filter',
  values = rs_list,
  mart = snp_mart
)

##################################################

gr <- GRanges(
  seqnames = snp_info$chr_name,
  ranges = IRanges(start = snp_info$chrom_start, width = 1),
  names = snp_info$refsnp_id
)
seqlevelsStyle(gr) <- "UCSC"

gr_dt <- as.data.table(gr)
gr_dt$start <- gr_dt$start - 1
gr_dt <- gr_dt[, .(seqnames, start, end, names)]

out_paper_beds <- sapply(unique(rs_dt_long$name), function(paper) {

  rs <- rs_dt_long[name == paper]$rs

  dt <- gr_dt[names %in% rs]

  order_key <- suppressWarnings(as.numeric(sub("^chr", "", dt$seqnames)))
  order_key[dt$seqnames %in% c("chrX")] <- 23
  order_key[dt$seqnames %in% c("chrY")] <- 24
  order_key[dt$seqnames %in% c("chrM", "chrMT")] <- 25
  order_key[is.na(order_key)] <- Inf
  dt[order(order_key, start)]

}, simplify = FALSE)

for (paper in names(out_paper_beds)) {
  bed_dt <- out_paper_beds[[paper]]

  write.table(
    bed_dt,
    file = file.path(out_bed_dir, paste0(paper, ".bed")),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )

}
