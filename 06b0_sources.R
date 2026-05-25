# load libraries
library(data.table)
library(ggplot2)
library(patchwork)
library(dplyr)
library(caret)
library(eulerr)
library(grid)
library(ggplotify)
library(R.utils)
library(future.apply)
library(ggpubr)
library(rtracklayer)
library(ggupset)
library(argparse)
library(ggrepel)
library(ggforce)
library(grImport2)

# Output directory for figures and tables
out_dir <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/figues_and_tables"

# /lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/06a1_make_pred_table.R
in_rna <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/RNA_results.csv"
in_wgs <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/WGS_results.csv"

# https://docs.google.com/spreadsheets/d/1h456SdQU5GM8cJQroTU4F0G3PEmkA77HBbEJbcBD1Tg/edit?gid=0#gid=0
in_sire <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/SIRE.tsv"

# /lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/06a3_GATK_vs_JAX.sh
in_RNA_csv <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/GATK_vs_JAX_RNA.csv"
in_WGS_csv <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/GATK_vs_JAX_WGS.csv"
in_RNA_unfiltered_csv <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/GATK_vs_JAX_RNA.unfiltered.csv"
in_WGS_unfiltered_csv <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/GATK_vs_JAX_WGS.unfiltered.csv"

# /lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/06a4_overlap_aims.sh
rna_tools_dir <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/aims_tools/results/RNA"
wgs_tools_dir <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/aims_tools/results/WGS"
rna_papers_dir <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/aims_tools/results/RNA_paper"
wgs_papers_dir <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/aims_tools/results/WGS_paper"
model_tools_dir <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/aims_tools/models"
model_papers_dir <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/aims_tools/papers/bed"

# /lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/06a2_Count_reads.sh
rna_reads_dir <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results_GATK/Data/RNA_GATK_reads.tsv"
wgs_reads_dir <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results_GATK/Data/WGS_GATK_reads.tsv"

# /lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/resources/get_resources.sh
usage_path <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/resources/usage.tsv"

# ?
rna_eigenvec_file <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/Admixture_RNA_GATK_0_0/pca_pruned/pruned_pca_results.eigenvec"
wgs_eigenvec_file <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/Admixture_WGS_GATK_0_0/pca_pruned/pruned_pca_results.eigenvec"

# ?
wgs_results_path <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/WGS_results.csv"

# ?
metadata_file <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Admixture/20130606_g1k_3202_samples_ped_population.txt"

# local
fig1a_path <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/figues_and_tables/fig1a.cairo.svg"

# Functions
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
process_gatk_vs_jax <- function(in_csv) {
  in_dt <- fread(in_csv)

  in_dt[Type == "shared_1", Type := "GATK Shared"]
  in_dt[Type == "shared_2", Type := "JAX Shared"]
  in_dt[Type == "unique_1", Type := "GATK Unique"]
  in_dt[Type == "unique_2", Type := "JAX Unique"]

  if (all(in_dt[Type == "JAX Shared"]$SNPs == in_dt[Type == "GATK Shared"]$SNPs)) {

    in_dt[Type == "JAX Shared", Type := "Shared"]
    in_dt[Type == "GATK Shared", Type := "Shared"]

    # Remove redundant rows (duplicates) based on SNPs and Type
    in_dt <- unique(in_dt, by = c("Type", "SNPs"))

  }

  process_samples(in_dt)

}

official_names <- c(
  Admixture = "Admixture",
  admixture = "Admixture",
  AEon = "AEon",
  Aeon = "AEon",
  aeon = "AEon",
  EthSEQ = "EthSEQ",
  EthSeq = "EthSEQ",
  Ethseq = "EthSEQ",
  ethseq = "EthSEQ",
  gnomAD = "gnomAD",
  gnomad = "gnomAD",
  RAIDS = "RAIDS",
  Raids = "RAIDS",
  raids = "RAIDS",
  `SNPWeights (JAX)` = "SNPWeights (JAX)",
  `JAX_SNPWeights` = "SNPWeights (JAX)",
  `JAX` = "SNPWeights (JAX)",
  `Kosoy et al.` = "Kosoy et al.",
  Kosoy = "Kosoy et al.",
  kosoy = "Kosoy et al.",
  `Nassir et al.` = "Nassir et al.",
  Nassir = "Nassir et al.",
  nassir = "Nassir et al.",
  `Torres et al.` = "Torres et al.",
  Torres = "Torres et al.",
  torres = "Torres et al."
)

# Colors

# https://tools.mysidewalk.com/palette-generator/index.html
main_colors <- c(
  Admixture = "#005ca4",
  AEon = "#845bb2",
  EthSEQ = "#d254a0",
  gnomAD = "#ff5e77",
  RAIDS = "#ff8843",
  `SNPWeights (JAX)` = "#F3BC00",
  `Kosoy et al.` = "#005ca4",
  `Nassir et al.` = "#ef558e",
  `Torres et al.` = "#F3BC00",
  `RNA GATK Unfiltered` = "#005ca4",
  `RNA JAX Unfiltered` = "#655db1",
  `RNA GATK Filtered` = "#a159af",
  `RNA JAX Filtered` = "#d254a0",
  `WGS GATK Unfiltered` = "#ff5e77",
  `WGS JAX Unfiltered` = "#ff7756",
  `WGS GATK Filtered` = "#ff9832",
  `WGS JAX Filtered` = "#F3BC00",
  AFR = "#005ca4",
  AMR = "#9a5ab0",
  EAS = "#ef558e",
  EUR = "#ff7b51",
  SAS = "#F3BC00"
)
