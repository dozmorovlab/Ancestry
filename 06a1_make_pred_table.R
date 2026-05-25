library(openxlsx)
library(data.table)

exclude_regex_list <- c(
  "VCU-BC-0(1|2)$",
  "VCU-BC-0(1|2)_",
  "VCU-BC-0(1|2)-SGR_",
  "VCU-BC-0(1|2)-SGR$",
  "^WHIM"
)
remove_samples <- function(samples, regex_list = exclude_regex_list) {
  for (pattern in regex_list) {
    samples <- samples[!grepl(pattern, samples, perl = TRUE)]
  }
  samples
}

aeon_RNA_path <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/Aeon_RNA_GATK/Aeon_RNA_GATK_ae.csv"
ethseq_RNA_path <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/EthSEQ_RNA_GATK"
gnomad_RNA_path <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/gnomAD_RNA_GATK/gnomAD_sample_pred.csv"
jax_snpweights_RNA_dir <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/JAX_RNA_GATK/results/snpweights_inferanc"
raids_RNA_path <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/RAIDS_RNA_GATK/merged_results.csv"
admixture_RNA_path <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/Admixture_RNA_GATK/predictions/admixture_results.xlsx"

aeon_WGS_path <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/Aeon_WGS_GATK/Aeon_WGS_GATK_ae.csv"
ethseq_WGS_path <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/EthSEQ_WGS_GATK"
gnomad_WGS_path <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/gnomAD_WGS_GATK/gnomAD_sample_pred.csv"
jax_snpweights_WGS_dir <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/JAX_WGS/results/snpweights_inferanc"
raids_WGS_path <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/RAIDS_WGS_GATK/merged_results.csv"
admixture_WGS_path <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/Admixture_WGS_GATK/predictions/admixture_results.xlsx"

# in_dir <- ethseq_RNA_path
read_ethseq <- function(in_dir) {

  out_df <- data.frame()

  for (n in list.files(in_dir)) {
    report_dir <- file.path(
      in_dir,
      n,
      "All",
      "Report.txt"
    )

    if (file.exists(report_dir)) {  # Check if the file exists
      report_df <- read.csv(file = report_dir, sep = "\t")
      # Append to the pop-specific data frame
      out_df <- rbind(out_df, report_df)
    } else {
      next
    }
  }
  if (nrow(out_df) == 0 || ncol(out_df) == 0) {
    return(NULL)
  }
  out_df

  # Convert to data.table
  dt <- copy(as.data.table(out_df))

  # Fill NA
  dt[is.na(contribution), contribution := paste0(
    unlist(lapply(strsplit(pop, "\\|", fixed = TRUE), 
                  function(p) paste0(p, "(100.00%)", collapse = "|")))
  )]

  # Step 1: Expand contribution column into rows
  dt_expanded <- dt[
    , .(pop_info = unlist(strsplit(contribution, "\\|"))),
    by = .(sample.id)
  ]

  # Step 2: Extract population and percentage
  dt_expanded[
    !is.na(pop_info) & grepl("\\(", pop_info),
    `:=`(
      pop = sub("\\(.*", "", pop_info),
      percent = as.numeric(gsub(".*\\(([^)%]+)%\\)", "\\1", pop_info)) / 100
    )
  ]

  # Step 3: Wide format: one column per population
  dt_wide <- dcast(
    dt_expanded,
    sample.id ~ pop,
    value.var = "percent",
    fill = 0,
    fun.aggregate = sum
  )

  # Step 4: Add column for top prediction
  pop_cols <- c("AFR", "AMR", "EAS", "EUR", "SAS")

  # Add missing columns with 0s
  missing_cols <- setdiff(pop_cols, names(dt_wide))
  for (col in missing_cols) {
    dt_wide[[col]] <- 0
  }

  # Compute top_prediction safely
  dt_wide[
    , top_prediction := pop_cols[max.col(.SD, ties.method = "first")],
    .SDcols = pop_cols
  ]

  setnames(dt_wide, c("sample.id"), c("Sample"), skip_absent = TRUE)

  return(copy(dt_wide))
}
read_raids <- function(in_path) {
  dt <- tryCatch({
    fread(in_path)
  }, error = function(e) {
    return(data.table(c("Sample" = NULL, "top_prediction" = NULL)))
  })

  if (is.null(dt)) return(NULL)  # Exit early if fread failed

  # Convert to data.table
  dt <- copy(as.data.table(dt))

  dt$sample.id <- gsub("\\..*", "", dt$sample.id)
  setnames(dt, "SuperPop", "top_prediction", skip_absent = TRUE)
  setnames(dt, "sample.id", "Sample", skip_absent = TRUE)

  return(copy(dt))
}
read_aeon <- function(in_path) {
  dt <- fread(in_path)
  dt_sum <- dt[, lapply(.SD, sum), by = Superpopulation, .SDcols = is.numeric]

  dt_transposed <- transpose(dt_sum, keep.names = "Sample", make.names = "Superpopulation")

  dt_transposed[, "top_prediction" := names(.SD)[max.col(.SD, ties.method = "first")], .SDcols = !c("Sample")]

  copy(dt_transposed)
}
read_gnomad <- function(in_path) {
  dt <- fread(in_path, skip = 1)
  out_dt <- dt[, .(V1, V3)]
  setnames(out_dt, c("V1", "V3"), c("Sample", "top_prediction"), skip_absent = TRUE)
  copy(out_dt)
}
read_jax_snpweights <- function(in_dir) {

  out_dt <- rbindlist(lapply(list.files(
    in_dir, pattern = "tsv", full.names = TRUE
  ), function(in_file) {
    dt <- fread(in_file)
    dt$Sample <- basename(in_file)
    dt[, "top_prediction" := names(.SD)[max.col(.SD, ties.method = "first")], .SDcols = !c("Sample", "nSites")]
    copy(dt)
  }))

  copy(out_dt)

}
read_admixture <- function(in_path, start_at = 32) {

  in_df <- tryCatch({
    read.xlsx(in_path, sheet = 1)
  }, error = function(e) {
    return(data.frame(Sample = NULL, top_prediction = NULL))
  })
  out_dt <- as.data.table(in_df[start_at:nrow(in_df), ])
  setnames(
    out_dt,
    c("SampleID", "AssignedSuperpopulation"), c("Sample", "top_prediction"),
    skip_absent = TRUE
  )
  copy(out_dt)
}

results_RNA <- list(
  EthSEQ = read_ethseq(ethseq_RNA_path),
  RAIDS = read_raids(raids_RNA_path),
  Aeon = read_aeon(aeon_RNA_path),
  gnomAD = read_gnomad(gnomad_RNA_path),
  `JAX SNPWeights` = read_jax_snpweights(jax_snpweights_RNA_dir),
  Admixture = read_admixture(admixture_RNA_path)
)

results_WGS <- list(
  EthSEQ = read_ethseq(ethseq_WGS_path),
  RAIDS = read_raids(raids_WGS_path),
  Aeon = read_aeon(aeon_WGS_path),
  gnomAD = read_gnomad(gnomad_WGS_path),
  `JAX SNPWeights` = read_jax_snpweights(jax_snpweights_WGS_dir),
  Admixture = read_admixture(admixture_WGS_path)
)

dt_list_RNA <- lapply(names(results_RNA), function(n) {
  in_dt <- results_RNA[[n]]
  try(in_dt <- in_dt[Sample %in% remove_samples(in_dt$Sample)], silent = TRUE)
  in_dt$Sample <- gsub("_.*|\\.ancestry\\.tsv", "", in_dt$Sample)

  cols_needed <- c("Sample", "top_prediction")
  cols_available <- intersect(cols_needed, names(in_dt))

  # Fill in missing columns with NA if needed
  for (col in setdiff(cols_needed, cols_available)) {
    in_dt[, (col) := NA]
  }

  # Now this is safe
  result <- in_dt[, ..cols_needed]

  out_dt <- in_dt[, .(Sample, top_prediction)]
  setnames(out_dt, c("top_prediction"), c(n), skip_absent = TRUE)
  copy(out_dt)
})

dt_list_WGS <- lapply(names(results_WGS), function(n) {
  in_dt <- results_WGS[[n]]
  try(in_dt <- in_dt[Sample %in% remove_samples(in_dt$Sample)], silent = TRUE)
  in_dt$Sample <- gsub("_.*|\\.ancestry\\.tsv", "", in_dt$Sample)

  cols_needed <- c("Sample", "top_prediction")
  cols_available <- intersect(cols_needed, names(in_dt))

  # Fill in missing columns with NA if needed
  for (col in setdiff(cols_needed, cols_available)) {
    in_dt[, (col) := NA]
  }

  # Now this is safe
  result <- in_dt[, ..cols_needed]


  out_dt <- in_dt[, .(Sample, top_prediction)]
  setnames(out_dt, c("top_prediction"), c(n), skip_absent = TRUE)
  copy(out_dt)
})

out_dt_RNA <- Reduce(function(x, y) merge(x, y, by = "Sample", all = TRUE), dt_list_RNA)

out_dt_WGS <- Reduce(function(x, y) merge(x, y, by = "Sample", all = TRUE), dt_list_WGS)

write.csv(out_dt_RNA, file = "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/RNA_results.csv", row.names = FALSE)
write.csv(out_dt_WGS, file = "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/WGS_results.csv", row.names = FALSE)
