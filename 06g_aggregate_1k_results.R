library(openxlsx)
library(data.table)
library(kgp)
library(caret)
library(ggplot2)
library(patchwork)

in_dir <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/overlap_chrs/predicted_RNA_10_select52"
out_dir <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/figues_and_tables"

dir.create(out_dir, recursive = TRUE)

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
  dt <- as.data.table(out_df)

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

  return(dt_wide)
}
read_raids <- function(in_path = NULL, in_dir = NULL) {
  if (is.null(in_path)) {
    in_path <- file.path(in_dir, "merged_results.csv")
  }
  dt <- tryCatch({
    fread(in_path)
  }, error = function(e) {
    return(data.table(c("Sample" = NULL, "top_prediction" = NULL)))
  })

  if (is.null(dt)) return(NULL)  # Exit early if fread failed

  # Convert to data.table
  dt <- as.data.table(dt)

  dt$sample.id <- gsub("\\..*", "", dt$sample.id)
  setnames(dt, "SuperPop", "top_prediction", skip_absent = TRUE)
  setnames(dt, "sample.id", "Sample", skip_absent = TRUE)

  return(dt)
}
read_aeon <- function(in_path = NULL, in_dir = NULL, aeon_name = "Aeon_ae.csv") {
  if (is.null(in_path)) {
    in_path <- file.path(in_dir, aeon_name)
  }
  dt <- fread(in_path)
  dt_sum <- dt[, lapply(.SD, sum), by = Superpopulation, .SDcols = is.numeric]

  dt_transposed <- transpose(dt_sum, keep.names = "Sample", make.names = "Superpopulation")

  dt_transposed[, "top_prediction" := names(.SD)[max.col(.SD, ties.method = "first")], .SDcols = !c("Sample")]

  dt_transposed
}
read_gnomad <- function(in_path = NULL, in_dir = NULL, gnomad_name = "gnomAD_sample_pred.rna_filtered.csv") {
  if (is.null(in_path)) {
    in_path <- file.path(in_dir, gnomad_name)
  }
  dt <- fread(in_path, skip = 1)
  out_dt <- dt[, .(V1, V3)]
  setnames(out_dt, c("V1", "V3"), c("Sample", "top_prediction"), skip_absent = TRUE)
  out_dt
}
read_jax_snpweights <- function(in_dir) {

  out_dt <- rbindlist(lapply(list.files(
    in_dir, pattern = "tsv", full.names = TRUE
  ), function(in_file) {
    dt <- fread(in_file)
    dt$Sample <- basename(in_file)
    dt[, "top_prediction" := names(.SD)[max.col(.SD, ties.method = "first")], .SDcols = !c("Sample", "nSites")]
    dt
  }))

  out_dt

}
read_admixture <- function(in_path = NULL, in_dir = NULL, start_at = 32) {
  if (is.null(in_path)) {
    in_path <- file.path(in_dir, "predictions", "pruned_merged_data.5_Proportions.csv")
  }
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
  out_dt
}

read_dict <- list(
  "Admixture" = read_admixture,
  "Aeon" = read_aeon,
  "EthSEQ" = read_ethseq,
  "gnomAD" = read_gnomad,
  "RAIDS" = read_raids
)

results <- rbindlist(
  lapply(list.files(file.path(in_dir), pattern = ".*_snps", full.names = TRUE),
    function(snps_dir) {
      rep_dirs <- list.files(snps_dir, full.names = TRUE)
      if (length(rep_dirs) == 0L) return(data.table())  # guard empty case

      rep_results <- rbindlist(
        lapply(rep_dirs, function(rep_dir) {
          pred_dirs <- list.files(rep_dir, full.names = TRUE)
          if (length(pred_dirs) == 0L) return(data.table())

          pred_results <- rbindlist(
            lapply(pred_dirs, function(pred_dir) {
              tryCatch(
                {
                  dt <- as.data.table(
                    read_dict[[basename(pred_dir)]](in_dir = pred_dir)
                  )
                  dt[, tool := basename(pred_dir)]
                  dt
                },
                error = function(e) data.table()
              )
            }),
            fill = TRUE
          )

          if (nrow(pred_results) == 0L) return(data.table())
          pred_results[, rep := basename(rep_dir)]
          pred_results
        }),
        fill = TRUE
      )

      if (nrow(rep_results) == 0L) return(data.table())
      rep_results[, n_snps := basename(snps_dir)]
      rep_results
    }
  ),
  fill = TRUE
)

results$top_prediction[results$top_prediction == "nfe"] <- "EUR"

results_meta <- merge(results, kgpe, by.x = "Sample", by.y = "id", all.x = TRUE, all.y = FALSE)
results_meta <- results_meta[! is.na(results_meta$Sample)]

get_confusion <- function(in_dt, in_pred_name, in_actual_name, target_pop) {

  in_pred <- toupper(as.list(in_dt[[in_pred_name]])) == target_pop
  in_actual <- toupper(as.list(in_dt[[in_actual_name]])) == target_pop

  in_df <- na.omit(data.frame(
    pred = in_pred,
    actual = in_actual
  ))
  pred <- factor(in_df$pred, levels = c(TRUE, FALSE))
  actual <- factor(in_df$actual, levels = c(TRUE, FALSE))
  c_matrix <- confusionMatrix(
    pred,
    actual,
    mode = "everything",
    positive = "TRUE"
  )

  conf_df <- as.data.frame(c_matrix$table)


  # Convert table
  TP <- conf_df[
    conf_df$Prediction == TRUE & conf_df$Reference == TRUE, "Freq"
  ]
  TN <- conf_df[
    conf_df$Prediction == FALSE & conf_df$Reference == FALSE, "Freq"
  ]
  FP <- conf_df[
    conf_df$Prediction == TRUE & conf_df$Reference == FALSE, "Freq"
  ]
  FN <- conf_df[
    conf_df$Prediction == FALSE & conf_df$Reference == TRUE, "Freq"
  ]
  table_flat <- c("TP" = TP, "TN" = TN, "FP" = FP, "FN" = FN)

  # Assemble output
  v <- c(table_flat, c_matrix$overall, c_matrix$byClass)
  out_dt <- as.data.table(as.list(v))
  out_dt$target_pop <- target_pop
  out_dt
}

confusion_results <- data.table()

for (in_pop in unique(results_meta$reg)) {
  for (in_tool in unique(results_meta$tool)) {
    for (in_rep in unique(results_meta$rep)) {
      for (in_n_snps in unique(results_meta$n_snps)) {

        next_row <- get_confusion(
          in_dt = results_meta[tool == in_tool][rep == in_rep][n_snps == in_n_snps],
          in_pred_name = "top_prediction",
          in_actual_name = "reg",
          target_pop = in_pop
        )
        next_row$tool <- in_tool
        next_row$rep <- in_rep
        next_row$n_snps <- in_n_snps
        confusion_results <- rbind(confusion_results, next_row)

      }
    }
  }
}

out_plot <- wrap_plots(sapply(c("10_snps", "100_snps", "1k_snps", "10k_snps", "100k_snps"), function(in_n_snps) {
  ggplot(
    confusion_results[n_snps == in_n_snps],
    aes(x = target_pop, y = F1, fill = tool)
  ) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = in_n_snps) +
    ylim(0, 1)

}, simplify = FALSE), nrow = 1)

ggsave(
  file.path(out_dir, "1k_F1_RNA.png"),
  out_plot,
  bg = "white",
  scale = 5,
  limitsize = FALSE
)

write.csv(
  results_meta,
  file = "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/1k_RNA_results.csv",
  row.names = FALSE
)

#############################################################

# Actual snps

snp_counts_1k_path <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/overlap_chrs/random_RNA_10_select52.counts.tsv"

snp_counts_1k_dt <- fread(snp_counts_1k_path)

# dir1	dir2	file	snp_count

plots

wrap_plots(sapply(c("10_snps", "100_snps", "1k_snps", "10k_snps", "100k_snps"), function(in_n_snps) {
  ggplot(
    snp_counts_1k_dt[dir1 == in_n_snps],
    aes(x = dir2, y = snp_count, fill = dir2)
  ) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = in_n_snps)

}, simplify = FALSE), nrow = 1)
