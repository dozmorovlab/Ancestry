source("/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/06b0_sources.R")

dir.create(out_dir, recursive = TRUE)

rna_pred_dt <- fread(in_rna)
wgs_pred_dt <- fread(in_wgs)
sire_dt <- fread(in_sire)
#####################################
# adjust gnomAD
gnomad_map <- c("afr" = "AFR", "nfe" = "EUR", "asj" = "EUR")
rna_pred_dt[
  , gnomAD := fifelse(
    gnomAD %in% names(gnomad_map),
    gnomad_map[gnomAD],
    toupper(gnomAD)
  )
]
wgs_pred_dt[
  , gnomAD := fifelse(
    gnomAD %in% names(gnomad_map),
    gnomad_map[gnomAD],
    toupper(gnomAD)
  )
]
##############################
pred_dt <- merge(
  rna_pred_dt,
  wgs_pred_dt,
  by = "Sample",
  all = TRUE,
  suffixes = c(" RNA", " WGS")
)

#######################

# # Sample adjustment
# pred_dt$Sample <- ifelse(
#   pred_dt$Sample %in% names(sample_map),
#   sample_map[pred_dt$Sample],
#   pred_dt$Sample
# )

# # Exclude
# pred_dt <- pred_dt[! Sample %in% exclude_list]
# sire_dt <- sire_dt[! Sample %in% exclude_list]

pred_dt <- process_samples(pred_dt)
sire_dt <- process_samples(sire_dt)
######################

pred_dt <- merge(
  pred_dt,
  sire_dt,
  by = "Sample",
  all = TRUE
)

table(sire_dt$SIRE)

pred_dt[, use_sample := SIRE %in% c("Black", "White")]
pred_dt[, positive := SIRE == "Black" & use_sample]

write.csv(pred_dt, file.path(out_dir, "pred_dt.csv"))

in_col <- "EthSEQ WGS"
in_use_pred_dt <- pred_dt[use_sample == TRUE]
results <- t(sapply(c(
  setdiff(names(pred_dt), c("Sample", "SIRE", "use_sample", "positive"))
), function(in_col) {

  in_df <- data.frame(
    pred = in_use_pred_dt[[in_col]] == "AFR",
    actual = in_use_pred_dt$positive
  )
  n_before_na <- nrow(in_df)
  in_df <- na.omit(in_df)
  n_after_na <- nrow(in_df)
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
  # OT <- n_before_na - n_after_na
  table_flat <- c("TP" = TP, "TN" = TN, "FP" = FP, "FN" = FN)#, "Other/Unused" = OT)

  # Assemble output
  c(table_flat, c_matrix$overall, c_matrix$byClass)
}))
results <- as.data.frame(results)
results$Name <- rownames(results)

# move Name to the front
results <- results[, c("Name", setdiff(names(results), "Name"))]

# add DataType
results$DataType <- unlist(
  regmatches(results$Name, gregexpr("RNA|WGS", results$Name))
)

results$Method <-  gsub(" .*", "", results$Name)

results <- as.data.table(results)
results[, Method := factor(official_names[Method], levels = unique(official_names[Method]))]

write.table(
  results,
  file.path(out_dir, "Black-AFR_vs_White.SIRE.csv"),
  sep = ",",
  row.names = FALSE
)
write.table(
  results[, .(Name, DataType, Method, TP, FP, TN, FN)],
  file.path(out_dir, "Black-AFR_vs_White.SIRE.table.csv"),
  sep = ",",
  row.names = FALSE
)

results[, color := main_colors[Method]]

wgs_acc_plot <- ggplot(
  results[DataType == "WGS"],
  aes(x = Method, y = `Balanced Accuracy`, fill = color)
) +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  labs(
    title = "WGS - SIRE",
    y = "Balanced Accuracy: Black \u2192 AFR",
    x = NULL
  ) +
  scale_fill_identity() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

wgs_f1_plot <- ggplot(
  results[DataType == "WGS"],
  aes(x = Method, y = F1, fill = color)
) +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  labs(
    title = "WGS - SIRE",
    y = "F1: Black \u2192 AFR",
    x = NULL
  ) +
  scale_fill_identity() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

####################################################

# ggsave(
#   file.path(out_dir, "WGS_Black-AFR_vs_White.accuracy.plot.png"),
#   wgs_acc_plot,
#   scale = 2,
#   height = 4,
#   width = 4,
#   bg = "white"
# )

# ggsave(
#   file.path(out_dir, "WGS_Black-AFR_vs_White.f1.plot.png"),
#   wgs_f1_plot,
#   scale = 2,
#   height = 4,
#   width = 4,
#   bg = "white"
# )

######################################################

# wgs_cols <- colnames(pred_dt)[unlist(
#   regmatches(colnames(pred_dt), gregexpr("RNA|WGS", colnames(pred_dt)))
# ) == "WGS"]



# consensus_dt$`gnomAD WGS` <- c("afr" = "AFR", "nfe" = "EUR", "asj" = "EUR")[
#   consensus_dt$`gnomAD WGS`
# ]
wgs_cols <- c(
  "Sample",
  colnames(pred_dt)[unlist(
    regmatches(colnames(pred_dt), gregexpr("RNA|WGS", colnames(pred_dt)))
  ) == "WGS"]
)
consensus_dt <- pred_dt[, ..wgs_cols]
consensus_dt[
  , consensus := as.character(
    apply(.SD, 1, function(x) {
      tab <- table(x)
      names(tab)[which.max(tab)]
    })
  ),
  .SDcols = setdiff(wgs_cols, c("Sample"))
]

# consensus_dt$Sample <- pred_dt$Sample

consensus_dt <- consensus_dt[, .(Sample, consensus)]

pred_dt[pred_dt$Sample %in% consensus_dt[consensus_dt$consensus == "NULL"]$Sample]

sample_id <- table(gsub(".*-", "", pred_dt$Sample))

########################################################################

# in_ex <- consensus_dt[["consensus"]]  %in% c("AFR", "EUR")
# in_ex[in_ex == FALSE] <- NA
# in_actual <- consensus_dt[["consensus"]][in_ex] == "AFR"

consensus_dt[, use_sample_consensus := consensus %in% c("AFR", "EUR")]
consensus_dt[, positive_consensus := consensus == "AFR" & use_sample_consensus]

write.csv(consensus_dt, file.path(out_dir, "consensus_dt.csv"))

in_col <- "EthSEQ WGS"
in_use_consensus_dt <- consensus_dt[use_sample_consensus == TRUE]

in_use_consensus_dt <- merge(in_use_consensus_dt, in_use_pred_dt, by = "Sample")

consensus_results <- t(sapply(c(
  setdiff(names(pred_dt), c("Sample", "SIRE", "use_sample", "positive"))
), function(in_col) {

  in_df <- data.frame(
    pred = in_use_consensus_dt[[in_col]] == "AFR",
    actual = in_use_consensus_dt$positive_consensus
  )
  n_before_na <- nrow(in_df)
  in_df <- na.omit(in_df)
  n_after_na <- nrow(in_df)
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
  # OT <- n_before_na - n_after_na
  table_flat <- c("TP" = TP, "TN" = TN, "FP" = FP, "FN" = FN)#, "Other/Unused" = OT)

  # Assemble output
  c(table_flat, c_matrix$overall, c_matrix$byClass)
}))

consensus_results <- as.data.frame(consensus_results)
consensus_results$Name <- rownames(consensus_results)

# move Name to the front
consensus_results <- consensus_results[, c("Name", setdiff(names(consensus_results), "Name"))]

# add DataType
consensus_results$DataType <- unlist(
  regmatches(consensus_results$Name, gregexpr("RNA|WGS", consensus_results$Name))
)

consensus_results$Method <-  gsub(" .*", "", consensus_results$Name)

consensus_results <- as.data.table(consensus_results)

consensus_results[, Method := factor(official_names[Method], levels = unique(official_names[Method]))]

write.table(
  consensus_results,
  file.path(out_dir, "AFR_RNA_WGS.consensus.csv"),
  sep = ",",
  row.names = FALSE
)
write.table(
  consensus_results[, .(Name, DataType, Method, TP, FP, TN, FN)],
  file.path(out_dir, "AFR_RNA_WGS.consensus.table.csv"),
  sep = ",",
  row.names = FALSE
)

consensus_results[, color := main_colors[Method]]

rna_acc_plot <- ggplot(
  consensus_results[DataType == "RNA"],
  aes(x = Method, y = `Balanced Accuracy`, fill = color)
) +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  labs(
    title = "RNA vs consensus (WGS)",
    y = "Balanced Accuracy - AFR",
    x = NULL
  ) +
  scale_fill_identity() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

rna_f1_plot <- ggplot(
  consensus_results[DataType == "RNA"],
  aes(x = Method, y = F1, fill = color)
) +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  labs(
    title = "RNA vs consensus (WGS)",
    y = "F1 - AFR",
    x = NULL
  ) +
  scale_fill_identity() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

####################################################

# ggsave(
#   file.path(out_dir, "AFR_RNA_WGS.consensus.accuracy.plot.png"),
#   rna_acc_plot,
#   scale = 2,
#   height = 4,
#   width = 4,
#   bg = "white"
# )

# ggsave(
#   file.path(out_dir, "AFR_RNA_WGS.consensus.f1.plot.png"),
#   rna_f1_plot,
#   scale = 2,
#   height = 4,
#   width = 4,
#   bg = "white"
# )

consensus_results$Method <- factor(
  consensus_results$Method,
  levels = consensus_results[DataType == "RNA"]$Method[order(consensus_results[DataType == "RNA"]$F1, decreasing = TRUE)]
)
results$Method <- factor(
  results$Method,
  levels = consensus_results[DataType == "RNA"]$Method[order(consensus_results[DataType == "RNA"]$F1, decreasing = TRUE)]
)
wgs_f1_plot <- ggplot(
  results[DataType == "WGS"],
  aes(
    x = Method,
    y = F1,
    fill = color
  )
) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  labs(
    title = "WGS",
    subtitle = "SIRE as standard",
    y = "F1: Black \u2192 AFR",
    x = NULL
  ) +
  scale_fill_identity() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

rna_f1_plot <- ggplot(
  consensus_results[DataType == "RNA"],
  aes(
    x = Method,
    y = F1,
    fill = color
  )
) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  labs(
    title = "RNA",
    subtitle = "Consensus of WGS as standard",
    y = "F1: AFR",
    x = NULL
  ) +
  scale_fill_identity() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

out_f1_plot <- wgs_f1_plot + rna_f1_plot

# ggsave(
#   file.path(out_dir, "f1_acc.plot.png"),
#   out_plot,
#   scale = 3,
#   height = 4,
#   width = 4,
#   bg = "white"
# )
