source("/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/06b0_sources.R")

process_inputs <- function(in_eigenvec_file, in_metadata_file) {

  # Read metadata
  metadata <- read.table(
    in_metadata_file,
    header = TRUE,
    stringsAsFactors = FALSE
  )

  # Read PCA results
  pca_data <- read.table(
    in_eigenvec_file,
    header = FALSE,
    comment.char = ""
  )

  # Assign column names for PCA data
  # colnames(pca_data) <- c("FamilyID", "SampleID", paste0("PC", 1:20))
  colnames(pca_data) <- c("FamilyID", "SampleID", paste0("PC", 1:(ncol(pca_data) - 2)))

  pca_data$FamilyID <- sapply(
    strsplit(as.character(pca_data$FamilyID), split = "_"), `[`, 1
  )

  # Join pca_data with metadata on SampleID
  pca_data_join <- pca_data %>%
    left_join(metadata, by = "SampleID") %>%
    mutate(FID = ifelse(!is.na(Superpopulation), Superpopulation, FamilyID.x)) %>%
    select(-Superpopulation)  # Remove the Superpopulation column

  pca_data_join$SID <- ifelse(
    pca_data_join$FID %in% c("AFR", "AMR", "EAS", "EUR", "SAS"),
    NA,
    pca_data_join$FID
  )
  pca_data_join$PID <- ifelse(
    pca_data_join$FID %in% c("AFR", "AMR", "EAS", "EUR", "SAS"),
    pca_data_join$FID,
    "Sample"
  )
  as.data.table(pca_data_join)
}

plot_pca <- function(
  in_data,
  in_title = "PCA of Genetic Data",
  flip_y = FALSE,
  flip_x = FALSE,
  colors_ = main_colors
) {
  if (flip_y) {
    in_data$PC2 <- in_data$PC2 * -1
  }
  if (flip_x) {
    in_data$PC1 <- in_data$PC1 * -1
  }

  ggplot(
    in_data,
    aes(x = PC1, y = PC2, color = PID)
  ) +
    geom_point(alpha = 0.5, size = 2) +
    ggrepel::geom_text_repel(
      aes(label = SID),
      color = "black",
      size = 3,
      segment.alpha = 0.5,
      max.overlaps = Inf,
      min.segment.length = 0,
      max.iter = Inf,
      max.time = 2
    ) +
    scale_color_manual(values = colors_) +
    theme_minimal() +
    labs(
      title = in_title,
      x = "Principal Component 1",
      y = "Principal Component 2",
      color = NULL
    ) +
    theme(legend.position = "bottom")
}
###################################################
wgs_results <- fread(wgs_results_path)

wgs_results$`gnomAD` <- c("afr" = "AFR", "nfe" = "EUR", "asj" = "EUR")[
  wgs_results$`gnomAD`
]

wgs_results[, consensus := as.character(apply(.SD, 1, function(x) {
  tab <- table(x)
  names(tab)[which.max(tab)]
}))]

wgs_results <- process_samples(wgs_results)

#################################################

rna_data <- process_samples(
  process_inputs(rna_eigenvec_file, metadata_file),
  sample_name = "FID"
)
wgs_data <- process_samples(
  process_inputs(wgs_eigenvec_file, metadata_file),
  sample_name = "FID"
)

rna_data <- merge(wgs_results, rna_data, by.x = "Sample", by.y = "FID", all = TRUE)
wgs_data <- merge(wgs_results, wgs_data, by.x = "Sample", by.y = "FID", all = TRUE)

# rna_data[PID == "Sample", PID := consensus]
# wgs_data[PID == "Sample", PID := consensus]

rna_pca_plot <- plot_pca(rna_data, "RNA")
wgs_pca_plot <- plot_pca(wgs_data, NULL, colors_ = main_colors)

###################################################################
# library(uwot)

# set.seed(42)  # for reproducibility

# get_umap <- function(
#   in_data,
#   in_metric = "euclidean",
#   in_n_neighbors = 15,
#   in_min_dist = 0.1
# ) {
#   umap_df <- in_data[, colnames(in_data)[grep("PC\\d*", colnames(in_data))]]
#   row.names(umap_df) <- in_data$SampleID
#   umap_result <- umap(
#     umap_df,
#     init = "spca",
#     n_neighbors = in_n_neighbors,
#     metric = in_metric,
#     min_dist = in_min_dist,
#     n_threads = parallel::detectCores(),
#     seed = 42
#   )
#   in_data$UMAP1 <- umap_result[, 1]
#   in_data$UMAP2 <- umap_result[, 2]
#   in_data
# }
# plot_umap <- function(in_data) {
#   ggplot(mapping = aes(x = UMAP1, y = UMAP2)) +
#     geom_point(
#       mapping = aes(color = PID, alpha = 0.5),
#       data = in_data[in_data$PID != "Sample", ],
#       size = 1
#     ) +
#     geom_point(
#       mapping = aes(alpha = NULL),
#       data = in_data[in_data$PID == "Sample", ],
#       size = 1
#     ) +
#     geom_text_repel(
#       mapping = aes(label = SID),
#       data = in_data,
#       size = 3,
#       segment.alpha = 0.25,
#       max.overlaps = Inf,
#       min.segment.length = 0,
#       max.iter = 100000
#     ) +
#     theme_minimal() +
#     labs(
#       x = "UMAP1",
#       y = "UMAP2",
#       color = "Ancestry"
#     ) +
#     guides(alpha = "none")
# }

# rna_data_u <- get_umap(
#   rna_data,
#   in_n_neighbors = 100,
#   in_min_dist = 0.001
# )
# wgs_data_u <- get_umap(
#   wgs_data,
#   in_n_neighbors = 100,
#   in_min_dist = 0.001
# )

# rna_plot <- plot_umap(rna_data_u) + labs(title = "RNA")
# wgs_plot <- plot_umap(wgs_data_u) + labs(title = "WGS")

# out_plot <- rna_plot + wgs_plot

# ggsave(
#   file.path(out_dir, "RNA_UMAP.png"),
#   rna_plot,
#   bg = "white",
#   height = 6,
#   width = 6,
#   scale = 2
# )
# ggsave(
#   file.path(out_dir, "WGS_UMAP.png"),
#   wgs_plot,
#   bg = "white",
#   height = 6,
#   width = 6,
#   scale = 2
# )

###################################################

# library(Rtsne)

# get_tsne <- function(
#   in_data
# ) {
#   tsne_df <- in_data[, colnames(in_data)[grep("PC\\d*", colnames(in_data))]]
#   row.names(tsne_df) <- in_data$SampleID
#   tsne_result <- Rtsne(
#     tsne_df,
#     dims = 2,
#     perplexity = 30,
#     verbose = FALSE,
#     max_iter = 500
#   )

#   in_data$TSNE1 <- tsne_result$Y[,1]
#   in_data$TSNE2 <- tsne_result$Y[,2]
#   in_data
# }
# plot_tsne <- function(in_data) {
#   ggplot(mapping = aes(x = TSNE1, y = TSNE2)) +
#     geom_point(
#       mapping = aes(color = PID, alpha = 0.5),
#       data = in_data[in_data$PID != "Sample", ],
#       size = 1
#     ) +
#     geom_point(
#       mapping = aes(alpha = NULL),
#       data = in_data[in_data$PID == "Sample", ],
#       size = 1
#     ) +
#     geom_text_repel(
#       mapping = aes(label = SID),
#       data = in_data,
#       size = 3,
#       segment.alpha = 0.25,
#       max.overlaps = Inf,
#       min.segment.length = 0,
#       max.iter = 100000
#     ) +
#     theme_minimal() +
#     labs(
#       x = "TSNE1",
#       y = "TSNE2",
#       color = "Ancestry"
#     ) +
#     guides(alpha = "none")
# }

# rna_data_t <- get_tsne(
#   rna_data
# )
# wgs_data_t <- get_tsne(
#   wgs_data
# )

# rna_plot <- plot_tsne(rna_data_t) + labs(title = "RNA")
# wgs_plot <- plot_tsne(wgs_data_t) + labs(title = "WGS")

# out_plot <- rna_plot + wgs_plot
