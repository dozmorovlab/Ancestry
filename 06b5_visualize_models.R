source("/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/06b0_sources.R")

out_dir <- "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/figues_and_tables"

# Count model SNPs
bed_dirs <- list(
  models = "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/aims_tools/models",
  papers = "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/aims_tools/papers/bed"
)

in_dir <- bed_dirs[[1]]
beds <- sapply(bed_dirs, function(in_dir) {

  in_beds <- list.files(in_dir, pattern = "\\.bed", full.names = TRUE)

  in_beds <- setNames(in_beds, gsub("\\.bed", "", lapply(in_beds, basename)))

  sapply(in_beds, function(in_bed) {
    import(in_bed)
  }, simplify = FALSE)

}, simplify = FALSE)

bed_list <- beds[[1]]
model_bar_plots <- sapply(beds, function(bed_list) {

  lengths <- sapply(bed_list, length, simplify = FALSE)
  l_dt <- data.table(
    name = official_names[names(lengths)],
    n = c(unlist(lengths)),
    color = main_colors[official_names[names(lengths)]]
  )

  ggplot(
    l_dt,
    aes(
      x = reorder(name, -n),
      y = n,
      fill = color
    )
  ) +
    geom_col(color = "black") +
    theme_minimal() +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = "Models",
      x = NULL,
      y = "Number of model SNPs"
    ) +
    scale_fill_identity() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

}, simplify = FALSE)

out_model_bar_plot <- wrap_plots(model_bar_plots)

model_table <- rbindlist(sapply(beds, function(bed_list) {

  lengths <- sapply(bed_list, length, simplify = FALSE)
  l_dt <- data.table(
    name = names(lengths),
    n = c(unlist(lengths))
  )

}, simplify = FALSE))

# ggsave(
#   file.path(out_dir, "model_bar_plot.png"),
#   out_model_bar_plot,
#   width = 6,
#   height = 3,
#   dpi = 600,
#   bg = "white",
#   scale = 2
# )

write.csv(
  model_table,
  file = file.path(out_dir, "model_counts.csv"),
  row.names = FALSE
)

####################################################################

# Tool BEDs
tool_dir <- bed_dirs$models
tool_paths <- list.files(tool_dir, full.names = TRUE, pattern = ".bed")
names(tool_paths) <- gsub("\\.bed$", "", basename(tool_paths))
tool_grs <- sapply(tool_paths, function(bed_path) {
  import(bed_path, format = "BED")
}, simplify = FALSE)

paper_dir <- bed_dirs$papers
paper_paths <- list.files(paper_dir, full.names = TRUE, pattern = ".bed")
names(paper_paths) <- gsub("\\.bed$", "", basename(paper_paths))
paper_grs <- sapply(paper_paths, function(bed_path) {
  import(bed_path, format = "BED")
}, simplify = FALSE)

tool_loci <- sapply(tool_grs, function(gr) {
  paste(as.character(seqnames(gr)), start(gr), sep = ":")
}, simplify = FALSE)

paper_loci <- sapply(paper_grs, function(gr) {
  paste(as.character(seqnames(gr)), start(gr), sep = ":")
}, simplify = FALSE)

###################################################################

# p1 <- as.ggplot(~eulerr:::plot.eulergram(plot(fit, quantities = TRUE)))

in_list <- paper_loci
get_pairwise_eulerr <- function(
  in_list,
  featured_row = names(in_list),
  exclude_col = c()
) {
  row_list <- names(in_list)[names(in_list) %in% featured_row]
  wrap_plots(sapply(row_list, function(x) {
    col_list <- names(in_list)[
      names(in_list) != x & ! names(in_list) %in% exclude_col
    ]
    wrap_plots(sapply(col_list, function(y) {
      e_list <- setNames(list(in_list[[x]], in_list[[y]]), c(x, y))
      fit <- euler(e_list)

      as.ggplot(function() {
        eulerr:::plot.eulergram(plot(
          fit,
          quantities = TRUE,
          fills = list(alpha = 0.6)
        ))
      })

    }, simplify = FALSE), nrow = 1)
  }, simplify = FALSE), nrow = length(in_list))
}
get_pairwise_bar <- function(
  in_list,
  featured_row = names(in_list),
  exclude_col = c()
) {
  row_list <- names(in_list)[names(in_list) %in% featured_row]
  wrap_plots(sapply(row_list, function(x) {
    col_list <- names(in_list)[
      names(in_list) != x & ! names(in_list) %in% exclude_col
    ]
    z <- sapply(col_list, function(y) {
      Count <- sum(duplicated(c(in_list[[x]], in_list[[y]])))
      Percent <- Count / length(in_list[[x]])
      c(Count = Count, Percent = Percent)
    }, simplify = FALSE)

    df <- data.frame(
      Sample = names(z),
      Count = sapply(z, `[[`, "Count"),
      Percent = sapply(z, `[[`, "Percent")
    )

    c_plot <- ggplot(df, aes(x = Sample, y = Count)) +
      geom_bar(stat = "identity") +
      labs(
        title = paste("Overlap with", x),
        y = "Count"
      ) +
      theme_minimal() +
      coord_flip()

    p_plot <- ggplot(df, aes(x = Sample, y = Percent)) +
      geom_bar(stat = "identity") +
      labs(
        title = paste("Overlap with", x),
        y = paste("Percent of", x)
      ) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 0.01)) +
      theme_minimal() +
      coord_flip()

    c_plot + p_plot

  }, simplify = FALSE), nrow = length(in_list))
}

# Tools vs tools
e_plots <- get_pairwise_eulerr(tool_loci)
#f_plots <- get_pairwise_bar(tool_loci)

# Tools vs papers
# f_plots <- get_pairwise_eulerr(
#   all_loci,
#   featured_row = c("Kossoy", "Nassir", "Torres"),
#   exclude_col = c("Kossoy", "Nassir", "Torres")
# )
# ggsave(
#   file.path(out_dir, "pairwise_tools.png"),
#   e_plots,
#   bg = "white",
#   scale = 4,
#   limitsize = FALSE
# )
#########################################################

get_combos <- function(
  in_list,
  in_names = names(in_list)
) {

  result <- lapply(seq_along(in_list), function(k) {
    combn(in_list, k, function(comb_list) {
      Reduce(intersect, comb_list)
    }, simplify = FALSE)
  })

  result <- lapply(seq_along(result), function(i) {
    #lapply(result[i], function(x) {
    #  unique(c(unlist(x), unlist(result[(i + 1): length(result)])))
    #})
    x <- unlist(result[i], recursive = FALSE)
    if (i == length(result)) {
      return(x)
    }
    y <- x[[1]]
    lapply(x, function(y) {
      y[! y %in% unlist(result[(i + 1): length(result)])]
    })
  })

  out_names <- lapply(seq_along(in_names), function(k) {
    combn(in_names, k, function(comb_names) {
      paste(comb_names, collapse = "-")
    }, simplify = FALSE)
  })

  out_result <- unlist(result, recursive = FALSE)

  names(out_result) <- unlist(out_names, recursive = FALSE)

  out_result
}

plots_by_order <- function(df) {
  df$order <- unlist(lapply(strsplit(df$Combo, "-"), length))
  sapply(unique(df$order), function(o) {
    df_subset <- df[df$order == o, ]
    df_subset$Combo <- factor(
      df_subset$Combo,
      levels = df_subset$Combo[order(df_subset$Count, decreasing = TRUE)]
    )
    ggplot(
      df_subset,
      aes(x = Combo, y = Count)
    ) +
      geom_bar(stat = "identity") +
      axis_combmatrix(
        sep = "-"
      ) +
      theme_minimal()
  }, simplify = FALSE)
}

tool_combo_loci <- get_combos(tool_loci)

tool_combo_counts <- sapply(tool_combo_loci, length, simplify = FALSE)
tool_combo_counts <- tool_combo_counts[tool_combo_counts != 0]

# Verify
sum(unlist(sapply(tool_combo_loci, length, simplify = FALSE))) ==
  length(unique(unlist(tool_loci)))

tool_combo_counts_plotting_df <- data.table(
  Combo = names(tool_combo_counts),
  Count = as.numeric(tool_combo_counts)
)

u_plot <- ggplot(
  tool_combo_counts_plotting_df,
  aes(x = Combo, y = Count)
) +
  geom_bar(stat = "identity") +
  axis_combmatrix(
    sep = "-"
  ) +
  theme_minimal()

u_plot_order <- wrap_plots(plots_by_order(tool_combo_counts_plotting_df), nrow = 1)

# tool_long_dt <- rbindlist(
#   lapply(names(tool_loci), function(nm) data.table(name = nm, value = tool_loci[[nm]]))
# )

# ggplot(
#   tool_long_dt,
#   aes(x = value)
# ) +
#   geom_bar() +
#   scale_x_upset() +
#   labs(title = NULL)

# ggsave(
#   file.path(out_dir, "upset_tools.png"),
#   u_plot_order,
#   bg = "white",
#   scale = 2,
#   limitsize = FALSE
# )

#####################################################################
paper_fit <- euler(list(
  `Kosoy et al.` = unique(beds$papers$Kosoy$name),
  `Nassir et al.` = unique(beds$papers$Nassir$name),
  `Torres et al.` = unique(beds$papers$Torres$name)
))

paper_euler <- as.ggplot(~eulerr:::plot.eulergram(plot(
  paper_fit,
  quantities = TRUE,
  fill = main_colors[names(paper_fit)]
)))
