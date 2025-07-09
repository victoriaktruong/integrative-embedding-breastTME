# Load libraries
library(dplyr)
library(readr)
library(tidyr)

# Load marker gene table from the concatenated embedding
df <- read.csv("scgpt_06-16_11clusters_markers_with_CTA")
df <- read.csv("mxbai_06-16-21clusters_markers_with_CTA.csv")
#df <- read.csv("concat_05-20_21clusters_markers_with_CTA.csv")

# Filter marker genes based on log2FC, adjusted p-value, and expression threshold
filtered_df <- df %>%
  filter(avg_log2FC > 1, p_val_adj < 0.05, pct.1 >= 0.25)

# Identify predicted cell types with no genes passing the filter
all_cell_types <- unique(df$predicted_cell_type)
missing_cell_types <- setdiff(all_cell_types, unique(filtered_df$predicted_cell_type))

# For each missing cell type, select top 100 genes by p-value, log2FC, and expression
top_genes_df <- df %>%
  filter(predicted_cell_type %in% missing_cell_types) %>%
  arrange(predicted_cell_type, p_val_adj, desc(avg_log2FC), desc(pct.1)) %>%
  group_by(predicted_cell_type) %>%
  slice_head(n = 100) %>%
  ungroup()

# Combine filtered genes with top-ranked supplements and remove unused columns
final_df <- bind_rows(filtered_df, top_genes_df) %>%
  select(-X, -cluster_size)

# Save filtered marker gene output
write.csv(final_df, "scgpt_06-16_11clusters_markers_with_CTA_filtered_ranked.csv", row.names = FALSE)
#write.csv(final_df, "mxbai_06-16_21clusters_markers_with_CTA_filtered_ranked.csv", row.names = FALSE)
#write.csv(final_df, "concat_05-20_21clusters_markers_with_CTA_filtered_ranked.csv", row.names = FALSE)

# Load filtered marker outputs from each embedding
df_scgpt <- read_csv("scgpt_06-16_11clusters_markers_with_CTA_filtered_ranked.csv")
df_mxbai <- read_csv("mxbai_06-16_21clusters_markers_with_CTA_filtered_ranked.csv")
df_concat <- read_csv("concat_05-20_21clusters_markers_with_CTA_filtered_ranked.csv")

# Count number of genes per cell type for each embedding
summary_concat <- df_concat %>%
  group_by(predicted_cell_type) %>%
  summarise(genes_concat = n(), .groups = "drop")

summary_mxbai <- df_mxbai %>%
  group_by(predicted_cell_type) %>%
  summarise(genes_mxbai = n(), .groups = "drop")

summary_scgpt <- df_scgpt %>%
  group_by(predicted_cell_type) %>%
  summarise(genes_scgpt = n(), .groups = "drop")

# Merge summaries across embeddings and fill missing values with zeros
summary_all <- full_join(summary_concat, summary_mxbai, by = "predicted_cell_type") %>%
  full_join(summary_scgpt, by = "predicted_cell_type") %>%
  replace_na(list(genes_concat = 0, genes_mxbai = 0, genes_scgpt = 0))

# Save the summary table for gene counts by cell type
write_csv(summary_all, "summary_genes_by_celltype.csv")