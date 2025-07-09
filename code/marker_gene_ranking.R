df <- read.csv("bassez_concat_06-28_15clusters_markers_with_CTA.csv")

library(dplyr)

# FURTHER FILTER 
filtered_df <- df %>%
  filter(avg_log2FC > 1, p_val_adj < 0.05, pct.1 >= 0.25)

# View the result
print(dim(filtered_df))
head(filtered_df)

# Get the set of all predicted cell types
all_cell_types <- unique(df$predicted_cell_type)

# Get cell types with no genes passing the filter
missing_cell_types <- setdiff(all_cell_types, unique(filtered_df$predicted_cell_type))

# For each missing cell type, select top 100 by p_val_adj, logFC, pct.1
top_genes_df <- df %>%
  filter(predicted_cell_type %in% missing_cell_types) %>%
  arrange(predicted_cell_type, p_val_adj, desc(avg_log2FC), desc(pct.1)) %>%
  group_by(predicted_cell_type) %>%
  slice_head(n = 100) %>%
  ungroup()

# Combine filtered + top selected genes
final_df <- bind_rows(filtered_df, top_genes_df)

final_df <- final_df %>% select(-X, -cluster_size)

# View result
dim(final_df)
head(final_df)

length(missing_cell_types)  
nrow(top_genes_df) 

# Save filtered output if desired
write.csv(final_df, "bassez_concat_06-28_15clusters_markers_with_CTA_filtered_100.csv", row.names = FALSE)


library(dplyr)
library(readr)
library(tidyr)


# Load the filtered marker datasets
df_scgpt<- read_csv("bassez_scgpt_06-28_9clusters_all_markers_with_CTA_filtered_100.csv")
df_mxbai <- read_csv("bassez_mxbai_06-28_14clusters_all_markers_with_CTA_filtered_100.csv")
df_concat <- read_csv("bassez_concat_06-28_15clusters_markers_with_CTA_filtered_100.csv")

# Count genes per predicted_cell_type
summary_concat <- df_concat %>%
  group_by(predicted_cell_type) %>%
  summarise(genes_concat = n(), .groups = "drop")

summary_mxbai <- df_mxbai %>%
  group_by(predicted_cell_type) %>%
  summarise(genes_mxbai = n(), .groups = "drop")

summary_scgpt <- df_scgpt %>%
  group_by(predicted_cell_type) %>%
  summarise(genes_scgpt = n(), .groups = "drop")

# Merge summaries by predicted_cell_type
summary_all <- full_join(summary_concat, summary_mxbai, by = "predicted_cell_type") %>%
  full_join(summary_scgpt, by = "predicted_cell_type") %>%
  replace_na(list(genes_concat = 0, genes_mxbai = 0, genes_scgpt = 0))


# View or export
print(summary_all, n=14)
write_csv(summary_all, "summary_genes_by_celltype_06-30.csv")



#--------------------------------------------------------------------------------------------------------

#EnrichR (Enrichment Analysis)

library(readr)
library(dplyr)
library(enrichR)

# Load your file
#df <- read_csv("concat_05-20_21clusters_markers_with_CTA_filtered_ranked.csv")
#df <- read_csv("mxbai_06-16_11clusters_markers_with_CTA_filtered_ranked.csv")
df <- read_csv("scgpt_06-16_11clusters_markers_with_CTA_filtered_100.csv")

# Get all unique genes (if you’re not subsetting by cell type)
#gene_list <- unique(df$gene)
cd8_genes <- unique(df$gene[df$predicted_cell_type == "CD8+ T-cells"])

# Select the enrichment database
dbs <- c("GO_Biological_Process_2025")

# Run enrichment
enrich_results <- enrichr(cd8_genes, dbs)

# Extract results
go_results <- enrich_results[["GO_Biological_Process_2025"]]

# View top 10 enriched terms sorted by Adjusted P-value
top_go <- go_results %>%
  arrange(Adjusted.P.value) %>%
  head(5)

print(top_go)

# Save results
write.csv(top_go, "scgpt_top5_GO_CD8+.csv", row.names = FALSE)


#---------------------------------------------------------
# COMBINE DATSSETS INTO 1
# Load necessary libraries
library(dplyr)
library(readr)

# Load each dataset
concat_df <- read_csv("concat_top5_GO_CD8+.csv")
mxbai_df <- read_csv("mxbai_top5_GO_CD8+.csv")
scgpt_df <- read_csv("scgpt_top5_GO_CD8+.csv")

# Add embedding labels
concat_df <- concat_df %>% mutate(Embedding = "scMIXE")
mxbai_df <- mxbai_df %>% mutate(Embedding = "Mixedbread")
scgpt_df <- scgpt_df %>% mutate(Embedding = "scGPT")

# Combine all into one dataframe
combined_df <- bind_rows(concat_df, mxbai_df, scgpt_df)

# Extract the gene count from the Overlap column (formatted as "x/y")
combined_df <- combined_df %>%
  mutate(Count = as.numeric(sub("/.*", "", Overlap)))

# Calculate -log10 adjusted p-value for plotting
combined_df <- combined_df %>%
  mutate(neg_log10_pval = -log10(`Adjusted.P.value`))

# Save the combined file for plotting
write_csv(combined_df, "combined_enrichr_results.csv")


# -----------------------------------------------------------

# DOT PLOT

# Load necessary libraries
library(ggplot2)
library(readr)
library(dplyr)

# Load the combined dataset
combined_df <- read_csv("combined_enrichr_results.csv")

# Set the desired embedding order
combined_df$Embedding <- factor(combined_df$Embedding, 
                                levels = c("scGPT", "Mixedbread", "scMIXE"))

# Clean GO terms by removing the GO ID
combined_df <- combined_df %>%
  mutate(Term = sub(" \\(GO:.*\\)", "", Term))

# Create the multi-panel dot plot
ggplot(combined_df, aes(x = Combined.Score, y = reorder(Term, Combined.Score), 
                        color = Adjusted.P.value, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue", name = "Adjusted P-value") +
  facet_wrap(~ Embedding, scales = "free_x") +
  labs(#title = "Top Enriched GO Terms Across Embedding Strategies",
       x = "Combined Score",
       y = "GO Term",
       size = "Gene Count") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9),
        strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))





