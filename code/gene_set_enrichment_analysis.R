# Load libraries
library(readr)
library(dplyr)
library(enrichR)
library(ggplot2)

# Load filtered marker genes from one embedding
df <- read_csv("scgpt_06-16_11clusters_markers_with_CTA_filtered_ranked.csv")
# df <- read_csv("mxbai_06-16_21clusters_markers_with_CTA_filtered_ranked.csv")
# df <- read_csv("concat_05-20_21clusters_markers_with_CTA_filtered_ranked.csv")

# Extract unique genes from CD8+ T-cells
cd8_genes <- unique(df$gene[df$predicted_cell_type == "CD8+ T-cells"])

# Run GO enrichment analysis using Enrichr
dbs <- c("GO_Biological_Process_2025")
enrich_results <- enrichr(cd8_genes, dbs)
go_results <- enrich_results[["GO_Biological_Process_2025"]]

# Save top 5 enriched GO terms sorted by adjusted p-value
top_go <- go_results %>%
  arrange(Adjusted.P.value) %>%
  head(5)

write.csv(top_go, "scgpt_top5_GO_CD8+.csv", row.names = FALSE)
# write.csv(top_go, "mxbai_top5_GO_CD8+.csv", row.names = FALSE)
# write.csv(top_go, "concat_top5_GO_CD8+.csv", row.names = FALSE)

# Load top GO terms for each embedding method and add labels
scgpt_df  <- read_csv("scgpt_top5_GO_CD8+.csv")  %>% mutate(Embedding = "scGPT")
mxbai_df  <- read_csv("mxbai_top5_GO_CD8+.csv")  %>% mutate(Embedding = "Mixedbread")
concat_df <- read_csv("concat_top5_GO_CD8+.csv") %>% mutate(Embedding = "scMIXE")

# Combine and process for visualization
combined_df <- bind_rows(concat_df, mxbai_df, scgpt_df) %>%
  mutate(
    Count = as.numeric(sub("/.*", "", Overlap)),  # Extract gene count from "x/y"
    neg_log10_pval = -log10(`Adjusted.P.value`)   # For coloring by significance
  )

# Save the combined enrichment results
write_csv(combined_df, "combined_enrichr_results.csv")


# Clean GO terms and set embedding order
combined_df$Embedding <- factor(combined_df$Embedding, levels = c("scGPT", "Mixedbread", "scMIXE"))
combined_df <- combined_df %>%
  mutate(Term = sub(" \\(GO:.*\\)", "", Term))  # Remove GO ID from term names

# Create multi-panel dot plot
ggplot(combined_df, aes(x = Combined.Score, y = reorder(Term, Combined.Score),
                        color = Adjusted.P.value, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue", name = "Adjusted P-value") +
  facet_wrap(~ Embedding, scales = "free_x") +
  labs(x = "Combined Score", y = "GO Term", size = "Gene Count") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 9),
    strip.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  )