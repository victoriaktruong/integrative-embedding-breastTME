# ---- Load Libraries ----
library(arrow)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(SingleR)
library(celldex)
library(SummarizedExperiment)
library(tidyr)

# ---- Load Embeddings and Create Seurat Object ----
# Load concatenated embeddings
embeddings <- read_feather("concat_emb_mxbai_scgpt_hvg2.feather")

# Convert embeddings to a matrix format suitable for Seurat
embeddings_matrix <- as.matrix(embeddings)

# Load cell barcode information
barcodes <- read.table("barcodes.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
rownames(embeddings_matrix) <- barcodes$V1
embeddings_sparse <- as(embeddings_matrix, "dgCMatrix")

# Create initial Seurat object from sparse matrix of embeddings
seurat_obj <- CreateSeuratObject(counts = t(embeddings_sparse), project = "BreastCancer", assay = "RNA")
rownames(embeddings_matrix) <- colnames(seurat_obj)
colnames(embeddings_matrix) <- 1:ncol(embeddings_matrix)
seurat_obj[["embeddings"]] <- CreateDimReducObject(embeddings = embeddings_matrix, key = "embed_", assay = "RNA")


# ---- Perform Clustering using Louvain Algorithm ----
options(future.globals.maxSize = 10 * 1024^3)  

# Find nearest neighbors using the embeddings
seurat_obj <- FindNeighbors(seurat_obj, reduction = "embeddings", dims = 1:1536)

# Cluster at multiple resolutions to determine the optimal resolution
resolutions <- seq(0.005, 0.1, by = 0.005)
for (res in resolutions) {
  seurat_obj <- FindClusters(seurat_obj, resolution = res)
  print(paste("Resolution:", res, "Number of clusters:", length(unique(seurat_obj@meta.data$seurat_clusters))))
}

# Perform cclustering at the optimal resolution
seurat_obj <- FindClusters(seurat_obj, resolution = 0.075)


# ---- Run UMAP and Visualize Clusters----
seurat_obj <- RunUMAP(seurat_obj, reduction = "embeddings", dims = 1:1536)

DimPlot(
  seurat_obj,
  reduction = "umap",
  label = TRUE,
  raster = FALSE
) + ggtitle("UMAP of Clustering Results (Resolution = 0.075)") + 
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size=12) 
  )

# Print cell count per clusters 
table(seurat_obj$seurat_clusters)


# ---- Load Raw Matrix and Assign Clusters ----
# Load raw count matrix
raw_matrix <- readMM("matrix.mtx")
genes <- read.table("genes.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
barcodes <- read.table("barcodes.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
rownames(raw_matrix) <- genes$V1
colnames(raw_matrix) <- barcodes$V1

# Get cluster assignments
cluster_assignments <- data.frame(Cell = colnames(seurat_obj), Cluster = Idents(seurat_obj))

# Create new Seurat object using raw counts
seurat_raw <- CreateSeuratObject(counts = raw_matrix)

# Add cluster assignments to metadata
rownames(cluster_assignments) <- cluster_assignments$Cell
seurat_raw$seurat_clusters <- cluster_assignments[colnames(seurat_raw), "Cluster"]


# ---- Normalize Data and Identify Marker Genes ----
# Normalize expression matrix
seurat_raw <- NormalizeData(seurat_raw)

# Find variable features
seurat_raw <- FindVariableFeatures(seurat_raw)

# Identify marker genes
markers <- FindAllMarkers(seurat_raw, only.pos = TRUE, group.by = "seurat_clusters")
write.csv(markers, "concat_04-18_22clusters_all_markers.csv", row.names = FALSE)


# ---- Create Pseudo-Bulk Matrix and Annotate with SingleR ----
# Load previously saved marker genes
all_markers <- read.csv("concat_05-20_21clusters_all_markers.csv")

# Load reference dataset
ref <- celldex::BlueprintEncodeData()

# Create pseudo-bulk matrix
pseudo_bulk <- all_markers %>%
  select(cluster, gene) %>%
  distinct() %>%
  mutate(pseudo_expr = 1) %>%
  pivot_wider(names_from = cluster, values_from = pseudo_expr, values_fill = list(pseudo_expr = 0))

pseudo_bulk_mat <- as.matrix(as.data.frame(pseudo_bulk[-1]))
rownames(pseudo_bulk_mat) <- pseudo_bulk$gene

# Annotate clusters with SingleR
singleR_result <- SingleR(test = pseudo_bulk_mat, ref = ref, labels = ref$label.main)
table(singleR_result$labels)
cluster_annotations <- data.frame(cluster = colnames(pseudo_bulk_mat), predicted_cell_type = singleR_result$labels)
write.csv(cluster_annotations, "concat_04-18_22clusters_all_markers_pred_celltypes.csv", row.names = FALSE)


# ---- Summarize Predicted Cell Types ----
predicted <- read.csv("concat_05-20_21clusters_all_markers_pred_celltypes.csv")
cluster_sizes <- table(seurat_raw$seurat_clusters)
predicted$cluster <- as.character(predicted$cluster)
predicted$cluster_size <- cluster_sizes[predicted$cluster]

# Summarize cell type predictions
predicted_summary <- predicted %>%
  group_by(predicted_cell_type) %>%
  summarise(predicted_cells = sum(cluster_size)) %>%
  arrange(desc(predicted_cells))
print(predicted_summary)


# ---- Add Annotations to Metadata and Run UMAP ----
seurat_raw$predicted_cell_type <- predicted$predicted_cell_type[match(seurat_raw$seurat_clusters, predicted$cluster)]
seurat_raw[["embeddings"]] <- seurat_obj[["embeddings"]]
seurat_raw <- RunUMAP(seurat_raw, reduction = "embeddings", dims = 1:1536)

# ---- UMAP Visualization by Predicted Cell Types ----
custom_colors <- c(
  "CD4+ T-cells"       = "#ea7f5a", 
  "CD8+ T-cells"       = "#5cbd4d",  
  "Epithelial cells"   = "#cc77ab",  
  "Fibroblasts"        = "#f2cf63", 
  "Keratinocytes"      = "#936bd6",  
  "Macrophages"        = "#ee9a34",  
  "Monocytes"          = "#265734",  
  "NK cells"           = "#9dbfd0",  
  "Smooth muscle"      = "#018bb2",  
  "Endothelial cells"  = "#806e49", 
  "B-cells"            = "#6c0297",
  "NA"                 = "#787878",
  "Smooth muscle"      = "#591834")

DimPlot(
  seurat_raw,
  reduction = "umap",
  group.by = "predicted_cell_type",
  cols = custom_colors,
  raster = FALSE
) + ggtitle("UMAP of Predicted Cell Types from 21 Clusters") +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )

# ---- Cluster Frequency Table ----
table(seurat_raw$predicted_cell_type, useNA = "ifany")