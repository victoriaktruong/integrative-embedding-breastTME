# ---- Load Required Libraries ----
library(arrow)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(SingleR)
library(celldex)
library(SummarizedExperiment)
library(tidyr)
library(readr)


# ---- Load Individual Embeddings and Create Seurat Object ----
# Load individual embedding matrix (either scGPT or Mixedbread, comment out the embeddings you do not want to use)
embeddings <- read_feather("cell_emb_HVG3000_768.feather")
embeddings <- read_feather("embeddings_mxbai.feather")
embeddings_matrix <- as.matrix(embeddings)

# Load cell barcodes and set rownames
barcodes <- read.table("barcodes.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
rownames(embeddings_matrix) <- barcodes$V1
embeddings_sparse <- as(embeddings_matrix, "dgCMatrix")

# Create Seurat object with embeddings as counts (transpose to match format)
seurat_obj <- CreateSeuratObject(counts = t(embeddings_sparse), project = "BreastCancer", assay = "RNA")
rownames(embeddings_matrix) <- colnames(seurat_obj)
colnames(embeddings_matrix) <- 1:ncol(embeddings_matrix)
seurat_obj[["embeddings"]] <- CreateDimReducObject(embeddings = embeddings_matrix, key = "embed_", assay = "RNA")


# ---- Perform Clustering using Louvain Algorithm ----
options(future.globals.maxSize = 10 * 1024^3)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "embeddings", dims = 1:768)

# Try a range of resolutions to determine the optimal resolution
resolutions <- seq(0.005, 0.1, by = 0.005)
for (res in resolutions) {
  seurat_obj <- FindClusters(seurat_obj, resolution = res)
  print(paste("Resolution:", res, "Number of clusters:", length(unique(seurat_obj$seurat_clusters))))
}

# Perform clustering at the optimal resolution
seurat_obj <- FindClusters(seurat_obj, resolution = 0.070)


# ---- Run UMAP and Visualize Clusters----
seurat_obj <- RunUMAP(seurat_obj, reduction = "embeddings", dims = 1:768)

DimPlot(
  seurat_obj,
  reduction = "umap",
  label = TRUE,
  raster = FALSE
) + ggtitle("UMAP of Clustering Results (Resolution = 0.070)") + 
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size=12) 
  )

# Print cell count per clusters 
table(seurat_obj$seurat_clusters)


# ---- Load Raw Matrix and Assign Clusters ----
# Load count matrix and corresponding gene/barcode info
raw_matrix <- readMM("matrix.mtx")
genes <- read.table("genes.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
barcodes <- read.table("barcodes.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
rownames(raw_matrix) <- genes$V1
colnames(raw_matrix) <- barcodes$V1


# Create cluster assignment metadata
cluster_assignments <- data.frame(Cell = colnames(seurat_obj), Cluster = Idents(seurat_obj))
write.csv(cluster_assignments, "cluster_assignments_mxbai.csv", row.names = FALSE)

# Create Seurat object from raw counts and add cluster info
seurat_raw <- CreateSeuratObject(counts = raw_matrix)
rownames(cluster_assignments) <- cluster_assignments$Cell
seurat_raw$seurat_clusters <- cluster_assignments[colnames(seurat_raw), "Cluster"]


# ---- Normalize Data and Identify Marker Genes ----
# Normalize expression matrix
seurat_raw <- NormalizeData(seurat_raw)

# Find variable features
seurat_raw <- FindVariableFeatures(seurat_raw)

# Identify marker genes
markers <- FindAllMarkers(seurat_raw, only.pos = TRUE, group.by = "seurat_clusters")
write.csv(markers, "mxbai_06-16_21clusters_all_markers.csv", row.names = FALSE)


# ---- Create Pseudo-Bulk Matrix and Annotate with SingleR ----
# Load previously saved marker genes
all_markers <- read.csv("mxbai_06-16_21clusters_all_markers.csv")

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
write.csv(cluster_annotations, "mxbai_06-16_21clusters_all_markers_pred_celltypes.csv", row.names = FALSE)


# ---- Summarize Predicted Cell Types ----
predicted <- read.csv("mxbai_06-16_21clusters_all_markers_pred_celltypes.csv")
cluster_sizes <- table(seurat_raw$seurat_clusters)
predicted$cluster <- as.character(predicted$cluster)
predicted$cluster_size <- cluster_sizes[predicted$cluster]

predicted_summary <- predicted %>%
  group_by(predicted_cell_type) %>%
  summarise(predicted_cells = sum(cluster_size)) %>%
  arrange(desc(predicted_cells))
print(predicted_summary)



# ---- Add Annotations and Visualize by Cell Type ----
seurat_raw$predicted_cell_type <- predicted$predicted_cell_type[match(seurat_raw$seurat_clusters, predicted$cluster)]
seurat_raw[["embeddings"]] <- seurat_obj[["embeddings"]]
seurat_raw <- RunUMAP(seurat_raw, reduction = "embeddings", dims = 1:768)

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
  label = TRUE,
  cols = custom_colors,
  raster = FALSE
) + ggtitle("UMAP of Predicted Cell Types from 10 Clusters") +
  theme(plot.title = element_text(size = 14, hjust = 0.5))


