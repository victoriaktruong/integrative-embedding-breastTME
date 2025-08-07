This repository contains all code, scripts, and results for the research project titled **Foundation Model-Based Framework for Cell Type Classification and Biomarker Discovery in Breast Cancer Single-Cell RNA Sequencing Data.**  

## Abstract  
_Background:_ Breast cancer is the most common cancer worldwide and a leading cause of mortality in women. Cellular heterogeneity within the tumour microenvironment (TME) complicates cell type identification and gene expression analysis, limiting the effectiveness of conventional methods. Foundation models offer a promising approach to overcome these analytical challenges.

_Objectives:_ This study compares three embedding strategies for single-cell RNA sequencing: scGPT (a single-cell foundation model), Mixedbread (a semantic language model), and scMIXE, which integrates both embeddings to improve cell type identification and marker gene discovery.

_Methods:_ Using data from 236,363 cells and 58,892 genes (Xu et al., 2017), embeddings from scGPT and Mixedbread were concatenated to form scMIXE. Clustering was performed using k-nearest neighbour graphs and Louvain modularity optimization. Marker genes were identified by differential expression, and clusters were annotated using the Blueprint/ENCODE reference. Gene set enrichment and ablation analyses compared the performance of all three embeddings.

_Results:_ scMIXE achieved the highest clustering resolution, most clusters and cell types, and fewest singletons. scMIXE and Mixedbread detected over 4,700 marker genes each, compared to 556 genes for scGPT. All cell types found by Mixedbread were captured by scMIXE, which uniquely detected adipocytes. CD8+ T cell enrichment showed scGPT prioritized specific immune pathways, while Mixedbread and scMIXE highlighted broader immune activity. Ablation analysis supported that scMIXE best captured cellular diversity.

_Conclusions:_ scMIXE embeddings improved clustering and annotation over scGPT and matched Mixedbread in gene detection, offering a more comprehensive view of the breast TME. These findings support the value of combining embeddings for advancing single-cell cancer analysis.

---

## Repository Structure

### `code/`  
Contains all scripts used to generate embeddings, perform clustering, identify markers, and annotate cell types. Scripts are listed in order of use:

- **`scGPT_embeddings.ipynb`**  
  Colab Python notebook used to generate scGPT embeddings and save them for downstream processing in R. 

- **`Mixedbread_embeddings.ipynb`**  
  Colab Python notebook used to generate Mixedbread embeddings.

- **`individual_embeddings_analysis.R`**  
   Pipeline for analyzing scGPT and Mixedbread individually to compare performance with the integrative approach.

- **`concatenated_embeddings_scMIXE.R`**  
  Complete analysis pipeline using concatenated scGPT and Mixedbread embeddings, including Seurat clustering, marker gene detection, and cell type annotation.

- **`marker_gene_filtering.R`**  
  Filters marker genes by log2FC, adjusted p-value, and expression threshold. Ensures each predicted cell type is represented by selecting top-ranked genes if none pass the filter.

- **`gene_set_enrichment_analysis.R`**  
  Performs GO enrichment analysis using Enrichr on filtered marker genes from each embedding strategy.

---

### `figures/`  
Visualizations including UMAP plots and project diagram.

---

### `results/`  
Stores all outputs from the analysis pipelines. Includes the following subfolders:

- **`scGPT/`** – Results from the scGPT embedding.
- **`Mixedbread/`** – Results from the Mixedbread embedding.
- **`scMIXE/`** – Outputs from the concatenated scGPT + Mixedbread embedding (scMIXE).
- **`final_filtered_markers/`** – Final set of ranked marker genes used for gene set enrichment analysis.

Additional files in this directory may include:
- Cluster assignments
- Marker gene tables 
- Predicted cell type labels

---

## Foundation Models Used
- [scGPT](https://github.com/bowang-lab/scGPT)
- [Mixedbread](https://huggingface.co/mxbai)

---

## Environment  
Code was developed and run in:
- R (version ≥ 4.2.0)
- Python 3.10+
- Google Colab Pro (L4 GPU for embedding extraction)

--- 

## Acknowledgements
Thank you to Dr. Pingzhao Hu for his guidance and support throughout this project. 
Special thanks to Emmett Peng for his Mixedbread contributions.

---

## Contact  
Victoria Truong 

LinkedIn: [@victoriaktruongg](https://www.linkedin.com/in/victoriaktruongg/)
GitHub: [@victoriaktruong](https://github.com/victoriaktruong)  
Email: victoria.truong@mail.utoronto.ca

