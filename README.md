This repository contains all code, scripts, and results for the research project titled **Foundation Model-Based Framework for Cell Type Classification and Biomarker Discovery in Breast Cancer Single-Cell RNA Sequencing Data.**  

## Abstract  
*To be added.*

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
Visual outputs including UMAP plots and project diagram.

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

