# Single-cell mapping of extrachromosomal circular DNA architecture and evolution in cutaneous squamous cell carcinoma

![Image text](https://github.com/Leelab-Kmmu/Single-cell-mapping-of-extrachromosomal-circular-DNA/blob/main/workflow.png)

## This project provides users with the following functionalities:
* Process Cell Ranger scATAC-seq results into an eccDNA object (Seurat object), followed by dimensionality reduction, clustering, gene activity calculation, and cell type annotation.
* Quantify, for each cell, the number of fragments with copy number amplification, adapter (breakpoint) events, and eccDNA-associated genes. Visualize the distribution of eccDNA genomic regions and eccDNA localization in specific cells.
* Generate a data frame recording the copy number information of each eccDNA-related gene in each cell across all samples, and perform clone trajectory analysis.
* Infer cell lineages from single eccDNA molecules and validate them using chromatin accessibility similarities.

## Script discription

* Single_cell_copy_number_matrix.R: calculates the number of amplicons in each single cell.
* Breakpoint_identification.R: identifies breakpoint information for each single cell.
* eccDNA_formation.R: detects eccDNA information for each single cell.
* Figure.R: generates figures; requires the above three scripts to be run beforehand.

