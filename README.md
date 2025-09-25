# Single-cell mapping of extrachromosomal circular DNA architecture and evolution in cutaneous squamous cell carcinoma

![Image text](https://github.com/Leelab-Kmmu/Single-cell-mapping-of-extrachromosomal-circular-DNA/blob/main/workflow.png)

## This project provides users with the following functionalities:
* Process Cell Ranger results of scATAC data into an eccDNA object(Seurat object), and complete the processes of dimensionality reduction, clustering, gene activity calculation, and cell annotation.
* Calculate the number of fragments with copy number amplification, adapter(breakpoint), eccDNA genes in each cell. Plot the distribution of eccDNA genomic regions and Visualizing eccDNA in specific cell.
* Create a data frame to record the copy number information of each eccDNA-related gene in each cell across all samples and perform clone trajectory analysis.
* Determin cell lineages by single eccDNA molecules and corroborat by chromatin accessibility similarities.

## Script discription

* Single_cell_copy_number_matrix.R provides the amplicons number of each single cell.
* Breakpoint_identification.R identifies the breakpoint information of each single cell.
* eccDNA_formation.R identifies the eccDNA information of each single cell.
* The above three scripts must be running before running Figure.R script.

