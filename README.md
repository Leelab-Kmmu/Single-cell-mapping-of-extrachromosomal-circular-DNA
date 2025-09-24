# Single-cell mapping of extrachromosomal circular DNA architecture and evolution in cutaneous squamous cell carcinoma

![Image text](https://github.com/Leelab-Kmmu/Single-cell-mapping-of-extrachromosomal-circular-DNA/main/workflow.png)

## This project provides users with the following functionalities:
Process Cell Ranger results of scATAC data into an eccDNA object(Seurat object), and complete the processes of dimensionality reduction, clustering, gene activity calculation, and cell annotation.
Calculate the number of fragments with copy number amplification, adapter(breakpoint), eccDNA genes in each cell. Plot the distribution of eccDNA genomic regions and Visualizing eccDNA in specific cell.
Create a data frame to record the copy number information of each eccDNA-related gene in each cell across all samples and perform clone trajectory analysis.
Determin cell lineages by single eccDNA molecules and corroborat by chromatin accessibility similarities.

## Dependiencies

* Python 3.0

* [BWA](https://github.com/lh3/bwa)

* [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html)
* [Cell Ranger ATAC v1.2](https://support.10xgenomics.com/single-cell-atac/software/downloads/1.2/)
* [Sinto](https://timoast.github.io/sinto/index.html)
* [Signac](https://github.com/stuart-lab/signac)
* [ArchR](https://github.com/GreenleafLab/ArchR)
