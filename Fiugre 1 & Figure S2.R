options(stringsAsFactors=F)
setwd('/data2/home/lijie/CSCC/result/all_cells')
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(future)
# read in peak sets
Nomal1 <- read.table(
  file = "/data2/home/lijie/CSCC/cell_ranger/Nomal1/filtered_peak_bc_matrix/peaks.bed",
  col.names = c("chr", "start", "end")
)
Nomal2 <- read.table(
  file = "/data2/home/lijie/CSCC/cell_ranger/Nomal2/filtered_peak_bc_matrix/peaks.bed",
  col.names = c("chr", "start", "end")
)
Nomal3 <- read.table(
  file = "/data2/home/lijie/CSCC/cell_ranger/Nomal3/filtered_peak_bc_matrix/peaks.bed",
  col.names = c("chr", "start", "end")
)
AK1 <- read.table(
  file = "/data2/home/lijie/CSCC/cell_ranger/AK1/outs/filtered_peak_bc_matrix/peaks.bed",
  col.names = c("chr", "start", "end")
)
AK2 <- read.table(
  file = "/data2/home/lijie/CSCC/cell_ranger/AK2/outs/filtered_peak_bc_matrix/peaks.bed",
  col.names = c("chr", "start", "end")
)
AK3 <- read.table(
  file = "/data2/home/lijie/CSCC/cell_ranger/AK3/outs/filtered_peak_bc_matrix/peaks.bed",
  col.names = c("chr", "start", "end")
)
CSCC1 <- read.table(
  file = "/data2/home/lijie/CSCC/cell_ranger/CSCC1/outs/filtered_peak_bc_matrix/peaks.bed",
  col.names = c("chr", "start", "end")
)
CSCC2 <- read.table(
  file = "/data2/home/lijie/CSCC/cell_ranger/CSCC2/outs/filtered_peak_bc_matrix/peaks.bed",
  col.names = c("chr", "start", "end")
)
CSCC3 <- read.table(
  file = "/data2/home/lijie/CSCC/cell_ranger/CSCC3/outs/filtered_peak_bc_matrix/peaks.bed",
  col.names = c("chr", "start", "end")
)

# convert to genomic ranges
gr.Nomal1 <- makeGRangesFromDataFrame(Nomal1)
gr.Nomal2 <- makeGRangesFromDataFrame(Nomal2)
gr.Nomal3 <- makeGRangesFromDataFrame(Nomal3)
gr.AK1 <- makeGRangesFromDataFrame(AK1)
gr.AK2 <- makeGRangesFromDataFrame(AK2)
gr.AK3 <- makeGRangesFromDataFrame(AK3)
gr.CSCC1 <- makeGRangesFromDataFrame(CSCC1)
gr.CSCC2 <- makeGRangesFromDataFrame(CSCC2)
gr.CSCC3 <- makeGRangesFromDataFrame(CSCC3)


# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.Nomal1, gr.Nomal2, gr.Nomal3, gr.AK1,gr.AK2,gr.AK3,gr.CSCC1,gr.CSCC2,gr.CSCC3))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
# load metadata
md.Nomal1 <- read.table(
  file = "/data2/home/lijie/CSCC/cell_ranger/Nomal1/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row
md.Nomal2 <- read.table(
  file = "/data2/home/lijie/CSCC/cell_ranger/Nomal2/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row
md.Nomal3 <- read.table(
  file = "/data2/home/lijie/CSCC/cell_ranger/Nomal3/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

md.AK1 <- read.table(
  file = "/data2/home/lijie/CSCC/cell_ranger/AK1/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row
md.AK2 <- read.table(
  file = "/data2/home/lijie/CSCC/cell_ranger/AK2/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row
md.AK3 <- read.table(
  file = "/data2/home/lijie/CSCC/cell_ranger/AK3/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row
md.CSCC1 <- read.table(
  file = "/data2/home/lijie/CSCC/cell_ranger/CSCC1/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row
md.CSCC2 <- read.table(
  file = "/data2/home/lijie/CSCC/cell_ranger/CSCC2/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row
md.CSCC3 <- read.table(
  file = "/data2/home/lijie/CSCC/cell_ranger/CSCC3/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row



# perform an initial filtering of low count cells
md.Nomal1 <- md.Nomal1[md.Nomal1$passed_filters > 500, ]
md.Nomal2 <- md.Nomal2[md.Nomal2$passed_filters > 500, ]
md.Nomal3 <- md.Nomal3[md.Nomal3$passed_filters > 500, ]
md.AK1 <- md.AK1[md.AK1$passed_filters > 500, ]
md.AK2 <- md.AK2[md.AK2$passed_filters > 500, ]
md.AK3 <- md.AK3[md.AK3$passed_filters > 500, ]
md.CSCC1 <- md.CSCC1[md.CSCC1$passed_filters > 500, ]
md.CSCC2 <- md.CSCC2[md.CSCC2$passed_filters > 500, ]
md.CSCC3 <- md.CSCC3[md.CSCC3$passed_filters > 500, ]

# create fragment objects
frags.Nomal1 <- CreateFragmentObject(
  path = "/data2/home/lijie/CSCC/cell_ranger/Nomal1/Nomal1_fragments.tsv.gz",
  cells = rownames(md.Nomal1)
)
frags.Nomal2 <- CreateFragmentObject(
  path = "/data2/home/lijie/CSCC/cell_ranger/Nomal2/Nomal2_fragments.tsv.gz",
  cells = rownames(md.Nomal2)
)
frags.Nomal3 <- CreateFragmentObject(
  path = "/data2/home/lijie/CSCC/cell_ranger/Nomal3/Nomal3_fragments.tsv.gz",
  cells = rownames(md.Nomal3)
)
frags.AK1 <- CreateFragmentObject(
  path = "/data2/home/lijie/CSCC/cell_ranger/AK1/outs/AK1_fragments.tsv.gz",
  cells = rownames(md.AK1)
)
frags.AK2 <- CreateFragmentObject(
  path = "/data2/home/lijie/CSCC/cell_ranger/AK2/outs/AK2_fragments.tsv.gz",
  cells = rownames(md.AK2)
)
frags.AK3 <- CreateFragmentObject(
  path = "/data2/home/lijie/CSCC/cell_ranger/AK3/outs/AK3_fragments.tsv.gz",
  cells = rownames(md.AK3)
)
frags.CSCC1 <- CreateFragmentObject(
  path = "/data2/home/lijie/CSCC/cell_ranger/CSCC1/outs/CSCC1_fragments.tsv.gz",
  cells = rownames(md.CSCC1)
)
frags.CSCC2 <- CreateFragmentObject(
  path = "/data2/home/lijie/CSCC/cell_ranger/CSCC2/outs/CSCC2_fragments.tsv.gz",
  cells = rownames(md.CSCC2)
)
frags.CSCC3 <- CreateFragmentObject(
  path = "/data2/home/lijie/CSCC/cell_ranger/CSCC3/outs/CSCC3_fragments.tsv.gz",
  cells = rownames(md.CSCC3)
)

Nomal1.counts <- FeatureMatrix(
  fragments = frags.Nomal1,
  features = combined.peaks,
  cells = rownames(md.Nomal1)
)
Nomal2.counts <- FeatureMatrix(
  fragments = frags.Nomal2,
  features = combined.peaks,
  cells = rownames(md.Nomal2)
)
Nomal3.counts <- FeatureMatrix(
  fragments = frags.Nomal3,
  features = combined.peaks,
  cells = rownames(md.Nomal3)
)
AK1.counts <- FeatureMatrix(
  fragments = frags.AK1,
  features = combined.peaks,
  cells = rownames(md.AK1)
)
AK2.counts <- FeatureMatrix(
  fragments = frags.AK2,
  features = combined.peaks,
  cells = rownames(md.AK2)
)
AK3.counts <- FeatureMatrix(
  fragments = frags.AK3,
  features = combined.peaks,
  cells = rownames(md.AK3)
)
CSCC1.counts <- FeatureMatrix(
  fragments = frags.CSCC1,
  features = combined.peaks,
  cells = rownames(md.CSCC1)
)
CSCC2.counts <- FeatureMatrix(
  fragments = frags.CSCC2,
  features = combined.peaks,
  cells = rownames(md.CSCC2)
)
CSCC3.counts <- FeatureMatrix(
  fragments = frags.CSCC3,
  features = combined.peaks,
  cells = rownames(md.CSCC3)
)
Nomal1_assay <- CreateChromatinAssay(Nomal1.counts, fragments = frags.Nomal1)
Nomal1 <- CreateSeuratObject(Nomal1_assay, assay = "ATAC", meta.data=md.Nomal1)

Nomal2_assay <- CreateChromatinAssay(Nomal2.counts, fragments = frags.Nomal2)
Nomal2 <- CreateSeuratObject(Nomal2_assay, assay = "ATAC", meta.data=md.Nomal2)

Nomal3_assay <- CreateChromatinAssay(Nomal3.counts, fragments = frags.Nomal3)
Nomal3 <- CreateSeuratObject(Nomal3_assay, assay = "ATAC", meta.data=md.Nomal3)

AK1_assay <- CreateChromatinAssay(AK1.counts, fragments = frags.AK1)
AK1 <- CreateSeuratObject(AK1_assay, assay = "ATAC", meta.data=md.AK1)

AK2_assay <- CreateChromatinAssay(AK2.counts, fragments = frags.AK2)
AK2 <- CreateSeuratObject(AK2_assay, assay = "ATAC", meta.data=md.AK2)

AK3_assay <- CreateChromatinAssay(AK3.counts, fragments = frags.AK3)
AK3 <- CreateSeuratObject(AK3_assay, assay = "ATAC", meta.data=md.AK3)

CSCC1_assay <- CreateChromatinAssay(CSCC1.counts, fragments = frags.CSCC1)
CSCC1 <- CreateSeuratObject(CSCC1_assay, assay = "ATAC", meta.data=md.CSCC1)

CSCC2_assay <- CreateChromatinAssay(CSCC2.counts, fragments = frags.CSCC2)
CSCC2 <- CreateSeuratObject(CSCC2_assay, assay = "ATAC", meta.data=md.CSCC2)

CSCC3_assay <- CreateChromatinAssay(CSCC3.counts, fragments = frags.CSCC3)
CSCC3 <- CreateSeuratObject(CSCC3_assay, assay = "ATAC", meta.data=md.CSCC3)
# add information to identify dataset of origin
Nomal1$dataset <- 'Nomal1'
Nomal2$dataset <- 'Nomal2'
Nomal3$dataset <- 'Nomal3'
AK1$dataset <- 'AK1'
AK2$dataset <- 'AK2'
AK3$dataset <- 'AK3'
CSCC1$dataset <- 'CSCC1'
CSCC2$dataset <- 'CSCC2'
CSCC3$dataset <- 'CSCC3'

# merge all datasets, adding a cell ID to make sure cell names are unique
combined_atac <- merge(
  x = Nomal1,
  y = list(Nomal2,Nomal3,AK1,AK2,AK3,CSCC1,CSCC2,CSCC3),
  add.cell.ids = c("Nomal1","Nomal2","Nomal3","AK1","AK2","AK3","CSCC1","CSCC2","CSCC3")
)
combined_atac[["ATAC"]]
save(combined_atac,file='combine_atac.Rdata')



######################################################################################################################################
# BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86'))
library(Seurat)
library(Signac)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)  # Human hg38
library(ggplot2)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(motifmatchr)

rm(list=ls())
options(stringsAsFactors = FALSE)
setwd('all_cells')
##
plan("multicore", workers = 1)
options(future.globals.maxSize = 40000 * 1024^2) # for 40 Gb RAM

load("combine_atac.Rdata")
View(combined_atac@meta.data)

######################################################################################################################################


# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg38
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(combined_atac) <- annotations

##########################################################################################################################################################

combined_atac <- NucleosomeSignal(object = combined_atac)
combined_atac <- TSSEnrichment(object = combined_atac, fast = FALSE)
combined_atac$pct_reads_in_peaks <- combined_atac$peak_region_fragments / combined_atac$passed_filters * 100
combined_atac$blacklist_ratio <- combined_atac$blacklist_region_fragments / combined_atac$peak_region_fragments

save(combined_atac, file = "combined_atac_preQC.Rdata")

# For Figure 1B & Figure S2A
DensityScatter(combined_atac, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

combined_atac$high.tss <- ifelse(combined_atac$TSS.enrichment > 3, 'TSS > 3', 'TSS < 3')
# TSS < 3    TSS > 3 
# 458637     99786 
TSSPlot(combined_atac, group.by = 'high.tss') + NoLegend()

combined_atac$nucleosome_group <- ifelse(combined_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = combined_atac, group.by = 'nucleosome_group')

VlnPlot(
  object = combined_atac,
  features = c('nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0,
  ncol = 5,
  group.by = "orig.ident"
)

summary(combined_atac$nCount_ATAC)
# Min. 1st Qu.  Median  Mea     3rd Qu.    Max. 
# 53     237     441    1461     713      275507

summary(combined_atac$TSS.enrichment)
# Min. 1st Qu.  Median    Mean    3rd Qu.    Max. 
# 0.000   1.182   1.638   2.095   2.464     30.170 

summary(combined_atac$blacklist_ratio) # 全是0

summary(combined_atac$nucleosome_signal)
# Min.   1st Qu.    Median      Mean      3rd Qu.      Max. 
# 0.04301   0.58883   0.87865   1.07810   1.39323     111.11741 

summary(combined_atac$pct_reads_in_peaks)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6969  3.9308  7.2262  9.7088 10.2415 82.5540 

###########################################################################

###after QC
summary(combined_atac$nCount_ATAC)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1001    3531    6861    9339   13426   29995 

summary(combined_atac$TSS.enrichment)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.000   4.327   4.939   5.051   5.654  15.105 

summary(combined_atac$blacklist_ratio) # 全是0

summary(combined_atac$nucleosome_signal)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0955  0.5645  0.6814  0.7681  0.8811  3.9815 

summary(combined_atac$pct_reads_in_peaks)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 15.00   29.73   48.05   46.44   62.54   82.55  

###########################################################################

rm(list=ls())
options(stringsAsFactors = FALSE)
library(Seurat)
library(Signac)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)  # Human hg38
library(ggplot2)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(motifmatchr)
library(ggseqlogo)
setwd('/data2/home/lijie/CSCC/result/all_cells')
load("combined_atac_preQC.Rdata")
combined_atac@assays[["ATAC"]]@fragments[[1]]@path<-"/data2/home/lijie/CSCC/cell_ranger/Nomal1/Nomal1_fragments.tsv.gz"
combined_atac@assays[["ATAC"]]@fragments[[2]]@path<-"/data2/home/lijie/CSCC/cell_ranger/Nomal2/Nomal2_fragments.tsv.gz"
combined_atac@assays[["ATAC"]]@fragments[[3]]@path<-"/data2/home/lijie/CSCC/cell_ranger/Nomal3/Nomal3_fragments.tsv.gz"
combined_atac@assays[["ATAC"]]@fragments[[4]]@path<-"/data2/home/lijie/CSCC/cell_ranger/AK1/outs/AK1_fragments.tsv.gz"
combined_atac@assays[["ATAC"]]@fragments[[5]]@path<-"/data2/home/lijie/CSCC/cell_ranger/AK2/outs/AK2_fragments.tsv.gz"
combined_atac@assays[["ATAC"]]@fragments[[6]]@path<-"/data2/home/lijie/CSCC/cell_ranger/AK3/outs/AK3_fragments.tsv.gz"
combined_atac@assays[["ATAC"]]@fragments[[7]]@path<-"/data2/home/lijie/CSCC/cell_ranger/CSCC1/outs/CSCC1_fragments.tsv.gz"
combined_atac@assays[["ATAC"]]@fragments[[8]]@path<-"/data2/home/lijie/CSCC/cell_ranger/CSCC2/outs/CSCC2_fragments.tsv.gz"
combined_atac@assays[["ATAC"]]@fragments[[9]]@path<-"/data2/home/lijie/CSCC/cell_ranger/CSCC3/outs/CSCC3_fragments.tsv.gz"



combined_atac <- subset(
  x = combined_atac,
  subset = nCount_ATAC > 1000 &  #
    nCount_ATAC < 30000 &  #
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3
)

table(combined_atac$dataset)
# AK1    AK2    AK3  CSCC1  CSCC2  CSCC3 Nomal1 Nomal2 Nomal3 
# 1189   6958   2816   1518   2441   3054   1314   4986   5526 

table(combined_atac@meta.data$ATAC_snn_res.0.3)
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17 
# 4272 3915 3287 2703 2648 2362 1924 1890 1187 1095 1045  901  850  522  452  368  219  162 

######################################################################################################################################
#
combined_atac <- RunTFIDF(combined_atac) 
combined_atac <- FindTopFeatures(combined_atac, min.cutoff = 10) 
combined_atac <- RunSVD(combined_atac)
DepthCor(combined_atac)
combined_atac <- RunUMAP(combined_atac, dims = 2:30, reduction = 'lsi')
DimPlot(combined_atac, group.by = 'dataset', pt.size = 0.1, raster=FALSE)

######################################################################################################################################
#For Figure S2B

DefaultAssay(combined_atac) <- 'ATAC'
combined_atac <- FindNeighbors(object = combined_atac, reduction = 'lsi', dims = 2:30)
combined_atac <- FindClusters(object = combined_atac, verbose = T, algorithm = 3, resolution = 0.3)
DimPlot(object = combined_atac, label = TRUE,raster=FALSE)
View(combined_atac@meta.data)
table(combined_atac@meta.data$ATAC_snn_res.0.3)

######################################################################################################################################
gene.activities <- GeneActivity(combined_atac)
# combined_atac@assays$peaks$counts[1:5,1:5]
# combined_atac@assays$RNA@counts [1:5,1:5]

# add the gene activity matrix to the Seurat object as a new assay and normalize it
combined_atac[['RNA']] <- CreateAssayObject(counts = gene.activities)
combined_atac <- NormalizeData(
  object = combined_atac,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined_atac$nCount_RNA)
)

######################################################################################################################################
#For Figure 1C &Fiure S2C, cell annotation
DefaultAssay(combined_atac) <- 'RNA'

celltype_marker<-c(                                                                                                         
  "COL17A1", "KRT14","KRT5" ,                                # basal                                        ### cluster1、2、4
  "SKAP1", "CD2", "CD3D" ,                                   # T cells                                      ### cluster7、
  "COL1A2", "DCN" ,"CCDC80", "COL3A1", "COL6A1", "COL1A1",   # mFib                                         ### cluster6、
  "KRT1", "KRT10" ,"KRTDAP",                                 # spinosum                                     ### cluster0、
  "S100A7", "S100A8", "S100A9",                              # granulosum                                   ### cluster3、
    "KRT17", "ADGRL3", "TENM2", "COL17A1", "KRT14","KRT5" ,    #                    ### cluster9、12
  "CD68", "AIF1" ,"C1QB" ,"C1QA" ,                           #     ### cluster8
  "AQP1", "ADGRL4", "VWF" ,                                  # Endothelial cells                            ### cluster5、
  "TAGLN", "ACTA2", "RGS5",                                  # CAFs                                         ### cluster10、17
  "COL6A1", "COL1A1", "COL3A1",  "COL1A2", "DCN" ,"CCDC80",  # mesFib                                       ###     
  "CD79A", "BANK1", "MS4A1", "CD19","CD79B",                 # B cells                    ### 
  "MS4A1", "IGKC", "IGHG1", "IGHA1",                         # Plasma cells1             ### 
  "IGHG1", "IGHGP", "IGLC2", "IGLC3", "IGHA1" ,              # Plasma cells2                    ### 
  "TPSB2", "TPSAB1", "CPA3" ,                                # Mast cells                                   ### 
  "DCT", "TYRP1", "PMEL"   )                                 # melanocyte                                   ### cluster14、


DotPlot(combined_atac, features = unique(celltype_marker),cluster.idents=T,dot.min = 0.1,)+  # col.min = 0,col.max = 1,
  coord_flip()

View(combined_atac@meta.data)
combined_atac$celltype_ATAC<-c(rep(0,nrow(combined_atac)))
for(i in 1:nrow(combined_atac@meta.data)){
  if(combined_atac@meta.data[i,33]%in%c("11","13","15","16")){combined_atac@meta.data[i,34]<-"undefined"}
  if(combined_atac@meta.data[i,33]%in%c("1","2","4")){combined_atac@meta.data[i,34]<-"basal"}
  if(combined_atac@meta.data[i,33]=="7"){combined_atac@meta.data[i,34]<-"T cells"}
  if(combined_atac@meta.data[i,33]=="6"){combined_atac@meta.data[i,34]<-"mFib"}
  if(combined_atac@meta.data[i,33]=="0"){combined_atac@meta.data[i,34]<-"spinosum"}
  if(combined_atac@meta.data[i,33]=="3"){combined_atac@meta.data[i,34]<-"granulosum"}
  if(combined_atac@meta.data[i,33]%in%c("9","12")){combined_atac@meta.data[i,34]<-"follicular"}
  if(combined_atac@meta.data[i,33]=="8"){combined_atac@meta.data[i,34]<-"APC"}
  if(combined_atac@meta.data[i,33]=="5"){combined_atac@meta.data[i,34]<-"Endothelial cells"}
  if(combined_atac@meta.data[i,33]%in%c("10","17")){combined_atac@meta.data[i,34]<-"CAFs"}
  if(combined_atac@meta.data[i,33]=="14"){combined_atac@meta.data[i,34]<-"melanocyte"}
}
DimPlot(combined_atac, group.by = 'celltype_ATAC', label = F,label.size = 4 ,pt.size = 0.5, raster=FALSE,cols=c('APC'='#1AB4B8','basal'='#E0925F','CAFs'='#E06AA4',
'Endothelial cells'='#DBD970','follicular'='#BF8CB3','granulosum'='#FC4F68','melanocyte'="#FF8C00",'mFib'="#7DC6AC",'spinosum'="#E9989A",'T cells'="#8B3A3A"))
table(combined_atac$celltype_ATAC)

######################################################################################################################################

saveRDS(combined_atac, file = file.path(files_dir,"output_QC3/combined_atac.rds"))
combined_atac <- readRDS(file.path(files_dir,"output_QC3/combined_atac.rds"))


combined_atac$dataset[which(combined_atac$dataset == "Nomal1")] <- "Normal1"
combined_atac$dataset[which(combined_atac$dataset == "Nomal2")] <- "Normal2"
combined_atac$dataset[which(combined_atac$dataset == "Nomal3")] <- "Normal3"
######################################################################################################################################

###For Figure S2D
combined_atac <- SetIdent(combined_atac, value = "celltype_ATAC")
current_idents <- unique(Idents(combined_atac))
combined_atac <- subset(combined_atac, idents = current_idents[current_idents != 'undefined'])

combined_atac$sampletype <- ifelse(combined_atac$dataset%in%c("AK1","AK2","AK3"),"AK",
                                   ifelse(combined_atac$dataset%in%c("CSCC1","CSCC2","CSCC3"),"CSCC","Normal"))

combined_atac$patient <- ifelse(combined_atac$dataset%in%c("AK1","CSCC1","Normal1"),"patient1",
                                ifelse(combined_atac$dataset%in%c("AK2","CSCC2","Normal2"),"patient2","patient3"))

# celltype_sample_distribution <- as.data.frame.matrix(table(combined_atac@meta.data$celltype_ATAC,combined_atac@meta.data$dataset))
celltype_sample_distribution <- as.data.frame(table(combined_atac@meta.data$celltype_ATAC,combined_atac@meta.data$dataset))
colnames(celltype_sample_distribution) <- c("celltype","sample","cellnumber")

celltype_sampletype_distribution <- as.data.frame(table(combined_atac@meta.data$celltype_ATAC,combined_atac@meta.data$sampletype))
colnames(celltype_sampletype_distribution) <- c("celltype","sampletype","cellnumber")

celltype_patient_distribution <- as.data.frame(table(combined_atac@meta.data$celltype_ATAC,combined_atac@meta.data$patient))
colnames(celltype_patient_distribution) <- c("celltype","patient","cellnumber")
  
windowsFonts(A=windowsFont("Times New Roman"),
             B=windowsFont("Arial"))

celltype_sample_distribution$sample <- factor(celltype_sample_distribution$sample,
            levels = c("Normal1","Normal2","Normal3","AK1","AK2","AK3","CSCC1","CSCC2","CSCC3"))

celltype_sampletype_distribution$sampletype <- factor(celltype_sampletype_distribution$sampletype,
                       levels = c("Normal","AK","CSCC"))

library(ggplot2)
p1 <- ggplot(celltype_sample_distribution, aes(celltype, cellnumber), position = "stack") +
      geom_bar(aes(fill = sample), stat = "identity", color = "black", size = 0.2,
               position = "fill", width = 0.65) +
      scale_fill_manual(values = c("#56B4E9", 'gray', '#CCEBC5', '#BC80BD', '#FCCDE5',
                                   '#BEBADA', '#FDB462',  '#FB8072','#8DD3C7',
                                   '#FFFFB3','#B3DE69' )) +
      # scale_fill_manual(values = c("#56B4E9", '#BC80BD', '#FDB462')) +
      # scale_fill_manual(values =  c("#D8B5A5", '#B3DE69', '#BEBADA')) +
      theme(
        legend.position = "right",
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.4),
        axis.text.x = element_text(angle = 45, hjust = 1,margin = margin(t = 3, r = 0, b = 0, l = -3),colour = "black"),  # 调整 x 轴标签与坐标轴之间的距离
        text = element_text(family = "A", size = 18),
        axis.ticks.length = unit(-0.25, "cm"),
        axis.text.y = element_text(margin = margin(t = 0, r = 3, b = 0, l = 2),colour = "black"),  # 调整y轴标签与坐标轴之间的距离
        legend.text = element_text(margin = margin(t = 0, r = 5, b = 0, l = -5),face = "plain"),  # 调整图例文本与颜色块之间的距离
        legend.margin = margin(t = 0, r = 0, b = 0, l = -3)  # 调整整个图例的外边距
      ) +
      labs(x = NULL, y = NULL)

p2 <- ggplot(celltype_sampletype_distribution, aes(celltype, cellnumber), position = "stack") +
  geom_bar(aes(fill = sampletype), stat = "identity", color = "black", size = 0.2,
           position = "fill", width = 0.65) +
  # scale_fill_manual(values = c("#56B4E9", 'gray', '#CCEBC5', '#BC80BD', '#FCCDE5',
  #                              '#BEBADA', '#FDB462',  '#FB8072','#8DD3C7',
  #                              '#FFFFB3','#B3DE69' )) +
  scale_fill_manual(values = c("#56B4E9", '#BC80BD', '#FDB462')) +
  # scale_fill_manual(values =  c("#D8B5A5", '#B3DE69', '#BEBADA')) +
  theme(
    legend.position = "right",
    panel.background = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.4),
    axis.text.x = element_text(angle = 45, hjust = 1,margin = margin(t = 3, r = 0, b = 0, l = -3),colour = "black"),  # 调整 x 轴标签与坐标轴之间的距离
    text = element_text(family = "A", size = 18),
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.y = element_text(margin = margin(t = 0, r = 3, b = 0, l = 2),colour = "black"),  # 调整y轴标签与坐标轴之间的距离
    legend.text = element_text(margin = margin(t = 0, r = 5, b = 0, l = -5),face = "plain"),  # 调整图例文本与颜色块之间的距离
    legend.margin = margin(t = 0, r = 0, b = 0, l = -3)  # 调整整个图例的外边距
  ) +
  labs(x = NULL, y = NULL)

p3 <- ggplot(celltype_patient_distribution, aes(celltype, cellnumber), position = "stack") +
  geom_bar(aes(fill = patient), stat = "identity", color = "black", size = 0.2,
           position = "fill", width = 0.65) +
  # scale_fill_manual(values = c("#56B4E9", 'gray', '#CCEBC5', '#BC80BD', '#FCCDE5',
  #                              '#BEBADA', '#FDB462',  '#FB8072','#8DD3C7',
  #                              '#FFFFB3','#B3DE69' )) +
  # scale_fill_manual(values = c("#56B4E9", '#BC80BD', '#FDB462')) +
  scale_fill_manual(values =  c("#D8B5A5", '#B3DE69', '#BEBADA')) +
  theme(
    legend.position = "right",
    panel.background = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.4),
    axis.text.x = element_text(angle = 45, hjust = 1,margin = margin(t = 3, r = 0, b = 0, l = -3),colour = "black"),  # 调整 x 轴标签与坐标轴之间的距离
    text = element_text(family = "A", size = 18),
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.y = element_text(margin = margin(t = 0, r = 3, b = 0, l = 2),colour = "black"),  # 调整y轴标签与坐标轴之间的距离
    legend.text = element_text(margin = margin(t = 0, r = 5, b = 0, l = -5),face = "plain"),  # 调整图例文本与颜色块之间的距离
    legend.margin = margin(t = 0, r = 0, b = 0, l = -3)  # 调整整个图例的外边距
  ) +
  labs(x = NULL, y = NULL)

p1+p2+p3



#For Figure 1D
#####################################CNV distribution, adapter distribution and eccDNA distribution
options(stringsAsFactors = FALSE)
library(Seurat)
library(Signac)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)  # Human hg38
library(ggplot2)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(motifmatchr)
library(ggseqlogo)
setwd('/data2/home/lijie/CSCC/result/all_cells')
load('/data2/home/lijie/scATAC_eccDNA/results/CSCC3/copy/cell_type.Rdata')
load("combined_atac_final.Rdata")

new_combined_atac <- RunUMAP(combined_atac, dims = 2:10, reduction = 'lsi')
DimPlot(new_combined_atac, group.by = 'dataset', label = F,label.size = 4 ,pt.size = 0.5, raster=FALSE,cols=c('Normal1'='#8B3A3A','Normal2'='#E0925F','Normal3'='#E06AA4',
'AK1'='#DBD970','AK2'='#BF8CB3','AK3'='#FC4F68','CSCC1'="#FF8C00",'CSCC2'="#1E90FF",'CSCC3'="#E9989A"))
DimPlot(new_combined_atac, group.by = 'celltype_ATAC', label = F,label.size = 4 ,pt.size = 0.5, raster=FALSE,cols=c('APC'='#1AB4B8','basal'='#E0925F','CAFs'='#E06AA4',
'Endothelial cells'='#DBD970','follicular'='#BF8CB3','granulosum'='#FC4F68','melanocyte'="#FF8C00",'mFib'="#7DC6AC",'spinosum'="#E9989A",'T cells'="#8B3A3A"))+NoLegend()


load('/data2/home/lijie/scATAC_eccDNA/results/eccDNA_features/CNV_statistic.Rdata')##CNV distribution
CNV_allcopy<-CNV_allcopy[names(cell_type)]
CNV_allcopy[is.na(CNV_allcopy)]<-0
new_combined_atac@meta.data[["CNV_copy"]]<-CNV_allcopy
length(CNV_allcopy[CNV_allcopy!=0])
summary(CNV_allcopy[CNV_allcopy!=0])
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   # 5.00   18.00   36.00   53.58   70.00  532.00
load('/data2/home/lijie/scATAC_eccDNA/results/eccDNA_features/adapter_num.Rdata')##adpter distribution
adapter_num<-adapter_num[names(cell_type)]
adapter_num[is.na(adapter_num)]<-0
new_combined_atac@meta.data[["adapter_num"]]<-adapter_num
length(adapter_num[adapter_num!=0])
summary(adapter_num[adapter_num!=0])
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 1.000   2.000   2.000   2.679   2.750  18.000
load('/data2/home/lijie/scATAC_eccDNA/results/eccDNA_features/eccDNA_gene_num.Rdata')
eccDNA_gene_num<-eccDNA_gene_num[names(cell_type)]
eccDNA_gene_num[is.na(eccDNA_gene_num)]<-0
new_combined_atac@meta.data[["eccDNA_gene_num"]]<-eccDNA_gene_num
length(eccDNA_gene_num[eccDNA_gene_num!=0])
summary(eccDNA_gene_num[eccDNA_gene_num!=0])
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 1.000   1.000   1.000   1.153   1.000   9.000

DimPlot(new_combined_atac, group.by = 'dataset', pt.size = 0.1, raster=FALSE)
FeaturePlot(new_combined_atac, features = 'CNV_copy',label.size = 4 ,pt.size = 0.5,max.cutoff=200,cols=c("grey", "red"),order =T, raster=FALSE )+ NoLegend()
FeaturePlot(new_combined_atac, features = 'adapter_num',label.size = 4 ,pt.size = 0.5,max.cutoff=3,cols=c("grey", "red"),order =T, raster=FALSE )+ NoLegend()
FeaturePlot(new_combined_atac, features = 'eccDNA_gene_num',label.size = 4 ,pt.size = 0.5,max.cutoff=1,cols=c("grey", "red"),order =T, raster=FALSE )+ NoLegend()

  

######CNV总量
celltype_CNV_distribution <- data.frame(cell_type=new_combined_atac@meta.data$celltype_ATAC,CNV=new_combined_atac@meta.data$CNV_copy)
celltype_CNV<-data.frame(cell_type=c('basal','granulosum','spinosum','follicular','melanocyte','mFib','Endothelial cells','T cells','APC','CAFs'),CNV=c(sum(celltype_CNV_distribution[celltype_CNV_distribution$cell_type%in%"basal",]$CNV),
sum(celltype_CNV_distribution[celltype_CNV_distribution$cell_type%in%"granulosum",]$CNV),sum(celltype_CNV_distribution[celltype_CNV_distribution$cell_type%in%"spinosum",]$CNV),
sum(celltype_CNV_distribution[celltype_CNV_distribution$cell_type%in%"follicular",]$CNV),sum(celltype_CNV_distribution[celltype_CNV_distribution$cell_type%in%"melanocyte",]$CNV),
sum(celltype_CNV_distribution[celltype_CNV_distribution$cell_type%in%"mFib",]$CNV),sum(celltype_CNV_distribution[celltype_CNV_distribution$cell_type%in%"Endothelial cells",]$CNV),
sum(celltype_CNV_distribution[celltype_CNV_distribution$cell_type%in%"T cells",]$CNV),sum(celltype_CNV_distribution[celltype_CNV_distribution$cell_type%in%"APC",]$CNV),
sum(celltype_CNV_distribution[celltype_CNV_distribution$cell_type%in%"CAFs",]$CNV)))
######adpter总量
celltype_adapter_distribution <- data.frame(cell_type=new_combined_atac@meta.data$celltype_ATAC,adapter_num=new_combined_atac@meta.data$adapter_num)
celltype_adapter<-data.frame(cell_type=c('basal','granulosum','spinosum','follicular','melanocyte','mFib','Endothelial cells','T cells','APC','CAFs'),adapter_num=c(sum(celltype_adapter_distribution[celltype_adapter_distribution$cell_type%in%"basal",]$adapter_num),
sum(celltype_adapter_distribution[celltype_adapter_distribution$cell_type%in%"granulosum",]$adapter_num),sum(celltype_adapter_distribution[celltype_adapter_distribution$cell_type%in%"spinosum",]$adapter_num),
sum(celltype_adapter_distribution[celltype_adapter_distribution$cell_type%in%"follicular",]$adapter_num),sum(celltype_adapter_distribution[celltype_adapter_distribution$cell_type%in%"melanocyte",]$adapter_num),
sum(celltype_adapter_distribution[celltype_adapter_distribution$cell_type%in%"mFib",]$adapter_num),sum(celltype_adapter_distribution[celltype_adapter_distribution$cell_type%in%"Endothelial cells",]$adapter_num),
sum(celltype_adapter_distribution[celltype_adapter_distribution$cell_type%in%"T cells",]$adapter_num),sum(celltype_adapter_distribution[celltype_adapter_distribution$cell_type%in%"APC",]$adapter_num),
sum(celltype_adapter_distribution[celltype_adapter_distribution$cell_type%in%"CAFs",]$adapter_num)))
######eccDNA总量
celltype_eccDNA_gene_num_distribution <- data.frame(cell_type=new_combined_atac@meta.data$celltype_ATAC,eccDNA_gene_num=new_combined_atac@meta.data$eccDNA_gene_num)
celltype_eccDNA_gene_num<-data.frame(cell_type=c('basal','granulosum','spinosum','follicular','melanocyte','mFib','Endothelial cells','T cells','APC','CAFs'),eccDNA_gene_num=c(sum(celltype_eccDNA_gene_num_distribution[celltype_eccDNA_gene_num_distribution$cell_type%in%"basal",]$eccDNA_gene_num),
sum(celltype_eccDNA_gene_num_distribution[celltype_eccDNA_gene_num_distribution$cell_type%in%"granulosum",]$eccDNA_gene_num),sum(celltype_eccDNA_gene_num_distribution[celltype_eccDNA_gene_num_distribution$cell_type%in%"spinosum",]$eccDNA_gene_num),
sum(celltype_eccDNA_gene_num_distribution[celltype_eccDNA_gene_num_distribution$cell_type%in%"follicular",]$eccDNA_gene_num),sum(celltype_eccDNA_gene_num_distribution[celltype_eccDNA_gene_num_distribution$cell_type%in%"melanocyte",]$eccDNA_gene_num),
sum(celltype_eccDNA_gene_num_distribution[celltype_eccDNA_gene_num_distribution$cell_type%in%"mFib",]$eccDNA_gene_num),sum(celltype_eccDNA_gene_num_distribution[celltype_eccDNA_gene_num_distribution$cell_type%in%"Endothelial cells",]$eccDNA_gene_num),
sum(celltype_eccDNA_gene_num_distribution[celltype_eccDNA_gene_num_distribution$cell_type%in%"T cells",]$eccDNA_gene_num),sum(celltype_eccDNA_gene_num_distribution[celltype_eccDNA_gene_num_distribution$cell_type%in%"APC",]$eccDNA_gene_num),
sum(celltype_eccDNA_gene_num_distribution[celltype_eccDNA_gene_num_distribution$cell_type%in%"CAFs",]$eccDNA_gene_num)))



library(ggpubr)
library(digest)
p<-ggbarplot(celltype_CNV, "cell_type", "CNV",color = "cell_type", size = 1,width=0.7, palette =c("#56B4E9", 'gray', '#CCEBC5', '#BC80BD', '#FCCDE5','#BEBADA', '#FDB462',  '#FB8072','#8DD3C7','#FFFFB3','#B3DE69'),add = "jitter")
p
p<-ggbarplot(celltype_adapter, "cell_type", "adapter_num",color = "cell_type", size = 1,width=0.7, palette =c("#56B4E9", 'gray', '#CCEBC5', '#BC80BD', '#FCCDE5','#BEBADA', '#FDB462',  '#FB8072','#8DD3C7','#FFFFB3','#B3DE69'),add = "jitter")
p
p<-ggbarplot(celltype_eccDNA_gene_num, "cell_type", "eccDNA_gene_num",color = "cell_type", size = 1,width=0.7, palette =c("#56B4E9", 'gray', '#CCEBC5', '#BC80BD', '#FCCDE5','#BEBADA', '#FDB462',  '#FB8072','#8DD3C7','#FFFFB3','#B3DE69'),add = "jitter")
p

