######For Figure 4A, 4B, S5A,S5B,S5C
library(Seurat)
library(ggplot2)
library(dplyr)
setwd('/data2/home/lijie/CSCC/cell_ranger/RNA-data')
load('/data2/home/lijie/scATAC_eccDNA/results/scRNA_eccDNA/eccDNA_gene.Rdata')
load('/data2/home/lijie/scATAC_eccDNA/results/eccDNA_clone/tissue_DEG.Rdata')
load('/data2/home/lijie/CSCC/cell_ranger/RNA-data/scRNA_harmony.modified.Rdata')
candidate_gene<-c('PTMA','NRP2','KLF6','MYC','PCAT1','PVT1')
celltype_marker<-c(                                                                                                          
  "COL17A1", "KRT14","KRT5" ,                                # basal                                        ### cluster1、2、4
  "CD3D" ,                                   # T cells                                      ### cluster7、
  "COL3A1", "COL1A1",   # mFib                                         ### cluster6、
  "KRT1", "KRT10" ,"KRTDAP",                                 # spinosum                                     ### cluster0、
  "S100A7", "S100A8", "S100A9",                              # granulosum                                   ### cluster3、
    "KRT17", "ADGRL3", "TENM2",    # follicular                   ### cluster9、12
  "CD68", "AIF1" ,"C1QB" ,"C1QA" ,                           #"HLA−DRA", "HLA−DPB1", "HLA−DPA1"    ### cluster8
  "VWF" ,                                  # Endothelial cells                            ### cluster5、
  "ACTA2", "RGS5",                                  # CAFs                                         ### cluster10、17
  "DCN" ,  # mesFib                                       ### marker    
  "CD79A", "BANK1", "MS4A1", "CD19","CD79B",                 # B cells                    
  "IGKC", "IGHG1", "IGHA1",                         # Plasma cells1             
  "TPSB2", "TPSAB1", "CPA3" ,                                # Mast cells                                   ### 
  "DCT", "TYRP1", "PMEL"   )                                 # melanocyte                                   ### cluster14、
DimPlot(object = scRNA_harmony.modified, reduction = "umap",raster=FALSE)
eccDNA_sample_pca<-scRNA_harmony.modified
set.seed(12345678)
eccDNA_sample_umap <- RunUMAP(eccDNA_sample_pca,dims = 1:20)
eccDNA_sample_umap@meta.data[["seurat_clusters"]][eccDNA_sample_umap@meta.data[["seurat_clusters"]]=='Plasma cells2']<-'Plasma cells1'
DimPlot(object = eccDNA_sample_umap, reduction = "umap",raster=FALSE)+NoLegend()
DimPlot(object = eccDNA_sample_umap, reduction = "umap",raster=FALSE,group.by = 'seurat_clusters')+NoLegend()
DimPlot(object = eccDNA_sample_umap, reduction = "umap",raster=FALSE,group.by = 'patient_id')+NoLegend()
DimPlot(object = eccDNA_sample_umap, reduction = "umap",raster=FALSE,group.by = 'tissue_type',cols = c("#2CAC3E", "#EB736A","#6692CC"))+NoLegend()
DotPlot(eccDNA_sample_umap,features=celltype_marker,cols=c('#BEBEBE','#FF0000'),scale=F)+coord_flip()
DotPlot(eccDNA_sample_umap,features=candidate_gene,cols=c('#BEBEBE','#FF0000'),scale=F)
FeaturePlot(eccDNA_sample_umap,features='PTMA',cols=c('#BEBEBE','#FF0000'),order=T,raster = F)

VlnPlot(scRNA_harmony.modified,features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = 'tissue_type',raster = F,pt.size = -1)
DotPlot(scRNA_harmony.modified,features=celltype_marker,cols=c('#BEBEBE','#FF0000'),scale=F)
DotPlot(scRNA_harmony.modified,features=candidate_gene,cols=c('#BEBEBE','#FF0000'),scale=F)

########
options(stringsAsFactors = FALSE)
setwd('/data2/home/lijie/scATAC_eccDNA/results/scRNA_eccDNA')
load('/data2/home/lijie/scATAC_eccDNA/results/eccDNA_clone/eccDNA_frame.Rdata')
load('/data2/home/lijie/scATAC_eccDNA/results/CSCC3/copy/cell_type.Rdata')
colnames(eccDNA_frame)<-gsub('N3','Nomal3',colnames(eccDNA_frame))
colnames(eccDNA_frame)<-gsub('N2','Nomal2',colnames(eccDNA_frame))
colnames(eccDNA_frame)<-gsub('N1','Nomal1',colnames(eccDNA_frame))
colnames(eccDNA_frame)<-gsub('A3','AK3',colnames(eccDNA_frame))
colnames(eccDNA_frame)<-gsub('A2','AK2',colnames(eccDNA_frame))
colnames(eccDNA_frame)<-gsub('A1','AK1',colnames(eccDNA_frame))
colnames(eccDNA_frame)<-gsub('C3','CSCC3',colnames(eccDNA_frame))
colnames(eccDNA_frame)<-gsub('C2','CSCC2',colnames(eccDNA_frame))
colnames(eccDNA_frame)<-gsub('C1','CSCC1',colnames(eccDNA_frame))


####Normal->AK->CSCC
Nomal3_gene_annotation<-eccDNA_frame[,colnames(eccDNA_frame)[grep('Nomal3',colnames(eccDNA_frame))]]
Nomal3_gene<-apply(Nomal3_gene_annotation,1,sum)
Nomal3_duplicated_genes<-names(sort(Nomal3_gene[Nomal3_gene!=0],decreasing = T))
Nomal2_gene_annotation<-eccDNA_frame[,colnames(eccDNA_frame)[grep('Nomal2',colnames(eccDNA_frame))]]
Nomal2_gene<-apply(Nomal2_gene_annotation,1,sum)
Nomal2_duplicated_genes<-names(sort(Nomal2_gene[Nomal2_gene!=0],decreasing = T))
Nomal1_gene_annotation<-eccDNA_frame[,colnames(eccDNA_frame)[grep('Nomal1',colnames(eccDNA_frame))]]
Nomal1_gene<-apply(Nomal1_gene_annotation,1,sum)
Nomal1_duplicated_genes<-names(sort(Nomal1_gene[Nomal1_gene!=0],decreasing = T))

AK3_gene_annotation<-eccDNA_frame[,colnames(eccDNA_frame)[grep('AK3',colnames(eccDNA_frame))]]
AK3_gene<-apply(AK3_gene_annotation,1,sum)
AK3_duplicated_genes<-names(sort(AK3_gene[AK3_gene!=0],decreasing = T))
AK2_gene_annotation<-eccDNA_frame[,colnames(eccDNA_frame)[grep('AK2',colnames(eccDNA_frame))]]
AK2_gene<-apply(AK2_gene_annotation,1,sum)
AK2_duplicated_genes<-names(sort(AK2_gene[AK2_gene!=0],decreasing = T))
AK1_gene_annotation<-eccDNA_frame[,colnames(eccDNA_frame)[grep('AK1',colnames(eccDNA_frame))]]
AK1_gene<-apply(AK1_gene_annotation,1,sum)
AK1_duplicated_genes<-names(sort(AK1_gene[AK1_gene!=0],decreasing = T))

CSCC3_gene_annotation<-eccDNA_frame[,colnames(eccDNA_frame)[grep('CSCC3',colnames(eccDNA_frame))]]
CSCC3_gene<-apply(CSCC3_gene_annotation,1,sum)
CSCC3_duplicated_genes<-names(sort(CSCC3_gene[CSCC3_gene!=0],decreasing = T))
CSCC2_gene_annotation<-eccDNA_frame[,colnames(eccDNA_frame)[grep('CSCC2',colnames(eccDNA_frame))]]
CSCC2_gene<-apply(CSCC2_gene_annotation,1,sum)
CSCC2_duplicated_genes<-names(sort(CSCC2_gene[CSCC2_gene!=0],decreasing = T))
CSCC1_gene_annotation<-eccDNA_frame[,colnames(eccDNA_frame)[grep('CSCC1',colnames(eccDNA_frame))]]
CSCC1_gene<-apply(CSCC1_gene_annotation,1,sum)
CSCC1_duplicated_genes<-names(sort(CSCC1_gene[CSCC1_gene!=0],decreasing = T))

library(VennDiagram)
granulosum_gene_annotation<-eccDNA_frame[,intersect(colnames(eccDNA_frame),names(cell_type[cell_type=='granulosum']))]
granulosum_gene<-apply(granulosum_gene_annotation,1,sum)
granulosum_gene<-names(sort(granulosum_gene[granulosum_gene!=0],decreasing = T))
basal_gene_annotation<-eccDNA_frame[,intersect(colnames(eccDNA_frame),names(cell_type[cell_type=='basal']))]
basal_gene<-apply(basal_gene_annotation,1,sum)
basal_gene<-names(sort(basal_gene[basal_gene!=0],decreasing = T))
spinosum_gene_annotation<-eccDNA_frame[,intersect(colnames(eccDNA_frame),names(cell_type[cell_type=='spinosum']))]
spinosum_gene<-apply(spinosum_gene_annotation,1,sum)
spinosum_gene<-names(sort(spinosum_gene[spinosum_gene!=0],decreasing = T))
follicular_gene_annotation<-eccDNA_frame[,intersect(colnames(eccDNA_frame),names(cell_type[cell_type=='follicular']))]
follicular_gene<-apply(follicular_gene_annotation,1,sum)
follicular_gene<-names(sort(follicular_gene[follicular_gene!=0],decreasing = T))
melanocyte_gene_annotation<-eccDNA_frame[,intersect(colnames(eccDNA_frame),names(cell_type[cell_type=='melanocyte']))]
melanocyte_gene<-apply(melanocyte_gene_annotation,1,sum)
melanocyte_gene<-names(sort(melanocyte_gene[melanocyte_gene!=0],decreasing = T))

####
APC_gene_annotation<-eccDNA_frame[,intersect(colnames(eccDNA_frame),names(cell_type[cell_type=='APC']))]
APC_gene<-apply(APC_gene_annotation,1,sum)
APC_gene<-names(sort(APC_gene[APC_gene!=0],decreasing = T))
CAFs_gene_annotation<-eccDNA_frame[,intersect(colnames(eccDNA_frame),names(cell_type[cell_type=='CAFs']))]
CAFs_gene<-apply(CAFs_gene_annotation,1,sum)
CAFs_gene<-names(sort(CAFs_gene[CAFs_gene!=0],decreasing = T))
Endothelial_gene_annotation<-eccDNA_frame[,intersect(colnames(eccDNA_frame),names(cell_type[cell_type=='Endothelial cells']))]
Endothelial_gene<-apply(Endothelial_gene_annotation,1,sum)
Endothelial_gene<-names(sort(Endothelial_gene[Endothelial_gene!=0],decreasing = T))
mFib_gene_annotation<-eccDNA_frame[,intersect(colnames(eccDNA_frame),names(cell_type[cell_type=='mFib']))]
mFib_gene<-apply(mFib_gene_annotation,1,sum)
mFib_gene<-names(sort(mFib_gene[mFib_gene!=0],decreasing = T))
T_gene_annotation<-eccDNA_frame[,intersect(colnames(eccDNA_frame),names(cell_type[cell_type=='T cells']))]
T_gene<-apply(T_gene_annotation,1,sum)
T_gene<-names(sort(T_gene[T_gene!=0],decreasing = T))

eccDNA_gene<-unique(c(granulosum_gene,basal_gene,spinosum_gene,melanocyte_gene,follicular_gene))
save(granulosum_gene,basal_gene,spinosum_gene,melanocyte_gene,follicular_gene,file='eccDNA_gene.Rdata')


basal_index<-scRNA_harmony.modified@meta.data[["cell_id"]][grep('basal',scRNA_harmony.modified@meta.data[["seurat_clusters"]])]
basal_cell<-subset(scRNA_harmony.modified,cells=basal_index)
spinosum_index<-scRNA_harmony.modified@meta.data[["cell_id"]][grep('spinosum',scRNA_harmony.modified@meta.data[["seurat_clusters"]])]
spinosum_cell<-subset(scRNA_harmony.modified,cells=spinosum_index)
granulosum_index<-scRNA_harmony.modified@meta.data[["cell_id"]][grep('granulosum',scRNA_harmony.modified@meta.data[["seurat_clusters"]])]
granulosum_cell<-subset(scRNA_harmony.modified,cells=granulosum_index)
follicular_index<-scRNA_harmony.modified@meta.data[["cell_id"]][grep('follicular',scRNA_harmony.modified@meta.data[["seurat_clusters"]])]
follicular_cell<-subset(scRNA_harmony.modified,cells=follicular_index)
melanocyte_index<-scRNA_harmony.modified@meta.data[["cell_id"]][grep('melanocyte',scRNA_harmony.modified@meta.data[["seurat_clusters"]])]
melanocyte_cell<-subset(scRNA_harmony.modified,cells=melanocyte_index)


#######Figure 4C,boxplot
setwd('/data2/home/lijie/scATAC_eccDNA/results/scRNA_eccDNA')
subset_cell<-melanocyte_cell
uni_eccDNA<-melanocyte_gene

cell_id<-subset_cell@meta.data[["cell_id"]]
scRNA_SCT<-as.matrix(subset_cell@assays[["SCT"]]@data)
uni_eccDNA<-intersect(uni_eccDNA,rownames(scRNA_SCT))

frame_eccDNA<-apply(scRNA_SCT[uni_eccDNA,],1,mean)
diff_frame<-setdiff(rownames(scRNA_SCT),uni_eccDNA)
diff_frame_eccDNA<-apply(scRNA_SCT[diff_frame,],1,mean)
wilcox.test(frame_eccDNA,diff_frame_eccDNA)
wilcox.test(frame_eccDNA,diff_frame_eccDNA,alternative = 'great')
mean(frame_eccDNA)/mean(diff_frame_eccDNA)

df<-data.frame(c(frame_eccDNA,diff_frame_eccDNA),c(rep("frame_eccDNA",length(frame_eccDNA)),rep("diff_frame_eccDNA",length(diff_frame_eccDNA))))
colnames(df)<-c("num","group")

library(ggpubr)

library(digest)
p<-ggboxplot(df, "group", "num",color = "group", size = 1,width=0.5, palette =c("#D45252","#6EA6D9"),add = "jitter", shape = "group")
p+scale_y_log10()

selected_markers<-uni_eccDNA[1:5]
FeaturePlot(eccDNA_sample_umap,features=selected_markers,raster=FALSE)

