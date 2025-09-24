######For Figure 5A
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
#######marker gene 细胞分布与表达情况
candidate_gene<-c('PTMA','NRP2','KLF6','MYC','PCAT1','PVT1')
candidate_frame<-eccDNA_frame[,]
index<-apply(candidate_frame, 1,function(x){return(x!=0)})
table(cell_type[rownames(index[index[,1]==TRUE,])])
feature_frame<-data.frame(basal=c(),)
library(pheatmap)
# clone_matrix[clone_matrix!=0]<-1
pheatmap(clone_matrix,cluster_cols = T,cluster_rows = F,
color = colorRampPalette(c("#D3D3D3","white","#C8504F"))(50),main=candidate_gene,annotation_col = annotation_col,show_rownames = F)

######For Figure 5B
######scRNA
library(Seurat)
library(ggplot2)
library(dplyr)
setwd('/data2/home/lijie/CSCC/cell_ranger/RNA-data')
load('/data2/home/lijie/scATAC_eccDNA/results/scRNA_eccDNA/eccDNA_gene.Rdata')
load('/data2/home/lijie/scATAC_eccDNA/results/eccDNA_clone/tissue_DEG.Rdata')
load('/data2/home/lijie/CSCC/cell_ranger/RNA-data/scRNA_harmony.modified.Rdata')
candidate_gene<-c('PTMA','NRP2','KLF6','MYC','PCAT1','PVT1')
celltype_marker<-c(                                                            ### cluster11、13、15、16没有明显分类                                                 
  "COL17A1", "KRT14","KRT5" ,                                # basal                                        ### cluster1、2、4
  "CD3D" ,                                   # T cells                                      ### cluster7、
  "COL3A1", "COL1A1",   # mFib                                         ### cluster6、
  "KRT1", "KRT10" ,"KRTDAP",                                 # spinosum                                     ### cluster0、
  "S100A7", "S100A8", "S100A9",                              # granulosum                                   ### cluster3、
    "KRT17", "ADGRL3", "TENM2",    # follicular没有23两个marker                   ### cluster9、12
  "CD68", "AIF1" ,"C1QB" ,"C1QA" ,                           # APC新加 "HLA−DRA", "HLA−DPB1", "HLA−DPA1"    ### cluster8
  "VWF" ,                                  # Endothelial cells                            ### cluster5、
  "ACTA2", "RGS5",                                  # CAFs                                         ### cluster10、17
  "DCN" ,  # mesFib                                       ### marker跟mFib一样，区分不开    
  "CD79A", "BANK1", "MS4A1", "CD19","CD79B",                 # B cells没有"BANK1"marker                     ### 没有
  "IGKC", "IGHG1", "IGHA1",                         # Plasma cells1 没有后3个marker                ### 没有
  "TPSB2", "TPSAB1", "CPA3" ,                                # Mast cells                                   ### 没有
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

#####For Figure 5C

########candidate_gene有扩增的细胞类型和没扩增细胞类型表达对比
basal_existence<-match(candidate_gene,basal_gene)
basal_existence[!is.na(basal_existence)]<-1
basal_existence[is.na(basal_existence)]<-0
spinosum_existence<-match(candidate_gene,spinosum_gene)
spinosum_existence[!is.na(spinosum_existence)]<-1
spinosum_existence[is.na(spinosum_existence)]<-0
granulosum_existence<-match(candidate_gene,granulosum_gene)
granulosum_existence[!is.na(granulosum_existence)]<-1
granulosum_existence[is.na(granulosum_existence)]<-0
follicular_existence<-match(candidate_gene,follicular_gene)
follicular_existence[!is.na(follicular_existence)]<-1
follicular_existence[is.na(follicular_existence)]<-0
melanocyte_existence<-match(candidate_gene,melanocyte_gene)
melanocyte_existence[!is.na(melanocyte_existence)]<-1
melanocyte_existence[is.na(melanocyte_existence)]<-0
candidate_gene_frame<-data.frame(basal_existence,spinosum_existence,granulosum_existence,follicular_existence,melanocyte_existence)
rownames(candidate_gene_frame)<-candidate_gene
library(pheatmap)
pheatmap(candidate_gene_frame,cluster_cols = F,cluster_rows = F,show_rownames = T)

######expression
setwd('/data2/home/lijie/scATAC_eccDNA/results/scRNA_eccDNA')
epi_index<-eccDNA_sample_umap@meta.data[["cell_id"]][c(grep('basal',eccDNA_sample_umap@meta.data[["seurat_clusters"]]),grep('spinosum',eccDNA_sample_umap@meta.data[["seurat_clusters"]]),grep('granulosum',eccDNA_sample_umap@meta.data[["seurat_clusters"]]),grep('follicular',eccDNA_sample_umap@meta.data[["seurat_clusters"]]),grep('melanocyte',eccDNA_sample_umap@meta.data[["seurat_clusters"]]))]
epi_cell<-subset(eccDNA_sample_umap,cells=epi_index)

DimPlot(object = epi_cell, reduction = "umap",raster=FALSE)+NoLegend()
DimPlot(object = epi_cell, reduction = "umap",raster=FALSE,group.by = 'seurat_clusters')+NoLegend()
DimPlot(object = epi_cell, reduction = "umap",raster=FALSE,group.by = 'patient_id')+NoLegend()
DotPlot(epi_cell,features=celltype_marker,cols=c('#BEBEBE','#FF0000'),scale=F)+coord_flip()
DotPlot(epi_cell,features=CSCC_eccDNA[1:10],cols=c('#BEBEBE','#FF0000'),scale=F,group.by = 'tissue_type')
FeaturePlot(epi_cell,features='PTMA',cols=c('#BEBEBE','#FF0000'),order=T,raster = F)

other_index<-eccDNA_sample_umap@meta.data[["cell_id"]][c(grep('granulosum',eccDNA_sample_umap@meta.data[["seurat_clusters"]]),grep('follicular',eccDNA_sample_umap@meta.data[["seurat_clusters"]]),grep('melanocyte',eccDNA_sample_umap@meta.data[["seurat_clusters"]]))]
other_cell<-subset(eccDNA_sample_umap,cells=other_index)
subset_index<-eccDNA_sample_umap@meta.data[["cell_id"]][c(grep('basal',eccDNA_sample_umap@meta.data[["seurat_clusters"]]),grep('spinosum',eccDNA_sample_umap@meta.data[["seurat_clusters"]]))]
subset_cell<-subset(eccDNA_sample_umap,cells=subset_index)


other_id<-other_cell@meta.data[["cell_id"]]
cell_id<-subset_cell@meta.data[["cell_id"]]
other_SCT<-as.matrix(other_cell@assays[["SCT"]]@data)
scRNA_SCT<-as.matrix(subset_cell@assays[["SCT"]]@data)

uni_eccDNA<-'NRP2'
frame_eccDNA<-other_cell@assays[["SCT"]]@data[uni_eccDNA,]
frame_eccDNA<-frame_eccDNA[frame_eccDNA>0]
diff_frame_eccDNA<-subset_cell@assays[["SCT"]]@data[uni_eccDNA,]
diff_frame_eccDNA<-diff_frame_eccDNA[diff_frame_eccDNA>0]
wilcox.test(frame_eccDNA,diff_frame_eccDNA)
wilcox.test(frame_eccDNA,diff_frame_eccDNA,alternative = 'great')
mean(frame_eccDNA)/mean(diff_frame_eccDNA)

df<-data.frame(c(frame_eccDNA,diff_frame_eccDNA),c(rep("frame_eccDNA",length(frame_eccDNA)),rep("diff_frame_eccDNA",length(diff_frame_eccDNA))))
colnames(df)<-c("num","group")

library(ggpubr)

library(digest)
p<-ggboxplot(df, "group", "num",color = "group", size = 1,width=0.5, palette =c("#D45252","#6EA6D9"),add = "jitter", shape = "group")
p

########Figure 5F
###采用seurat对象来研究细胞间克隆的关系
options(stringsAsFactors = FALSE)
options(scipen = 999)
setwd('/data2/home/lijie/scATAC_eccDNA/results/eccDNA_clone')
load('/data2/home/lijie/scATAC_eccDNA/results/CSCC3/copy/cell_type.Rdata')
load('eccDNA_frame.Rdata')
cell_eccDNA_name<-colnames(eccDNA_frame)
cell_eccDNA_name<-gsub('N3','Nomal3',cell_eccDNA_name)
cell_eccDNA_name<-gsub('N2','Nomal2',cell_eccDNA_name)
cell_eccDNA_name<-gsub('N1','Nomal1',cell_eccDNA_name)
cell_eccDNA_name<-gsub('A3','AK3',cell_eccDNA_name)
cell_eccDNA_name<-gsub('A2','AK2',cell_eccDNA_name)
cell_eccDNA_name<-gsub('A1','AK1',cell_eccDNA_name)
cell_eccDNA_name<-gsub('C3','CSCC3',cell_eccDNA_name)
cell_eccDNA_name<-gsub('C2','CSCC2',cell_eccDNA_name)
cell_eccDNA_name<-gsub('C1','CSCC1',cell_eccDNA_name)
colnames(eccDNA_frame)<-cell_eccDNA_name
post_eccDNA_frame<-eccDNA_frame[,intersect(colnames(eccDNA_frame),names(cell_type))]
post_eccDNA_frame<-post_eccDNA_frame[rowMeans(post_eccDNA_frame)!=0,]
save(post_eccDNA_frame,file='post_eccDNA_frame.Rdata')

library(dplyr)
library(Seurat)
library('harmony')
eccDNA_sample <- CreateSeuratObject(counts = post_eccDNA_frame, project = "eccDNA_sample",min.cells = 2,min.features = 2)
eccDNA_sample_normalized <- NormalizeData(eccDNA_sample,verbose = F)
eccDNA_sample_feature <- FindVariableFeatures(eccDNA_sample_normalized,selection.method = 'vst',nfeatures = length(eccDNA_sample@assays[["RNA"]]@features),verbose = F)

#######################################################
eccDNA_sample_scale <- ScaleData(eccDNA_sample_feature,verbose = F)
eccDNA_sample_pca <- RunPCA(eccDNA_sample_scale,verbose = F)

eccDNA_sample_pca <- RunHarmony(object = eccDNA_sample_pca, group.by.vars = "orig.ident")#theta：整合强度参数，数值越大表示批次校正越强（默认 2，常调 1–3）

set.seed(12345678)
ElbowPlot(eccDNA_sample_pca)
eccDNA_sample_pca <- FindNeighbors(eccDNA_sample_pca,dims = 1:20, reduction = "harmony")
eccDNA_sample_pca <- FindClusters(eccDNA_sample_pca,resolution = 3)
eccDNA_sample_umap <- RunUMAP(eccDNA_sample_pca,dims = 2:20,verbose = F,reduction = "harmony")
eccDNA_sample_umap@meta.data[["cell_type"]]<-cell_type[names(eccDNA_sample@active.ident)]
sample_type<-c()
sample_type[grep('Nomal',eccDNA_sample_umap@meta.data[["orig.ident"]])]<-'Nomal'
sample_type[grep('AK',eccDNA_sample_umap@meta.data[["orig.ident"]])]<-'AK'
sample_type[grep('CSCC',eccDNA_sample_umap@meta.data[["orig.ident"]])]<-'CSCC'
eccDNA_sample_umap@meta.data[["sample_type"]]<-factor(sample_type,levels=c('CSCC','AK','Nomal'))

DimPlot(object = eccDNA_sample_umap, reduction = "umap")
DimPlot(object = eccDNA_sample_umap, reduction = "umap",group.by = 'cell_type')
DimPlot(object = eccDNA_sample_umap, reduction = "umap",group.by = 'orig.ident')
DimPlot(object = eccDNA_sample_umap, reduction = "umap",group.by = 'sample_type')

#####Figure 5G
######各个sample type的marker gene
eccDNA_sample_umap <- SetIdent(eccDNA_sample_umap, value = "sample_type")
sample_type_eccDNA_markers <- FindAllMarkers(eccDNA_sample_umap, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1)
save(sample_type_eccDNA_markers,file='sample_type_eccDNA_markers.Rdata')
write.table(sample_type_eccDNA_markers,file='sample_type_eccDNA_markers.txt',sep='\t',row.names = F,col.names = T,quote=F)
load('sample_type_eccDNA_markers.Rdata')
top10 <- sample_type_eccDNA_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
library(ggplot2)
DoHeatmap(eccDNA_sample_umap, features = top10$gene,raster =F)+scale_fill_gradient2(low = "#BEBEBE",high = "#FF0000", midpoint = 0)
DotPlot(eccDNA_sample_umap,features=top10$gene,cols=c('#BEBEBE','#FF0000'),scale=F)
######Figure S6A,B,C
###功能富集
CSCC_eccDNA<-sample_type_eccDNA_markers[sample_type_eccDNA_markers$cluster%in%'CSCC',]$gene
AK_eccDNA<-sample_type_eccDNA_markers[sample_type_eccDNA_markers$cluster%in%'AK',]$gene
Nomal_eccDNA<-sample_type_eccDNA_markers[sample_type_eccDNA_markers$cluster%in%'Nomal',]$gene
save(CSCC_eccDNA,AK_eccDNA,Nomal_eccDNA,file='tissue_DEG.Rdata')
library(org.Hs.eg.db)
k=keys(org.Hs.eg.db,keytype = "ENSEMBL")
list=select(org.Hs.eg.db,keys=k,columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
pro_entrezid=list[match(unique(Nomal_eccDNA),list[,"SYMBOL"]),][,2]

library(clusterProfiler)
KEGG_gene<-enrichKEGG(pro_entrezid,organism = "hsa", pvalueCutoff = 0.05)
KEGG_gene<-setReadable(KEGG_gene,OrgDb='org.Hs.eg.db',keyType='ENTREZID')
GO_gene<-enrichGO(pro_entrezid,'org.Hs.eg.db',ont = "BP", pvalueCutoff = 0.05)
GO_gene<-setReadable(GO_gene,OrgDb='org.Hs.eg.db',keyType='ENTREZID')
View(GO_gene@result)
View(KEGG_gene@result)
write.table(KEGG_gene@result,file='Nomal_KEGG_tss.txt',sep='\t',row.names = F,col.names = T,quote = F)
write.table(GO_gene@result,file='Nomal_GO_tss.txt',sep='\t',row.names = F,col.names = T,quote = F)


KEGG_tss<-read.table(file='Nomal_KEGG_tss.txt',sep='\t',header=T,quote='')
GO_tss<-read.table(file='Nomal_GO_tss.txt',sep='\t',header=T,quote='')
observed_KEGG<-c('Cell cycle','DNA replication')
observed_GO<-c('G1/S transition of mitotic cell cycle','cell cycle G1/S phase transition','translational elongation','double-strand break repair via break-induced replication','transforming growth factor beta activation','mesenchymal cell differentiation','regulation of epithelial to mesenchymal transition',
'mesenchyme development','pigmentation','axon extension','neuron projection extension','axonogenesis')
all_path<-GO_tss
path_index<-match(observed_GO,all_path[,2])
new_kegg<-all_path[path_index,c(2,5,9)]
all_kegg<-new_kegg

final_frame<-all_kegg
final_frame<-na.omit(final_frame)
#final_frame$Description<-factor(final_frame$Description,levels=rev(observed_GO))
final_frame$Description<-factor(final_frame$Description,levels=final_frame[order(final_frame$Count),]$Description)
library(ggplot2)
p = ggplot(final_frame,aes(Count,Description))
p = ggplot(final_frame,aes(Count,Description))
p=p + geom_point()  
# 修稿点的大小
p=p + geom_point(aes(size=Count))
# 展示三维数据
pbubble = p+ geom_point(aes(size=Count,color=pvalue))
# 设置渐变色
pr = pbubble+scale_color_gradient(low="red",high = "yellow")
# 绘制p气泡图
pr = pr+labs(color=expression(pvalue),size="Count",  
                           x="GeneRatio",y="Pathway name",title="Pathway enrichment")
pr + theme_bw()


####Figure 5H,I
############组织特异的eccDNA表达情况
CSCC_index<-scRNA_harmony.modified@meta.data[["cell_id"]][grep('CSCC',scRNA_harmony.modified@meta.data[["tissue_type"]])]
CSCC_cell<-subset(scRNA_harmony.modified,cells=CSCC_index)
AK_index<-scRNA_harmony.modified@meta.data[["cell_id"]][grep('AK',scRNA_harmony.modified@meta.data[["tissue_type"]])]
AK_cell<-subset(scRNA_harmony.modified,cells=AK_index)
Nomal_index<-scRNA_harmony.modified@meta.data[["cell_id"]][grep('Normal',scRNA_harmony.modified@meta.data[["tissue_type"]])]
Nomal_cell<-subset(scRNA_harmony.modified,cells=Nomal_index)


patient_id<-granulosum_cell@meta.data[["patient_id"]]


#######eccDNA与非eccDNA在各种细胞中表达情况,boxplot
setwd('/data2/home/lijie/scATAC_eccDNA/results/scRNA_eccDNA')
subset_cell<-Nomal_cell
uni_eccDNA<-Nomal_eccDNA

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


