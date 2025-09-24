#For Figure 2A
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

  

######CNV
celltype_CNV_distribution <- data.frame(cell_type=new_combined_atac@meta.data$celltype_ATAC,CNV=new_combined_atac@meta.data$CNV_copy)
celltype_CNV<-data.frame(cell_type=c('basal','granulosum','spinosum','follicular','melanocyte','mFib','Endothelial cells','T cells','APC','CAFs'),CNV=c(sum(celltype_CNV_distribution[celltype_CNV_distribution$cell_type%in%"basal",]$CNV),
sum(celltype_CNV_distribution[celltype_CNV_distribution$cell_type%in%"granulosum",]$CNV),sum(celltype_CNV_distribution[celltype_CNV_distribution$cell_type%in%"spinosum",]$CNV),
sum(celltype_CNV_distribution[celltype_CNV_distribution$cell_type%in%"follicular",]$CNV),sum(celltype_CNV_distribution[celltype_CNV_distribution$cell_type%in%"melanocyte",]$CNV),
sum(celltype_CNV_distribution[celltype_CNV_distribution$cell_type%in%"mFib",]$CNV),sum(celltype_CNV_distribution[celltype_CNV_distribution$cell_type%in%"Endothelial cells",]$CNV),
sum(celltype_CNV_distribution[celltype_CNV_distribution$cell_type%in%"T cells",]$CNV),sum(celltype_CNV_distribution[celltype_CNV_distribution$cell_type%in%"APC",]$CNV),
sum(celltype_CNV_distribution[celltype_CNV_distribution$cell_type%in%"CAFs",]$CNV)))
######adpter
celltype_adapter_distribution <- data.frame(cell_type=new_combined_atac@meta.data$celltype_ATAC,adapter_num=new_combined_atac@meta.data$adapter_num)
celltype_adapter<-data.frame(cell_type=c('basal','granulosum','spinosum','follicular','melanocyte','mFib','Endothelial cells','T cells','APC','CAFs'),adapter_num=c(sum(celltype_adapter_distribution[celltype_adapter_distribution$cell_type%in%"basal",]$adapter_num),
sum(celltype_adapter_distribution[celltype_adapter_distribution$cell_type%in%"granulosum",]$adapter_num),sum(celltype_adapter_distribution[celltype_adapter_distribution$cell_type%in%"spinosum",]$adapter_num),
sum(celltype_adapter_distribution[celltype_adapter_distribution$cell_type%in%"follicular",]$adapter_num),sum(celltype_adapter_distribution[celltype_adapter_distribution$cell_type%in%"melanocyte",]$adapter_num),
sum(celltype_adapter_distribution[celltype_adapter_distribution$cell_type%in%"mFib",]$adapter_num),sum(celltype_adapter_distribution[celltype_adapter_distribution$cell_type%in%"Endothelial cells",]$adapter_num),
sum(celltype_adapter_distribution[celltype_adapter_distribution$cell_type%in%"T cells",]$adapter_num),sum(celltype_adapter_distribution[celltype_adapter_distribution$cell_type%in%"APC",]$adapter_num),
sum(celltype_adapter_distribution[celltype_adapter_distribution$cell_type%in%"CAFs",]$adapter_num)))
######eccDNA
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

#For Figure 2B & Figure S2B

#########eccDNA clone
options(stringsAsFactors = FALSE)
options(scipen = 999)
setwd('/data2/home/lijie/scATAC_eccDNA/results/eccDNA_clone')
load('/data2/home/lijie/scATAC_eccDNA/results/CSCC3/copy/cell_type.Rdata')
load('/data2/home/lijie/scATAC_eccDNA/results/scRNA_eccDNA/eccDNA_gene.Rdata')
load('/data2/home/lijie/scATAC_eccDNA/results/Nomal3/copy/annotation_combine.Rdata')
Nomal3_duplicated_genes<-duplicated_genes
Nomal3_new_annotation<-annotation_combine
load('/data2/home/lijie/scATAC_eccDNA/results/Nomal2/copy/annotation_combine.Rdata')
Nomal2_duplicated_genes<-duplicated_genes
Nomal2_new_annotation<-annotation_combine
load('/data2/home/lijie/scATAC_eccDNA/results/Nomal1/copy/annotation_combine.Rdata')
Nomal1_duplicated_genes<-duplicated_genes
Nomal1_new_annotation<-annotation_combine
load('/data2/home/lijie/scATAC_eccDNA/results/AK3/copy/annotation_combine.Rdata')
AK3_duplicated_genes<-duplicated_genes
AK3_new_annotation<-annotation_combine
load('/data2/home/lijie/scATAC_eccDNA/results/AK2/copy/annotation_combine.Rdata')
AK2_duplicated_genes<-duplicated_genes
AK2_new_annotation<-annotation_combine
load('/data2/home/lijie/scATAC_eccDNA/results/AK1/copy/annotation_combine.Rdata')
AK1_duplicated_genes<-duplicated_genes
AK1_new_annotation<-annotation_combine
load('/data2/home/lijie/scATAC_eccDNA/results/CSCC3/copy/annotation_combine.Rdata')
CSCC3_duplicated_genes<-duplicated_genes
CSCC3_new_annotation<-annotation_combine
load('/data2/home/lijie/scATAC_eccDNA/results/CSCC2/copy/annotation_combine.Rdata')
CSCC2_duplicated_genes<-duplicated_genes
CSCC2_new_annotation<-annotation_combine
load('/data2/home/lijie/scATAC_eccDNA/results/CSCC1/copy/annotation_combine.Rdata')
CSCC1_duplicated_genes<-duplicated_genes
CSCC1_new_annotation<-annotation_combine

eccDNA_name<-unique(c(Nomal3_duplicated_genes,Nomal2_duplicated_genes,Nomal1_duplicated_genes,AK3_duplicated_genes,AK2_duplicated_genes,AK1_duplicated_genes,CSCC3_duplicated_genes,CSCC2_duplicated_genes,CSCC1_duplicated_genes))



eccDNA_frame<-data.frame()
for(i in 1:length(eccDNA_name)){
temp_frame_N3<-Nomal3_new_annotation[Nomal3_new_annotation$SYMBOL%in%eccDNA_name[i],]
temp_frame_N3<-temp_frame_N3[!duplicated(temp_frame_N3[,1:3]),]
temp_frame_N2<-Nomal2_new_annotation[Nomal2_new_annotation$SYMBOL%in%eccDNA_name[i],]
temp_frame_N2<-temp_frame_N2[!duplicated(temp_frame_N2[,1:3]),]
temp_frame_N1<-Nomal1_new_annotation[Nomal1_new_annotation$SYMBOL%in%eccDNA_name[i],]
temp_frame_N1<-temp_frame_N1[!duplicated(temp_frame_N1[,1:3]),]
temp_frame_A3<-AK3_new_annotation[AK3_new_annotation$SYMBOL%in%eccDNA_name[i],]
temp_frame_A3<-temp_frame_A3[!duplicated(temp_frame_A3[,1:3]),]
temp_frame_A2<-AK2_new_annotation[AK2_new_annotation$SYMBOL%in%eccDNA_name[i],]
temp_frame_A2<-temp_frame_A2[!duplicated(temp_frame_A2[,1:3]),]
temp_frame_A1<-AK1_new_annotation[AK1_new_annotation$SYMBOL%in%eccDNA_name[i],]
temp_frame_A1<-temp_frame_A1[!duplicated(temp_frame_A1[,1:3]),]
temp_frame_C3<-CSCC3_new_annotation[CSCC3_new_annotation$SYMBOL%in%eccDNA_name[i],]
temp_frame_C3<-temp_frame_C3[!duplicated(temp_frame_C3[,1:3]),]
temp_frame_C2<-CSCC2_new_annotation[CSCC2_new_annotation$SYMBOL%in%eccDNA_name[i],]
temp_frame_C2<-temp_frame_C2[!duplicated(temp_frame_C2[,1:3]),]
temp_frame_C1<-CSCC1_new_annotation[CSCC1_new_annotation$SYMBOL%in%eccDNA_name[i],]
temp_frame_C1<-temp_frame_C1[!duplicated(temp_frame_C1[,1:3]),]

copy_number_N3<-apply(temp_frame_N3[,28:dim(temp_frame_N3)[2]],2,function(x){sum(as.numeric(x))})
copy_number_N2<-apply(temp_frame_N2[,28:dim(temp_frame_N2)[2]],2,function(x){sum(as.numeric(x))})
copy_number_N1<-apply(temp_frame_N1[,28:dim(temp_frame_N1)[2]],2,function(x){sum(as.numeric(x))})
copy_number_A3<-apply(temp_frame_A3[,28:dim(temp_frame_A3)[2]],2,function(x){sum(as.numeric(x))})
copy_number_A2<-apply(temp_frame_A2[,28:dim(temp_frame_A2)[2]],2,function(x){sum(as.numeric(x))})
copy_number_A1<-apply(temp_frame_A1[,28:dim(temp_frame_A1)[2]],2,function(x){sum(as.numeric(x))})
copy_number_C3<-apply(temp_frame_C3[,28:dim(temp_frame_C3)[2]],2,function(x){sum(as.numeric(x))})
copy_number_C2<-apply(temp_frame_C2[,28:dim(temp_frame_C2)[2]],2,function(x){sum(as.numeric(x))})
copy_number_C1<-apply(temp_frame_C1[,28:dim(temp_frame_C1)[2]],2,function(x){sum(as.numeric(x))})

eccDNA_frame<-rbind(eccDNA_frame,t(as.data.frame(c(copy_number_N3,copy_number_N2,copy_number_N1,copy_number_A3,copy_number_A2,copy_number_A1,copy_number_C3,copy_number_C2,copy_number_C1))))
}
rownames(eccDNA_frame)<-eccDNA_name
colnames(eccDNA_frame)<-c(paste0('N3_',names(copy_number_N3)),paste0('N2_',names(copy_number_N2)),paste0('N1_',names(copy_number_N1)),paste0('A3_',names(copy_number_A3)),paste0('A2_',names(copy_number_A2)),paste0('A1_',names(copy_number_A1)),
paste0('C3_',names(copy_number_C3)),paste0('C2_',names(copy_number_C2)),paste0('C1_',names(copy_number_C1)))

sum_num<-apply(eccDNA_frame,2,sum)
eccDNA_frame<-eccDNA_frame[,sum_num!=0]
save(eccDNA_frame,file='eccDNA_frame.Rdata')

######eccDNA amplification fraction
options(stringsAsFactors = FALSE)
options(scipen = 999)
load('/data2/home/lijie/scATAC_eccDNA/results/eccDNA_clone/eccDNA_frame.Rdata')
setwd('/data2/home/lijie/scATAC_eccDNA/results/CSCC3/copy')
load('annotation_results.Rdata')
new_annotation$sample<-rep('CSCC3',dim(new_annotation)[1])
all_annotation<-new_annotation
setwd('/data2/home/lijie/scATAC_eccDNA/results/CSCC2/copy')
load('annotation_results.Rdata')
new_annotation$sample<-rep('CSCC2',dim(new_annotation)[1])
all_annotation<-rbind(all_annotation,new_annotation)
setwd('/data2/home/lijie/scATAC_eccDNA/results/CSCC1/copy')
load('annotation_results.Rdata')
new_annotation$sample<-rep('CSCC1',dim(new_annotation)[1])
all_annotation<-rbind(all_annotation,new_annotation)
setwd('/data2/home/lijie/scATAC_eccDNA/results/AK3/copy')
load('annotation_results.Rdata')
new_annotation$sample<-rep('AK3',dim(new_annotation)[1])
all_annotation<-rbind(all_annotation,new_annotation)
setwd('/data2/home/lijie/scATAC_eccDNA/results/AK2/copy')
load('annotation_results.Rdata')
new_annotation$sample<-rep('AK2',dim(new_annotation)[1])
all_annotation<-rbind(all_annotation,new_annotation)
setwd('/data2/home/lijie/scATAC_eccDNA/results/AK1/copy')
load('annotation_results.Rdata')
new_annotation$sample<-rep('AK1',dim(new_annotation)[1])
all_annotation<-rbind(all_annotation,new_annotation)
setwd('/data2/home/lijie/scATAC_eccDNA/results/Nomal3/copy')
load('annotation_results.Rdata')
new_annotation$sample<-rep('Nomal3',dim(new_annotation)[1])
all_annotation<-rbind(all_annotation,new_annotation)
setwd('/data2/home/lijie/scATAC_eccDNA/results/Nomal2/copy')
load('annotation_results.Rdata')
new_annotation$sample<-rep('Nomal2',dim(new_annotation)[1])
all_annotation<-rbind(all_annotation,new_annotation)
setwd('/data2/home/lijie/scATAC_eccDNA/results/Nomal1/copy')
load('annotation_results.Rdata')
new_annotation$sample<-rep('Nomal1',dim(new_annotation)[1])
all_annotation<-rbind(all_annotation,new_annotation)
setwd('/data2/home/lijie/scATAC_eccDNA/results/eccDNA_features')
new_exist_cell<-c()
for(i in 1:dim(all_annotation)[1]){
temp_data<-all_annotation[i,12:15]
if(length(grep('-1',temp_data))!=0){
temp_index<-strsplit(as.character(temp_data[grep('-1',temp_data)]),':')
new_exist_cell[i]<-temp_index[[1]][[3]]
}else{
new_exist_cell[i]<-0
}
}
all_annotation$exist_cell<-new_exist_cell
load('/data2/home/lijie/scATAC_eccDNA/results/CSCC3/copy/cell_type.Rdata')
cell_name<-unlist(lapply(strsplit(names(cell_type),'_'), function(x){x[[2]]}))
names(cell_type)<-cell_name
all_annotation<-all_annotation[all_annotation$exist_cell%in%cell_name,]

new_annotation<-all_annotation
write.table(all_annotation,file='all_annotation.txt',sep='\t',row.names = F,col.names = T,quote = F)
eccDNA_copy<-apply(eccDNA_frame,1,sum)
duplicated_genes<-names(sort(eccDNA_copy,decreasing = T))
save(duplicated_genes,file='duplicated_genes.Rdata')
candidate_gene<-duplicated_genes[1:20]
# candidate_gene<-c('PTMA','NRP2','KLF6','MYC','PCAT1','PVT1')
break_point<-c()
CNV_num<-c()
cell_num<-c()
for(i in 1:length(candidate_gene)){
candidate_frame<-new_annotation[new_annotation$SYMBOL%in%candidate_gene[i],]
break_point[i]<-dim(candidate_frame)[1]
CNV_num[i]<-eccDNA_copy[candidate_gene[i]]
cell_name_id<-c()
for(j in 1:dim(candidate_frame)[1]){
temp_data<-as.character(candidate_frame[j,12:14])
if(length(grep('-1',temp_data))!=0){
temp_index<-strsplit(temp_data[grep('-1',temp_data)],':')
cell_name_id[j]<-temp_index[[1]][[3]]
}else{
cell_name_id[j]<-0
}
}
cell_num[i]<-length(unique(cell_name_id))
}
names(break_point)<-candidate_gene
names(CNV_num)<-candidate_gene
names(cell_num)<-candidate_gene

final_frame<-data.frame(candidate_gene=candidate_gene,break_point=break_point,CNV_num=CNV_num,cell_num=cell_num)
final_frame$candidate_gene<-factor(final_frame$candidate_gene,levels=rev(final_frame$candidate_gene))

library(ggplot2)
p = ggplot(final_frame,aes(CNV_num,candidate_gene))
pbubble = p+ geom_point(aes(size=break_point,color=cell_num))
#
pr = pbubble+scale_color_gradient(low="grey",high = "red")
#
pr = pr+labs(color=expression(cell_num),size="break_point",  
                           x="CNV_num",y="candidate_gene",title="selected eccDNA")
pr + theme_bw()
######barplot
ggplot(final_frame, aes(x = candidate_gene, y = break_point, fill = CNV_num)) +
  geom_col(width = 0.8) +
  scale_fill_gradient(low = "grey", high = "red") +  
  labs(x = NULL, y = "eccDNA_num", fill = "CNV_num") +
  theme_classic()

###Figure 2D&2E
#eccDNA genome distribution
library(ChIPseeker)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GeneSummary)
library(dplyr)
library(ggplot2)

#read data
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(GenomicRanges)
# copyPeaks_GR <- GRanges(seqnames = one_cell[,"chr"], IRanges(as.numeric(one_cell[,"start"]), as.numeric(one_cell[,"end"])))
# mcols(copyPeaks_GR) <- one_cell[,4:12]

#add gene symbol
copyPeaks_GR <- GRanges(seqnames = all_annotation[,1], IRanges(as.numeric(all_annotation[,2]), as.numeric(all_annotation[,3])))
annotation_results = annotatePeak(peak=copyPeaks_GR, tssRegion=c(-3000, 3000), TxDb=txdb,annoDb = "org.Hs.eg.db",overlap = "all")
plotAnnoPie(annotation_results)
plotDistToTSS(annotation_results)

#scATAC-seq TSS region fraction
promoter_distribution<-data.frame(sample_type=c('CSCC3','CSCC2','CSCC1','AK3','AK2','AK1','Normal3','Normal2','Normal1'),Fraction_of_high_quality_fragments_overlapping_TSS=c(27.5,10.9,11.6,31.5,21.4,10.4,34,24,12.2))

library(ggpubr)
library(digest)
p<-ggbarplot(promoter_distribution, "sample_type", "Fraction_of_high_quality_fragments_overlapping_TSS",color = "sample_type", size = 1,width=0.7, palette =c("#56B4E9", 'gray', '#CCEBC5', '#BC80BD', '#FCCDE5','#BEBADA', '#FDB462',  '#FB8072','#8DD3C7','#FFFFB3','#B3DE69'),add = "jitter")
p





