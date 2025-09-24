######For Figure 6A &Figure S7A
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

eccDNA_sample_pca <- RunHarmony(object = eccDNA_sample_pca, group.by.vars = "orig.ident")

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

######For Figure 6D &Figure S7B
FeaturePlot(object = eccDNA_sample_umap,cols=c('#BEBEBE','#FF0000'),order=T,raster = F,features = 'nCount_RNA',max.cutoff=50)
FeaturePlot(object = eccDNA_sample_umap,cols=c('#BEBEBE','#FF0000'),order=T,raster = F,features = 'nFeature_RNA',max.cutoff=15)

######For Figure 6B
#########0,1,3 trajectory
cell_index<-c(names(eccDNA_sample_umap@active.ident[eccDNA_sample_umap@active.ident==0]),names(eccDNA_sample_umap@active.ident[eccDNA_sample_umap@active.ident==1]),names(eccDNA_sample_umap@active.ident[eccDNA_sample_umap@active.ident==3]))
trajectory_umap<-subset(eccDNA_sample_umap,cells=cell_index)
DimPlot(object = trajectory_umap, reduction = "umap")

######
cluster_eccDNA_markers <- FindAllMarkers(eccDNA_sample_umap, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1)
save(cluster_eccDNA_markers,file='cluster_eccDNA_markers.Rdata')
load('cluster_eccDNA_markers.Rdata')
top1 <- cluster_eccDNA_markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)

central_cell<-c()
for(i in 1:length(table(top1$cluster))){
cluster_0<-names(eccDNA_sample_umap@active.ident[eccDNA_sample_umap@active.ident==c(i-1)])
temp_eccDNA<-eccDNA_frame[top1[top1$cluster==c(i-1),]$gene[1],cluster_0]
final_eccDNA<-as.numeric(temp_eccDNA)
names(final_eccDNA)<-colnames(temp_eccDNA)
central_cell[i]<-names(final_eccDNA[final_eccDNA==max(final_eccDNA)])[1]
}
central_location<-as.data.frame(eccDNA_sample_umap@reductions[["umap"]]@cell.embeddings[na.omit(central_cell),])
central_location$cluster<-as.character(c(match(rownames(central_location),central_cell)-1))

#########
cluster_eccDNA<-data.frame()
for(i in 1:length(table(top1$cluster))){
cluster_0<-names(eccDNA_sample_umap@active.ident[eccDNA_sample_umap@active.ident==c(i-1)])
temp_eccDNA<-eccDNA_frame[,cluster_0]
final_eccDNA<-apply(temp_eccDNA,1,sum)
cluster_eccDNA<-rbind(cluster_eccDNA,t(final_eccDNA))
}
cluster_eccDNA<-t(cluster_eccDNA)
colnames(cluster_eccDNA)<-0:c(dim(cluster_eccDNA)[2]-1)
gene_count<-log10(apply(cluster_eccDNA,1,sum))
cluster_count<-apply(cluster_eccDNA,1,function(x){length(which(x!=0))})
#####For Figure 6C
library(ggplot2)
library(ggrepel)
dfm<-data.frame(gene_count=gene_count,cluster_count=cluster_count,name=names(cluster_count))
data1 = dfm[dfm$gene_count>=2&dfm$cluster_count>=6,]
p<-ggplot(dfm,aes(x=gene_count,y=cluster_count))+geom_point(size=0.5,col='grey')+
geom_text_repel(
    data = dfm[dfm$gene_count>=2&dfm$cluster_count>=6,],
    aes(label = name),
    size = 3,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+geom_point(data=data1,aes(x=gene_count,y=cluster_count),size=0.5,col='red')
p+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

eccDNA_count<-apply(cluster_eccDNA,2,sum)
central_location$eccDNA_count<-eccDNA_count[match(central_location$cluster,names(eccDNA_count))]
# cluster_eccDNA[cluster_eccDNA!=0]<-1
pearson_cluster<-cor(cluster_eccDNA)
high_cor <- which(abs(pearson_cluster) > 0.1 & abs(pearson_cluster) < 1, arr.ind = TRUE)


result <- data.frame(
  row = rownames(pearson_cluster)[high_cor[, 1]],
  col = colnames(pearson_cluster)[high_cor[, 2]],
  correlation = pearson_cluster[high_cor]
)
result <-result[duplicated(result[,3]),]
write.table(result,file='cluster_eccDNA_clone.txt',sep='\t',row.names = F,col.names = T,quote = F)

link_pot<-data.frame()
for(i in 1:dim(result)[1]){
link_temp<-rbind(central_location[central_location$cluster==result[i,1],],central_location[central_location$cluster==result[i,2],])
link_pot<-rbind(link_pot,link_temp)
}
link_pot$correlation<-rep(result$correlation,each=2)


umap_plot <- DimPlot(object = eccDNA_sample_umap, reduction = "umap")


library(ggplot2)
library(igraph)
library(ggraph)
# 
first_link_pot<-link_pot[1:2,]
first_link_pot<-first_link_pot[order(first_link_pot$eccDNA_count),]
if(first_link_pot$umap_1[1]<first_link_pot$umap_1[2]){derection<-'last'}else{derection<-'first'}
final_plot <- umap_plot +
geom_point(data = central_location, aes(x = umap_1, y = umap_2), color = "red", size = log10(central_location$eccDNA_count)) +  # 添加红色点
geom_line(data = first_link_pot, aes(x = umap_1, y = umap_2), color = "grey", size = 0.5,arrow = arrow(type = "open", length = unit(0.1, "inches"), angle = 35,ends = derection)) +  # 添加折线
  geom_text(data = central_location,aes(x = umap_1, y = umap_2,label = cluster), vjust = -1, hjust = 0.5, color = "black", size = 3)+
  labs(title = "UMAP with Custom Trajectory")
final_plot
  
for(i in 2:c(dim(link_pot)[1]/2)){
j<-2*i
m<-j-1
order_link_pot<-link_pot[m:j,]
order_link_pot<-order_link_pot[order(order_link_pot$eccDNA_count),]
if(order_link_pot$umap_1[1]<order_link_pot$umap_1[2]){derection<-'last'}else{derection<-'first'}
final_plot<-final_plot+geom_line(data = order_link_pot, aes(x = umap_1, y = umap_2), color = "grey", size = 0.5,arrow = arrow(type = "open", length = unit(0.1, "inches"), angle = 20,ends = derection))
}
final_plot +theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

#####For Figure 6E
k<-c(3,1,0)
group_eccDNA<-rownames(cluster_eccDNA)[which(cluster_eccDNA[,c(k[2]+1)]-cluster_eccDNA[,c(k[1]+1)]>1&cluster_eccDNA[,c(k[3]+1)]-cluster_eccDNA[,c(k[2]+1)]>1)]
cell_group<-c()
annotation_col<-c()
for(i in 1:length(k)){
temp_group<-names(eccDNA_sample_umap@active.ident[eccDNA_sample_umap@active.ident==k[i]])
annotation_col<-c(annotation_col,rep(k[i],length(eccDNA_sample_umap@active.ident[eccDNA_sample_umap@active.ident==k[i]])))
cell_group<-c(cell_group,temp_group)
}

heat_frame<-eccDNA_frame[group_eccDNA,cell_group]
heat_frame[heat_frame!=0]<-1

library(pheatmap)
annotation_col<- data.frame(Type = factor(annotation_col,levels=k))
rownames(annotation_col)<- colnames(heat_frame)
pheatmap(heat_frame,cluster_cols = F,cluster_rows = T,show_colnames = F,show_rownames = T,color = colorRampPalette(c("white","#D45252"))(50),annotation_col=annotation_col)

#####For Figure 6G
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

eccDNA_gene<-unique(c(granulosum_gene,basal_gene,spinosum_gene,melanocyte_gene,follicular_gene))
save(granulosum_gene,basal_gene,spinosum_gene,melanocyte_gene,follicular_gene,file='eccDNA_gene.Rdata')

a<-venn.diagram(list(basal_gene=basal_gene,spinosum_gene=spinosum_gene,granulosum_gene=granulosum_gene,follicular_gene=follicular_gene,melanocyte_gene=melanocyte_gene),
fill=c("red","green","blue",'yellow','grey'), alpha=c(0.5,0.5,0.5,0.5,0.5), cex=2, cat.fontface=4, fontfamily=3, filename=NULL)
grid.draw(a)

#####For Figure 6H
library(pheatmap)
R_frame<-data.frame(basal_eccDNA<-c(length(intersect(basal_gene,basal_gene)),length(intersect(basal_gene,spinosum_gene)),length(intersect(basal_gene,granulosum_gene)),length(intersect(basal_gene,follicular_gene)),length(intersect(basal_gene,melanocyte_gene))),
spinosum_eccDNA<-c(length(intersect(spinosum_gene,basal_gene)),length(intersect(spinosum_gene,spinosum_gene)),length(intersect(spinosum_gene,granulosum_gene)),length(intersect(spinosum_gene,follicular_gene)),length(intersect(spinosum_gene,melanocyte_gene))),
granulosum_eccDNA<-c(length(intersect(granulosum_gene,basal_gene)),length(intersect(granulosum_gene,spinosum_gene)),length(intersect(granulosum_gene,granulosum_gene)),length(intersect(granulosum_gene,follicular_gene)),length(intersect(granulosum_gene,melanocyte_gene))),
follicular_eccDNA<-c(length(intersect(follicular_gene,basal_gene)),length(intersect(follicular_gene,spinosum_gene)),length(intersect(follicular_gene,granulosum_gene)),length(intersect(follicular_gene,follicular_gene)),length(intersect(follicular_gene,melanocyte_gene))),
melanocyte_eccDNA<-c(length(intersect(melanocyte_gene,basal_gene)),length(intersect(melanocyte_gene,spinosum_gene)),length(intersect(melanocyte_gene,granulosum_gene)),length(intersect(melanocyte_gene,follicular_gene)),length(intersect(melanocyte_gene,melanocyte_gene))))


pheatmap(R_frame,cluster_cols = F,cluster_rows = F,
color = colorRampPalette(c('#FFFAFA',"#D45252"))(50),show_colnames = F,show_rownames = F,display_numbers = T)


