######link CNV and breakpoint,use CSCC3 as a example
options(stringsAsFactors = FALSE)
options(scipen = 999)
setwd('CSCC3/copy')

differ_chr<-read.table(file='possorted_bam_mate_differ_chr.txt',sep='\t',header=F,fill=T)
differ_chr_1<-data.frame(differ_chr$V1,differ_chr$V2,differ_chr$V2+1,differ_chr$V3,differ_chr[,c(4:5,7:10)])
colnames(differ_chr_1)<-c('chr','start','end','mapq','mate_chr','mate_start','V7','V8','V9','V10')
differ_chr_2<-data.frame(differ_chr$V4,differ_chr$V5-1,differ_chr$V5,differ_chr$V3,differ_chr$V1,differ_chr$V2,differ_chr[,7:10])
colnames(differ_chr_2)<-c('chr','start','end','mapq','mate_chr','mate_start','V7','V8','V9','V10')
differ_chr_frame<-rbind(differ_chr_1,differ_chr_2)
write.table(differ_chr_frame,file='differ_chr_frame.bed',sep='\t',row.names = F,col.names = F,quote = F)

mate_longer<-read.table(file='possorted_bam_mate_longer.txt',sep='\t',header=F,fill=T)
mate_longer_1<-data.frame(mate_longer$V1,mate_longer$V2,mate_longer$V2+1,mate_longer$V3,mate_longer$V1,mate_longer$V5,mate_longer[,7:10])
colnames(mate_longer_1)<-c('chr','start','end','mapq','mate_chr','mate_start','V7','V8','V9','V10')
mate_longer_2<-data.frame(mate_longer$V1,mate_longer$V5-1,mate_longer$V5,mate_longer$V3,mate_longer$V1,mate_longer$V2,mate_longer[,7:10])
colnames(mate_longer_2)<-c('chr','start','end','mapq','mate_chr','mate_start','V7','V8','V9','V10')
mate_longer_frame<-rbind(mate_longer_1,mate_longer_2)
write.table(mate_longer_frame,file='mate_longer_frame.bed',sep='\t',row.names = F,col.names = F,quote = F)

load('copy_matrix_ecDNA.Rdata')
chr_name<-unlist(lapply(strsplit(rownames(copy_matrix_ecDNA),'_'),function(x){return(x[[1]])}))
start_name<-as.numeric(unlist(lapply(strsplit(rownames(copy_matrix_ecDNA),'_'),function(x){return(x[[2]])})))
end_name<-as.numeric(unlist(lapply(strsplit(rownames(copy_matrix_ecDNA),'_'),function(x){return(x[[3]])})))
copy_mat<-data.frame(chr_name,start_name,end_name)
colnames(copy_mat)<-c('chr_name','start_name','end_name')

copy_mat$start_name<-copy_mat$start_name-2000
copy_mat$end_name<-copy_mat$end_name+2000
write.table(copy_mat,file='copy.bed',sep='\t',row.names = F,col.names = F,quote = F)
save(copy_mat,file='copy_mat.Rdata')
####
cd /data2/home/lijie/scATAC_ecDNA/results/AK1/copy
bedtools intersect -a copy.bed -b differ_chr_frame.bed -wa -wb > differ_chr_copy.bed
bedtools intersect -a copy.bed -b mate_longer_frame.bed -wa -wb > mate_longer_copy.bed


#####breakpoint located in two different chromosome.
options(stringsAsFactors = FALSE)
options(scipen = 999)
setwd('CSCC3/copy')
differ_chr_copy<-read.table(file='differ_chr_copy.bed',sep='\t',header=F)
load('copy_matrix_ecDNA.Rdata')
load('copy_mat.Rdata')
rownames(copy_matrix_ecDNA)<-paste0(copy_mat[,1],'_',copy_mat[,2],'_',copy_mat[,3])
temp_copy<-copy_matrix_ecDNA[unique(paste0(differ_chr_copy[,1],'_',differ_chr_copy[,2],'_',differ_chr_copy[,3])),]
all_copy<-cbind(differ_chr_copy,temp_copy[match(paste0(differ_chr_copy[,1],'_',differ_chr_copy[,2],'_',differ_chr_copy[,3]),rownames(temp_copy)),])
copy_name<-apply(all_copy[,14:as.numeric(dim(all_copy)[2])],1,function(x){names(x[x>=5])})#####cell copy number>=5
temp_differ_chr_cell<-c()
for(i in 1:dim(all_copy)[1]){
temp_data<-all_copy[i,10:13]
temp_index<-unlist(apply(temp_data,2,function(x){y<-strsplit(x,':');return(y[[1]][3])}))
if(length(grep('-1',temp_index))!=0){
temp_differ_chr_cell[i]<-temp_index[grep('-1',temp_index)]
}else{
temp_differ_chr_cell[i]<-0
}
}

cell_copy<-list()
for(i in 1:length(copy_name)){
cell_copy[[i]]<-na.omit(match(as.character(temp_differ_chr_cell[i]),unlist(copy_name[[i]])))
}
names(cell_copy)<-names(copy_name)
cell_info<-as.character(temp_differ_chr_cell[which(cell_copy!=0)])
cell_info<-unique(cell_info[grep('-1',cell_info)])
final_copy_differ_chr<-all_copy[which(cell_copy!=0),c(colnames(all_copy)[1:13],cell_info)]
colnames(final_copy_differ_chr)[1:13]<-c('chr','start','end','break_chr','break_start','break_end','mapq','mate_chr','mate_start','break_cell_1','break_cell_2','break_cell_3','break_cell_4')
save(final_copy_differ_chr,file='final_copy_differ_chr.Rdata')

col_sum<-apply(final_copy_differ_chr[,14:dim(final_copy_differ_chr)[2]],2,function(x){length(which(x!=0))})
head(sort(col_sum,decreasing = T),100)
row_sum<-apply(final_copy_differ_chr[,14:dim(final_copy_differ_chr)[2]],1,function(x){length(which(x>=5))})
head(sort(row_sum,decreasing = T),100)


#####break point located in the same chromosome,span>500bp
options(stringsAsFactors = FALSE)
options(scipen = 999)
setwd('CSCC3/copy')
mate_longer_copy<-read.table(file='mate_longer_copy.bed',sep='\t',header=F)
load('copy_matrix_ecDNA.Rdata')
load('copy_mat.Rdata')
rownames(copy_matrix_ecDNA)<-paste0(copy_mat[,1],'_',copy_mat[,2],'_',copy_mat[,3])
temp_copy<-copy_matrix_ecDNA[unique(paste0(mate_longer_copy[,1],'_',mate_longer_copy[,2],'_',mate_longer_copy[,3])),]
all_copy<-cbind(mate_longer_copy,temp_copy[match(paste0(mate_longer_copy[,1],'_',mate_longer_copy[,2],'_',mate_longer_copy[,3]),rownames(temp_copy)),])
copy_name<-apply(all_copy[,14:as.numeric(dim(all_copy)[2])],1,function(x){names(x[x>=5])})#####cell copy number>=5
temp_mate_longer_cell<-c()
for(i in 1:dim(all_copy)[1]){
temp_data<-all_copy[i,10:13]
temp_index<-unlist(apply(temp_data,2,function(x){y<-strsplit(x,':');return(y[[1]][3])}))
if(length(grep('-1',temp_index))!=0){
temp_mate_longer_cell[i]<-temp_index[grep('-1',temp_index)]
}else{
temp_mate_longer_cell[i]<-0
}
}

cell_copy<-list()
for(i in 1:length(copy_name)){
cell_copy[[i]]<-na.omit(match(as.character(temp_mate_longer_cell[i]),unlist(copy_name[[i]])))
}

cell_info<-as.character(temp_mate_longer_cell[which(cell_copy!=0)])
cell_info<-unique(cell_info[grep('-1',cell_info)])
final_copy_mate_longer<-all_copy[which(cell_copy!=0),c(colnames(all_copy)[1:13],cell_info)]
final_copy_mate_longer_matrix<-final_copy_mate_longer
colnames(final_copy_mate_longer_matrix)[1:13]<-c('chr','start','end','break_chr','break_start','break_end','mapq','mate_chr','mate_start','break_cell_1','break_cell_2','break_cell_3','break_cell_4')
save(final_copy_mate_longer_matrix,file='final_copy_mate_longer_matrix.Rdata')


col_sum<-apply(final_copy_mate_longer_matrix[,14:dim(final_copy_mate_longer_matrix)[2]],2,function(x){length(which(x!=0))})
head(sort(col_sum,decreasing = T),100)
row_sum<-apply(final_copy_mate_longer_matrix[,14:dim(final_copy_mate_longer_matrix)[2]],1,function(x){length(which(x>=5))})
head(sort(row_sum,decreasing = T),100)



#####form a eccDNA cycle
options(stringsAsFactors = FALSE)
options(scipen = 999)
setwd('CSCC3/copy')
load('final_copy_mate_longer_matrix.Rdata')
load('final_copy_differ_chr.Rdata')
name_col<-unique(c(colnames(final_copy_mate_longer_matrix),colnames(final_copy_differ_chr)))
name_row<-unique(c(rownames(final_copy_mate_longer_matrix),rownames(final_copy_differ_chr)))
final_copy<-matrix(0,nrow=length(name_row),ncol=length(name_col))
rownames(final_copy)<-name_row
colnames(final_copy)<-name_col
final_copy[rownames(final_copy_mate_longer_matrix),colnames(final_copy_mate_longer_matrix)]<-as.matrix(final_copy_mate_longer_matrix)
final_copy[rownames(final_copy_differ_chr),colnames(final_copy_differ_chr)]<-as.matrix(final_copy_differ_chr)


exist_cell<-c()
for(i in 1:dim(final_copy)[1]){
temp_data<-final_copy[i,10:13]
if(length(grep('-1',temp_data))!=0){
temp_index<-strsplit(temp_data[grep('-1',temp_data)],':')
exist_cell[i]<-temp_index[[1]][[3]]
}else{
exist_cell[i]<-0
}
}
break_cell_point<-sort(table(exist_cell),decreasing = T)
head(break_cell_point,20)

retain_index<-c()
for(i in 1:length(exist_cell)){
if(final_copy[i,exist_cell[i]]=='0'){
retain_index[i]<-FALSE
}else{
retain_index[i]<-TRUE
}
}
table(retain_index)


one_cell_copy<-sort(apply(final_copy[,14:dim(final_copy)[2]],2,function(x){max(as.numeric(x))}),decreasing=T)
head(one_cell_copy,15)
copy_cells<-sort(apply(final_copy[,14:dim(final_copy)[2]],1,function(x){max(as.numeric(x))}),decreasing=T)
head(copy_cells,15)
colnames(final_copy)[1:13]<-c('chr','start','end','break_chr','break_start','break_end','mapq','mate_chr','mate_start','break_cell_1','break_cell_2','break_cell_3','break_cell_4')
save(final_copy,file='final_copy.Rdata')

#######break point statistic
options(stringsAsFactors = FALSE)
options(scipen = 999)
break_cell<-function(temp_copy){
exist_cell<-c()
for(i in 1:dim(temp_copy)[1]){
temp_data<-temp_copy[i,10:13]
if(length(grep('-1',temp_data))!=0){
temp_index<-strsplit(as.character(temp_data[grep('-1',temp_data)]),':')
exist_cell[i]<-temp_index[[1]][[3]]
}else{
exist_cell[i]<-0
}
}
return(exist_cell)
}
setwd('/data2/home/lijie/scATAC_ecDNA/results/CSCC3/copy')
load('final_copy.Rdata')
temp_name<-paste0(final_copy[,4],'_',final_copy[,5],'_',final_copy[,9])
temp_copy<-final_copy[match(unique(temp_name),temp_name),]
CSCC3_breakpoint<-break_cell(temp_copy)
CSCC3_breakpoint<-table(CSCC3_breakpoint)
names(CSCC3_breakpoint)<-paste0('CSCC3_',names(CSCC3_breakpoint))
setwd('/data2/home/lijie/scATAC_ecDNA/results/CSCC2/copy')
load('final_copy.Rdata')
temp_name<-paste0(final_copy[,4],'_',final_copy[,5],'_',final_copy[,9])
temp_copy<-final_copy[match(unique(temp_name),temp_name),]
CSCC2_breakpoint<-break_cell(temp_copy)
CSCC2_breakpoint<-table(CSCC2_breakpoint)
names(CSCC2_breakpoint)<-paste0('CSCC2_',names(CSCC2_breakpoint))
setwd('/data2/home/lijie/scATAC_ecDNA/results/CSCC1/copy')
load('final_copy.Rdata')
temp_name<-paste0(final_copy[,4],'_',final_copy[,5],'_',final_copy[,9])
temp_copy<-final_copy[match(unique(temp_name),temp_name),]
CSCC1_breakpoint<-break_cell(temp_copy)
CSCC1_breakpoint<-table(CSCC1_breakpoint)
names(CSCC1_breakpoint)<-paste0('CSCC1_',names(CSCC1_breakpoint))
setwd('/data2/home/lijie/scATAC_ecDNA/results/AK3/copy')
load('final_copy.Rdata')
temp_name<-paste0(final_copy[,4],'_',final_copy[,5],'_',final_copy[,9])
temp_copy<-final_copy[match(unique(temp_name),temp_name),]
AK3_breakpoint<-break_cell(temp_copy)
AK3_breakpoint<-table(AK3_breakpoint)
names(AK3_breakpoint)<-paste0('AK3_',names(AK3_breakpoint))
setwd('/data2/home/lijie/scATAC_ecDNA/results/AK2/copy')
load('final_copy.Rdata')
temp_name<-paste0(final_copy[,4],'_',final_copy[,5],'_',final_copy[,9])
temp_copy<-final_copy[match(unique(temp_name),temp_name),]
AK2_breakpoint<-break_cell(temp_copy)
AK2_breakpoint<-table(AK2_breakpoint)
names(AK2_breakpoint)<-paste0('AK2_',names(AK2_breakpoint))
setwd('/data2/home/lijie/scATAC_ecDNA/results/AK1/copy')
load('final_copy.Rdata')
temp_name<-paste0(final_copy[,4],'_',final_copy[,5],'_',final_copy[,9])
temp_copy<-final_copy[match(unique(temp_name),temp_name),]
AK1_breakpoint<-break_cell(temp_copy)
AK1_breakpoint<-table(AK1_breakpoint)
names(AK1_breakpoint)<-paste0('AK1_',names(AK1_breakpoint))
setwd('/data2/home/lijie/scATAC_ecDNA/results/Nomal3/copy')
load('final_copy.Rdata')
temp_name<-paste0(final_copy[,4],'_',final_copy[,5],'_',final_copy[,9])
temp_copy<-final_copy[match(unique(temp_name),temp_name),]
Nomal3_breakpoint<-break_cell(temp_copy)
Nomal3_breakpoint<-table(Nomal3_breakpoint)
names(Nomal3_breakpoint)<-paste0('Nomal3_',names(Nomal3_breakpoint))
setwd('/data2/home/lijie/scATAC_ecDNA/results/Nomal2/copy')
load('final_copy.Rdata')
temp_name<-paste0(final_copy[,4],'_',final_copy[,5],'_',final_copy[,9])
temp_copy<-final_copy[match(unique(temp_name),temp_name),]
Nomal2_breakpoint<-break_cell(temp_copy)
Nomal2_breakpoint<-table(Nomal2_breakpoint)
names(Nomal2_breakpoint)<-paste0('Nomal2_',names(Nomal2_breakpoint))
setwd('/data2/home/lijie/scATAC_ecDNA/results/Nomal1/copy')
load('final_copy.Rdata')
temp_name<-paste0(final_copy[,4],'_',final_copy[,5],'_',final_copy[,9])
temp_copy<-final_copy[match(unique(temp_name),temp_name),]
Nomal1_breakpoint<-break_cell(temp_copy)
Nomal1_breakpoint<-table(Nomal1_breakpoint)
names(Nomal1_breakpoint)<-paste0('Nomal1_',names(Nomal1_breakpoint))

adapter_num<-c(CSCC3_breakpoint,CSCC2_breakpoint,CSCC1_breakpoint,AK3_breakpoint,AK2_breakpoint,AK1_breakpoint,Nomal3_breakpoint,Nomal2_breakpoint,Nomal1_breakpoint)
save(adapter_num,file='/data2/home/lijie/scATAC_ecDNA/results/ecDNA_features/adapter_num.Rdata')


##############ecDNA number statistic
options(stringsAsFactors = FALSE)
options(scipen = 999)
break_cell<-function(new_annotation){
exist_cell<-c()
for(i in 1:dim(new_annotation)[1]){
temp_data<-new_annotation[i,12:15]
if(length(grep('-1',temp_data))!=0){
temp_index<-strsplit(as.character(temp_data[grep('-1',temp_data)]),':')
exist_cell[i]<-temp_index[[1]][[3]]
}else{
exist_cell[i]<-0
}
}
return(exist_cell)
}
setwd('/data2/home/lijie/scATAC_ecDNA/results/CSCC3/copy')
load('annotation_results.Rdata')
CSCC3_breakpoint<-break_cell(new_annotation)
new_annotation$break_cell<-CSCC3_breakpoint
cell_id<-unique(new_annotation$break_cell)
ecDNA_gene_num<-c()
for(i in 1:length(cell_id)){
temp_frame<-new_annotation[new_annotation$break_cell%in%cell_id[i],]
ecDNA_gene_num[i]<-length(unique(temp_frame$SYMBOL))
}
CSCC3_ecDNA_gene_num<-ecDNA_gene_num
names(CSCC3_ecDNA_gene_num)<-paste0('CSCC3_',cell_id)

setwd('/data2/home/lijie/scATAC_ecDNA/results/CSCC2/copy')
load('annotation_results.Rdata')
CSCC2_breakpoint<-break_cell(new_annotation)
new_annotation$break_cell<-CSCC2_breakpoint
cell_id<-unique(new_annotation$break_cell)
ecDNA_gene_num<-c()
for(i in 1:length(cell_id)){
temp_frame<-new_annotation[new_annotation$break_cell%in%cell_id[i],]
ecDNA_gene_num[i]<-length(unique(temp_frame$SYMBOL))
}
CSCC2_ecDNA_gene_num<-ecDNA_gene_num
names(CSCC2_ecDNA_gene_num)<-paste0('CSCC2_',cell_id)

setwd('/data2/home/lijie/scATAC_ecDNA/results/CSCC1/copy')
load('annotation_results.Rdata')
CSCC1_breakpoint<-break_cell(new_annotation)
new_annotation$break_cell<-CSCC1_breakpoint
cell_id<-unique(new_annotation$break_cell)
ecDNA_gene_num<-c()
for(i in 1:length(cell_id)){
temp_frame<-new_annotation[new_annotation$break_cell%in%cell_id[i],]
ecDNA_gene_num[i]<-length(unique(temp_frame$SYMBOL))
}
CSCC1_ecDNA_gene_num<-ecDNA_gene_num
names(CSCC1_ecDNA_gene_num)<-paste0('CSCC1_',cell_id)

setwd('/data2/home/lijie/scATAC_ecDNA/results/AK3/copy')
load('annotation_results.Rdata')
AK3_breakpoint<-break_cell(new_annotation)
new_annotation$break_cell<-AK3_breakpoint
cell_id<-unique(new_annotation$break_cell)
ecDNA_gene_num<-c()
for(i in 1:length(cell_id)){
temp_frame<-new_annotation[new_annotation$break_cell%in%cell_id[i],]
ecDNA_gene_num[i]<-length(unique(temp_frame$SYMBOL))
}
AK3_ecDNA_gene_num<-ecDNA_gene_num
names(AK3_ecDNA_gene_num)<-paste0('AK3_',cell_id)

setwd('/data2/home/lijie/scATAC_ecDNA/results/AK2/copy')
load('annotation_results.Rdata')
AK2_breakpoint<-break_cell(new_annotation)
new_annotation$break_cell<-AK2_breakpoint
cell_id<-unique(new_annotation$break_cell)
ecDNA_gene_num<-c()
for(i in 1:length(cell_id)){
temp_frame<-new_annotation[new_annotation$break_cell%in%cell_id[i],]
ecDNA_gene_num[i]<-length(unique(temp_frame$SYMBOL))
}
AK2_ecDNA_gene_num<-ecDNA_gene_num
names(AK2_ecDNA_gene_num)<-paste0('AK2_',cell_id)

setwd('/data2/home/lijie/scATAC_ecDNA/results/AK1/copy')
load('annotation_results.Rdata')
AK1_breakpoint<-break_cell(new_annotation)
new_annotation$break_cell<-AK1_breakpoint
cell_id<-unique(new_annotation$break_cell)
ecDNA_gene_num<-c()
for(i in 1:length(cell_id)){
temp_frame<-new_annotation[new_annotation$break_cell%in%cell_id[i],]
ecDNA_gene_num[i]<-length(unique(temp_frame$SYMBOL))
}
AK1_ecDNA_gene_num<-ecDNA_gene_num
names(AK1_ecDNA_gene_num)<-paste0('AK1_',cell_id)

setwd('/data2/home/lijie/scATAC_ecDNA/results/Nomal3/copy')
load('annotation_results.Rdata')
Nomal3_breakpoint<-break_cell(new_annotation)
new_annotation$break_cell<-Nomal3_breakpoint
cell_id<-unique(new_annotation$break_cell)
ecDNA_gene_num<-c()
for(i in 1:length(cell_id)){
temp_frame<-new_annotation[new_annotation$break_cell%in%cell_id[i],]
ecDNA_gene_num[i]<-length(unique(temp_frame$SYMBOL))
}
Nomal3_ecDNA_gene_num<-ecDNA_gene_num
names(Nomal3_ecDNA_gene_num)<-paste0('Nomal3_',cell_id)

setwd('/data2/home/lijie/scATAC_ecDNA/results/Nomal2/copy')
load('annotation_results.Rdata')
Nomal2_breakpoint<-break_cell(new_annotation)
new_annotation$break_cell<-Nomal2_breakpoint
cell_id<-unique(new_annotation$break_cell)
ecDNA_gene_num<-c()
for(i in 1:length(cell_id)){
temp_frame<-new_annotation[new_annotation$break_cell%in%cell_id[i],]
ecDNA_gene_num[i]<-length(unique(temp_frame$SYMBOL))
}
Nomal2_ecDNA_gene_num<-ecDNA_gene_num
names(Nomal2_ecDNA_gene_num)<-paste0('Nomal2_',cell_id)

setwd('/data2/home/lijie/scATAC_ecDNA/results/Nomal1/copy')
load('annotation_results.Rdata')
Nomal1_breakpoint<-break_cell(new_annotation)
new_annotation$break_cell<-Nomal1_breakpoint
cell_id<-unique(new_annotation$break_cell)
ecDNA_gene_num<-c()
for(i in 1:length(cell_id)){
temp_frame<-new_annotation[new_annotation$break_cell%in%cell_id[i],]
ecDNA_gene_num[i]<-length(unique(temp_frame$SYMBOL))
}
Nomal1_ecDNA_gene_num<-ecDNA_gene_num
names(Nomal1_ecDNA_gene_num)<-paste0('Nomal1_',cell_id)

ecDNA_gene_num<-c(CSCC3_ecDNA_gene_num,CSCC2_ecDNA_gene_num,CSCC1_ecDNA_gene_num,AK3_ecDNA_gene_num,AK2_ecDNA_gene_num,AK1_ecDNA_gene_num,Nomal3_ecDNA_gene_num,Nomal2_ecDNA_gene_num,Nomal1_ecDNA_gene_num)
save(ecDNA_gene_num,file='/data2/home/lijie/scATAC_ecDNA/results/ecDNA_features/ecDNA_gene_num.Rdata')


