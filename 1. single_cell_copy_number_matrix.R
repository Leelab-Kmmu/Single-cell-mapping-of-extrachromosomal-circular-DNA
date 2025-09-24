#!/bin/bash
#user

export PATH=cellranger-atac-2.1.0:$PATH

echo "******************************BEGIN******************************"
begin=$(date "+%Y-%m-%d %H:%M:%S")
echo "The script starts at ${begin}."

samples=('FKDL202607050-1a-SI-NA-C2')
dirs=('new_1836591_FKDN202559715-1A/')

for ((i=0; i<1; i++))
    do
        echo "--------------------------${samples[$i]} START--------------------------"
        starttime=$(date "+%Y-%m-%d %H:%M:%S")
        echo "Start at ${starttime}"
        cellranger-atac count --id=${samples[$i]} \
                              --reference=refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                              --fastqs=${dirs[$i]} \
                              --sample=${samples[$i]}
        finishtime=$(date "+%Y-%m-%d %H:%M:%S")
        echo "End at ${finishtime}"
        echo "---------------------------${samples[$i]} END---------------------------"
    done
doen=$(date "+%Y-%m-%d %H:%M:%S")
echo "The script ends at ${doen}."
echo "*****************************FINISH******************************"

########
########use AMULET to overlap count matrix,remove doublets, use CSCC3 as example

#####doublets calculation in CSCC,AK,Nomal
bash AMULET-v1.1/AMULET.sh possorted_bam.bam singlecell.csv AMULET-v1.1/human_autosomes.txt \
blacklist_repeats_segdups_rmsk_hg38.bed \
CSCC3 AMULET-v1.1 \
--maxinsertsize 30000000000 --mapqthresh 5


##overlap counts
conda activate pyscenic
cd CSCC3/outs
java -jar AMULET-v1.1/snATACOverlapCounter.jar possorted_bam.bam singlecell.csv AMULET-v1.1/human_autosomes.txt CSCC3 \
--maxinsertsize 300000000 --mapqthresh 5
cd CSCC3
python3 AMULET-v1.1/AMULET.py --rfilter blacklist_repeats_segdups_rmsk_hg38.bed Overlaps.txt OverlapSummary.txt CSCC3


#######make peak-cell matrix
options(stringsAsFactors = FALSE)
setwd('CSCC3')
overlap_count<-read.table(file='Overlaps.txt',sep='\t',header=T)

copy_nun_max<-list()
copy_nun_min<-list()
unique_cell_id<-unique(overlap_count$cell.id)
for(i in 1:length(unique_cell_id)){
copy_nun_max[[i]]<-as.numeric(overlap_count[overlap_count$cell.id%in%unique_cell_id[i],]$Max.Overlap.Count)
copy_nun_min[[i]]<-as.numeric(overlap_count[overlap_count$cell.id%in%unique_cell_id[i],]$Min.Overlap.Count)
}
names(copy_nun_max)<-unique_cell_id
names(copy_nun_min)<-unique_cell_id
max_mean<-sort(unlist(lapply(copy_nun_max,mean)),decreasing = T)
min_mean<-sort(unlist(lapply(copy_nun_min,mean)),decreasing = T)
save(max_mean,min_mean,file='copy/cell_mean_copy.Rdata')
save(overlap_count,file='Overlaps.Rdata')
######specific region sites
specific_sites<-data.frame()
chr_list<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18',
'chr19','chr20','chr21','chr22')
for(j in 1:22){
chr_overlap_count<-overlap_count[overlap_count$chr%in%chr_list[j],]
chr_overlap_count<-chr_overlap_count[order(chr_overlap_count$start),]
rownames(chr_overlap_count)<-1:dim(chr_overlap_count)[1]


line_index<-1
i<-1
while(i<=dim(chr_overlap_count)[1]){
site_start<-max(which(chr_overlap_count$start<=chr_overlap_count$end[i]))+1
line_index<-c(line_index,site_start)
i=site_start
}
line_index<-line_index[line_index<=dim(chr_overlap_count)[1]]
specific_sites<-rbind(specific_sites,chr_overlap_count[line_index,c('chr','start','end')])
}
write.table(specific_sites,file='copy/specific_sites.bed',sep='\t',row.names = F,col.names = F,quote=F)

black_region<-read.table(file='blacklist_repeats_segdups_rmsk_hg38.bed',sep='\t',header=F)
colnames(black_region)<-c('chr','start','end')########remove blacklist repeats

cd CSCC3/copy
bedtools subtract -A -a specific_sites.bed -b blacklist_repeats_segdups_rmsk_hg38.bed > specific_sites_filter.bed
######
devtools::install_github("PhanstielLab/bedtoolsr")
library(bedtoolsr)

#####copy num matrix
specific_sites_filter<-read.table(file='copy/specific_sites_filter.bed',sep='\t',header=F)
colnames(specific_sites_filter)<-c('chr','start','end')

copy_matrix_frame<-data.frame()
chr_list<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18',
'chr19','chr20','chr21','chr22')
for(j in 1:22){
temp_specific_sites_filter<-specific_sites_filter[specific_sites_filter$chr%in%chr_list[j],]
chr_overlap_count<-overlap_count[overlap_count$chr%in%chr_list[j],]
chr_overlap_count<-chr_overlap_count[order(chr_overlap_count$start),]
rownames(chr_overlap_count)<-1:dim(chr_overlap_count)[1]

filter_line_index<-match(temp_specific_sites_filter$start,chr_overlap_count$start)
end_index<-c()
for(k in 1:length(temp_specific_sites_filter$end)){
end_index[k]<-max(which(temp_specific_sites_filter$end[k]>=chr_overlap_count$start))
}

copy_matrix<-matrix(0,nrow=length(filter_line_index),ncol = length(unique_cell_id))
copy_matrix<-as.data.frame(copy_matrix)
colnames(copy_matrix)<-unique_cell_id
rownames(copy_matrix)<-paste0(chr_overlap_count[filter_line_index,]$chr,'_',chr_overlap_count[filter_line_index,]$start,'_',chr_overlap_count[filter_line_index,]$end)

for(i in 1:length(filter_line_index)){
if(i<length(filter_line_index)){
temp_mat<-chr_overlap_count[filter_line_index[i]:end_index[i],]
temp_ob<-with(temp_mat,tapply(Max.Overlap.Count,list (cell.id),sum))
copy_matrix[i,names(temp_ob)]<-temp_ob
}else{
temp_mat<-chr_overlap_count[filter_line_index[i]:end_index[i],]
temp_ob<-with(temp_mat,tapply(Max.Overlap.Count,list (cell.id),sum))
copy_matrix[i,names(temp_ob)]<-temp_ob
}
}
copy_matrix_frame<-rbind(copy_matrix_frame,copy_matrix)
}
copy_matrix<-copy_matrix_frame
#######remove doublets
doublet<-read.table(file='MultipletBarcodes_01.txt',sep='\t',header=F)
copy_matrix<-copy_matrix[,!colnames(copy_matrix)%in%c(as.character(doublet[,1]))]
col_sum<-apply(copy_matrix,2,function(x){length(which(x!=0))})
head(sort(col_sum,decreasing = T),100)
row_sum<-apply(copy_matrix,1,function(x){length(which(x>=3))})
head(sort(row_sum,decreasing = T),100)
save(copy_matrix,file='copy/copy_matrix.Rdata')

###CNV>=5 as eccDNA cutoff
options(stringsAsFactors = FALSE)
options(scipen = 999)
setwd('CSCC3/copy')
load('copy_matrix.Rdata')
df_sub <- copy_matrix[apply(copy_matrix, 1, function(x) any(x >= 5)), ]      
copy_matrix_eccDNA<-df_sub[, apply(df_sub, 2, function(x) any(x >= 5))]
save(copy_matrix_eccDNA,file='copy_matrix_eccDNA.Rdata')

######CNV statistic
options(stringsAsFactors = FALSE)
options(scipen = 999)
setwd('CSCC3/copy')
load('copy_matrix_eccDNA.Rdata')
CSCC3_copy<-unlist(apply(copy_matrix_eccDNA,2,function(x){length(which(x!=0))}))
CSCC3_allcopy<-unlist(apply(copy_matrix_eccDNA,2,function(x){sum(x)}))
names(CSCC3_copy)<-paste0('CSCC3_',colnames(copy_matrix_eccDNA))
names(CSCC3_allcopy)<-paste0('CSCC3_',colnames(copy_matrix_eccDNA))
setwd('CSCC2/copy')
load('copy_matrix_eccDNA.Rdata')
CSCC2_copy<-unlist(apply(copy_matrix_eccDNA,2,function(x){length(which(x!=0))}))
CSCC2_allcopy<-unlist(apply(copy_matrix_eccDNA,2,function(x){sum(x)}))
names(CSCC2_copy)<-paste0('CSCC2_',colnames(copy_matrix_eccDNA))
names(CSCC2_allcopy)<-paste0('CSCC2_',colnames(copy_matrix_eccDNA))
setwd('CSCC1/copy')
load('copy_matrix_eccDNA.Rdata')
CSCC1_copy<-unlist(apply(copy_matrix_eccDNA,2,function(x){length(which(x!=0))}))
CSCC1_allcopy<-unlist(apply(copy_matrix_eccDNA,2,function(x){sum(x)}))
names(CSCC1_copy)<-paste0('CSCC1_',colnames(copy_matrix_eccDNA))
names(CSCC1_allcopy)<-paste0('CSCC1_',colnames(copy_matrix_eccDNA))
setwd('Nomal3/copy')
load('copy_matrix_eccDNA.Rdata')
Normal3_copy<-unlist(apply(copy_matrix_eccDNA,2,function(x){length(which(x!=0))}))
Normal3_allcopy<-unlist(apply(copy_matrix_eccDNA,2,function(x){sum(x)}))
names(Normal3_copy)<-paste0('Nomal3_',colnames(copy_matrix_eccDNA))
names(Normal3_allcopy)<-paste0('Nomal3_',colnames(copy_matrix_eccDNA))
setwd('Nomal2/copy')
load('copy_matrix_eccDNA.Rdata')
Normal2_copy<-unlist(apply(copy_matrix_eccDNA,2,function(x){length(which(x!=0))}))
Normal2_allcopy<-unlist(apply(copy_matrix_eccDNA,2,function(x){sum(x)}))
names(Normal2_copy)<-paste0('Nomal2_',colnames(copy_matrix_eccDNA))
names(Normal2_allcopy)<-paste0('Nomal2_',colnames(copy_matrix_eccDNA))
setwd('Nomal1/copy')
load('copy_matrix_eccDNA.Rdata')
Normal1_copy<-unlist(apply(copy_matrix_eccDNA,2,function(x){length(which(x!=0))}))
Normal1_allcopy<-unlist(apply(copy_matrix_eccDNA,2,function(x){sum(x)}))
names(Normal1_copy)<-paste0('Nomal1_',colnames(copy_matrix_eccDNA))
names(Normal1_allcopy)<-paste0('Nomal1_',colnames(copy_matrix_eccDNA))
setwd('AK3/copy')
load('copy_matrix_eccDNA.Rdata')
AK3_copy<-unlist(apply(copy_matrix_eccDNA,2,function(x){length(which(x!=0))}))
AK3_allcopy<-unlist(apply(copy_matrix_eccDNA,2,function(x){sum(x)}))
names(AK3_copy)<-paste0('AK3_',colnames(copy_matrix_eccDNA))
names(AK3_allcopy)<-paste0('AK3_',colnames(copy_matrix_eccDNA))
setwd('AK2/copy')
load('copy_matrix_eccDNA.Rdata')
AK2_copy<-unlist(apply(copy_matrix_eccDNA,2,function(x){length(which(x!=0))}))
AK2_allcopy<-unlist(apply(copy_matrix_eccDNA,2,function(x){sum(x)}))
names(AK2_copy)<-paste0('AK2_',colnames(copy_matrix_eccDNA))
names(AK2_allcopy)<-paste0('AK2_',colnames(copy_matrix_eccDNA))
setwd('AK1/copy')
load('copy_matrix_eccDNA.Rdata')
AK1_copy<-unlist(apply(copy_matrix_eccDNA,2,function(x){length(which(x!=0))}))
AK1_allcopy<-unlist(apply(copy_matrix_eccDNA,2,function(x){sum(x)}))
names(AK1_copy)<-paste0('AK1_',colnames(copy_matrix_eccDNA))
names(AK1_allcopy)<-paste0('AK1_',colnames(copy_matrix_eccDNA))

CNV_copy<-c(CSCC3_copy,CSCC2_copy,CSCC1_copy,Normal3_copy,Normal2_copy,Normal1_copy,AK3_copy,AK2_copy,AK1_copy)
CNV_allcopy<-c(CSCC3_allcopy,CSCC2_allcopy,CSCC1_allcopy,Normal3_allcopy,Normal2_allcopy,Normal1_allcopy,AK3_allcopy,AK2_allcopy,AK1_allcopy)
save(CNV_copy,CNV_allcopy,file='CNV_statistic.Rdata')
