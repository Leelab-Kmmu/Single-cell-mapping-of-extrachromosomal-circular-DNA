#######For Figure 3A &3B
options(stringsAsFactors = FALSE)
options(scipen = 999)
setwd('/data2/home/lijie/scATAC_eccDNA/results/CSCC3/copy')
load('final_copy.Rdata')
load('/data2/home/lijie/scATAC_eccDNA/results/CSCC3/copy/cell_type.Rdata')
names(cell_type)[grep('TGGCAATCAGGCATCC-1',names(cell_type))]#PTMA
names(cell_type)[grep('TGATGCAGTCTAACCA-1',names(cell_type))]#NRP2
names(cell_type)[grep('CGCGCAACAGAGCCAA-1',names(cell_type))]#KLF7
names(cell_type)[grep('ACAAAGATCTCTGAGA-1',names(cell_type))]#MYC
names(cell_type)[grep('TATTGCTGTATCATGC-1',names(cell_type))]#PCAT1
names(cell_type)[grep('CGCGCAAGTTTAGACC-1',names(cell_type))]#PVT1



cell_name<-unlist(lapply(strsplit(names(cell_type),'_'), function(x){x[[2]]}))
names(cell_type)<-cell_name
table(cell_type[colnames(final_copy[,14:dim(final_copy)[2]])])



signac_cell<-intersect(colnames(final_copy)[14:dim(final_copy)[2]],cell_name)######signac QC
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
exist_row<-which(exist_cell%in%signac_cell==T)
exist_final_copy<-cbind(final_copy[exist_row,c(1:13)],final_copy[exist_row,signac_cell])
new_final_copy<-exist_final_copy
dim(new_final_copy)
table(cell_type[colnames(new_final_copy[,14:dim(new_final_copy)[2]])])

final_copy<-new_final_copy
###
new_exist_cell<-c()
for(i in 1:dim(new_final_copy)[1]){
temp_data<-new_final_copy[i,10:13]
if(length(grep('-1',temp_data))!=0){
temp_index<-strsplit(temp_data[grep('-1',temp_data)],':')
new_exist_cell[i]<-temp_index[[1]][[3]]
}else{
new_exist_cell[i]<-0
}
}
new_index_row<-c()
for(i in 1:dim(new_final_copy)[1]){
if(new_final_copy[i,new_exist_cell[i]]>=3){
new_index_row[i]<-TRUE
}else{
new_index_row[i]<-FALSE
}
}
new_frame<-cbind(new_final_copy[new_index_row,1:13],new_final_copy[new_index_row,unique(new_exist_cell[new_index_row])])
final_copy<-new_frame
dim(final_copy)


break_cell_point<-sort(table(new_exist_cell),decreasing = T)
head(break_cell_point,20)
table(cell_type[names(break_cell_point)])
cell_type[names(break_cell_point[1:100])]


break_cell_barcode<-'TGGCAATCAGGCATCC-1'#PTMA
break_cell_barcode<-'TGATGCAGTCTAACCA-1'#NRP2
break_cell_barcode<-'CGCGCAACAGAGCCAA-1'#KLF7
break_cell_barcode<-'ACAAAGATCTCTGAGA-1'#MYC
break_cell_barcode<-'TATTGCTGTATCATGC-1'#PCAT1
break_cell_barcode<-'CGCGCAAGTTTAGACC-1'#PVT1


one_cell<-cbind(final_copy[,c(1:13)],final_copy[,break_cell_barcode])
one_cell<-one_cell[one_cell[,14]!='0',]
colnames(one_cell)<-c('chr','start','end','break_chr','break_start','break_end','mapq','mate_chr','mate_start','break_cell_1','break_cell_2','break_cell_3','break_cell_4','copy_num')

one_longer_cell<-c()
for(i in 1:dim(one_cell)[1]){
temp_index<-one_cell[i,10:13]
if(length(grep(break_cell_barcode,temp_index))!=0){
one_longer_cell[i]<-temp_index[grep(break_cell_barcode,temp_index)]
}else{
one_longer_cell[i]<-0
}
}
one_cell<-one_cell[which(one_longer_cell!=0),]
one_cell<-as.data.frame(one_cell)
write.table(one_cell,file='MYC_one_cell.bed',sep='\t',row.names = F,col.names = T,quote = F)
one_cell<-read.table(file='MYC_one_cell.bed',sep='\t',header=T)

#chipseeker
setwd('/data2/home/lijie/scATAC_eccDNA/results/eccDNA_features')
library(ChIPseeker)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GeneSummary)
library(dplyr)
library(ggplot2)

#read data
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(GenomicRanges)
copyPeaks_GR <- GRanges(seqnames = one_cell[,"chr"], IRanges(as.numeric(one_cell[,"start"]), as.numeric(one_cell[,"end"])))
mcols(copyPeaks_GR) <- one_cell[,4:14]

#add gene symbol
# copyPeaks_GR <- GRanges(seqnames = final_copy[,1], IRanges(as.numeric(final_copy[,2]), as.numeric(final_copy[,3])))
# mcols(copyPeaks_GR) <- final_copy[,4:13]
annotation_results = annotatePeak(peak=copyPeaks_GR, tssRegion=c(-3000, 3000), TxDb=txdb,annoDb = "org.Hs.eg.db",overlap = "all")
annotation_results=as.data.frame(annotation_results)

#####OmicCircos
###############
library(OmicCircos)
library(circlize)
one_cell$SYMBOL<-annotation_results$SYMBOL
one_cell$mate_end<-one_cell$mate_start+1
circos_frame<-data.frame(one_cell[,c(1:3,11,14,15)])
circos_frame<-circos_frame[!duplicated(circos_frame),]

break_start<-one_cell[,c(4,5,6,11)]
break_start$copy_num<-rep(2,dim(break_start)[1])
break_start$SYMBOL<-rep('break_start',dim(break_start)[1])
colnames(break_start)<-colnames(circos_frame)
break_start<-break_start[!duplicated(break_start),]
mate_start<-one_cell[,c(8,9,16,11)]
mate_start$copy_num<-rep(2,dim(mate_start)[1])
mate_start$SYMBOL<-rep('mate_start',dim(mate_start)[1])
colnames(mate_start)<-colnames(circos_frame)
mate_start<-mate_start[!duplicated(mate_start),]

circos_frame<-rbind(circos_frame,break_start,mate_start)
break_frame<-data.frame(one_cell[,c(4,5,8,9,15)])
break_frame<-break_frame[!duplicated(break_frame),]

CNV_frame<-circos_frame[,-5]
colnames(CNV_frame)[1:3]<-c("seg.name","seg.Start","seg.End")
CNV_num<-circos_frame[!circos_frame$SYMBOL%in%'break_start'&!circos_frame$SYMBOL%in%'mate_start',c(1,2,5)]
CNV_num<-rbind(CNV_num,c('chrM',1,1))
colnames(CNV_num)<-c("seg.name","seg.po","name1")
link_data<-break_frame[,c(1,2,5,3,4,5,5)]
colnames(link_data)<-c("seg1","po1","name1","seg2","po2","name2","name3")
lable_frame<-circos_frame[!circos_frame$SYMBOL%in%'break_start'&!circos_frame$SYMBOL%in%'mate_start',c(1,2,6)]

# seg.f=sim.out$seg.frame
# seg.v=sim.out$seg.mapping
# link.v=sim.out$seg.link
# link.pg.v=sim.out$seg.link.pg
seg.num=length(unique(CNV_frame[,1]))

#segment(option)
seg.name=unique(CNV_frame$seg.name)
seg.Start=CNV_frame$seg.Start
db=segAnglePo(CNV_frame,seg=seg.name)
#
colors=rainbow(seg.num,alpha=0.5)
par(mar=c(2,2,2,2))
plot(c(1,800),c(1,800),type="n",axes=FALSE,xlab="",ylab="",main="")
#
circos(R=400,cir=db,type="chr",col=colors,print.chr.lab=TRUE,W=4,scale=TRUE)
##circos
circos(R=220,cir=db,W=100,mapping=CNV_num,col.v=3,type="ls",B=FALSE,col=colors
[1],lwd=2,scale=TRUE)
circos(R=150,cir=db,W=100,mapping=link_data,type="link",lwd=2,col=colors[c(1,7)])
circos(R=410,cir=db,W=50,mapping=lable_frame,type="label",side="out",col="blue",cex=0.4)

##
one_cell<-one_cell[1:2,]
one_cell$mate_end<-one_cell$mate_start+1
circos_frame<-data.frame(one_cell[,c(1:3,11,14,15)])
circos_frame<-circos_frame[!duplicated(circos_frame),]
break_start<-one_cell[,c(4,5,6,11)]
break_start$copy_num<-rep(2,dim(break_start)[1])
break_start$SYMBOL<-rep('break_start',dim(break_start)[1])
colnames(break_start)<-colnames(circos_frame)
break_start<-break_start[!duplicated(break_start),]
mate_start<-one_cell[,c(8,9,16,11)]
mate_start$copy_num<-rep(2,dim(mate_start)[1])
mate_start$SYMBOL<-rep('mate_start',dim(mate_start)[1])
colnames(mate_start)<-colnames(circos_frame)
mate_start<-mate_start[!duplicated(mate_start),]

circos_frame<-rbind(circos_frame,break_start,mate_start)
break_frame<-data.frame(one_cell[,c(4,5,8,9,15)])
break_frame<-break_frame[!duplicated(break_frame),]

CNV_frame<-circos_frame[,-5]
colnames(CNV_frame)[1:3]<-c("seg.name","seg.Start","seg.End")
CNV_num<-circos_frame[!circos_frame$SYMBOL%in%'break_start'&!circos_frame$SYMBOL%in%'mate_start',c(1,2,5)]
CNV_num_2<-circos_frame[!circos_frame$SYMBOL%in%'break_start'&!circos_frame$SYMBOL%in%'mate_start',c(1,3,5)]
colnames(CNV_num_2)<-c("seg.name","seg.po","name1")
CNV_num<-rbind(CNV_num,c('chrM',1,1))
colnames(CNV_num)<-c("seg.name","seg.po","name1")
CNV_num<-rbind(CNV_num,CNV_num_2)
link_data<-break_frame[,c(1,2,5,3,4,5,5)]
colnames(link_data)<-c("seg1","po1","name1","seg2","po2","name2","name3")
lable_frame<-circos_frame[!circos_frame$SYMBOL%in%'break_start'&!circos_frame$SYMBOL%in%'mate_start',c(1,2,6)]
break_location_frame<-circos_frame[circos_frame$SYMBOL%in%'break_start'|circos_frame$SYMBOL%in%'mate_start',c(1,2,3)]
chr_start_frame<-data.frame(chr=rep(circos_frame[1,1],2),starts=as.numeric(circos_frame[1,2:3]),ends=as.numeric(circos_frame[1,2:3]))

# seg.f=sim.out$seg.frame
# seg.v=sim.out$seg.mapping
# link.v=sim.out$seg.link
# link.pg.v=sim.out$seg.link.pg
seg.num=length(unique(CNV_frame[,1]))

#命名segment(option)
seg.name=unique(CNV_frame$seg.name)
db=segAnglePo(CNV_frame,seg=seg.name)
#设置颜色，设置图形范围
colors=rainbow(seg.num,alpha=0.5)
par(mar=c(2,2,2,2))
plot(c(1,800),c(1,800),type="n",axes=FALSE,xlab="",ylab="",main="")
#绘制外层坐标
circos(R=400,cir=db,type="chr",col=colors,print.chr.lab=TRUE,W=4,scale=TRUE)
##逐层绘制circos
circos(R=220,cir=db,W=100,mapping=CNV_num,col.v=3,type="ls",B=FALSE,col=colors
[1],lwd=2,scale=TRUE)
circos(R=150,cir=db,W=100,mapping=link_data,type="link",lwd=2,col=colors[c(1,7)])
circos(R=410,cir=db,W=50,mapping=lable_frame,type="label",side="out",col="blue",cex=0.4)
circos(R=150,cir=db,W=50,mapping=break_location_frame,type="label",side="out",col="blue",cex=0.4)
circos(R=210,cir=db,W=50,mapping=chr_start_frame,type="label",side="out",col="blue",cex=0.4)
