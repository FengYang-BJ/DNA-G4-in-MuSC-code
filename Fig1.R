#######Fig 1B
library("ggplot2")
data<-read.table("data_merge.bed",sep="\t",header=T)
data$Type<-factor(data$Type, levels = c("FISC", "ASC-24h", "ASC-48h","DSC"))
pdf("peak_merge.pdf")
ggplot(data,aes(Type,Number,fill=Type))+geom_bar(position = "dodge",stat = "identity",width = 0.5)+theme_classic()
#theme_bw() + theme(panel.grid=element_blank())
dev.off()

########Fig 1C
grep chr ASC_rep123.bed > chr_ASC_rep123.bed
bedtools getfasta -fi mm10.fa -bed chr_ASC_rep123.bed -fo chr_ASC_rep123.fa
meme-chip -meme-p 20 -oc 123_Denoval -dreme-m 20 -meme-nmotifs 20 chr_ASC_rep123.fa -meme-maxw 20 -meme-minw 6

########Fig 1D
bedtools intersect -a chr_FISC_rep123.bed -b chr_ASC_rep123.bed -wo |cut -f 1-3|sort|uniq|wc -l

########Fig 1H
library(ggplot2)
library(ChIPseeker)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
x <- readPeakFile("chr_ASC_rep123.bed")
peakAnno <- annotatePeak(x,tssRegion=c(-2000, 2000),TxDb=txdb)

print(peakAnno)
pdf("pie.pdf")
names = c("Promoter", "Exon","Intron","Distal Intergenic")
info<-c(peakAnno@annoStat[1,2]+peakAnno@annoStat[2,2],peakAnno@annoStat[3,2]+peakAnno@annoStat[4,2]+peakAnno@annoStat[5,2]+peakAnno@annoStat[6,2],peakAnno@annoStat[7,2]+peakAnno@annoStat[8,2],peakAnno@annoStat[9,2]+peakAnno@annoStat[10,2])
cols<-c("#A6CEE3","#E31A1C","#FF7F00","#6A3D9A")
print(info)
pie(info, labels=names,col=cols)
dev.off()

########Fig 1J,K,M,N,O analysis pipeline is similar with mouse, using hg38 to map reads and annotation. 