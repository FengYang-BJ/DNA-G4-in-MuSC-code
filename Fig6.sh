#################Fig6B,bash
##############MAX example
bedtools inetrect -a chr_ASC_rep123.bed -b MAX.bed -wo |cut -f 1-3|sort|uniq|wc -l
#########################
#########################

#################Fig6C,bash
computeMatrix reference-point --referencePoint center -R 1.bed -S Merged_ASC48.mm10.fltd.uniq.RPGCnorm.bw MAX.bw YY1.mm10.pe.uniq.bw USF1.bw CTCF.bw MYOD1.bw CEBPB.bw FOSL1.bw  -b 2000 -a 2000 -o result_C2C12.gz -p 10
plotHeatmap -m result_C2C12.gz --colorMap OrRd -out result_G4_signal.pdf --heatmapHeight 14
############################

#################Fig6D,R
library(ggplot2)
library(ChIPseeker)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
x <- readPeakFile("MAX.bed")
peakAnno <- annotatePeak(x,tssRegion=c(-2000, 2000),TxDb=txdb)

print(peakAnno)
plotAnnoPie(peakAnno)
dev.off()
#############################

#################Fig6F,bash
computeMatrix reference-point --referencePoint center -R 1.bed -S Merged_ASC48.mm10.fltd.uniq.RPGCnorm.bw ASC_MAX.bw -b 2000 -a 2000 -o result_ASC.gz -p 10
plotHeatmap -m result_ASC.gz --colorMap OrRd -out result_G4_MAX_signal.pdf --heatmapHeight 14
#############################

#################Fig6G,bash
ngs.plot.r -G mm10 -R bed -C config.txt -O ./MAX_signal -N 0.5
#############################

#################Fig6H,bash
python3 bw_correlation.py loop_anchor.bed Merged_ASC48.mm10.fltd.uniq.RPGCnorm.bw ASC_MAX.bw signal_cor.pdf
#############################

#################Fig6I,bash
computeMatrix reference-point --referencePoint center -R chr_ASC_rep123.bed -S DMSO_MAX.bw PDS_MAX.bw -b 2000 -a 2000 -o result_MAX.gz -p 10
plotHeatmap -m result_MAX.gz --colorMap Blues -out result_MAX_DMSO_PDS_signal.pdf --heatmapHeight 14
#############################

#################Fig6J,bash
ngs.plot.r -G mm10 -R bed -C config_MAX.txt -O ./MAX_DMSO_PDS_signal -N 0.5
#############################


#################Fig6M
#################bash
bedtools makewindows -g mm10.genome -w 5000 > windows.bed
bedtools intersect -a windows.bed -b chr_ASC_rep123.bed -wo |cut -f 1-3|sort|uniq > G4_windows.bed
bedtools intersect -a G4_windows.bed -b chr_ASC48_MAX.narrowPeak -wo|cut -f 1-3|sort|uniq > G4_MAX_windows.bed
#####################

#################R
library(BiocParallel)
mm10_blacklist<-read.table("mm10.blacklist.bed")
colnames(mm10_blacklist)<-c("chromosome","start","end")
UT242_DMSO <- read.table("NC/matrix_5kb.txt", header = FALSE)
UT242_PDS <- read.table("siMAX/matrix_5kb.txt", header = FALSE)

# load("cnv.Rdata")
# combine cnv excluded regions with blacklist regions
# exclude <- rbind(cnv, mm10_blacklist)
exclude<-mm10_blacklist

chr <- paste0('chr', c(1:19, 'X')) # First chromosome, or use c(1:22, 'X') for all 

# Read data
chr_list    <- list()
for (i in 1:length(chr)) {
        print(i)
        UT242_DMSO_temp<-UT242_DMSO[UT242_DMSO$V1 == chr[i] & UT242_DMSO$V4 == chr[i], ]
        sparse_DMSO <- cooler2sparse(UT242_DMSO_temp)
        UT242_DPS_temp<-UT242_PDS[UT242_PDS$V1 == chr[i] & UT242_PDS$V4 == chr[i], ]
        sparse_PDS <- cooler2sparse(UT242_DPS_temp)
        exclude_temp<-exclude[exclude$chromosome==chr[i],]
        chr.table <- create.hic.table(sparse_DMSO, sparse_PDS, chr = chr[i],scale=FALSE, exclude.regions = exclude_temp, exclude.overlap = 0.2)
        chr_list[[i]] <- chr.table
}

chr_list <- total_sum(chr_list)
register(MulticoreParam(workers = 10), default = TRUE)
hic.list <- chr_list
hic.list <- hic_loess(hic.list, Plot=TRUE, parallel=TRUE)
save(hic.list,file='removed_hic.list.Rdata')
load("removed_hic.list.Rdata")

hic.list <- hic_compare(hic.list, A.min = 15, adjust.dist = TRUE, p.method = 'fdr', Plot = TRUE, parallel=TRUE)
hic.list <- do.call(rbind, hic.list)

G4_window_data<-read.table("G4_MAX_windows.bed",sep="\t")
result=data.frame()
for (i in 1:length(G4_window_data$V1)) {
        chr=G4_window_data$V1[i]
        start=G4_window_data$V2[i]
        end=G4_window_data$V3[i]
        temp=hic.list[((hic.list$chr1==chr)&(hic.list$start1==start)&(hic.list$end1==end))|((hic.list$chr2==chr)&(hic.list$start2==start)&(hic.list$end2==end))]
        if(length(temp$chr1)<1){
                next
        }
        else{
                result=rbind(result,temp)
        }
}

write.table(result,file="G4_MAX_window_NC_vs_siMAX_all.bed",row.names=FALSE,quote=FALSE,sep="\t")

result<-read.csv("G4_MAX_window_NC_vs_siMAX_all.bed",sep="\t")
result<-result[!duplicated(result), ]
result$type="NoDiff"
result$type[(result$adj.M>0.5)&(result$p.adj<0.1)]="Up"
result$type[(result$adj.M< -0.5)&(result$p.adj<0.1)]="Down"
result$type<-factor(result$type,levels<-c("Up","NoDiff","Down"))

valcano <- ggplot(data=result, aes(x=adj.M, y=-log10(p.adj),fill=type,color=type)) +
    geom_point(alpha=0.8, size=1) +
    theme_bw(base_size=15) +
    theme(
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()
    ) +
    ggtitle("DESeq2 Valcano") +
    geom_vline(xintercept=c(-0.5, 0.5), lty=2, col="gray", lwd=0.5) +
    geom_hline(yintercept=-log10(0.1), lty=1, col="gray", lwd=0.5)
# pdf("result_select_valcano15.pdf")
# valcano
# dev.off()
png("result_valcano15.png",height=1000,width=1000)
        valcano
dev.off()
write.table(result[result$type!="NoDiff",],file="G4_MAX_window_NC_vs_siMAX_all_diff.bed",row.names=FALSE,quote=FALSE,sep="\t")
###############################

#################Fig6M
#################bash
grep Down G4_MAX_window_NC_vs_siMAX_all_diff.bed |awk '{print $1"\t"$2"\t"$3"\n"$4"\t"$5"\t"$6"\n"}' - |sort|uniq > decreased_G4_MAX_bin.bed
bedtools intersect -a mm10_promoter_2k_2k.bed -b decreased_G4_MAX_bin.bed -wo |cut -f 1-6|sort|uniq > promoter_decreased_G4_MAX_bin.bed
######################
##################R
library(DESeq2)
sampleNames <- c("siNC1","siMax1","siNC2","siMax2")
data <- read.table("out.bed", header=TRUE, quote="\t", skip=1)
names(data)[7:10] <- sampleNames
countData <- as.matrix(data[7:10])
rownames(countData) <- data$Geneid
database <- data.frame(name=sampleNames, condition=c("siNC","siMax","siNC","siMax"))
rownames(database) <- sampleNames
dds <- DESeqDataSetFromMatrix(countData, colData=database, design= ~ condition)

dds <- DESeq(dds)
library(ggplot2)
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca<-ggplot(pcaData, aes(PC1, PC2, color=condition)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) +
geom_text(aes(label=name),vjust=2)+
xlim(-8,8)+
theme_bw(base_size=15) + 
theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
)
pdf("PCA.pdf")
pca
dev.off()

res <- lfcShrink(dds,coef="condition_siNC_vs_siMax", type="apeglm")
pdf("res_shrink.pdf")
plotMA(res, main = "Shrinkage by ashr", alpha=0.05, ylim=c(-4,4))
dev.off()
write.csv(res, "res_shrink_des_output.csv")

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
promoter<-read.table("mm10_promoter_1k_1k.bed")
resdata$genename="None"
for(i in 1:length(resdata$Row.names)){
    temp_gene=promoter$V5[promoter$V4==resdata$Row.names[i]]
    if(length(temp_gene)>0){
        resdata$genename[i]=temp_gene
    }
}
write.csv(resdata, "all_des_output.csv", row.names=FALSE)

resdata<-resdata[complete.cases(resdata),]
resdata$log2FoldChange<- -resdata$log2FoldChange
library(ggplot2)
resdata$change <- as.factor(
    ifelse(
        resdata$padj<0.05 & abs(resdata$log2FoldChange)>0.5,
        ifelse(resdata$log2FoldChange>0.5, "Up", "Down"),
        "NoDiff"
    )
)


decreased_MAX_G4_anchor_promoter<-read.table("MicroC/promoter_decreased_G4_MAX_bin.bed")

resdata_select<-resdata[resdata$Row.names %in% decreased_MAX_G4_anchor_promoter$V4,]
valcano_select <- ggplot(data=resdata_select, aes(x=log2FoldChange, y=-log10(padj), color=change)) + 
    geom_point(alpha=0.8, size=1) + 
    theme_bw(base_size=15) + 
    theme(
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()
    ) + 
    ggtitle("DESeq2 Valcano") + 
    scale_color_manual(name="", values=c("#F8766D", "#619CFF", "#00BA38"), limits=c("Up", "Down", "NoDiff")) + 
    geom_vline(xintercept=c(-0.5, 0.5), lty=2, col="gray", lwd=0.5) + 
    geom_hline(yintercept=-log10(0.05), lty=1, col="gray", lwd=0.5)+xlim(-5,5)
pdf("valcano_select_MAX_G4.pdf")
valcano_select
dev.off()
write.csv(resdata_select,"Gene_decreased_MAX_G4_anchor_promoter.csv",row.names =FALSE)
