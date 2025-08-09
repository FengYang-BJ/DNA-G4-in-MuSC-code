################Fig4I,bash
###########call loop
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ../DMSO/new_loop/5k_loop.bed > DMSO_loop_6.bed
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ../PDS/new_loop/5k_loop.bed > PDS_loop_6.bed
sort DMSO_6.bed PDS_6.bed PDS_6.bed | uniq -u|wc -l
################

################Fig4J,bash
grep G4 ../DMSO/new_loop/5k_loop.bed|awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' - > DMSO_G4_6.bed
coolpup.py ../DMSO/matrix_5kb.cool DMSO_G4_6.bed --nshifts 10 --unbalanced --coverage_norm --mindist 100000 --outname DMSO_G4_loop.txt
coolpup.py ../PDS/matrix_5kb.cool DMSO_G4_6.bed --nshifts 10 --unbalanced --coverage_norm --mindist 100000 --outname PDS_G4_loop.txt
############################


######################Fig4K
######################bash
bedtools makewindows -g mm10.genome -w 5000 > windows.bed
bedtools intersect -a windows.bed -b chr_ASC_rep123.bed -wo |cut -f 1-3|sort|uniq > G4_windows.bed
##########################
######################R
library(HiCcompare)
library(ggplot2)
library(QDNAseq)
library(BiocParallel)
mm10_blacklist<-read.table("/lustre/home/fengyang/Ref/blacklist/mm10.blacklist.bed")
colnames(mm10_blacklist)<-c("chromosome","start","end")
UT242_DMSO <- read.table("DMSO/matrix_5kb.txt", header = FALSE)
UT242_PDS <- read.table("PDS/matrix_5kb.txt", header = FALSE)

cnv <- get_CNV(path2bam = 'bam', out.file = 'outfile',
               bin.size = 5, genome = 'mm10', CNV.level = 2)
save(cnv,file='cnv.Rdata')
load("cnv.Rdata")
# combine cnv excluded regions with blacklist regions
exclude <- rbind(cnv, mm10_blacklist)

chr <- paste0('chr', c(1:19, 'X')) # First chromosome, or use c(1:22, 'X') for all chromosomes

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
save(hic.list,file='hic.list.Rdata')
load("hic.list.Rdata")

hic.list <- hic_compare(hic.list, A.min = 15, adjust.dist = TRUE, p.method = 'fdr', Plot = TRUE, parallel=TRUE)
hic.list <- do.call(rbind, hic.list)

G4_window_data<-read.table("G4_windows.bed",sep="\t")
result=data.frame()
for (i in 1:length(G4_window_data$V1)) {
        chr=G4_window_data$V1[i]
        start=G4_window_data$V2[i]
        end=G4_window_data$V3[i]
        temp=hic.list[((hic.list$chr1==chr)&(hic.list$start1==start)&(hic.list$end1==end)&(hic.list$p.adj<0.1))|((hic.list$chr2==chr)&(hic.list$start2==start)&(hic.list$end2==end)&(hic.list$p.adj<0.1))]
        if(length(temp$chr1)<1){
                next
        }
        else{
                result=rbind(result,temp)
        }
}

result$type="NoDiff"
result$type[(result$adj.M>0.5)&(result$p.adj<0.1)]="Up"
result$type[(result$adj.M< -0.5)&(result$p.adj<0.1)]="Down"
write.table(result[result$type!="NoDiff"],file="G4_window_DNSO_vs_PDS_diff.bed",row.names=FALSE,quote=FALSE,sep="\t")
valcano <- ggplot(data=result, aes(x=adj.M, y=-log10(p.adj),fill=type,color=type)) +
    geom_point(alpha=0.8, size=1) +
    theme_bw(base_size=15) +
    theme(
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()
    ) +
    ggtitle("DESeq2 Valcano") +
    # scale_color_manual(name="", values=c("#F8766D", "#619CFF", "#00BA38"), limits=c("Up", "Down", "NoDiff")) +
    geom_vline(xintercept=c(-0.5, 0.5), lty=2, col="gray", lwd=0.5) +
    geom_hline(yintercept=-log10(0.1), lty=1, col="gray", lwd=0.5)
pdf("resultvalcano15.pdf")
valcano
dev.off()
#################################

######################Fig4L,bash
ngs.plot.r -G mm10 -R bed -C config.txt -O H3K4me3_signal_across_G4_loops -N 0.5
###############################

######################Fig4M,N
######################Bash
awk 'NR>1 {print $1"\t"$2"\t"$3"\n"$4"\t"$5"\t"$6}' G4_window_DNSO_vs_PDS_diff.bed |sort|uniq > G4_window_DNSO_vs_PDS_diff_anchor.bed
bedtools intersect -a mm10_promoter_2k_05k_pcg.bed -b G4_window_DNSO_vs_PDS_diff_anchor.bed -wo |sort|uniq > 5k_loop_anchor_promoter.bed
bedtools intersect -a G4_window_DNSO_vs_PDS_diff_anchor.bed -b ASC_H3K27ac.bed -wo |cut -f 1-3|sort|uniq > enhancer_5k_loop_anchor.bed
bedtools intersect -a G4_window_DNSO_vs_PDS_diff_anchor.bed -b mm10_promoter_2k_05k_pcg.bed -wo |cut -f 1-3|sort|uniq > promoter_5k_loop_anchor.bed

python3 loop_G4_type.py G4_window_DNSO_vs_PDS_diff.bed G4_5k_loop_anchor.bed 5k_loop_anchor_promoter.bed promoter_5k_loop_anchor.bed enhancer_5k_loop_anchor.bed G4_window_DNSO_vs_PDS_diff_loop.bed

grep Down G4_window_DNSO_vs_PDS_diff.bed|awk '{print $1"\t"$2"\t"$6}' > G4_window_DNSO_vs_PDS_down_interaction.bed
############################
######################R
library(ggplot2)
data<-read.table("G4_window_DNSO_vs_PDS_diff_loop.bed",sep="\t")
gene_Up_list=c()
gene_Down_list=c()
for(i in 1:length(data$V1)){
        if(data$V19[i]=="Up"){
                if(data$V28[i]!=""){
                        gene_1=strsplit(data$V28[i],",")[[1]]
                        for (j in gene_1) {
                                if(j %in% gene_Up_list){
                                        next
                                }
                                else{
                                        gene_Up_list<-append(gene_Up_list,j)
                                }
                        }
                }
                if(data$V29[i]!=""){
                        gene_1=strsplit(data$V29[i],",")[[1]]
                        for (j in gene_1) {
                                if(j %in% gene_Up_list){
                                        next
                                }
                                else{
                                        gene_Up_list<-append(gene_Up_list,j)
                                }
                        }
                }
        }
        if(data$V19[i]=="Down"){
                if(data$V28[i]!=""){
                        gene_1=strsplit(data$V28[i],",")[[1]]
                        for (j in gene_1) {
                                if(j %in% gene_Down_list){
                                        next
                                }
                                else{
                                        gene_Down_list<-append(gene_Down_list,j)
                                }
                        }
                }
                if(data$V29[i]!=""){
                        gene_1=strsplit(data$V29[i],",")[[1]]
                        for (j in gene_1) {
                                if(j %in% gene_Down_list){
                                        next
                                }
                                else{
                                        gene_Down_list<-append(gene_Down_list,j)
                                }
                        }
                }
        }
}
gene_mix_list<-intersect(gene_Down_list,gene_Up_list)
gene_Down_list_only<-setdiff(gene_Down_list,gene_mix_list)
gene_Up_list_only<-setdiff(gene_Up_list,gene_mix_list)

resdata<-read.csv("RNA/DEG/all_des_output.csv")
Allgene <- read.table("RNA/DEG/out.bed", header=TRUE, quote="\t", skip=1)
sampleNames<-c("DMSO_1","DMSO_2","PDS_1","PDS_2")
library(DGEobj.utils)
rownames(Allgene) <- Allgene$Geneid
names(Allgene)[7:10] <- sampleNames
allgeneexpression <- Allgene[,7:10]

geneLength <- Allgene$Length

tpmMatrix <- convertCounts(as.matrix(allgeneexpression),unit = "TPM",geneLength = geneLength,log = F,

              normalize = "none",prior.count = NULL)
tpmMatrix1 <- as.data.frame(tpmMatrix)

for(i in 1:length(resdata$Row.names)){
    name=resdata$Row.names[i]
    resdata$DMSO_1_RNA[i]=tpmMatrix1$DMSO_1[rownames(tpmMatrix1)==name]
    resdata$DMSO_2_RNA[i]=tpmMatrix1$DMSO_2[rownames(tpmMatrix1)==name]
    resdata$PDS_1_RNA[i]=tpmMatrix1$PDS_1[rownames(tpmMatrix1)==name]
    resdata$PDS_2_RNA[i]=tpmMatrix1$PDS_2[rownames(tpmMatrix1)==name]
}

library(DESeq2)
library(pheatmap)
library(tidyr)

resdata$change <- as.factor(
    ifelse(
        resdata$padj<0.05 & abs(resdata$log2FoldChange)>0.5,
        ifelse(resdata$log2FoldChange>0.5, "Up", "Down"),
        "NoDiff"
    )
)
resdata<-resdata[complete.cases(resdata),]

data_frame=resdata[(resdata$genename %in% gene_Down_list_only)&(resdata$change!="NoDiff"),7:10]
p<-pheatmap(data_frame, scale="row",cluster_rows=T,cluster_cols=F,show_rownames=F,treeheight_col = 0, treeheight_row = 0,clustering_method = "ward")
pdf("PDS_interaction_down_diff_promoter.pdf")
p
dev.off()


number=0
activate_genes<-resdata$genename[(resdata$genename %in% gene_Down_list_only)&(resdata$change=="Down")]
for(i in 1:length(data$V1)){
        if(data$V19[i]=="Down"){
                log=0
                if(data$V28[i]!=""){
                        gene_1=strsplit(data$V28[i],",")[[1]]
                        for (j in gene_1) {
                                if(j %in% activate_genes){
                                        number = number+1
                                        log=1
                                        break
                                }
                        }
                }
                if((log==0)&(data$V29[i]!="")){
                        gene_1=strsplit(data$V29[i],",")[[1]]
                        for (j in gene_1) {
                                if(j %in% activate_genes){
                                        number = number+1
                                        log=1
                                        break
                                }
                        }
                }
        }
}
print(number)
sum(resdata$change[(resdata$genename %in% gene_Down_list_only)&(resdata$change!="NoDiff")]=="Up")
# 32
sum(resdata$change[(resdata$genename %in% gene_Down_list_only)&(resdata$change!="NoDiff")]=="Down")
# 507
write.csv(resdata[(resdata$genename %in% gene_Down_list_only),],file="down_G4_interaction_gene.csv")
sum(resdata$change[(resdata$genename %in% gene_Up_list_only)&(resdata$change!="NoDiff")]=="Up")
# 2
sum(resdata$change[(resdata$genename %in% gene_Up_list_only)&(resdata$change!="NoDiff")]=="Down")
# 52

library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(DOSE)
library(topGO)
library(clusterProfiler)
library(pathview)
library(ggplot2)

ensLookup <- resdata$genename[(resdata$change=="Up")&(resdata$genename %in% gene_Down_list_only)]
ego <- enrichGO(gene          = ensLookup,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                # minGSSize = 10, maxGSSize = 1000
                )
head(summary(ego))

ensLookup <- resdata$genename[(resdata$change=="Down")&(resdata$genename %in% gene_Down_list_only)]
ego <- enrichGO(gene          = ensLookup,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                # minGSSize = 10, maxGSSize = 1000
                )
head(summary(ego))
p<-dotplot(ego, showCategory=5,font.size=5,title="Enrichment GO")
pdf("PDS_loop_G4_RNA_down.pdf")
p
dev.off()
write.csv(summary(ego),"enrich-Down.csv",row.names =FALSE)
##################################