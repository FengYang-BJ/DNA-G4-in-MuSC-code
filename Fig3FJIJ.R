library(DESeq2)

sampleNames <- c("DMSO_1", "DMSO_2","PDS_1","PDS_2")
data <- read.table("out.bed", header=TRUE, quote="\t", skip=1)


names(data)[7:10] <- sampleNames
countData <- as.matrix(data[7:10])
rownames(countData) <- data$Geneid
database <- data.frame(name=sampleNames, condition=c("DMSO", "DMSO","PDS","PDS"))
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
theme_bw(base_size=15) +
theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
)
pdf("PCA.pdf")
pca
dev.off()

res <- lfcShrink(dds,coef="condition_PDS_vs_DMSO", type="apeglm")
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

resdata$change <- as.factor(
    ifelse(
        resdata$padj<0.05 & abs(resdata$log2FoldChange)>0.5,
        ifelse(resdata$log2FoldChange>0.5, "Up", "Down"),
        "NoDiff"
    )
)
resdata<-resdata[complete.cases(resdata),]

valcano <- ggplot(data=resdata, aes(x=log2FoldChange, y=-log10(padj), color=change)) +
    geom_point(alpha=0.8, size=1) +
    theme_bw(base_size=15) +
    theme(
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()
    ) +
    ggtitle("DESeq2 Valcano") +
    scale_color_manual(name="", values=c("#F8766D", "#619CFF", "#00BA38"), limits=c("Up", "Down", "NoDiff")) +
    geom_vline(xintercept=c(-0.5, 0.5), lty=2, col="gray", lwd=0.5) +
    geom_hline(yintercept=-log10(0.05), lty=1, col="gray", lwd=0.5)+xlim(-10,10)
pdf("valcano.pdf")
valcano
dev.off()

sum(resdata$change=="Up")
sum(resdata$change=="Down")


library(tidyr)
Allgene <- read.table("out.bed", header=TRUE, quote="\t", skip=1)
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

library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(DOSE)
library(topGO)
library(clusterProfiler)
library(ggplot2)
ensLookup <- resdata$genename[resdata$change=="Up"]
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
pdf("PDS_RNA_up.pdf")
dotplot(ego, showCategory=10,font.size=10,title="Enrichment GO")+theme_bw()
dev.off()
write.csv(summary(ego),"enrich-Up.csv",row.names =FALSE)

ensLookup <- resdata$genename[resdata$change=="Down"]
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
p=dotplot(ego, showCategory=10,font.size=10,title="Enrichment GO")+theme_bw()
graph2ppt(x=p,file="PDS_RNA_down.ppt")
write.csv(summary(ego),"enrich-down.csv",row.names =FALSE)

promoter<-read.table("mm10_promoter_1k_1k.bed",sep="\t")
DEG_promoter<-resdata[resdata$change!="NoDiff",]
write.table(DEG_promoter,"DEG_promoter.bed",quote=FALSE,row.names=FALSE)


ASC_G4_DEG_promoter_up=read.table("G4_promoter_up.bed",sep="\t")
ASC_G4_DEG_promoter_down=read.table("G4_promoter_down.bed",sep="\t")

ASC_G4_DEG_promoter<-union(ASC_G4_DEG_promoter_up$V5,ASC_G4_DEG_promoter_down$V5)

G4_promoter_up=resdata$genename[(resdata$change=="Up")&(resdata$genename %in% ASC_G4_DEG_promoter)]
G4_promoter_down=resdata$genename[(resdata$change=="Down")&(resdata$genename %in% ASC_G4_DEG_promoter)]
G4_promoter_deg<-union(G4_promoter_up,G4_promoter_down)
G4_promoter_deg_df<-resdata[resdata$genename %in% G4_promoter_deg,]

ASC48_G4_promoter_data=resdata[resdata$genename %in% ASC_G4_DEG_promoter_up$V5,]
sum(ASC48_G4_promoter_data$change=="Up")
#1
sum(ASC48_G4_promoter_data$change=="Down")
#310

ASC48_G4_up_PDS_down_data=promoter[promoter$V4 %in% ASC48_G4_promoter_data$Row.names[ASC48_G4_promoter_data$change=="Down"],]
write.table(ASC48_G4_up_PDS_down_data,"ASC_G4_up_PDS_down_promoter_df.bed",sep="\t",row.names =FALSE,quote=FALSE,col=FALSE)

valcano <- ggplot(data=ASC48_G4_promoter_data, aes(x=log2FoldChange, y=-log10(padj), color=change)) +
    geom_point(alpha=0.8, size=1) +
    theme_bw(base_size=15) +
    theme(
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()
    ) +
    ggtitle("DESeq2 Valcano") +
    scale_color_manual(name="", values=c("#F8766D", "#619CFF", "#00BA38"), limits=c("Up", "Down", "NoDiff")) +
    geom_vline(xintercept=c(-0.5, 0.5), lty=2, col="gray", lwd=0.5) +
    geom_hline(yintercept=-log10(0.05), lty=1, col="gray", lwd=0.5)+xlim(-2.5,2.5)
################Fig 3F
pdf("ASC_G4_up_G4_valcano.pdf")
valcano
dev.off()
################
write.table(ASC48_G4_promoter_data,"ASC_G4_DEG_promoter_up_data.bed",sep="\t",row.names =FALSE,quote=FALSE)
ASC48_G4_promoter_data_G4_down=ASC48_G4_promoter_data$genename[ASC48_G4_promoter_data$change=="Down"]
ego <- enrichGO(gene          = ASC48_G4_promoter_data_G4_down,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                # minGSSize = 10, maxGSSize = 1000
                )
head(summary(ego))
################Fig 3G
pdf("PDS_G4_RNA_down_ASCG4_up.pdf")
dotplot(ego, showCategory=5,font.size=10,title="Enrichment GO")
dev.off()
################
write.csv(summary(ego),"PDS_G4_RNA_down_ASCG4_up.csv",row.names =FALSE)
ASC48_G4_promoter_PDS_down<-promoter[(promoter$V5 %in% ASC48_G4_promoter_data$genename[ASC48_G4_promoter_data$change=="Down"]),]
write.table(ASC48_G4_promoter_PDS_down,"ASC_G4_up_PDS_down_promoter.bed",sep="\t",row.names =FALSE,quote=FALSE)

ASC48_G4_promoter_data=resdata[resdata$genename %in% ASC_G4_DEG_promoter_down$V5,]
sum(ASC48_G4_promoter_data$change=="Up")
#69
sum(ASC48_G4_promoter_data$change=="Down")
#72

ASC48_G4_down_PDS_up_data=promoter[promoter$V4 %in% ASC48_G4_promoter_data$Row.names[ASC48_G4_promoter_data$change=="Up"],]
write.table(ASC48_G4_down_PDS_up_data,"ASC_G4_down_PDS_up_promoter_df.bed",sep="\t",row.names =FALSE,quote=FALSE,col=FALSE)

valcano <- ggplot(data=ASC48_G4_promoter_data, aes(x=log2FoldChange, y=-log10(padj), color=change)) +
    geom_point(alpha=0.8, size=1) +
    theme_bw(base_size=15) +
    theme(
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()
    ) +
    ggtitle("DESeq2 Valcano") +
    scale_color_manual(name="", values=c("#F8766D", "#619CFF", "#00BA38"), limits=c("Up", "Down", "NoDiff")) +
    geom_vline(xintercept=c(-0.5, 0.5), lty=2, col="gray", lwd=0.5) +
    geom_hline(yintercept=-log10(0.05), lty=1, col="gray", lwd=0.5)+xlim(-2.5,2.5)
################Fig 3I
pdf("ASC_G4_down_G4_valcano.pdf")
valcano
dev.off()
################
write.table(ASC48_G4_promoter_data,"ASC_G4_DEG_promoter_down_data.bed",sep="\t",row.names =FALSE,quote=FALSE)

ASC48_G4_promoter_data_G4_down=ASC48_G4_promoter_data$genename[ASC48_G4_promoter_data$change=="Up"]
ego <- enrichGO(gene          = ASC48_G4_promoter_data_G4_down,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                # minGSSize = 10, maxGSSize = 1000
                )
head(summary(ego))
################Fig 3J
pdf("PDS_G4_RNA_up_ASCG4_down.pdf")
dotplot(ego, showCategory=5,font.size=10,title="Enrichment GO")
dev.off()
################
write.csv(summary(ego),"PDS_G4_RNA_up_ASCG4_down.csv",row.names =FALSE)

