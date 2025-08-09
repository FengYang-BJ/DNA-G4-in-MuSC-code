#######Fig 3CDE
library(DESeq2)
sampleNames <- c("FISC_1", "FISC_2", "ASC_1","ASC_2")
data <- read.table("out.bed", header=TRUE, quote="\t", skip=1)
names(data)[7:10] <- sampleNames
countData <- as.matrix(data[7:10])
rownames(countData) <- data$Geneid
database <- data.frame(name=sampleNames, condition=c("FISC", "FISC","ASC", "ASC"))
rownames(database) <- sampleNames
dds <- DESeqDataSetFromMatrix(countData, colData=database, design= ~ condition)
dds <- DESeq(dds)

res <- lfcShrink(dds,coef="condition_FISC_vs_ASC", type="apeglm")
pdf("res_shrink.pdf")
plotMA(res, main = "Shrinkage by ashr", alpha=0.05, ylim=c(-4,4))
dev.off()
write.csv(res, "res_shrink_des_output.csv")

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
promoter<-read.table("mm10_promoter_2k_2k.bed")
resdata$genename="None"
for(i in 1:length(resdata$Row.names)){
        temp_gene=promoter$V5[promoter$V4==resdata$Row.names[i]]
        if(length(temp_gene)>0){
                resdata$genename[i]=temp_gene
        }
}
write.csv(resdata, "all_des_output.csv", row.names=FALSE)

resdata$log2FoldChange<- 0-resdata$log2FoldChange
library(ggplot2)
resdata$change <- as.factor(
        ifelse(
                resdata$padj<0.05 & abs(resdata$log2FoldChange)>0.5,
                ifelse(resdata$log2FoldChange>0.5, "Up", "Down"),
                "NoDiff"
        )
)
resdata<-resdata[complete.cases(resdata),]

resdata$change<-factor(resdata$change,levels=c("Up","NoDiff","Down"))
valcano <- ggplot(data=resdata, aes(x=log2FoldChange, y=-log10(padj), color=change)) +
        geom_point(alpha=0.8, size=1) +
        theme_bw(base_size=15) +
        theme(
                panel.grid.minor=element_blank(),
                panel.grid.major=element_blank(),
        ) +
        # scale_color_manual(name="", values=c("red", "green", "black"), limits=c("Up", "Down", "NoDiff")) +
        geom_vline(xintercept=c(-0.5, 0.5), lty=2, col="gray", lwd=0.5) +
        geom_hline(yintercept=-log10(0.05), lty=1, col="gray", lwd=0.5)+xlim(-10,10)
# pdf("valcano.pdf")
png("valcano.png")
valcano
dev.off()

resdata$change<-factor(resdata$change,levels=c("Up","NoDiff","Down"))
valcano <- ggplot(data=resdata, aes(x=log10(baseMean), y=log2FoldChange, color=change)) +
        geom_point(alpha=0.8, size=1) +
        theme_bw(base_size=15) +
        theme(
                panel.grid.minor=element_blank(),
                panel.grid.major=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                legend.position="none"
        )
pdf("MA.pdf")
# png("MA.png")
valcano
dev.off()

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

library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(DOSE)
library(topGO)
library(clusterProfiler)
library(ggplot2)
ensLookup <- resdata$genename[resdata$change=="Down"]
table <- resdata[resdata$change=="Down",]
write.table(table,"Down_genes.bed",quote=FALSE,row.names=FALSE,sep="\t")
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
plot<-dotplot(ego, showCategory=5,font.size=10,title="Enrichment GO")
pdf("ASC_RNA_down.pdf")
plot
dev.off()
write.csv(summary(ego),"enrich-down.csv",row.names =FALSE)

ensLookup <- resdata$genename[resdata$change=="Up"]
table <- resdata[resdata$change=="Up",]
write.table(table,"Up_genes.bed",quote=FALSE,row.names=FALSE,sep="\t")
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
plot<-dotplot(ego, showCategory=5,font.size=10,title="Enrichment GO")
pdf("ASC_RNA_up.pdf")
plot
dev.off()
write.csv(summary(ego),"enrich-up.csv",row.names =FALSE)


Allgene <- read.table("out.bed", header=TRUE, quote="\t", skip=1)
library(DGEobj.utils)
rownames(Allgene) <- Allgene$Geneid
names(Allgene)[7:10] <- sampleNames
allgeneexpression <- Allgene[,7:10]

geneLength <- Allgene$Length

tpmMatrix <- convertCounts(as.matrix(allgeneexpression),unit = "TPM",geneLength = geneLength,log = F,

              normalize = "none",prior.count = NULL)
tpmMatrix1 <- as.data.frame(tpmMatrix)
resdata$FISC_1_RNA=0
resdata$FISC_2_RNA=0
resdata$ASC_1_RNA=0
resdata$ASC_2_RNA=0
for(i in 1:length(resdata$Row.names)){
        name=resdata$Row.names[i]
        resdata$FISC_1_RNA[i]=tpmMatrix1$FISC_1[rownames(tpmMatrix1)==name]
        resdata$FISC_2_RNA[i]=tpmMatrix1$FISC_2[rownames(tpmMatrix1)==name]
        resdata$ASC_1_RNA[i]=tpmMatrix1$ASC_1[rownames(tpmMatrix1)==name]
        resdata$ASC_2_RNA[i]=tpmMatrix1$ASC_2[rownames(tpmMatrix1)==name]
}

G4_FISC_promoter<-read.table("promoter_chr_FISC_rep123.bed",sep="\t")
G4_ASC_promoter=read.table("promoter_chr_ASC_rep123.bed",sep="\t")

ASC_unique_promoter<-setdiff(G4_ASC_promoter$V5,G4_FISC_promoter$V5)
ASC_unique_promoter_data=resdata[resdata$genename %in% ASC_unique_promoter,]
write.csv(ASC_unique_promoter_data,"ASC_unique_promoter_data.csv",row.names =FALSE)
sum(resdata$change[resdata$genename %in% ASC_unique_promoter]=="NoDiff")
sum(resdata$change[resdata$genename %in% ASC_unique_promoter]=="Down")
sum(resdata$change[resdata$genename %in% ASC_unique_promoter]=="Up")

resdata_new=resdata
resdata_new$G4="Non_ASC_unique_G4"
resdata_new$G4[resdata$genename %in% ASC_unique_promoter]="ASC_unique_G4"
write.csv(resdata_new,"ASC_FISC_RNA_data.csv",row.names =FALSE)


valcano <- ggplot(data=ASC_unique_promoter_data, aes(x=log10(baseMean), y=log2FoldChange, color=change))+
        geom_point(alpha=0.8, size=1) +
        theme_bw(base_size=15) +
        theme(
                panel.grid.minor=element_blank(),
                panel.grid.major=element_blank(),
        ) +
        ylim(-5,5)+
        ggtitle("") + xlab("")+ylab("")
########Fig 3C
pdf("G4_MA.pdf")
valcano
dev.off()
##################

G4_promoter_up=resdata$genename[(resdata$change=="Up")&(resdata$genename %in% ASC_unique_promoter)]
G4_promoter_down=resdata$genename[(resdata$change=="Down")&(resdata$genename %in% ASC_unique_promoter)]

write.table(G4_ASC_promoter[G4_ASC_promoter$V5 %in% G4_promoter_down,],"G4_promoter_down.bed",sep="\t",row.names=F,quote=FALSE,col.names=F)
write.table(G4_ASC_promoter[G4_ASC_promoter$V5 %in% G4_promoter_up,],"G4_promoter_up.bed",sep="\t",row.names=F,quote=FALSE,col.names=F)

G4_promoter_up_df=resdata[(resdata$change=="Up")&(resdata$genename %in% ASC_unique_promoter),]
G4_promoter_down_df=resdata[(resdata$change=="Down")&(resdata$genename %in% ASC_unique_promoter),]
write.table(G4_promoter_up_df,"G4_promoter_up_df.bed",sep="\t",row.names=F,quote=FALSE,)
write.table(G4_promoter_down_df,"G4_promoter_down_df.bed",sep="\t",row.names=F,quote=FALSE)

ego <- enrichGO(gene          = G4_promoter_up,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                # minGSSize = 10, maxGSSize = 1000
                )
head(summary(ego))
# plot<-dotplot(ego, showCategory=10,font.size=10,title="Enrichment GO")
plot<-dotplot(ego, showCategory=5,font.size=10,title="Enrichment GO")
########Fig 3D
pdf("ASC_G4_RNA_up.pdf")
plot
dev.off()
write.csv(summary(ego),"enrich-G4_up.csv",row.names =FALSE)

ego <- enrichGO(gene          = G4_promoter_down,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                # minGSSize = 10, maxGSSize = 1000
                )
head(summary(ego))
# plot<-dotplot(ego, showCategory=10,font.size=10,title="Enrichment GO")
plot<-dotplot(ego, showCategory=5,font.size=10,title="Enrichment GO")
########Fig 3E
pdf("ASC_G4_RNA_down.pdf")
plot
dev.off()
write.csv(summary(ego),"enrich-G4_down.csv",row.names =FALSE)
