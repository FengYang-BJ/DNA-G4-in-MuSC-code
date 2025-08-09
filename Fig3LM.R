#######Fig 3L,R
#######get ATAC,histone signal across promoter region from bw files
################H3K4ME3
# python3 bw2signal.py ASC_G4_up_PDS_down_promoter_df.bed QSC_y_H3K4me3_1.sorted.rmdup.RPGCnorm.bw H3K4me3_G4_promoter_ASCup_PDSdown_FISC_signal.bed
# python3 bw2signal.py ASC_G4_up_PDS_down_promoter_df.bed ASC48h_y_H3K4me3_1.sorted.rmdup.RPGCnorm.bw H3K4me3_G4_promoter_ASCup_PDSdown_ASC_signal.bed

# python3 bw2signal.py ASC_G4_down_PDS_up_promoter_df.bed QSC_y_H3K4me3_1.sorted.rmdup.RPGCnorm.bw H3K4me3_G4_promoter_ASCdown_PDSup_FISC_signal.bed
# python3 bw2signal.py ASC_G4_down_PDS_up_promoter_df.bed ASC48h_y_H3K4me3_1.sorted.rmdup.RPGCnorm.bw H3K4me3_G4_promoter_ASCdown_PDSup_ASC_signal.bed

library(ggplot2)
library(export)
data_G4_up_FISC<-read.table("H3K4me3_G4_promoter_ASCup_PDSdown_FISC_signal.bed")
data_G4_up_ASC<-read.table("H3K4me3_G4_promoter_ASCup_PDSdown_ASC_signal.bed")
data_G4_down_FISC<-read.table("H3K4me3_G4_promoter_ASCdown_PDSup_FISC_signal.bed")
data_G4_down_ASC<-read.table("H3K4me3_G4_promoter_ASCdown_PDSup_ASC_signal.bed")
data_G4_up<-cbind(data_G4_up_FISC$V1,data_G4_up_ASC$V1)
data_G4_up<-as.data.frame(data_G4_up)
data_NonG4_up<-cbind(data_G4_down_FISC$V1,data_G4_down_ASC$V1)
data_NonG4_up<-as.data.frame(data_NonG4_up)
data_G4_up$type="G4_ASCup_PDSdown"
data_NonG4_up$type="G4_ASCdown_PDSup"
new_data<-rbind(data_G4_up,data_NonG4_up)

p<-ggplot(data=new_data, aes(x=type, y=log2((V2/V1)),fill=type)) +
  geom_violin(trim=FALSE,draw_quantiles = c(0.25, 0.75), linetype = "dashed")+
  geom_violin(fill="transparent",draw_quantiles = 0.5)+
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")+
  labs(title="",x="", y = "log2(ASC/FISC)")+
  ylim(-4,4)+
  theme_classic()
pdf("result_H3K4me3_bg_signal.pdf")
p
dev.off()
shapiro.test(log2((data_G4_up_ASC$V1)/(data_G4_up_FISC$V1)))
wilcox.test(log2((data_G4_up_ASC$V1)/(data_G4_up_FISC$V1)),log2((data_NonG4_up_ASC$V1)/(data_NonG4_up_FISC$V1)))

#######Fig 3M bash
#######H3K4ME3
ngs.plot.r -G mm10 -R bed -C config.txt -O ./H3K4ME3_signal -L 2000