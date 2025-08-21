################Fig4B,bash
###########call loop
java -jar juicertools.jar hiccups --cpu -k KR -r 5000 -t 20 --ignore-sparsity contact_map.hic new_loop
cd new_loop
awk 'NR>2 {if(($3-$2)==5000) print $0}' merged_loops.bedpe > 5k_merged_loops.bedpe
awk '{print $1"\t"$2"\t"$6}' 5k_merged_loops.bedpe> ASC_loop.bed
ngs.plot.r -G mm10 -R bed -C config.txt -O G4_signal_across_loops -N 0.5
################

################Fig4C,D
######################Bash
awk '{print $1"\t"$2"\t"$3"\n"$4"\t"$5"\t"$6}' 5k_merged_loops.bedpe |sort|uniq> 5k_loop_anchor.bed
bedtools intersect -a 5k_loop_anchor.bed -b chr_ASC_rep123.bed -wo |cut -f 1-3|sort|uniq > G4_5k_loop_anchor.bed
grep "protein_coding" mm10_promoter_2k_05k.bed > mm10_promoter_2k_05k_pcg.bed
awk '{if(($2+$3)/2-1000<0)print $1"\t"0"\t"int(($2+$3)/2)+1000; else print $1"\t"int(($2+$3)/2)-1000"\t"int(($2+$3)/2)+1000}' GM_H3K27ac_1_2.bed| bedtools intersect -a - -b mm10_promoter_2k_05k_pcg.bed -v > GM_enhancer.bed
bedtools intersect -a mm10_promoter_2k_05k_pcg.bed -b 5k_loop_anchor.bed -wo |sort|uniq > 5k_loop_anchor_promoter.bed
bedtools intersect -a 5k_loop_anchor.bed -b GM_H3K27ac_1_2.bed -wo |cut -f 1-3|sort|uniq > enhancer_5k_loop_anchor.bed
bedtools intersect -a 5k_loop_anchor.bed -b mm10_promoter_2k_05k_pcg.bed -wo |cut -f 1-3|sort|uniq > promoter_5k_loop_anchor.bed
python3 loop_G4_h3k27me3.py 5k_merged_loops.bedpe G4_5k_loop_anchor.bed G4_5k_loop_anchor.bed 5k_loop_anchor_promoter.bed promoter_5k_loop_anchor.bed enhancer_5k_loop_anchor.bed 5k_loop.bed
##########################
######################R
library(ggplot2)
data<-read.table("5k_loop.bed",sep="\t")
data$type="None"
data$type[data$V25=="G4"&data$V26=="G4"]="Both"
data$type[(data$V25=="G4"&data$V26!="G4")|(data$V26=="G4"&data$V25!="G4")]="Either"
data$loop_type="None"
data$loop_type[(data$V29=="P")&(data$V30=="P")]="PP"
data$loop_type[((data$V29=="E"&data$V32=="P")|(data$V30=="P"&data$V31=="E"))&(data$loop_type!="PP")]="EP"
data$loop_type[(data$V31=="E"&data$V32=="E")&(data$loop_type!="PP")&(data$loop_type!="EP")]="EE"
write.csv(data,"5k_loop_classification.csv")
pdf("loop_type_number.pdf")
pie(c(sum(data$type=="Both"), sum(data$type=="Either"), sum(data$type=="None")), labels = c("Both","Either","None"))
dev.off()

shapiro.test(log2(data$V12[data$type=="None"]))
shapiro.test(log2(data$V12[data$type=="Either"]))
shapiro.test(log2(data$V12[data$type=="Both"]))
wilcox.test(log2(data$V12[data$type=="None"]),log2(data$V12[data$type=="Both"]))
wilcox.test(log2(data$V12[data$type=="None"]),log2(data$V12[data$type=="Either"]))
wilcox.test(log2(data$V12[data$type=="Either"]),log2(data$V12[data$type=="Both"]))
p<-ggplot(data, aes(x=type, y=log2(V12),fill=type)) +
  geom_boxplot()+
  ylim(0,12)+
  # scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_bw()+theme(panel.grid=element_blank())+
  labs(title="",x="loop type", y = "Interaction")
pdf("Interaction.pdf")
p
dev.off()
#########################

################Fig4E,F
######################R
RNA_data<-read.table("ASC_48h_1/02quantification/genes.fpkm_tracking",header=T)
gene_OneG4_list=c()
gene_TwoG4_list=c()
gene_NoneG4_list=c()
for(i in 1:length(data$V1)){
        if((data$V25[i]!="G4")&(data$V26[i]!="G4")){
                if(data$V33[i]!=""){
                        gene_1=strsplit(data$V33[i],",")[[1]]
                        for (j in gene_1) {
                                if(j %in% gene_NoneG4_list){
                                        next
                                }
                                else{
                                        gene_NoneG4_list<-append(gene_NoneG4_list,j)
                                }
                        }
                }
                if(data$V34[i]!=""){
                        gene_1=strsplit(data$V34[i],",")[[1]]
                        for (j in gene_1) {
                                if(j %in% gene_NoneG4_list){
                                        next
                                }
                                else{
                                        gene_NoneG4_list<-append(gene_NoneG4_list,j)
                                }
                        }
                }
        }
        if((data$V25[i]=="G4")&(data$V26[i]=="G4")){
                if(data$V33[i]!=""){
                        gene_1=strsplit(data$V33[i],",")[[1]]
                        for (j in gene_1) {
                                if(j %in% gene_TwoG4_list){
                                        next
                                }
                                else{
                                        gene_TwoG4_list<-append(gene_TwoG4_list,j)
                                }
                        }
                }
                if(data$V34[i]!=""){
                        gene_1=strsplit(data$V34[i],",")[[1]]
                        for (j in gene_1) {
                                if(j %in% gene_TwoG4_list){
                                        next
                                }
                                else{
                                        gene_TwoG4_list<-append(gene_TwoG4_list,j)
                                }
                        }
                }
        }
        if(((data$V25[i]=="G4")&(data$V26[i]!="G4"))|((data$V25[i]!="G4")&(data$V26[i]=="G4"))){
                if(data$V33[i]!=""){
                        gene_1=strsplit(data$V33[i],",")[[1]]
                        for (j in gene_1) {
                                if(j %in% gene_OneG4_list){
                                        next
                                }
                                else{
                                        gene_OneG4_list<-append(gene_OneG4_list,j)
                                }
                        }
                }
                if(data$V34[i]!=""){
                        gene_1=strsplit(data$V34[i],",")[[1]]
                        for (j in gene_1) {
                                if(j %in% gene_OneG4_list){
                                        next
                                }
                                else{
                                        gene_OneG4_list<-append(gene_OneG4_list,j)
                                }
                        }
                }
        }
}
RNA_data$G4type="0"
RNA_data$G4type[RNA_data$gene_short_name %in% gene_NoneG4_list]="None"
RNA_data$G4type[RNA_data$gene_short_name %in% gene_OneG4_list]="Either"
RNA_data$G4type[RNA_data$gene_short_name %in% gene_TwoG4_list]="Both"
data_frame=RNA_data[RNA_data$G4type!="0",]
pdf("promoter_type_number.pdf")
pie(c(sum(RNA_data$G4type=="Both"), sum(RNA_data$G4type=="Either"), sum(RNA_data$G4type=="None")), labels = c("Both","Either","None"))
dev.off()
p<-ggplot(data_frame, aes(x=G4type, y=log(FPKM+1),fill=G4type)) +
  geom_boxplot()+
  ylim(0,15)+
  # scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_bw()+theme(panel.grid=element_blank())+
  labs(title="",x="G4 type", y = "Gene expression")
pdf("RNA.pdf")
p
dev.off()
shapiro.test(log2(data_frame$FPKM[data_frame$G4type=="Either"]+1))
wilcox.test(log2(data_frame$FPKM[data_frame$G4type=="Either"]+1),log2(data_frame$FPKM[data_frame$G4type=="None"]+1))
wilcox.test(log2(data_frame$FPKM[data_frame$G4type=="Both"]+1),log2(data_frame$FPKM[data_frame$G4type=="None"]+1))
wilcox.test(log2(data_frame$FPKM[data_frame$G4type=="Both"]+1),log2(data_frame$FPKM[data_frame$G4type=="Either"]+1))
######################################