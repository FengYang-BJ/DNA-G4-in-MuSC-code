#######trim reads
fastp -g --detect_adapter_for_pe -i FISC_R1.fastq.gz -I FISC_R2.fastq.gz -o FISC_1.paired.fastq.gz -O  FISC_2.paired.fastq.gz --thread 8
##################


#######For mouse, using hg38 for human
filebase=$1

R1=$filebase"_R1.paired.fastq.gz"
R2=$filebase"_R2.paired.fastq.gz"

#######alignment
hisat2 --dta-cufflinks -p 20 -x mm10/hisat2_v25m/hisat2_index --downstream-transcriptome-assembly -1 $R1 -2 $R2 -S accepted_hits.sam 1>$filebase".hisat2.log" 2>$filebase".hisat2.err" && touch 01align.finished
01sam2bam.finished: 01align.finished
samtools sort -@ 4 -o accepted_hits.bam accepted_hits.sam && touch 01sam2bam.finished
01linkbam.finished: 01sam2bam.finished
ln -s accepted_hits.bam $filebase".bam" && samtools flagstat $filebase".bam" && touch 01linkbam.finished
rm accepted_hits.sam
########quantification
cufflinks -o 02quantification -p 20 -G mm10.annotation.gtf $filebase".bam"  >"02quantification/"$filebase".cufflinks.log" 2>"02quantification/"$filebase".cufflinks.err" && touch 02quantification.finished
########################################