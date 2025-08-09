#trim reads
fastp -g --detect_adapter_for_pe -i G4_R1.fastq.gz -I G4_R2.fastq.gz -o G4_1.paired.fastq.gz -O  G4_2.paired.fastq.gz --thread 8

#######For mouse, using hg38 for human
filebase=$1
#indexpath="mm10_indexpath"

R1=$filebase"_1.paired.fastq.gz"
R2=$filebase"_2.paired.fastq.gz"

#alignment
bowtie2 -p 20 -x $indexpath -q -1 $R1 -2 $R2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 2> $filebase".mm10.all.bowtie.log" | samtools view -bS -F 4  - | samtools sort -o $filebase".mm10.all.bam" -
samtools index $filebase".mm10.all.bam"
samtools flagstat $filebase".mm10.all.bam" > $filebase".mm10.all.flagstat.txt"
samtools idxstats $filebase".mm10.all.bam" > $filebase".mm10.all.idxstats.txt"
java -jar picard.jar MarkDuplicates REMOVE_DUPLICATES=TRUE I=$filebase".mm10.all.bam" O=$filebase".mm10.dedup.bam" M=$filebase".mm10.all.MarkDuplicates.txt"

samtools index $filebase".mm10.dedup.bam"
bedtools intersect -v -abam $filebase".mm10.dedup.bam" -b mm10.blacklist.bed > $filebase".mm10.fltd.bam"

rm -f $filebase".mm10.dedup.bam"
rm -f $filebase".mm10.dedup.bam.bai"

samtools index $filebase".mm10.fltd.bam"
samtools flagstat $filebase".mm10.fltd.bam" > $filebase".mm10.fltd.flagstat.txt"
samtools idxstats $filebase".mm10.fltd.bam" > $filebase".mm10.fltd.idxstats.txt"

samtools view -bh -q 20 $filebase".mm10.fltd.bam" > $filebase".mm10.fltd.uniq.bam"

samtools index $filebase".mm10.fltd.uniq.bam"
samtools flagstat $filebase".mm10.fltd.uniq.bam" > $filebase".mm10.fltd.uniq.flagstat.txt"
samtools idxstats $filebase".mm10.fltd.uniq.bam" > $filebase".mm10.fltd.uniq.idxstats.txt"

java -jar picard.jar CollectInsertSizeMetrics I=$filebase".mm10.fltd.bam" O=$filebase".mm10.insert_sizes.txt" H=$filebase".mm10.insert_sizes.pdf" M=0.5 &

igvtools count --pairs $filebase".mm10.fltd.bam" $filebase".mm10.fltd.tdf" mm10.chrom.sizes &
wait

bamCoverage -p 2 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 -bs 5 -e --minMappingQuality 20 -b $filebase".mm10.fltd.bam" -o $filebase".mm10.pe.uniq.bw" &
bamCoverage -p 2 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 -bs 5 -e --minMappingQuality 20 --minFragmentLength 50 --maxFragmentLength 300 -b $filebase".mm10.fltd.bam" -o $filebase".mm10.pe.50-300.uniq.bw" &
bamCoverage -p 2 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 -bs 5 -b $filebase".mm10.fltd.bam" -o $filebase".mm10.bw" &

plotFingerprint -b $filebase".mm10.fltd.bam" --minMappingQuality 30 --skipZeros --region 19 --numberOfSamples 50000 --plotFile $filebase".mm10.fingerprint.png"

mkdir macs2
cd macs2
macs2 callpeak -t ../$filebase".mm10.fltd.uniq.bam" -n $filebase -g mm
macs2 callpeak -t ../$filebase".mm10.fltd.uniq.bam" -n $filebase"_p001" -g mm -p 0.001
bedtools getfasta -fi mm10.fa -bed $filebase"_p001_peaks.narrowPeak" -fo $filebase"_p001.fa"
meme-chip -meme-p 20 -oc $filebase"_p001_Denoval" -dreme-m 20 -meme-nmotifs 20 $filebase"_p001.fa" -meme-maxw 20 -meme-minw 6