bwa mem -5SP -T0 -t16 GRCm38.p6.genome.fa MicroC_R1.fastq.gz MicroC_R1_R2.fastq.gz -o aligned.sam
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path mm10.genome aligned.sam >  parsed.pairsam
pairtools sort --nproc 16 --tmpdir=./  parsed.pairsam > sorted.pairsam
pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats stats.txt --output dedup.pairsam sorted.pairsam
pairtools split --nproc-in 8 --nproc-out 8 --output-pairs mapped.pairs --output-sam unsorted.bam dedup.pairsam
samtools sort -@16 -T temp.bam -o mapped.PT.bam unsorted.bam
samtools index mapped.PT.bam
python3 get_qc.py -p stats.txt
java -Xmx48000m  -Djava.awt.headless=true -jar juicertools.jar pre --threads 16 mapped.pairs contact_map.hic mm10.genome
bgzip mapped.pairs
pairix mapped.pairs.gz
cooler cload pairix -p 16 mm10.genome:5000 mapped.pairs.gz matrix_5kb.cool
