#!/bin/sh
#SBATCH -J counting
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=/home/user_25/scripts/log/6.read_counts_CD0.out
#SBATCH --error=/home/user_25/scripts/log/6.read_counts_CD0.err
export PATH=/data/mo/envsn:$PATH

sample=('CD0_1' 'CD0_2' 'CD0_3')

echo start CD0 read_counts `date`
for i in $(seq 1 ${#sample[@]});do
echo start ${sample[$i-1]} `date`
cd /data/mo/PartI.RNA-seq/
mkdir -p 5.read_counts/${sample[$i-1]}
cd 5.read_counts/${sample[$i-1]}

featureCounts \
-T 8 \
-s 0 \
-p -t CDS \
-g gene_id \
-a /data/TA_QUIZ_RNA_regulation/data/ATH/GTF/Arabidopsis_thaliana.TAIR10.34.gtf \
-o /data/mo/PartI.RNA-seq/5.read_counts/${sample[$i-1]}/${sample[$i-1]}.featurecounts.txt \
/data/mo/PartI.RNA-seq/4.mapping/${sample[$i-1]}/${sample[$i-1]}Aligned.sortedByCoord.out.bam

featureCounts \
-T 8 \
-s 0 \
-p -t exon \
-g gene_id \
-a /data/TA_QUIZ_RNA_regulation/data/ATH/GTF/Arabidopsis_thaliana.TAIR10.34.gtf \
-o /data/mo/PartI.RNA-seq/5.read_counts/${sample[$i-1]}/${sample[$i-1]}.featurecounts.all.txt \
/data/mo/PartI.RNA-seq/4.mapping/${sample[$i-1]}/${sample[$i-1]}Aligned.sortedByCoord.out.bam

echo read_counts success

cd /data/mo/PartI.RNA-seq
mkdir -p 5.read_counts/result/${sample[$i-1]}
cd 5.read_counts/result/${sample[$i-1]}

echo -e "gene_id	${sample[$i-1]}" >/data/mo/PartI.RNA-seq/5.read_counts/result/${sample[$i-1]}.txt
cat /data/mo/PartI.RNA-seq/5.read_counts/${sample[$i-1]}/${sample[$i-1]}.featurecounts.txt| grep -v '#' | grep -v 'Geneid' | cut -f 1,7 >> /data/mo/PartI.RNA-seq/5.read_counts/result/${sample[$i-1]}.txt

echo -e "gene_id	${sample[$i-1]}" >/data/mo/PartI.RNA-seq/5.read_counts/result/${sample[$i-1]}.all.txt
cat /data/mo/PartI.RNA-seq/5.read_counts/${sample[$i-1]}/${sample[$i-1]}.featurecounts.all.txt| grep -v '#' | grep -v 'Geneid' | cut -f 1,7 >> /data/mo/PartI.RNA-seq/5.read_counts/result/${sample[$i-1]}.all.txt

echo finish ${sample[$i-1]} `date`
done 

echo success CD0 `date`
