
#Nanopolish script, mapping with minimap2 to transcriptome includes nanopolishcomp script for nanocompore
#important - output of this nanopolish is not compatible with xpore or m6anet for those run "nanopolish for xpore and m6anet after this"

#note for the genome file we used gencode and shortened the headings, only keeping the transcript ids with version duie to downstream incomaptibilities with the pipes in the the headers - gencode.v33.transcripts_shorthead.fa
########################
#A56_1
########################

#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 10
#SBATCH -t 0
nanopolish index -s /home/samirwatson/faststorage/METTL3/DRS_data/A56_1_HQcalls/sequencing_summary.txt -d /home/samirwatson/faststorage/METTL3/DRS_data/DRS_METTL3_KO/A56_1/20200907_1535_MN26644_FAO43307_5892db6b/fast5 /home/samirwatson/faststorage/METTL3/DRS_data/basecalled/basecalled_A56_1.fastq
minimap2 -ax map-ont -t 10 -L /home/samirwatson/faststorage/NAT10/gencode.v33.transcripts_shorthead.fa /home/samirwatson/faststorage/METTL3/DRS_data/basecalled/basecalled_A56_1.fastq | samtools view -bh -F 2324 -q 10 | samtools sort -O bam > /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_A56_1.aligned.bam
samtools index /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_A56_1.aligned.bam
nanopolish eventalign --reads /home/samirwatson/faststorage/METTL3/DRS_data/basecalled/basecalled_A56_1.fastq --bam /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_A56_1.aligned.bam --genome /home/samirwatson/faststorage/NAT10/gencode.v33.transcripts_shorthead.fa --samples --print-read-names --scale-events --samples > /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_A56_1_eventalign_reads.tsv
NanopolishComp Eventalign_collapse -t 10 -i /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_A56_1_eventalign_reads.tsv -o /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_A56_1_eventalign_collapsed_reads.tsv

########################
#A56_2
########################
#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 7.5G
#SBATCH -c 64
#SBATCH -t 0
nanopolish index -s /home/samirwatson/faststorage/METTL3/DRS_data/A56_2_HQcalls/sequencing_summary.txt -d /home/samirwatson/faststorage/METTL3/DRS_data/DRS_METTL3_KO/A56_2/20200903_1348_MN30066_FAO43269_b5ee73c7/fast5 /home/samirwatson/faststorage/METTL3/DRS_data/basecalled/basecalled_A56_2.fastq
minimap2 -ax map-ont -t 64 -L  /home/samirwatson/faststorage/NAT10/gencode.v33.transcripts_shorthead.fa /home/samirwatson/faststorage/METTL3/DRS_data/basecalled/basecalled_A56_2.fastq | samtools view -bh -F 2324 -q 10 | samtools sort -O bam > /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_A56_2.aligned.bam
samtools index /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_A56_2.aligned.bam
nanopolish eventalign -t 64 --reads /home/samirwatson/faststorage/METTL3/DRS_data/basecalled/basecalled_A56_2.fastq --bam /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_A56_2.aligned.bam --genome /home/samirwatson/faststorage/NAT10/gencode.v33.transcripts_shorthead.fa --samples --print-read-names --scale-events --samples > /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_A56_2_eventalign_reads.tsv
NanopolishComp Eventalign_collapse -t 64 -i /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_A56_2_eventalign_reads.tsv -o /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_A56_2_eventalign_collapsed_reads.tsv


########################
#A56_3
########################
#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 10G
#SBATCH -c 48
#SBATCH -t 0
nanopolish index -s /home/samirwatson/faststorage/METTL3/DRS_data/A56_3_HQcalls/sequencing_summary.txt -d /home/samirwatson/faststorage/METTL3/DRS_data/DRS_METTL3_KO/A56_3/20201021_2153_MN26644_FAO44326_cb50bba3/fast5 /home/samirwatson/faststorage/METTL3/DRS_data/basecalled/basecalled_A56_3.fastq
minimap2 -ax map-ont -t 48 -L  /home/samirwatson/faststorage/NAT10/gencode.v33.transcripts_shorthead.fa /home/samirwatson/faststorage/METTL3/DRS_data/basecalled/basecalled_A56_3.fastq | samtools view -@48 -bh -F 2324 -q 10 | samtools sort -@48 -O bam > /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_A56_3.aligned.bam
samtools index -@48 /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_A56_3.aligned.bam
nanopolish eventalign -t 48 --reads /home/samirwatson/faststorage/METTL3/DRS_data/basecalled/basecalled_A56_3.fastq --bam /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_A56_3.aligned.bam --genome /home/samirwatson/faststorage/NAT10/gencode.v33.transcripts_shorthead.fa --samples --print-read-names --scale-events --samples > /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_A56_3_eventalign_reads.tsv
NanopolishComp Eventalign_collapse -t 48 -i /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_A56_3_eventalign_reads.tsv -o /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_A56_3_eventalign_collapsed_reads.tsv


########################
#WT_1
########################

#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 10
#SBATCH -t 0
nanopolish index -s /home/samirwatson/faststorage/METTL3/DRS_data/WT_1_HQcalls/sequencing_summary.txt -d /home/samirwatson/faststorage/METTL3/DRS_data/DRS_METTL3_WT/WT1/20200818_1356_MN30066_FAM96325_0c489307/fast5_pass /home/samirwatson/faststorage/METTL3/DRS_data/basecalled/basecalled_WT_1.fastq
minimap2 -ax map-ont -t 10 -L  /home/samirwatson/faststorage/NAT10/gencode.v33.transcripts_shorthead.fa /home/samirwatson/faststorage/METTL3/DRS_data/basecalled/basecalled_WT_1.fastq | samtools view -bh -F 2324 -q 10 | samtools sort -O bam > /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_WT_1.aligned.bam
samtools index /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_WT_1.aligned.bam
nanopolish eventalign --reads /home/samirwatson/faststorage/METTL3/DRS_data/basecalled/basecalled_WT_1.fastq --bam /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_WT_1.aligned.bam --genome /home/samirwatson/faststorage/NAT10/gencode.v33.transcripts_shorthead.fa --samples --print-read-names --scale-events --samples > /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_WT_1_eventalign_reads.tsv
NanopolishComp Eventalign_collapse -t 10 -i /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_WT_1_eventalign_reads.tsv -o /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_WT_1_eventalign_collapsed_reads.tsv



########################
#WT_2
########################

#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 10
#SBATCH -t 0
nanopolish index -s /home/samirwatson/faststorage/METTL3/DRS_data/WT_2_HQcalls/sequencing_summary.txt -d /home/samirwatson/faststorage/METTL3/DRS_data/DRS_METTL3_WT/WT2/20200907_1543_MN30066_FAO43597_6297976a/fast5 /home/samirwatson/faststorage/METTL3/DRS_data/basecalled/basecalled_WT_2.fastq
minimap2 -ax map-ont -t 10 -L /home/samirwatson/faststorage/NAT10/gencode.v33.transcripts_shorthead.fa /home/samirwatson/faststorage/METTL3/DRS_data/basecalled/basecalled_WT_2.fastq | samtools view -bh -F 2324 -q 10 | samtools sort -O bam > /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_WT_2.aligned.bam
samtools index /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_WT_2.aligned.bam
nanopolish eventalign --reads /home/samirwatson/faststorage/METTL3/DRS_data/basecalled/basecalled_WT_2.fastq --bam /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_WT_2.aligned.bam --genome /home/samirwatson/faststorage/NAT10/gencode.v33.transcripts_shorthead.fa --samples --print-read-names --scale-events --samples > /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_WT_2_eventalign_reads.tsv
NanopolishComp Eventalign_collapse -t 10 -i /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_WT_2_eventalign_reads.tsv -o /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_WT_2_eventalign_collapsed_reads.tsv



########################
#WT_3
########################

#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 15G
#SBATCH -c 32
#SBATCH -t 6-06:00:00
nanopolish index -s /home/samirwatson/faststorage/METTL3/DRS_data/WT_3_HQcalls/sequencing_summary.txt -d /home/samirwatson/faststorage/METTL3/DRS_data/DRS_METTL3_WT/WT3/20201030_1817_MN30066_FAO83012_8a18447a/fast5 /home/samirwatson/faststorage/METTL3/DRS_data/basecalled/basecalled_WT_3.fastq
minimap2 -ax map-ont -t 32 -L /home/samirwatson/faststorage/NAT10/gencode.v33.transcripts_shorthead.fa /home/samirwatson/faststorage/METTL3/DRS_data/basecalled/basecalled_WT_3.fastq | samtools view -@ 32 -bh -F 2324 -q 10 | samtools sort -@ 32 -O bam > /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_WT_3.aligned.bam
samtools index -@ 32 /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_WT_3.aligned.bam
nanopolish eventalign -t 32 --reads /home/samirwatson/faststorage/METTL3/DRS_data/basecalled/basecalled_WT_3.fastq --bam /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_WT_3.aligned.bam --genome /home/samirwatson/faststorage/NAT10/gencode.v33.transcripts_shorthead.fa --samples --print-read-names --scale-events --samples > /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_WT_3_eventalign_reads.tsv
NanopolishComp Eventalign_collapse -t 32 -i /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_WT_3_eventalign_reads.tsv -o /home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_WT_3_eventalign_collapsed_reads.tsv
