#must first run "nanopolish for xpore and m6anet" script

#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem 24G
#SBATCH --account nanopore_drs
#SBATCH -c 20
#SBATCH -t 10:00:00

#must do this for each replicate for WT and KO

xpore dataprep \
--eventalign /home/samirwatson/faststorage/METTL3/DRS_data/m6anet/preprocess/A56_1.m6anet.eventalign.txt \
--gtf_or_gff /home/samirwatson/faststorage/NAT10/gencode.v33.annotation.gtf \
--transcript_fasta /home/samirwatson/faststorage/NAT10/gencode.v33.transcripts_shorthead.fa \
--out_dir /home/samirwatson/faststorage/METTL3/DRS_data/xpore/preprocess/A56.1/ \
--n_processes 20
