#script for baseecalling DRS data using Guppy - change paths for each sample

#!/bin/bash
#SBATCH --partition gpu
#SBATCH -c 8
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu 20G
#SBATCH --gres=gpu:2

singularity exec --nv guppy4-0-11.simg nvidia-smi
singularity exec --nv guppy4-0-11.simg guppy_basecaller -i /home/samirwatson/faststorage/METTL3/DRS_data/DRS_METTL3_KO/A56_3/20201021_2153_MN26644_FAO44326_cb50bba3/fast5 \
-s /home/samirwatson/faststorage/METTL3/DRS_data/A56_3_HQcalls/ \
-x cuda:all:100% \
--flowcell FLO-MIN106 \
--kit SQK-RNA002 \
-r \
--calib_detect \
--enable_trimming true \
--trim_strategy rna \
--reverse_sequence true \

cat /home/samirwatson/faststorage/METTL3/DRS_data/A56_3_HQcalls/*.fastq > /home/samirwatson/faststorage/METTL3/DRS_data/basecalled/basecalled_A56_3.fastq
