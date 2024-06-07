# change the code for each sample
#must first run "nanopolish for xpore and m6anet" script

#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem 32G
#SBATCH --account nanopore_drs
#SBATCH -c 30
#SBATCH -t 4:00:00
xpore diffmod --config /home/samirwatson/faststorage/METTL3/DRS_data/xpore/diff_mod/A56/A56.yaml \
--n_processes 30

# now postpocessing
#can run  with minimal resources, can just use srun. can pipe these into the diffmod scripts

xpore postprocessing --diffmod_dir /home/samirwatson/faststorage/METTL3/DRS_data/xpore/diff_mod/A56
