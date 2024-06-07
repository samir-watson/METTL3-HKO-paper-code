#must first run "map_align_nanoplishcomp"
#!/usr/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 24
#SBATCH -t 1-00:00:00
python
from nanocompore.SampComp import SampComp

s = SampComp (
    eventalign_fn_dict ={
        'WT':{
        'rep1' : '/home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_WT_1_eventalign_collapsed_reads_tsv/out_eventalign_collapse.tsv',
        'rep2':'/home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_WT_2_eventalign_collapsed_reads_tsv/out_eventalign_collapse.tsv',
        'rep3':'/home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_WT_3_eventalign_collapsed_reads_tsv/out_eventalign_collapse.tsv'
        },
        'KO':{
        'rep1' : '/home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_A56_2_eventalign_collapsed_reads_tsv/out_eventalign_collapse.tsv',
        'rep2':'/home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_A56_3_eventalign_collapsed_reads_tsv/out_eventalign_collapse.tsv',
        'rep3':'/home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/DRS_METTL3_A56_1_eventalign_collapsed_reads_tsv/out_eventalign_collapse.tsv'
        }},
    outpath = "/home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/nanocompore_results/3rd_analysis",
    fasta_fn = "/home/samirwatson/faststorage/NAT10/gencode.v33.transcripts_shorthead.fa",
    bed_fn = "/home/samirwatson/faststorage/NAT10/gencode.v33.annotation.bed",
    allow_warnings = True,
    min_coverage = 30,
    nthreads = 24,
    logit = "True",
    comparison_methods=["GMM", "MW", "KS", "TT"],
    sequence_context=2,
    sequence_context_weights="harmonic",
    overwrite = True)
db = s()





from nanocompore.SampCompDB import SampCompDB, jhelp
db = SampCompDB (
  db_fn = "/home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/nanocompore_results/3rd_analysis/out_SampComp.db",
  fasta_fn = "/home/samirwatson/faststorage/NAT10/gencode.v33.transcripts_shorthead.fa",
  bed_fn = "/home/samirwatson/faststorage/NAT10/gencode.v33.annotation.bed",
  log_level = "warning"
)

db.save_all ("/home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/nanocompore_results/3rd_analysis/")
