#first run "nanopolish for xpore and m6anet script". run this script for each sample

#run m6anets dataprep - run for each sample
m6anet dataprep --eventalign /home/samirwatson/faststorage/METTL3/DRS_data/m6anet/preprocess/A56_1.m6anet.eventalign.txt \
--out_dir /home/samirwatson/faststorage/METTL3/DRS_data/m6anet/preprocess/a56.1 --n_processes 18

#run m6anet inference - run for each sample
m6anet inference --input_dir /home/samirwatson/faststorage/METTL3/DRS_data/m6anet/preprocess/a56.1 /home/samirwatson/faststorage/METTL3/DRS_data/m6anet/preprocess/a56.2 /home/samirwatson/faststorage/METTL3/DRS_data/m6anet/preprocess/a56.3 \
--out_dir /home/samirwatson/faststorage/METTL3/DRS_data/m6anet/inference/a56 \
--n_processes 12 \
--num_iterations 1000
