library(ggplot2)
library(Guitar)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

basedir <- ("/faststorage/home/samirwatson/METTL3/DRS_data/multi_guitarplots/a56_final_final_plots/combined_plots/")
setwd(basedir)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

stBedFiles <- list(file.path(basedir, "sig_sites_KS_intensity_pvalue_context_2_thr_0.01.bed"),
                   file.path(basedir, "A56_filtered_majority_direction_kmer_diffmod_with_gene_names_sorted_filtered.table_new_transcript_positions.bed"),
                   file.path(basedir, "final.wt.data.site_proba.csv_new_transcript_positions.filter.prob.0.9.bed"))

print(stBedFiles)

g <- GuitarPlot(txTxdb = txdb,
           headOrtail = FALSE,
           stBedFiles = stBedFiles,
           enableCI = FALSE,
           pltTxType = c("mrna"),
           stGroupName = c("nanocompore KS intensity HKO-1","xpore DRACH HKO-1","m6anet WT"),
           #miscOutFilePrefix = "nanocompore_combined-gmm",
)
# Customize the plot
g <- g + ggtitle(label = "method comparison Distribution on mRNA") +
  theme_bw() +  # Use a white background
  theme(text = element_text(size = 20),  # Change the font size
        legend.title = element_text(size = 13),  # Change the legend title size
        legend.text = element_text(size = 12))  # Change the legend text size


# Open a PNG device
png(filename = "metagene method comparison all KS intensity.png", width = 2000, height = 1600, res = 300)

# Plot
print(g)

# Close the device
dev.off()

# Open a PDF device
pdf(file = "metagene method comparison all KS intensity.pdf", width = 9, height = 7)

# Plot
print(g)

# Close the device
dev.off()

gc()# to clear memory use after R runs! otherwise it uses more and more RAM each time you run it

################# for nanocompore comparison of different statitical tests used

library(Guitar)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)
library(cowplot)
library(lemon)
library(ggrepel)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

basedir <- ("/faststorage/home/samirwatson/METTL3/DRS_data/multi_guitarplots/a56_final_final_plots/combined_plots/")
setwd(basedir)
#dir()


# Get list of all .bed files in the current directory
files <- list.files(pattern = "\\pvalue_context_2_thr_0.01.bed")
#files <- list.files(pattern = "^all_KO_filtered_majority_direction_kmer_diffmod")
#files <- list.files(pattern = "csv_transcript_positions.bed")
print(files)

diffSites <- list()

# Loop over each file
for (file in files) {
  # Process the file and create a plot
  ncmp_results_sel <- pipe(paste0("tail -n +2 ", file)) %>% read_tsv(, col_names=c("chr", "start", "end", "name", "score", "strand"))
  # Use all differential sites
  sites05 <- ncmp_results_sel %>% with(., GRanges(seqnames=chr, ranges=IRanges(start=start, end=end), strand=strand))
  sites001 <- filter(ncmp_results_sel, score>2) %>% with(., GRanges(seqnames=chr, ranges=IRanges(start=start, end=end), strand=strand))
  
  # Add the sites to the list
  diffSites[[file]] <- sites001
}

# Generate guitar plot with different groups
p <- GuitarPlot(stGRangeLists = diffSites, 
                txTxdb=txdb, pltTxType=c("mrna"), 
                headOrtail=FALSE, 
                enableCI=FALSE, 
                stGroupName = names(diffSites))
