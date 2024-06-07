library(dplyr)
library(readr)
library ('tidyverse')

# Define file paths
input_file <- "/home/samirwatson/nanocompore_results.tsv"
output_file <- "/home/samirwatson/3rd_analysis/filtered_data/new_sylamer_like.in.paper/TT_intensity_2.full_sites_ext.fa"
ids_file <- "/home/samirwatson/3rd_analysis/filtered_data/new_sylamer_like.in.paper/TT_intensity_2.full_sites_ext.fa.ids.txt"

# Read the input file
data <- read_tsv(input_file)

# Filter and transform the data
filtered_data <- data %>%
  mutate(TT_intensity_pvalue_context_2 = as.numeric(TT_intensity_pvalue_context_2)) %>% # Convert KS_dwell_pvalue_context_2 to numeric
  replace_na(list(TT_intensity_pvalue_context_2 = 1)) %>% # Replace NA values with 1 (or any other value that makes sense in your context)
  filter(TT_intensity_pvalue_context_2 < 0.5) %>%
  mutate(TT_intensity_pvalue_context_2 = -log10(TT_intensity_pvalue_context_2)) %>%
  arrange(chr, genomicPos, desc(TT_intensity_pvalue_context_2)) %>%
  distinct(chr, genomicPos, .keep_all = TRUE) %>%
  arrange(desc(TT_intensity_pvalue_context_2))

# Write the output file
writeLines(paste0(">", filtered_data$chr, "_", filtered_data$genomicPos, "_", filtered_data$TT_intensity_pvalue_context_2, "\n", filtered_data$ref_kmer), output_file)

# Extract IDs and write to the ids file
ids <- gsub(">", "", readLines(output_file))
writeLines(ids, ids_file)

sylamer -fasta TT_intensity_2.full_sites_ext.fa -universe TT_intensity_2.full_sites_ext.fa.ids.txt -k 5 -grow 100 --over  -o sylamer_TT_intensity_2.full.out
sylamer -fasta TT_intensity_2.full_sites_ext.fa -universe TT_intensity_2.full_sites_ext.fa.ids.txt -k 5 -grow 100 --logfold  -o sylamer_TT_intensity_2.lfc.out

library ('tidyverse')

### BASEDIR=paste0("/home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/nanocompore_results/3rd_analysis/filtered_data/new_sylamer_like.in.paper")
#RESULTS=paste0("/home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/nanocompore_results/3rd_analysis/filtered_data/new_sylamer_like.in.paper/results")

BASEDIR=paste0("/home/samirwatson/3rd_analysis/filtered_data/new_sylamer_like.in.paper")
RESULTS=paste0("/home/samirwatson/3rd_analysis/filtered_data/new_sylamer_like.in.paper/results")

syl <- read_tsv(paste0(BASEDIR, "/sylamer_TT_intensity_2.full.out")) %>% reshape2::melt() %>% mutate(variable=as.numeric(as.character(variable)))
scores <- read_delim(paste0(BASEDIR, "/TT_intensity_2.full_sites_ext.fa.ids.txt"), delim="_", col_names=F)
score_thr <- (filter(scores, X3>-log10(0.01)) %>% nrow )

top_motifs <- filter(syl, value>10) %>% pull(upper) %>% unique
top_auc <- group_by(syl, upper) %>% summarise(AUC=sum(value)) %>% top_n(1, AUC) %>% pull(upper)
top_100_auc <- group_by(syl, upper) %>% summarise(AUC=sum(value)) %>% top_n(100, AUC) %>% pull(upper)

syl <- syl %>% mutate(label=case_when(upper%in%top_auc~upper, T~"Other")) %>% mutate(label=gsub("T", "U", label))
syl[syl$label=="Other", 'label'] <- NA

#ggplot(syl, aes(x=variable, y=value, group=upper, colour=label)) + geom_line()

pdf(paste0(RESULTS, "/sylamer_TT_intensity_2.full.pdf"), width=12)
filter(syl, upper%in%top_100_auc) %>% ggplot(aes(x=variable, y=value, group=upper, colour=label)) + geom_line(size=0.3) +geom_vline(xintercept=score_thr, size=0.5, linetype=2) + theme_classic(24) + xlab("Ranked sequences") + ylab("Hypergeometric p-value\n(-log10)") + labs(colour="Motif") + theme(axis.line = element_line(colour = 'black', size = 0.5))
dev.off()

syl_fc <- read_tsv(paste0(BASEDIR, "/sylamer_TT_intensity_2.lfc.out")) %>% reshape2::melt() %>% mutate(variable=as.numeric(as.character(variable)))
filter(syl, upper%in% top_auc, variable<=score_thr)
filter(syl_fc, upper%in% top_auc, variable<=score_thr) %>% mutate(FC=2^value)



########## script for top 10 AUC kmers ########

library('tidyverse')

# Define file paths
BASEDIR=paste0("/home/samirwatson/3rd_analysis/filtered_data/new_sylamer_like.in.paper")
RESULTS=paste0("/home/samirwatson/3rd_analysis/filtered_data/new_sylamer_like.in.paper/results")

# Read the data
syl <- read_tsv(paste0(BASEDIR, "/sylamer_TT_intensity_2.full.out")) %>% reshape2::melt() %>% mutate(variable=as.numeric(as.character(variable)))
scores <- read_delim(paste0(BASEDIR, "/TT_intensity_2.full_sites_ext.fa.ids.txt"), delim="_", col_names=F)
score_thr <- (filter(scores, X3>-log10(0.01)) %>% nrow )

# Get the top motifs and AUC
top_motifs <- filter(syl, value>10) %>% pull(upper) %>% unique
top_auc <- group_by(syl, upper) %>% summarise(AUC=sum(value)) %>% top_n(10, AUC) %>% pull(upper) # Changed to top 10
top_100_auc <- group_by(syl, upper) %>% summarise(AUC=sum(value)) %>% top_n(100, AUC) %>% pull(upper)

# Modify the data
syl <- syl %>% mutate(label=case_when(upper%in%top_auc~upper, T~"Other")) %>% mutate(label=gsub("T", "U", label))
syl[syl$label=="Other", 'label'] <- NA

# Reorder the labels based on the AUC
syl <- syl %>% mutate(label = reorder(label, value, sum))

# Plot the data
pdf(paste0(RESULTS, "/AUC_top10sylamer_TT_intensity_2.full.pdf"), width=12)
filter(syl, upper%in%top_100_auc) %>% ggplot(aes(x=variable, y=value, group=upper, colour=label)) + geom_line(size=0.3) +geom_vline(xintercept=score_thr, size=0.5, linetype=2) + theme_classic(24) + xlab("Ranked sequences") + ylab("Hypergeometric p-value\n(-log10)") + labs(colour="Motif") + theme(axis.line = element_line(colour = 'black', size = 0.5))
dev.off()




############# here is for top pvalue kmers (highest peak)##################

library('tidyverse')

# Define file paths
BASEDIR=paste0("/home/samirwatson/3rd_analysis/filtered_data/new_sylamer_like.in.paper")
RESULTS=paste0("/home/samirwatson/3rd_analysis/filtered_data/new_sylamer_like.in.paper/results")

# Read the data
syl <- read_tsv(paste0(BASEDIR, "/sylamer_GMM_logit_2.full.out")) %>% reshape2::melt() %>% mutate(variable=as.numeric(as.character(variable)))
scores <- read_delim(paste0(BASEDIR, "/GMM_logit_2.full_sites_ext.fa.ids.txt"), delim="_", col_names=F)
score_thr <- (filter(scores, X3>-log10(0.01)) %>% nrow )

# Get the top motifs and highest value
top_motifs <- filter(syl, value>10) %>% pull(upper) %>% unique
top_values <- group_by(syl, upper) %>% summarise(max_value=max(value)) %>% top_n(10, max_value) %>% pull(upper) # Changed to top 10
top_100_values <- group_by(syl, upper) %>% summarise(max_value=max(value)) %>% top_n(100, max_value) %>% pull(upper)

# Modify the data
syl <- syl %>% mutate(label=case_when(upper%in%top_values~upper, T~"Other")) %>% mutate(label=gsub("T", "U", label))
syl[syl$label=="Other", 'label'] <- NA

# Reorder the labels based on the maximum value
syl <- syl %>% mutate(label = reorder(label, value, max))

# Define the PDF file name
pdf_file_name <- paste0(RESULTS, "/top_10_peaks_pval_", gsub(".tsv", ".pdf", basename(paste0(BASEDIR, "/sylamer_GMM_logit_2.full.pdf"))))

# Plot the data
pdf(pdf_file_name, width=12)
filter(syl, upper%in%top_100_values) %>% ggplot(aes(x=variable, y=value, group=upper, colour=label)) + geom_line(size=0.3) +geom_vline(xintercept=score_thr, size=0.5, linetype=2) + theme_classic(24) + xlab("Ranked sequences") + ylab("Hypergeometric p-value\n(-log10)") + labs(colour="Motif") + theme(axis.line = element_line(colour = 'black', size = 0.5))
dev.off()



################### CODE TO SHOW ONLY DRACH MOTIFS in SYLAMER PLOT   ###############################################################

library('tidyverse')

# Define file paths
BASEDIR=paste0("/home/samirwatson/3rd_analysis/filtered_data/new_sylamer_like.in.paper")
RESULTS=paste0("/home/samirwatson/3rd_analysis/filtered_data/new_sylamer_like.in.paper/results")

# Read the data
syl <- read_tsv(paste0(BASEDIR, "/sylamer_MW_dwell_2.full.out")) %>% reshape2::melt() %>% mutate(variable=as.numeric(as.character(variable)))
scores <- read_delim(paste0(BASEDIR, "/MW_dwell_2.full_sites_ext.fa.ids.txt"), delim="_", col_names=F)
score_thr <- (filter(scores, X3>-log10(0.01)) %>% nrow )

# Define the DRACH motifs
D <- c("A", "G", "T")
R <- c("A", "G")
A <- "A"
C <- "C"
H <- c("A", "C", "T")

# Generate all combinations of DRACH
DRACH_motifs <- expand.grid(D=D, R=R, A=A, C=C, H=H, stringsAsFactors = FALSE) %>%
  unite("DRACH", D:H, sep = "", remove = TRUE)

# Convert the DRACH motifs to a vector
DRACH_motifs <- DRACH_motifs$DRACH

# Get the top motifs and highest value
top_motifs <- filter(syl, value>10) %>% pull(upper) %>% unique
top_values <- group_by(syl, upper) %>% summarise(max_value=max(value)) %>% top_n(1, max_value) %>% pull(upper) # Changed to top 15
top_100_values <- group_by(syl, upper) %>% summarise(max_value=max(value)) %>% top_n(100, max_value) %>% pull(upper)

# Modify the data
syl <- syl %>% mutate(label=case_when(upper%in%top_values~upper, T~"Other")) %>% mutate(label=gsub("T", "U", label))
syl[syl$label=="Other", 'label'] <- NA

# Highlight the DRACH motifs
syl <- syl %>% mutate(highlight = ifelse(upper %in% DRACH_motifs, upper, "Other"))

# Reorder the labels based on the maximum value
syl <- syl %>% mutate(label = reorder(label, value, max))

# Define the PDF file name
pdf_file_name <- paste0(RESULTS, "/DRACH_peaks_", gsub(".tsv", ".pdf", basename(paste0(BASEDIR, "/sylamer_MW_dwell_2.full.pdf"))))

# Define colors for each DRACH motif present in the data
unique_motifs <- unique(syl$highlight)
colors <- setNames(rainbow(length(unique_motifs)), unique_motifs)
colors["Other"] <- "grey"  # Set the color for "Other" to grey

# Plot the data
pdf(pdf_file_name, width=12)
filter(syl, upper%in%top_100_values) %>% ggplot(aes(x=variable, y=value, group=upper, colour=highlight)) + geom_line(size=0.3) +geom_vline(xintercept=score_thr, size=0.5, linetype=2) + theme_classic(24) + xlab("Ranked sequences") + ylab("Hypergeometric p-value\n(-log10)") + labs(colour="Motif") + theme(axis.line = element_line(colour = 'black', size = 0.5)) + scale_color_manual(values = colors)
dev.off()
