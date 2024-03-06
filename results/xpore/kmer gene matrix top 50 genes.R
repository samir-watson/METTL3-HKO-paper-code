library(tidyverse)
library(viridis)
library(pheatmap)


basedir <- ("C:/Users/au656873/OneDrive - Aarhus Universitet/Desktop/METTL3/xpore/all_KO")
setwd(basedir)
dir()

# Read the data from the CSV file
data <- read_csv("all_KO_filtered_majority_direction_kmer_diffmod_with_gene_names_sorted_filtered.table")
head(data)

# Count the occurrences of each k-mer
kmer_counts <- data %>%
  group_by(kmer) %>%
  summarise(total_count = n()) %>%
  arrange(desc(total_count))

# Identify the top 5 kmers
top_kmers <- kmer_counts %>%
  slice_head(n = 5) %>%
  pull(kmer)

# Add a new column to the original data where non-top kmers are labeled as "others"
data <- data %>%
  mutate(kmer = ifelse(kmer %in% top_kmers, kmer, "others"))

# Now filter the top 25 genes
top_genes <- data %>%
  arrange(desc(diff_mod_rate_WT_vs_xKO)) %>%
  slice_head(n = 50)

# Get the order of gene names
gene_order <- top_genes$gene_name

# Count the occurrences of each k-mer for these genes
kmer_counts <- top_genes %>%
  group_by(gene_name, kmer) %>%
  summarise(count = n()) %>%
  ungroup()

# Convert gene_name to a factor and specify the levels
kmer_counts$gene_name <- factor(kmer_counts$gene_name, levels = unique(gene_order))

# Create a matrix-like plot # change x and y axis to change plot direction
p <- ggplot(kmer_counts, aes(y = kmer, x = gene_name, fill = count)) +
  geom_tile(color = "grey") +
  geom_text(aes(label = count), color = "black", size = 4) +
  scale_fill_gradient(low = "powderblue", high = "steelblue") +
  theme_set(theme_bw() + theme(panel.grid = element_blank())) +
  labs(y = "Gene", x = "K-mer", fill = "Count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

# Save the plot as a PNG file with 300 dpi
ggsave("kmer_counts_per_gene_xpore1.png", dpi = 300, bg = "white", width = 11, height =3)

