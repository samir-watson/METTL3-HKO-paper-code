library(tidyverse)
library(viridis)
library(pheatmap)

basedir <- ("C:/Users/au656873/OneDrive - Aarhus Universitet/Desktop/METTL3/xpore/A56")
setwd(basedir)
dir()

# Read the data from the CSV file
data <- read_csv("A56_filtered_majority_direction_kmer_diffmod_with_gene_names_sorted_filtered.table")
head(data)

# Count the number of times each gene appears
gene_counts <- data %>%
  group_by(gene_name) %>%
  summarise(count = n())

# Set the theme to remove gridlines and increase text size
theme_set(theme_bw() + theme(panel.grid = element_blank(), text = element_text(size=20)))

# Create a histogram of these counts
p <- ggplot(gene_counts, aes(x=count)) +
  geom_histogram(binwidth=1, fill="steelblue", color="black") +
  labs(x="Number of Differentially \nModified DRACH Motifs", y="Number of Genes") 
p

# Save the plot as a PNG file with 300 dpi
ggsave("A56 modified DRACH per gene xpore mid dataset.png", plot = p, dpi = 300)

# Save the plot as a PDF file
ggsave("A56 modified DRACH per gene xpore mid dataset.pdf", plot = p, dpi = 300)
