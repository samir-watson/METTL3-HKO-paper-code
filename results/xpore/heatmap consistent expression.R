library(tidyverse)
library(viridis)
library(pheatmap)


basedir <- ("C:/Users/au656873/OneDrive - Aarhus Universitet/Desktop/METTL3/xpore/a56")
setwd(basedir)
dir()

# Read the data from the CSV file
data <- read_csv("a56_filtered_majority_direction_kmer_diffmod_with_gene_names_sorted_filtered.strict.table")

# Specify the columns you're interested in
cols_to_plot <- c("gene_name", "mod_rate_WT-REP1", "mod_rate_WT-REP2", "mod_rate_WT-REP3", 
                  "mod_rate_xKO-REP1", "mod_rate_xKO-REP2", "mod_rate_xKO-REP3")

# Subset the data to include only the columns of interest
data_to_plot <- data[cols_to_plot]

# Extract the gene names
gene_names <- data_to_plot$gene_name

# Remove the gene names from the data to be plotted
data_to_plot <- data_to_plot[-1]

# Rename the columns
data_to_plot <- rename(data_to_plot, 
                       "WT_REP1" = "mod_rate_WT-REP1",
                       "WT_REP2" = "mod_rate_WT-REP2",
                       "WT_REP3" = "mod_rate_WT-REP3",
                       "HKO_1_REP1" = "mod_rate_xKO-REP1",
                       "HKO_1_REP2" = "mod_rate_xKO-REP2",
                       "HKO_1_REP3" = "mod_rate_xKO-REP3")

# Start the PNG device
png("final heatmap.strict.no.border.png", res = 300, width = 2000, height = 2000)

# Generate the heatmap
pheatmap(data_to_plot, 
         color = viridis(100), # Use the viridis color palette
         scale = "row", # Scale the data by row for better visualization
         labels_row = gene_names,  # Add gene names
         fontsize = 6,
         border_color =NA)

# Close the PNG device
dev.off()
