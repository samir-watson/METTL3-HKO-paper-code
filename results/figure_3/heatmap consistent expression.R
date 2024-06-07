library(tidyverse)
library(biomaRt)
library(viridis)
library (pheatmap)

basedir <- ("C:/Users/au656873/OneDrive - Aarhus Universitet/Desktop/METTL3/xpore/A56")
setwd(basedir)
dir()

data <- read_csv("A56_filtered_majority_direction_kmer_diffmod.table")
head(data)

# Convert transcript IDs to character type
data$id <- as.character(data$id)

# Remove version numbers from Gencode IDs
data$id_no_version <- sub("\\.[0-9]*", "", data$id)

#listensembl
listEnsembl(mirror ="asia")
# Connect to the Ensembl BioMart

#mart <- useEnsembl(biomart="genes", mirror = "www")
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror ="asia")
#or can use
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# Get gene names corresponding to the transcript IDs
gene_names <- getBM(
  filters = "ensembl_transcript_id",
  attributes = c("ensembl_transcript_id", "hgnc_symbol"),
  values = data$id_no_version,
  mart = mart
)

# Join the original data frame with the gene names
final_df <- left_join(data, gene_names, by = c("id_no_version" = "ensembl_transcript_id"))
head(final_df)
# Reorder the columns to put the gene name as the second column
final_df <- final_df[, c("id", "id_no_version", "hgnc_symbol", setdiff(names(final_df), c("id", "id_no_version", "hgnc_symbol")))]
head(final_df)
# Rename the 'hgnc_symbol' column to 'gene_name'
names(final_df)[names(final_df) == "hgnc_symbol"] <- "gene_name"
head(final_df)

# Sort the data frame by 'diff_mod_rate_WT_vs_xKO' in descending order
final_df <- arrange(final_df, desc(diff_mod_rate_WT_vs_xKO.A56))

# Filter the rows where at least 4 of the columns starting with 'coverage' have values greater than 0. the relaxed set would be at least 2 columns have counts greater than 0 and strict would have all 6
final_df2 <- final_df %>%
  rowwise() %>%
  mutate(coverage_count = sum(c_across(starts_with("coverage")) > 0, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::filter(coverage_count >= 6) %>%
  dplyr::select(-coverage_count)

head(final_df2)
# Write the final data frame to a new CSV file
write_tsv(final_df2, "A56_filtered_majority_direction_kmer_diffmod_with_gene_names_sorted_filtered.strict.tsv")
write_csv(final_df2, "A56_filtered_majority_direction_kmer_diffmod_with_gene_names_sorted_filtered.strict.table")
#might have to remove all NAs


######### script to make a heatmap
library(tidyverse)
library(viridis)
library(pheatmap)


basedir <- ("C:/Users/au656873/OneDrive - Aarhus Universitet/Desktop/METTL3/xpore/A56")
setwd(basedir)
dir()

# Read the data from the CSV file
data <- read_csv("A56_filtered_majority_direction_kmer_diffmod_with_gene_names_sorted_filtered.strict.table")

# Specify the columns you're interested in
cols_to_plot <- c("gene_name", "mod_rate_WT-REP1", "mod_rate_WT-REP2", "mod_rate_WT-REP3", 
                  "mod_rate_xKO.A56-REP1", "mod_rate_xKO.A56-REP2", "mod_rate_xKO.A56-REP3")

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
                       "HKO_REP1" = "mod_rate_xKO.A56-REP1",
                       "HKO_REP2" = "mod_rate_xKO.A56-REP2",
                       "HKO_REP3" = "mod_rate_xKO.A56-REP3")

# Start the PNG device
png("A56 final heatmap.strict.no.border.png", res = 300, width = 3000, height = 3200)

# Generate the heatmap
pheatmap(data_to_plot, 
         color = viridis(100), # Use the viridis color palette
         scale = "row", # Scale the data by row for better visualization
         labels_row = gene_names,  # Add gene names
         fontsize = 5,
         border_color =NA)
#now with bigger x axis labels and with bigger cell heights
pheatmap(data_to_plot, 
         color = viridis(100), # Use the viridis color palette
         scale = "row", # Scale the data by row for better visualization
         labels_row = gene_names,  # Add gene names
         fontsize = 5,
         fontsize_col = 10,  # Increase the font size of the column names
         border_color = NA,
         cellheight = 4.5)

# Close the PNG device
dev.off()
