library(tidyverse)
library(wesanderson)


basedir <- ("C:/Users/au656873/OneDrive - Aarhus Universitet/Desktop/METTL3/xpore/A56")
setwd(basedir)
dir()

D <- c("A", "G", "T")
R <- c("A", "G")
A <- "A"
C <- "C"
H <- c("A", "C", "T")

# Generate all combinations of DRACH
DRACH_motifs <- expand.grid(D=D, R=R, A=A, C=C, H=H, stringsAsFactors = FALSE) %>%
  unite("DRACH", D:H, sep = "", remove = TRUE)

# Read the TSV file
data <- read_csv("majority_direction_kmer_diffmod.table", col_names=TRUE)

# Filter the data
filtered_data <- data %>%
  filter(pval_WT_vs_xKO.A56 < 0.05, z_score_WT_vs_xKO.A56 > 0) %>%
  filter(kmer %in% DRACH_motifs$DRACH)

# Save the filtered data to a new TSV filec
write_csv(filtered_data, "A56_filtered_majority_direction_kmer_diffmod.table")

# Print a message to let the user know the file has been saved
print("The filtered data has been saved to 'all_KO_filtered_majority_direction_kmer_diffmod.table'.")


# Read the TSV file
data3 <- read_csv("A56_filtered_majority_direction_kmer_diffmod.table")

data_char <- data3 %>%
  mutate_all(as.character)

# Replace NA values with blank spaces
data_no_na <- data_char %>%
  mutate_all(~replace_na(., ""))

# Write the modified data to a new CSV file
write_csv(data_no_na, "A56_filtered_majority_direction_kmer_diffmod.table")

# Print a message to let the user know the file has been saved
print("The modified data has been saved to 'modified_file.csv'.")
