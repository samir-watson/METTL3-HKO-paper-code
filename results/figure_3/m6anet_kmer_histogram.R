library(tidyverse)
library(ggplot2)
library(wesanderson)


basedir <- ("C:/Users/au656873/OneDrive - Aarhus Universitet/Desktop/METTL3 nanocompore results/m6anet/")
setwd(basedir)
dir()


# Load the data
wt <- read_csv("wt.data.site_proba.csv")
a56 <- read_csv("a56.data.site_proba.csv")

# Filter the data
a56_filtered <- a56 %>% filter(mod_ratio < wt$mod_ratio - 0.05)

# Calculate the average and standard deviation of mod_ratio for each kmer
wt_summary <- wt %>% group_by(kmer) %>% summarise(avg = mean(mod_ratio), sd = sd(mod_ratio), se = sd/sqrt(n()))
a56_summary <- a56_filtered %>% group_by(kmer) %>% summarise(avg = mean(mod_ratio), sd = sd(mod_ratio), se = sd/sqrt(n()))

# Merge the two summaries
summary <- full_join(wt_summary, a56_summary, by = "kmer", suffix = c("_wt", "_a56"))

# Reorder kmer by average mod_ratio in wt
summary <- summary %>% mutate(kmer = reorder(kmer, -avg_wt))

# Plot the average values with error bars
p <- ggplot(summary, aes(x = kmer)) +
  geom_bar(aes(y = avg_wt), stat = "identity", fill = "blue", alpha = 0.5, position = position_dodge()) +
  geom_errorbar(aes(ymin = avg_wt - se_wt, ymax = avg_wt + se_wt), width = 0.2, position = position_dodge(0.9)) +
  geom_bar(aes(y = avg_a56), stat = "identity", fill = "red", alpha = 0.5, position = position_dodge()) +
  geom_errorbar(aes(ymin = avg_a56 - se_a56, ymax = avg_a56 + se_a56), width = 0.2, position = position_dodge(0.9)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Kmer", y = "Average mod_ratio", title = "Average mod_ratio for each kmer in wt and a56")
p
# Save the plot as a PNG file
ggsave("a56 avg mod ratio.png", plot = p, width = 10, height = 10, dpi = 300)
