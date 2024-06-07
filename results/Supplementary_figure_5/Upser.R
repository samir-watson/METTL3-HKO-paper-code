library("UpSetR")
library(tidyverse)
library(wesanderson)


basedir <- ("C:/Users/au656873/OneDrive - Aarhus Universitet/Desktop/METTL3/3rd analysis/filtered/GO")
setwd(basedir)
dir()

data <- read_csv("significant gene lists all stats tests.UPSETR.csv", col_names=TRUE)
head(data)


# Remove the second row(headers are true so its the first row)
data <- slice(data, -1)
head(data)

list_of_vectors <- lapply(data, as.vector)

# Create the UpSet plot
p1 <- upset(fromList(list_of_vectors), order.by = "freq")


############## exclude last column on table
# Exclude the last column
data_without_last_column <- data[, -ncol(data)]

# Convert each column to a vector and store them in a list
list_of_vectors2 <- lapply(data_without_last_column, as.vector)
png("Upset.plot.nanocompore.stats.png", width = 10, height = 10, units = "in", res = 300)
upset(fromList(list_of_vectors2), order.by = "freq")
dev.off()
