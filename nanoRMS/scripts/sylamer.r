################################
## CONFIGURATION AND SETTINGS ##
################################
# Name of sylamer executable (can use full path)
sylamer <- "sylamer"
# FASTA file with sequences
utrFile <- "/home/samirwatson/faststorage/NAT10/gencode.v33.transcripts_shorthead.fa"
# Sorted list of sequence identifiers
sorted  <- "/home/samirwatson/faststorage/METTL3/DRS_data/nanocompore/sylamer/nanocompore_sig_intensity_sorted_table.csv"
# Word size, has to be the same length as those in the "words" file
kSylamer <- 5
# Size of a smaller word to be used for correcting composition biases
kMarkov  <- 4
# How many sequences should be added in each consecutive window
winSize  <- 100
# Name of sylamer output
sylOutput <- "test_sylamer.output.tab"
# Name of image output
imageOut  <- "test_image.output.pdf"
# How many best words should I highlight in the image?
topOligos <- 10
lowOligos <- 10
# Any particular words that you want highlighted
chosenOligos <- c()
# Extra parameters
#extras <- "-a 10 -aa 5 --funny-ok"
#################
## RUN SYLAMER ##
#################
# Generate the command line
command <- paste(sylamer, "-fasta", utrFile, "-grow", winSize, "-k", kSylamer, "-m", kMarkov, "-universe", sorted, "-o", sylOutput)
print(command)
# Run the command
error <- system(command)
if (error != 0) {
   stop("Sylamer command generated an error!", call.=FALSE)
}
# Check to see if the file exists
if (!file.exists(sylOutput)) {
   stop(paste("Can't find the expected output file: ",sylOutput,sep=""), call.=FALSE)
}
############################
## PROCESS SYLAMER OUTPUT ##
############################
# hash out the 4 lines below if running the full script
topOligos <- 10
lowOligos <- 0
kSylamer <- 5
chosenOligos <- c()
imageOut  <- "TT_intensity_2.sig_fig_like_paper.pdf"
# Read in and process the result
sylTable <- read.table("sylamer_TT_intensity_2.full.out", sep="\t", row.names=1, header=T, check.names=F)
sylTable <- cbind("0"=0,sylTable)   # To add an initial column of 0s
pdf(file=imageOut, width=10, height=8)
xVals <- as.numeric(colnames(sylTable))
# Generate the plot area
yMin   <- min(sylTable)
yMax   <- max(sylTable)
yRange <- yMax - yMin
plot(NULL, xlab="Sorted sequences", ylab="log10(enrichment P-value)", axes=T, 
     main=paste("Sylamer landscape using words of length: ",kSylamer,sep=""),
     ylim=c(round(yMin-yRange/10),round(yMax+yRange/10)), xlim=range(xVals))
     #ylim=c(-3, 3), xlim=range(xVals))
# It can save time to plot no more than ~1,000 lines (particularly for all words of length 7 or 8)
maxPlot <- 1000
if (maxPlot > 0) {
   sylTable <- sylTable[order(apply(sylTable,1,function(x) max(abs(x))),decreasing=TRUE),]
   if (maxPlot > nrow(sylTable)) {
      maxPlot <- nrow(sylTable)
   }
} else {
   maxPlot <- nrow(sylTable)
}
# Plot the background lines
for (i in 1:maxPlot) {
   lines(xVals,sylTable[i,], col='grey')
}
# Draw a reference line at 0
abline(h=0)
# Up/Down best words
oligosUp <- c()
oligosDown <- c()
oligosChosen <- c()
if (topOligos > 0) { # Only if I really want these plots
   oligosUp <- names((sort(apply(sylTable,1, function(x) {max(x[is.finite(x)])}),decreasing=TRUE))[1:topOligos])
}
if (lowOligos > 0) {
   oligosDown <- names((sort(apply(sylTable,1,function(x) {min(x[is.finite(x)])}),decreasing=FALSE))[1:lowOligos])
}
if (length(chosenOligos) > 0) {
   oligosChosen <- chosenOligos[chosenOligos %in% rownames(sylTable)]
}
oligosAll <- unique(c(oligosUp,oligosDown,oligosChosen))
oligosAll <- names(sort(apply(sylTable,1, function(x) {max(abs(x[is.finite(x)]))})[oligosAll],decreasing=TRUE))
if (length(oligosAll) > 0) {
   colors   <- rainbow(length(oligosAll))
   names(colors) <- oligosAll
   for (i in rev(seq_along(oligosDown))) {
      lines(xVals, sylTable[oligosDown[i],], col=colors[oligosDown[i]], lwd=2)
   }
   for (i in rev(seq_along(oligosUp))) {
      lines(xVals, sylTable[oligosUp[i],], col=colors[oligosUp[i]], lwd=2)
   }
   for (i in rev(seq_along(oligosChosen))) {
      lines(xVals, sylTable[oligosChosen[i],], col=colors[oligosChosen[i]], lwd=2)
   }
   if (topOligos >0) {
      legend('topleft', inset=c(0.01,0.01), legend=oligosUp, lwd=2, lty=1, horiz=TRUE, col=colors[oligosUp],
             cex=0.6, bg='white', title="Words with highest peak")
   }
   if (lowOligos >0) {
      legend('bottomleft', inset=c(0.01,0.01), legend=oligosDown, lwd=2, lty=1, horiz=TRUE, col=colors[oligosDown],
             cex=0.6, bg='white', title="Words with lowest peak")
   }
   if (length(oligosChosen) >0) {
      legend('topright', inset=c(0.01,0.01), legend=oligosChosen, lwd=2, lty=1, horiz=F, ncol=1, col=colors[oligosChosen],
             cex=0.6, bg='white', title="Selected words")
   }
}
dev.off()