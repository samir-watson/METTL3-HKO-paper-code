# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
setwd("C:/Users/samir/Desktop/Nanocompore");
# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0002181","cytoplasmic translation",0.3974763229665376,35.45222529461218,0.6465054722887054,0,"cytoplasmic translation"),
c("GO:0006412","translation",4.38869169324396,23.844663962534938,0.533852427753572,0.6792736,"cytoplasmic translation"),
c("GO:0009058","biosynthetic process",29.004480601854056,6.377785977033705,0.8055178255746622,0.23252945,"cytoplasmic translation"),
c("GO:0009059","macromolecule biosynthetic process",16.15219652221219,12.8153085691824,0.6022056119624007,0.48020862,"cytoplasmic translation"),
c("GO:0010467","gene expression",12.663260691525247,4.868654058124139,0.5944558458870809,0.48245041,"cytoplasmic translation"),
c("GO:0016072","rRNA metabolic process",1.653535024007062,2.4738862006489364,0.693270269248906,0.30509215,"cytoplasmic translation"),
c("GO:0019538","protein metabolic process",14.600987650961638,4.81412934475893,0.6563739962375984,0.48568148,"cytoplasmic translation"),
c("GO:0034641","cellular nitrogen compound metabolic process",24.77172315999013,5.269610219703399,0.6676700428749966,0.32496334,"cytoplasmic translation"),
c("GO:0043170","macromolecule metabolic process",29.459939937239827,2.423037603678227,0.7703872886515689,0.23434422,"cytoplasmic translation"),
c("GO:0043603","amide metabolic process",6.707376287050344,17.899629454882437,0.781024057201327,0.13247309,"cytoplasmic translation"),
c("GO:0044237","cellular metabolic process",46.10331480933213,6.903089986991944,0.7768960249499441,0.13458491,"cytoplasmic translation"),
c("GO:0044249","cellular biosynthetic process",27.441482878386665,6.815308569182402,0.599566135664864,0.56396741,"cytoplasmic translation"),
c("GO:0044271","cellular nitrogen compound biosynthetic process",12.957025677761646,9.380906669373257,0.590475549042635,0.48531185,"cytoplasmic translation"),
c("GO:0046034","ATP metabolic process",0.836742315457756,8.63264407897398,0.7364816345062493,0.20362055,"cytoplasmic translation"),
c("GO:1901564","organonitrogen compound metabolic process",27.31959408292786,3.675129278736515,0.6978020890381297,0.40483001,"cytoplasmic translation"),
c("GO:1901566","organonitrogen compound biosynthetic process",14.093783560518295,17.459670525209127,0.6132486460950014,0.27972648,"cytoplasmic translation"),
c("GO:1901576","organic substance biosynthetic process",28.21764434528959,6.616184634019569,0.6256090359169569,0.68695064,"cytoplasmic translation"),
c("GO:0006091","generation of precursor metabolites and energy",2.5115061465197326,6.1979107421182675,0.8193467688889731,0.09060152,"generation of precursor metabolites and energy"),
c("GO:0008152","metabolic process",57.597931274565454,4.868022244988635,1,-0,"metabolic process"),
c("GO:0042254","ribosome biogenesis",2.136224808753416,3.100892785411888,0.8946079356997922,0.01288502,"ribosome biogenesis"),
c("GO:0022613","ribonucleoprotein complex biogenesis",2.5375929564430133,1.7814709278509522,0.9083096277748428,0.65931649,"ribosome biogenesis"),
c("GO:0042773","ATP synthesis coupled electron transport",0.3502371129335071,8.476253533188435,0.5402244279559913,0.07463964,"ATP synthesis coupled electron transport"),
c("GO:0022900","electron transport chain",0.8207017864535318,6.97061622231479,0.573446440097228,0.68695204,"ATP synthesis coupled electron transport"),
c("GO:1902600","proton transmembrane transport",1.312853708699458,5.92769862260329,0.9802707750797949,0.0106876,"proton transmembrane transport"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "Revigo TreeMap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

