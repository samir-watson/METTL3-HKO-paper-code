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
revigo.data <- rbind(c("GO:0002181","cytoplasmic translation",0.3974763229665376,53.96657624451305,0.7058766770557127,0,"cytoplasmic translation"),
c("GO:0006412","translation",4.38869169324396,37.08512818245995,0.6186095441207669,0.6792736,"cytoplasmic translation"),
c("GO:0006413","translational initiation",0.4742529791560868,2.269971936283751,0.7008291853837132,0.6934247,"cytoplasmic translation"),
c("GO:0009058","biosynthetic process",29.004480601854056,1.8854533641642286,0.8741751335812674,0.23252945,"cytoplasmic translation"),
c("GO:0009059","macromolecule biosynthetic process",16.15219652221219,15.536107011014092,0.7085771536965311,0.48020862,"cytoplasmic translation"),
c("GO:0009141","nucleoside triphosphate metabolic process",1.2420663711147437,3.7822923859327062,0.6379873374987007,0.65568623,"cytoplasmic translation"),
c("GO:0016072","rRNA metabolic process",1.653535024007062,4.632353403999483,0.73547291196913,0.30509215,"cytoplasmic translation"),
c("GO:0019538","protein metabolic process",14.600987650961638,5.747856005535273,0.756068435215616,0.48568148,"cytoplasmic translation"),
c("GO:0034641","cellular nitrogen compound metabolic process",24.77172315999013,3.2631534990451883,0.7627482465513585,0.40483001,"cytoplasmic translation"),
c("GO:0043170","macromolecule metabolic process",29.459939937239827,2.3959157099408683,0.8448373985629563,0.23434422,"cytoplasmic translation"),
c("GO:0043603","amide metabolic process",6.707376287050344,26.54211810326601,0.8358636799472879,0.13247309,"cytoplasmic translation"),
c("GO:0044237","cellular metabolic process",46.10331480933213,6.432973633840939,0.8523910642413096,0.13458491,"cytoplasmic translation"),
c("GO:0044249","cellular biosynthetic process",27.441482878386665,2.220510075461862,0.7047416389300166,0.56396741,"cytoplasmic translation"),
c("GO:0044271","cellular nitrogen compound biosynthetic process",12.957025677761646,6.209011524911184,0.6800924890956268,0.48531185,"cytoplasmic translation"),
c("GO:0046034","ATP metabolic process",0.836742315457756,22.826813731587727,0.5626208060409879,0.20362055,"cytoplasmic translation"),
c("GO:1901564","organonitrogen compound metabolic process",27.31959408292786,4.379765016304792,0.7831809193679614,0.33473133,"cytoplasmic translation"),
c("GO:1901566","organonitrogen compound biosynthetic process",14.093783560518295,27.609064892896622,0.7188351724587875,0.27972648,"cytoplasmic translation"),
c("GO:1901576","organic substance biosynthetic process",28.21764434528959,1.9331754628586473,0.7268881858244092,0.68695064,"cytoplasmic translation"),
c("GO:0006091","generation of precursor metabolites and energy",2.5115061465197326,17.18641901143181,0.8674790044037475,0.09281265,"generation of precursor metabolites and energy"),
c("GO:0006119","oxidative phosphorylation",0.4950208797998473,24.97469413473523,0.7089327867244508,0.07702283,"oxidative phosphorylation"),
c("GO:0008152","metabolic process",57.597931274565454,3.0936940558222883,1,-0,"metabolic process"),
c("GO:0042254","ribosome biogenesis",2.136224808753416,9.657577319177793,0.7656995481005959,0.01288502,"ribosome biogenesis"),
c("GO:0000027","ribosomal large subunit assembly",0.012821578346646382,1.724802659844416,0.8019393322658878,0.55914193,"ribosome biogenesis"),
c("GO:0007005","mitochondrion organization",0.8276597479438399,3.9890219878252577,0.8460557874088768,0.38245709,"ribosome biogenesis"),
c("GO:0010257","NADH dehydrogenase complex assembly",0.19259972609835477,5.524373841197004,0.8055352691157451,0.58873399,"ribosome biogenesis"),
c("GO:0022613","ribonucleoprotein complex biogenesis",2.5375929564430133,7.931814138253839,0.8005165634107279,0.65931649,"ribosome biogenesis"),
c("GO:0032981","mitochondrial respiratory chain complex I assembly",0.18765792206432438,5.524373841197004,0.7987640221560887,0.48971776,"ribosome biogenesis"),
c("GO:0071826","protein-RNA complex organization",0.6518695497816881,2.817364261783578,0.8228991046640499,0.55363498,"ribosome biogenesis"),
c("GO:0072593","reactive oxygen species metabolic process",0.3349261470185961,1.7316652636896754,0.8949475018806783,0.07434245,"reactive oxygen species metabolic process"),
c("GO:0098754","detoxification",0.9926322015147898,1.6629717726635977,1,-0,"detoxification"),
c("GO:1902600","proton transmembrane transport",1.312853708699458,13.388276691992658,0.9488436458548791,0.0106876,"proton transmembrane transport"),
c("GO:0006839","mitochondrial transport",0.470038273471203,3.086140072475195,0.9795706288281498,0.28337751,"proton transmembrane transport"),
c("GO:0051238","sequestering of metal ion",0.05868299862654819,1.9487484696972115,0.982864159442404,0.22103554,"proton transmembrane transport"),
c("GO:2001233","regulation of apoptotic signaling pathway",0.13887085520666995,3.1360018315840135,0.9403817710007568,-0,"regulation of apoptotic signaling pathway"),
c("GO:0006417","regulation of translation",1.3923070727099334,2.610049172462242,0.9723619374868295,0.18968422,"regulation of apoptotic signaling pathway"),
c("GO:0034248","regulation of amide metabolic process",1.4389695335439798,1.646953921282559,0.9734523084825246,0.40157068,"regulation of apoptotic signaling pathway"),
c("GO:1902255","positive regulation of intrinsic apoptotic signaling pathway by p53 class mediator",0.0004658358914871522,2.155854729666864,0.9430339545446346,0.53786774,"regulation of apoptotic signaling pathway"));

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

