library("tidyverse")
library("rtracklayer")
library("biomaRt")
library("ggrepel")

basedir <- ("C:/Users/au656873/OneDrive - Aarhus Universitet/Desktop/METTL3/3rd analysis/filtered data/")
setwd(basedir)
dir()
# Load nanocompore results
nanocompore <- read_tsv(paste0(basedir,"GMM_logit_pvalue_context_2_filter_nanocompore_results.tsv"), col_types="cciciccddddddddddddddddcicdd")
#col_types="icicccddddddddcicdc")
#asdf <- unique(nanocompore$ref_id)
#d <- nanocompore%>% group_by(ref_id)%>%summarise(n=n())
#d<- read_tsv("nanocompore.results.GMM_logit_pvalue.0.01_context_2.sorted.tsv")
#d <- ("nanocompore.results.GMM_logit_pvalue.0.01_context_2.sorted.tsv")
#ncol(d)
#d <- (nanocompore$ref_id)
# Annotate gene names
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl", host="mar99.ensembl.org")
id2name <- getBM(attributes=c("ensembl_transcript_id_version", "hgnc_symbol"), 
                 filters = "ensembl_transcript_id_version", 
                 values = unique(nanocompore$ref_id), 
                 mart = ensembl)
nanocompore <- left_join(nanocompore, id2name, by=c("ref_id"="ensembl_transcript_id_version"))

dplyr::select(nanocompore, pos, genomicPos, ref_id, ref_kmer, GMM_logit_pvalue, hgnc_symbol) %>%
  write_tsv(paste0(basedir,"gmm.logit.annotated_sites.txt"))

gmm.logit.volcano <- nanocompore %>% mutate(LOR=case_when(Logit_LOR=="NC"~0, TRUE~as.numeric(Logit_LOR))) %>% { ggplot(., aes(x=LOR, y=-log10(GMM_logit_pvalue))) + 
    geom_point(size=0.9, alpha=0.8) + 
    ggrepel::geom_text_repel(data=top_n(., 20, -log10(GMM_logit_pvalue)), size=4, aes(label=paste0(hgnc_symbol, "\n", ref_kmer, " (", pos, ")"))) +
    xlab("Logistic regression\nodds ratio") +
    ylab("Nanocompore GMM Logit p-value (-log10)") +
    theme_bw(22)
}


gmm.logit.volcano_abs <- nanocompore %>% mutate(LOR=case_when(Logit_LOR=="NC"~0, TRUE~as.numeric(Logit_LOR))) %>% filter(!is.na(GMM_logit_pvalue), GMM_logit_pvalue<0.1) %>% { ggplot(., aes(x=abs(LOR), y=-log10(GMM_logit_pvalue))) + 
    geom_point(alpha=0.8) + 
    ggrepel::geom_text_repel(data=top_n(., 15, -log10(GMM_logit_pvalue)), size=6, aes(label=paste0(hgnc_symbol, "\n", ref_kmer, " (", pos, ")"))) +
    xlab("Logistic regression\nodds ratio") +
    ylab("Nanocompore p-value (-log10)") +
    theme_bw(22)
}


pdf(paste0(basedir,"/gmm_logit_volcano_plot.pdf"), height=10, width=12)
print(volcano)
dev.off()

pdf(paste0(basedir,"/gmm_logit_volcano_abs_lor_plot.pdf"), height=10, width=12)
print(volcano_abs)
dev.off()
"......................................................................................................................................................................."
#ks intensity

library("tidyverse")
library("rtracklayer")
library("biomaRt")
library("ggrepel")

basedir <- ("C:/Users/au656873/OneDrive - Aarhus Universitet/Desktop/METTL3/3rd analysis/filtered data/")
setwd(basedir)
dir()
# Load nanocompore results
nanocompore <- read_tsv(paste0(basedir,"KS_intensity_pvalue_context_2_filter_nanocompore_results.tsv"), col_types="cciciccddddddddddddddddcicdd")
#col_types="icicccddddddddcicdc")
#asdf <- unique(nanocompore$ref_id)
#d <- nanocompore%>% group_by(ref_id)%>%summarise(n=n())
#d<- read_tsv("nanocompore.results.GMM_logit_pvalue.0.01_context_2.sorted.tsv")
#d <- ("nanocompore.results.GMM_logit_pvalue.0.01_context_2.sorted.tsv")
#ncol(d)
#d <- (nanocompore$ref_id)

# Annotate gene names
ensembl1 = useEnsembl("ensembl",dataset="hsapiens_gene_ensembl", mirror="useast")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
id2name <- getBM(attributes=c("ensembl_transcript_id_version", "hgnc_symbol"), 
                 filters = "ensembl_transcript_id_version", 
                 values = unique(nanocompore$ref_id), 
                 mart = ensembl)
nanocompore <- left_join(nanocompore, id2name, by=c("ref_id"="ensembl_transcript_id_version"))

dplyr::select(nanocompore, pos, genomicPos, ref_id, ref_kmer, KS_intensity_pvalue, hgnc_symbol) %>%
  write_tsv(paste0(basedir,"ks.intensity.annotated_sites.txt"))

volcano <- nanocompore %>% mutate(LOR=case_when(Logit_LOR=="NC"~0, TRUE~as.numeric(Logit_LOR))) %>% { ggplot(., aes(x=LOR, y=-log10(KS_intensity_pvalue))) + 
    geom_point(size=0.9, alpha=0.8) + 
    ggrepel::geom_text_repel(data=top_n(., 20, -log10(KS_intensity_pvalue)), size=4, aes(label=paste0(hgnc_symbol, "\n", ref_kmer, " (", pos, ")"))) +
    xlab("Logistic regression\nodds ratio") +
    ylab("Nanocompore ks intensity p-value (-log10)") +
    theme_bw(22)
}


volcano_abs <- nanocompore %>% mutate(LOR=case_when(Logit_LOR=="NC"~0, TRUE~as.numeric(Logit_LOR))) %>% filter(!is.na(KS_intensity_pvalue), KS_intensity_pvalue<0.1) %>% { ggplot(., aes(x=abs(LOR), y=-log10(KS_intensity_pvalue))) + 
    geom_point(alpha=0.8) + 
    ggrepel::geom_text_repel(data=top_n(., 15, -log10(KS_intensity_pvalue)), size=6, aes(label=paste0(hgnc_symbol, "\n", ref_kmer, " (", pos, ")"))) +
    xlab("Logistic regression\nodds ratio") +
    ylab("Nanocompore ks intensity p-value (-log10)") +
    theme_bw(22)
}


pdf(paste0(basedir,"/KS_intensity_volcano_plot.pdf"), height=10, width=12)
print(volcano)
dev.off()

pdf(paste0(basedir,"/KS_intensity_volcano_abs_lor_plot.pdf"), height=10, width=12)
print(volcano_abs)
dev.off()


############################ new code for paper figures ###########################################
library("tidyverse")
library("rtracklayer")
library("biomaRt")
library("ggrepel")

basedir <- ("C:/Users/au656873/OneDrive - Aarhus Universitet/Desktop/METTL3/3rd analysis/filtered data/")
setwd(basedir)
dir()

nanocompore <- read_tsv(paste0(basedir,"GMM_logit_pvalue_context_2_filter_nanocompore_results.tsv"), col_types="cciciccddddddddddddddddcicdd")

# Use the archive URL for release 99
listMarts(host="https://jan2020.archive.ensembl.org/")
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="https://jan2020.archive.ensembl.org")

id2name <- getBM(attributes=c("ensembl_transcript_id_version", "hgnc_symbol"), 
                 filters = "ensembl_transcript_id_version", 
                 values = unique(nanocompore$ref_id), 
                 mart = mart)
nanocompore <- left_join(nanocompore, id2name, by=c("ref_id"="ensembl_transcript_id_version"))

dplyr::select(nanocompore, pos, genomicPos, ref_id, ref_kmer, GMM_logit_pvalue_context_2, hgnc_symbol) %>%
  write_tsv(paste0(basedir,"gmm.logit.annotated_sites.txt"))

nanocompore$GMM_logit_pvalue_context_2 <- as.numeric(nanocompore$GMM_logit_pvalue_context_2)
gmm.logit.volcano <- nanocompore %>% mutate(LOR=case_when(Logit_LOR=="NC"~0, TRUE~as.numeric(Logit_LOR))) %>% { ggplot(., aes(x=LOR, y=-log10(GMM_logit_pvalue_context_2))) + 
    geom_point(size=0.9, alpha=0.8) + 
    ggrepel::geom_text_repel(data=top_n(., 20, -log10(GMM_logit_pvalue_context_2)), size=4, aes(label=paste0(hgnc_symbol, "\n", ref_kmer, " (", pos, ")"))) +
    xlab("Logistic regression\nodds ratio") +
    ylab("Nanocompore GMM Logit p-value (-log10)") +
    theme_bw(22)
}

gmm.logit.volcano_abs <- nanocompore %>% mutate(LOR=case_when(Logit_LOR=="NC"~0, TRUE~as.numeric(Logit_LOR))) %>% filter(!is.na(GMM_logit_pvalue_context_2), GMM_logit_pvalue_context_2<0.1) %>% { ggplot(., aes(x=abs(LOR), y=-log10(GMM_logit_pvalue_context_2))) + 
    geom_point(alpha=0.8) + 
    ggrepel::geom_text_repel(data=top_n(., 15, -log10(GMM_logit_pvalue_context_2)), size=6, aes(label=paste0(hgnc_symbol, "\n", ref_kmer, " (", pos, ")"))) +
    xlab("Logistic regression\nodds ratio") +
    ylab("Nanocompore p-value (-log10)") +
    theme_bw(22)
}

pdf(paste0(basedir,"gmm_logit_volcano_plot.pdf"), height=10, width=12)
print(gmm.logit.volcano)
dev.off()

pdf(paste0(basedir,"gmm_logit_volcano_abs_lor_plot.pdf"), height=10, width=12)
print(gmm.logit.volcano_abs)
dev.off()

"......................................................................................................................................................................."
#ks intensity

library("tidyverse")
library("rtracklayer")
library("biomaRt")
library("ggrepel")

basedir <- ("C:/Users/au656873/OneDrive - Aarhus Universitet/Desktop/METTL3/3rd analysis/filtered data/")
setwd(basedir)
dir()
# Load nanocompore results
nanocompore <- read_tsv(paste0(basedir,"KS_intensity_pvalue_context_2_filter_nanocompore_results.tsv"), col_types="cciciccddddddddddddddddcicdd")
#col_types="icicccddddddddcicdc")
#asdf <- unique(nanocompore$ref_id)
#d <- nanocompore%>% group_by(ref_id)%>%summarise(n=n())
#d<- read_tsv("nanocompore.results.GMM_logit_pvalue.0.01_context_2.sorted.tsv")
#d <- ("nanocompore.results.GMM_logit_pvalue.0.01_context_2.sorted.tsv")
#ncol(d)
#d <- (nanocompore$ref_id)
# Annotate gene names
listMarts(host="https://jan2020.archive.ensembl.org/")
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="https://jan2020.archive.ensembl.org")

id2name <- getBM(attributes=c("ensembl_transcript_id_version", "hgnc_symbol"), 
                 filters = "ensembl_transcript_id_version", 
                 values = unique(nanocompore$ref_id), 
                 mart = mart)
nanocompore <- left_join(nanocompore, id2name, by=c("ref_id"="ensembl_transcript_id_version"))

dplyr::select(nanocompore, pos, genomicPos, ref_id, ref_kmer, KS_intensity_pvalue_context_2, hgnc_symbol) %>%
  write_tsv(paste0(basedir, "ks.intensity.annotated_sites.txt"))

ks.logit.volcano <- nanocompore %>% mutate(LOR=case_when(Logit_LOR=="NC"~0, TRUE~as.numeric(Logit_LOR))) %>% { ggplot(., aes(x=LOR, y=-log10(KS_intensity_pvalue_context_2))) + 
    geom_point(size=0.9, alpha=0.8) + 
    ggrepel::geom_text_repel(data=top_n(., 20, -log10(KS_intensity_pvalue_context_2)), size=4, aes(label=paste0(hgnc_symbol, "\n", ref_kmer, " (", pos, ")"))) +
    xlab("Logistic regression\nodds ratio") +
    ylab("Nanocompore KS context 2 p-value (-log10)") +
    theme_bw(22)
}


ks.logit.volcano_abs <- nanocompore %>% mutate(LOR=case_when(Logit_LOR=="NC"~0, TRUE~as.numeric(Logit_LOR))) %>% filter(!is.na(KS_intensity_pvalue_context_2), KS_intensity_pvalue_context_2<0.1) %>% { ggplot(., aes(x=abs(LOR), y=-log10(KS_intensity_pvalue_context_2))) + 
    geom_point(alpha=0.8) + 
    ggrepel::geom_text_repel(data=top_n(., 15, -log10(KS_intensity_pvalue_context_2)), size=6, aes(label=paste0(hgnc_symbol, "\n", ref_kmer, " (", pos, ")"))) +
    xlab("Logistic regression\nodds ratio") +
    ylab("Nanocompore KS context 2 p-value (-log10)") +
    theme_bw(22)
}


pdf(paste0(basedir,"/KS_intensity_volcano_plot.pdf"), height=10, width=12)
print(ks.logit.volcano)
dev.off()

pdf(paste0(basedir,"KS_intensity_volcano_abs_lor_plot.pdf"), height=10, width=12)
print(ks.logit.volcano_abs)
dev.off()
