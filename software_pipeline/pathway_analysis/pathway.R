################Ananlysis of Reh samples#######################
################ Go analysis of DE genes#######################
################Go analysis of DE genes #######################
# author Rinki Ratnapriya ### 03.27.17

# to restsrt r
#.rs.restartR()

# Ensure clean R environment
rm(list=ls())

#install all the R packages
list.of.packages <- c("limma","biomaRt","edgeR", "apcluster", "rgl", "xlsx”, “dplyr", "ggbiplot", "ggplot2", "org.Hs.eg.db", "GO.db")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


# Change global default setting so every data frame created will not auto-convert to factors unless
# explicitly instructed
closeAllConnections()
options(stringsAsFactors = FALSE) 


library(limma)
library(biomaRt)
library(edgeR)
library(apcluster)
library(rgl)
library(xlsx)
library(dplyr)
library(ggbiplot)
library(ggplot2)
library(org.Hs.eg.db)
library(GO.db)


# Session Information
#output to xxx.txt
#sink("GOpathway/session_info_GO_human_fetal_DE.txt")
#sessionInfo()
#sink()

######perform enrichment analysis
#use listAtrribute(mart = ens) --> the name change from entrezgene --> entrezgene_id
# Obtain Entrez IDs as well as Ensembl IDs so that GOANA can be performed. You might want to add external_gene_name to the attributes you're grabbing from BioMart if you want to be able to look back at the gene names.
ens <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene.entrez <- data.frame(getBM(attributes = c("ensembl_gene_id","entrezgene_id"), mart = ens))
genes <- gene.entrez$ensembl_gene_id
length(genes) # 64,101 -->67,405
genes <- genes[genes != ""] # 64,101
length(genes)
length(unique(genes)) # 63,305 -->67,140
entrez <- gene.entrez$entrezgene_id
entrez <- entrez[entrez != ""]
length(entrez)
length(unique(entrez)) # 25,052....note the severe redundancy in entrez IDs  -->19,939


# These are the Ensembl IDs of interest that you want to run GOANA on.
#make sure read the correct directory
gene_ids <- read.csv("~/D2K_BCM_DATASET/ranking.csv")
entrez_list <- gene.entrez$entrezgene[gene.entrez$ensembl_gene_id %in% gene_ids$rf.smote]
length(entrez_list) # 298...note that there are more entrez IDs than Ensembl IDs....part of redundancy

# Remove non-terminal pathway nodes....BP = biological Pathway, MF = Molecular function, CC= cellular component
goana_bp <- as.list(GOBPCHILDREN)
bp_nochildren <- goana_bp[is.na(goana_bp)]
bp_nochildren <- names(bp_nochildren)

goana_cc <- as.list(GOCCCHILDREN)
cc_nochildren <- goana_cc[is.na(goana_cc)]
cc_nochildren <- names(cc_nochildren)

goana_mf <- as.list(GOMFCHILDREN)
mf_nochildren <- goana_mf[is.na(goana_mf)]
mf_nochildren <- names(mf_nochildren)



# Obtain all the GOANA results, then filter only those for biological process ontology (most meaningful to interpret), then remove all
# children
goana_results_bp <- goana(de = entrez_list, species = "Hs")
goana_results_bp <- goana_results_bp[goana_results_bp$Ont == "BP" & rownames(goana_results_bp) %in% bp_nochildren, ]
goana_results_bp <- goana_results_bp[order(goana_results_bp$P.DE), ]
rank.goana_results_bp <- 1:dim(goana_results_bp)[1]
dim(goana_results_bp) 



#Calculate FDR
goana_results_bp$FDR <- goana_results_bp$P.DE*dim(goana_results_bp[goana_results_bp$DE!=0,])[1]/rank.goana_results_bp
goana_results_bp$Ratio <- goana_results_bp$DE/goana_results_bp$N

# order the results by FDR
goana_results_bp <- goana_results_bp[order(goana_results_bp$FDR),]
write.csv(goana_results_bp, file = "~/ls_1.csv", row.names = FALSE)


#plotting the results
goana.df <- goana_results_bp
goana.df.filt = goana_results_bp[1:20, ] # get tp 20 enrighmed pathways

#have NaN
goana.df.filt$V1 <- factor(goana.df.filt$V1, levels = goana.df.filt$V1)

#store the plots
pdf("~/goanaBP_ls1.pdf", width = 8, height = 6)
p <- ggplot(goana.df.filt, aes(x = Ont, y = Term)) 
p + geom_point(aes(colour = P.DE, size = N))
dev.off()


pdf("~/goanaBPn_ls1.pdf", width = 8, height = 6)
p <- ggplot(goana.df.filt, aes(x = Ont, y = Term, color = FDR, size = Ratio)) 
p + geom_point() + scale_colour_gradient(low = "red", high = "blue")
dev.off()


