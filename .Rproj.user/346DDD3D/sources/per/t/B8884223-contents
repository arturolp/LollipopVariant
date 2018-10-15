##########################################################
# 
# CLINGEN project: PTEN analysis
#
# Author: Arturo Pineda <arturolp@stanford.edu>
#
# Date: Sep 27, 2018
#
##########################################################

rm(list=ls())


#--------------------------------------
# 0. Load Library
#--------------------------------------

library(rentrez)
library(ggplot2)
library(stringr)
library(ggpubr)
library(MASS)
library(viridis)

#--------------------------------------
# 1. Query Clinvar for automatic retrieval
#--------------------------------------

# Gene PTEN, ID: 5728
# geneID = 5728


# all_the_links <- entrez_link(dbfrom='gene', id=geneID, db='all')
# clinvar_ids <- all_the_links$links$gene_clinvar
# snp_ids <- all_the_links$links$gene_snp
# 
# clinvar_summ <- entrez_summary(db="clinvar", id=clinvar_ids[1])
# 
# 
# gene_clinvar <- all_the_links$links$gene_clinvar
# 
# 
# snp_ids<- entrez_search(db="snp", term="PTEN[gene]")$ids
# 
# snp_links <- entrez_link(dbfrom="clinvar", db="snp", id=clinvar_ids[1])
# 
# snp_summ <- entrez_summary(db="snp", id=1491555086)
# 
# #Obtain Summaries
# alleles = list()
# for(variant in gene_clinvar[1:50]){
#   clinvar_summ <- entrez_summary(db="clinvar", id=variant)
#   chr = clinvar_summ$chr_sort
#   source = unlist(clinvar_summ$genes)[4]
#   snp_summ <- entrez_link(dbfrom="snp", id="536571[alleleid]")
#   pos = str_split(snp_summ$chrpos, ":", simplify=TRUE)[2]
#   
#   alleles = rbind(alleles, c(variant, chr, pos, allele, source))
# }
# 
# 
# 
# res <- entrez_search(db = "clinvar", term = "PTEN")
# cv <- entrez_summary(db="clinvar", id=res$ids)
# 
# 
# # Search for PTEN entries in Clinvar
# res <- entrez_search(db = "clinvar", term = "PTEN[GENE]", retmax=100)



#--------------------------------------
# 1. Read Data
#--------------------------------------

pten.clinvar <- read.csv("data/20181004_pten_clinvar_mutations.csv", header=TRUE)
variants <- pten.clinvar[, c("Name", "Clinical.significance..Last.reviewed.")]
colnames(variants) = c("name", "significance")

#---------------
#Alleles ALL

#Extract locations
locations = str_match(variants[,1], "(p\\.[a-zA-Z]*)([0-9]*)")
#Extract significance
sig = str_split(variants[,2], "\\(", simplify = TRUE)[,1]

#Create Variant-Significance table
alleles = data.frame(allele=as.numeric(locations[,3]), significance=sig)

#Count
allele.count = table(locations[,3])
alleles.all = data.frame(x=as.numeric(names(allele.count)), y=as.vector(allele.count))


#---------------
#Categories

pathogenic = c("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic")
vus = c("Uncertain significance")
benign = c("Likely benign", "Benign", "Benign/Likely benign")
conflicting = c("Conflicting interpretations of pathogenicity")

#---------------
#Pathogenic
index = which(sig %in% pathogenic)
allele.count = table(locations[index,3])
alleles.pathogenic = data.frame(x=as.numeric(names(allele.count)), y=as.vector(allele.count))


#---------------
#VUS
index = which(sig %in% vus)
allele.count = table(locations[index,3])
alleles.vus = data.frame(x=as.numeric(names(allele.count)), y=as.vector(allele.count))


#---------------
#Benign
index = which(sig %in% benign)
allele.count = table(locations[index,3])
alleles.benign = data.frame(x=as.numeric(names(allele.count)), y=as.vector(allele.count))


#---------------
#Conflicting
index = which(sig %in% conflicting)
allele.count = table(locations[index,3])
alleles.conflicting = data.frame(x=as.numeric(names(allele.count)), y=as.vector(allele.count))


#--------------------------------------
# 2. Create Lollipop plot 
#--------------------------------------

# Specific data for the PROTEIN and DOMAINS
protein = data.frame(symbol="PTEN", start=1, end=403)
domains = data.frame(start=c(1,14,190,351,401), end=c(13,185,350,400,403), name=c("PBD", "P-Loop", "C2", "C-tail", "PDZ-PD"))


# Lollipop plot for ALL
getLollipopPlot <- function(protein, domains, mydata, maxV, title, col="gray33"){
  
  #Visualization options
  voff = 0
  seg_size = c(6, 10)
  
  plot <- ggplot(data=mydata) +
    #Add points
    geom_segment(aes(x=x, xend=x, y=voff, yend=y), color="gray88", alpha=0.8) +
    geom_point(aes(x=x, y=y), fill=col, colour="gray88", size=3, shape=21, alpha=0.8) +
    #Plot protein
    geom_segment(data=protein, aes(x=start, xend=end, y=voff, yend=voff), color="gray88", size=seg_size[1]) +
    #Add domain
    geom_segment(data=domains, aes(x=start, xend=end, y=voff, yend=voff, color=name), size=seg_size[2]) +
    geom_text(data=domains, aes(x=(start+(end-start)/2), y=voff, label=name), color="white", show.legend = FALSE) +
    #Labels
    ylab("# Mutations") +
    labs(subtitle = paste("N = ", sum(mydata$y), sep="")) +
    #Scaling
    scale_color_manual("Domains", 
                       values = c("#66c2a5", "#fc8d62", "#8da0cb", "#a6d854", "#e78ac3"),
                       breaks = domains$name) +
    scale_x_continuous(breaks=c(seq(0, protein$end-5, 25), protein$end)) +
    scale_y_continuous(breaks=seq(0, maxV, 1), limits=c(0, maxV)) +
    #Theme
    theme_light() +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.title.x=element_blank(),
          #plot.title = element_text(vjust=-3, face="bold", color=col),
          plot.subtitle = element_text(hjust=1, face="bold")
    )
  
  return(plot)
  
}


# Lollipop plot for All variants
maxV <- max(alleles.all$y)
gAll <- getLollipopPlot(protein, domains, alleles.all, maxV, "All variants", "black")

# Lollipop plot for Pathogenic
g1 <- getLollipopPlot(protein, domains, alleles.pathogenic,  maxV, "Pathogenic", "firebrick4")

# Lollipop plot for VUS
g2 <- getLollipopPlot(protein, domains, alleles.vus,  maxV, "Uncertain significance", "darkslategrey")

# Lollipop plot for Benign
g3 <- getLollipopPlot(protein, domains, alleles.benign,  maxV, "Benign", "darkgreen")

# Lollipop plot for Conflicting
g4 <- getLollipopPlot(protein, domains, alleles.conflicting,  maxV, "Conflicting", "dodgerblue4")

#Save the plot
gpanels <- ggarrange(gAll, g1, g2, g3, g4,
                     nrow = 5, ncol = 1,
                     labels = c("All variants", "Pathogenic", "Uncertain", "Benign", "Conflicting"),
                     legend = "bottom", common.legend = TRUE)
gpanels <- annotate_figure(gpanels,
                           top = text_grob("PTEN variants in ClinVar", face = "bold", size = 14))
ggexport(gpanels, filename="figures/pten-clinvar-significance.png", width = 4000, height = 4000, res=300)

