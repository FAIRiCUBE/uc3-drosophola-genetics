##analyse permutations

library(dplyr)
library(pegas)
library(ggplot2)
library(raster)
library(ggrepel)
library(LEA)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
library(ggpubr)
library(vegan)
library(ggVennDiagram)
library(cowplot)
library(corrplot)
library(robust)
#BiocManager::install("qvalue")
library(qvalue)
library(ggpubr)
library(geodata)
library('corrr')
#install.packages("ggcorrplot")
library(ggcorrplot)
#install.packages("FactoMineR")
library(FactoMineR)
#install.packages("factoextra")
library(factoextra)
#library(dplyr)
library(tidyverse)
library(UpSetR)


folder <- "/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/RDA_ResearchPlan_FilterMethod1/permutations_geo_control"
permutationFiles <- list.files("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/RDA_ResearchPlan_FilterMethod1/permutations_geo_control", pattern = "pvalues.csv$",recursive = TRUE,  full.names = TRUE)

permuted_pvalues <- list()

for (file in permutationFiles) {
  data <- read.csv(file)
  number <- gsub("\\D", "", basename(dirname(file)))
  permuted_pvalues[[number]] <- data$x
}

pv <- as.data.frame(permuted_pvalues)
#meanPvals <- rowMeans(pv)
pvm <- as.matrix(pv)
hist(as.numeric(pvm[1,]))

trueRDA <- readRDS("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/RDA_ResearchPlan_FilterMethod1/AssociationAnalysis/RDA_env.rds")
source("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA-landscape-genomics/src/rdadapt.R")
rdadapt_env<-rdadapt(trueRDA, 2)

SNPs <- colnames(trueRDA$Ybar)

write.csv(rdadapt_env$p.values, paste0(outdir, "/permut", permut,"/rdadapt_pvalues.csv"), row.names = FALSE)

MeanPermPval<-rowMeans(pv)
SortMeanPermPval<-sort(MeanPermPval)

##original
#empty<-c()
#SortObs<-sort(rdadapt_env$p.values)
#LEN<-length(SortObs)
#for (i in seq(1,LEN,1)){
#  empty<-rbind(empty,length(SortMeanPermPval[SortMeanPermPval <= SortObs[i]])/i)
#               }
##### naming repeat
SO <- as.data.frame(rdadapt_env$p.values)
colnames(SO) <- "pvalues"
rownames(SO) <- SNPs
SO_SNPs <- SO[order(SO$pvalues), , drop = FALSE]

empty2 <- c()
LEN<-nrow(SO_SNPs)
for (i in seq(1,LEN,1)){
  empty2<-rbind(empty2,length(SortMeanPermPval[SortMeanPermPval <= SO_SNPs[i,]])/i)
}

rownames(empty2) <- rownames(SO_SNPs)
d <- empty2[empty2 <= 0.05,]

# Apply strsplit to each element
split_vec <- strsplit(names(d), "\\.")
tsv_SNPs <- as.data.frame(do.call(rbind, split_vec))
write.table(names(d),"/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/RDA_ResearchPlan_FilterMethod1/AssociationAnalysis/PermutationOutliers.csv", col.names = FALSE, row.names = FALSE)
write.table(tsv_SNPs,"/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/RDA_ResearchPlan_FilterMethod1/AssociationAnalysis/PermutationOutliers.tsv", col.names = FALSE, row.names = FALSE, quote = FALSE)

### approach 

long_column <- unlist(pv, use.names = FALSE)
hist(long_column)

hist(rdadapt_env$p.values)

length(long_column[long_column <= sort(rdadapt_env$p.values)[2] ])/200


outliers <- as.data.frame(names(d))
colnames(outliers) <- "Loci"
ann="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/data/annotationdata/annotations.txt"
annotations <- read.table(ann, na.strings = NaN)
annotations <- unique(annotations)

colnames(annotations) <- c("Chr", "Pos", "Gene")
rownames(annotations) <- paste0(annotations$Chr,".", annotations$Pos)


locus_scores <- scores(trueRDA, choices=c(1:2), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = rownames(locus_scores), locus_scores)
TAB_loci$chromosome <- sub("\\..*", "", TAB_loci$names)
TAB_loci$type <- "All Loci"
TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "Outlier"
#TAB_loci$type[TAB_loci$names%in%top_1p_chromosome_outliers] <- "Top1"
TAB_loci$type <- ifelse(TAB_loci$type == "Top1", 
                        paste0("Chr ", TAB_loci$chromosome), 
                        TAB_loci$type)
TAB_loci$type <- factor(TAB_loci$type, levels = c(unique(TAB_loci$type)))
TAB_loci <- TAB_loci[order(TAB_loci$type),]
TAB_var <- as.data.frame(scores(RDA_env, choices=c(1,2), display="bp")) # pull the biplot scores

TAB_loci$Gene <- ifelse(
  rownames(TAB_loci) %in% rownames(annotations),
  annotations$Gene[match(rownames(TAB_loci), rownames(annotations))],
  NA
)

TAB_var <- as.data.frame(scores(trueRDA, choices=c(1,2), display="bp")) # pull the biplot scores


PLOT1 <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color=gray(0.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color=gray(0.80), size=0.6) +
  
  # Conditionally color points by chromosome if type != "neutral"
  geom_point(
    data = TAB_loci,
    aes(
      x = RDA1 * 20,
      y = RDA2 * 20,
      colour = ifelse(type == "All Loci", "All Loci", chromosome)
    ),
    size = 2
  ) +
  
  geom_text_repel(
    data = TAB_loci[TAB_loci$type != "All Loci", ],  # Only label colored genes
    aes(x = RDA1 * 20, y = RDA2 * 20, label = Gene),
    size = 2.5,
    vjust = 1.5,
    hjust = 0.5,
    family = "Times" 
    #max.overlaps = 80
  ) +
  
  scale_color_manual(
    values = c("All Loci" = "gray90", "2L"= "#F9A242FF", "3R"="#44ad61", "3L"="#7761d0", "2R"= "#4ab09c", "X"="#967dca", "Y"="#588dcc")
  ) +
  
  geom_segment(
    data = TAB_var,
    aes(xend = RDA1, yend = RDA2, x = 0, y = 0),
    colour = "black",
    size = 0.15,
    linetype = 1,
    arrow = arrow(length = unit(0.02, "npc"))
  ) +
  
  geom_text(
    data = TAB_var,
    aes(x = 1.1 * RDA1, y = 1.1 * RDA2, label = row.names(TAB_var)),
    size = 2.5,
    family = "Times"
  ) +
  
  xlab("RDA 1") + 
  ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  guides(color = guide_legend(title = "Locus type")) +
  theme_bw(base_size = 11, base_family = "Times")

ggsave("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/RDA_ResearchPlan_FilterMethod1/AssociationAnalysis/permuted_RDA_SNPs.png", plot = PLOT1, width = 8, height = 6, dpi = 300)


