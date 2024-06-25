#install.packages('tidyverse')
library(tidyverse)
library(gridExtra)
library(ggplot2)
library(readr)
library(maps)
args <- commandArgs(TRUE)

## set working directory
setwd(args[1])
arm=args[4]
directory_path=paste("results/",arm, "/GM", sep="")
dir.create(directory_path, recursive = TRUE)
## read AlleleFrequency Matrix
DATA=read.table(args[2],
  header=3)
AF=DATA[,3:ncol(DATA)]
meta <- read.table(args[3], header=1, sep=",",  dec = ".")
pvalcsv <- data.frame(DATA$Chr, DATA$Pos)
pvalcsv2 <- data.frame(DATA$Chr, DATA$Pos)
## read Metadata
#meta=read.csv("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k/data/metadata.csv",header=T)
#DATA=read.table("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k/data/Subsampled_Europe_50k.recode.vcf.gz.af",header=3)
cleaned_row_names <- gsub("\\.", "-", colnames(AF))
cleaned_row_names_DATA <- gsub("\\.", "-", colnames(DATA))

colnames(AF) <- cleaned_row_names
colnames(DATA) <- cleaned_row_names_DATA


actsamp<-intersect(cleaned_row_names, meta$sampleId)
actsamp2<-intersect(cleaned_row_names_DATA, meta$sampleId)
actsamp2 <- c("Chr", "Pos", actsamp2)
                    

#meta <- read_csv(file = commandArgs(trailingOnly = TRUE)[3])
#meta <- read.table("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/testNSAT/data/NSAT.csv", header=1, sep=",",  dec = ".")



AF <- AF[, actsamp, drop = FALSE]
DATA <- DATA[, actsamp2, drop = FALSE]
#AF2 <- AF[meta$sampleId %in% unlist(colnames(AF)),]
meta <- meta[meta$sampleId %in% colnames(AF),]

setwd(directory_path)

## make suuuuper simple function that fits a linear regression model and returns the p-values
lm.func <- function(x) {
  summary(lm(unlist(x)~y))[[4]][,4][2]
}

print(colnames(meta))
#pvalcsv <- data.frame(DATA$Chr, DATA$Pos)
facs <- c()
facs2 <- c()
for ( i in colnames(meta)){
  y=meta[[i]]
  print(y)
  if (is.numeric(y) && length(table(y)) > 2){
    p.val<-apply(AF,1,lm.func)
    p.val.arcsin<-apply(asin(sqrt(AF**2)),1,lm.func)
    ID=paste0(i,".pval")
    ID2=paste0(i,".pval.arcsin")
    DATA[[ID]]<-p.val
    pvalcsv[[ID]]<- p.val 
    facs <- append(facs,ID)
    DATA[[ID2]]<-p.val.arcsin
    pvalcsv2[[ID2]]<- p.val.arcsin
    facs2 <- append(facs2,ID2)
    write.csv(pvalcsv, paste(i,"_pvalues.csv", sep=""))
    write.csv(pvalcsv2, paste(i,"_pvalues.arcsin.csv", sep=""))
    pdf(paste0("Histograms_P_Values_",i,".pdf"),
    width=15,
    height=5)
    hist(as.numeric(p.val),breaks=100, main=paste0("Histogram of p-values", i), xlab="p-value")
    hist(as.numeric(p.val.arcsin),breaks=100, main=paste0("Histogram of p-values", i), xlab="p-value")
  } else {
    next
  }
}



#export DATA to directory
#pvalpath=paste(directoy_path, "/P_Values.csv")
write.csv(pvalcsv, "P_Values.csv")
write.csv(pvalcsv2, "P_Values.arcsin.csv")

## multiple testing problem??
y=runif(ncol(AF))
D=data.frame("X.1"=runif(nrow(AF)))
for (i in seq(2,ncol(AF),1)){
  D[[paste0("X.",i)]]<-runif(nrow(AF))
}

Bonf=0.05/nrow(DATA)
Test.p<-apply(D,1,lm.func)
Test.p[Test.p<Bonf]


pdf("Multtest_control.pdf",
    width=15,
    height=5)
hist(Test.p,breaks=100, main="Histogram of p-values overall", xlab="p-value")
#dev.off()

outliersBonf <- pvalcsv2[which(pvalcsv2$lat.pval.arcsin < Bonf),]
OlBonf <- as.data.frame(paste0(outliersBonf$DATA.Chr,".", outliersBonf$DATA.Pos))
write.csv(OlBonf, "LinearOutliers_ArcSin.csv")

## Boferroni-correctd p-value threshold
#Bonf=-log10(0.05/nrow(DATA))

##### plot with ggplot
#PLOT<-facs %>%
#  map(function(z){
#    PLOT.df<-data.frame(Pos = DATA$Pos/1000000, P.val = -log10(DATA[[z]]))
#    pl <- ggplot(PLOT.df, aes(x = Pos, y = P.val)) +
#      geom_point(alpha=0.3) +
#      xlab("Genomic Position on 3R [Mbp]") +
#      ylab("-log10(P.value)")+
#      geom_hline(yintercept=Bonf,
#                 col="blue",
#                 lty=2)+
#      geom_hline(yintercept=-log10(0.05),
#                 col="red",
#                 lty=2)+
#      ggtitle(z)+
#      theme_bw()
#    return(pl)
#  }) %>%
#  arrangeGrob(grobs = ., nrow=3) %>% #add ncol/nrow argument here
#  grid.arrange()
#
#ggsave("3R_Pvalues.png",
#       PLOT,
#       width=15,
#       height=7)

## Create a ggplot function for mapping
#ggplot_func <- function(PLOT.df, Bonf, title) {
#  ggplot(PLOT.df, aes(x = Pos, y = P.val)) +
#    geom_point(alpha = 0.3) +
#    xlab("Genomic Position [Mbp]") +
#    ylab("-log10(P.value)") +
#    geom_hline(yintercept = Bonf, col = "blue", lty = 2) +
#    geom_hline(yintercept = -log10(0.05), col = "red", lty = 2) +
#    ggtitle(title) +
#    theme_bw()
#}

## Apply ggplot function to each factor
#PLOT_list <- lapply(facs, function(z) {
#  PLOT.df <- data.frame(Pos = DATA$Pos / 1000000, P.val = -log10(DATA[[z]]))
#  ggplot_func(PLOT.df, Bonf, z)
#})
#
## Arrange and save the plots
#pdf("3R_Pvalues.pdf", width = 12, height = 7)
#grid.arrange(grobs = PLOT_list, nrow = 3)
#dev.off()

#for (i in seq_along(facs)) {
#  z <- facs[i]
#  PLOT.df <- data.frame(Pos = DATA$Pos / 100000, P.val = -log10(DATA[[z]]))
#  pl <- ggplot_func(PLOT.df, Bonf, z)
#  
#  # Export as PNG
#  png(filename = paste0(z, "_Pvalues.png"), width = 800, height = 600)
#  print(pl)
#  dev.off()
#}


