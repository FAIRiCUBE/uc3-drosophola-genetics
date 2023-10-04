
#install.packages('tidyverse')
library(tidyverse)
library(gridExtra)
library(ggplot2)
library(readr)

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

## read Metadata
#meta=read.csv(args[3],header=T)

meta <- read.table(args[3], header=1, sep=",",  dec = ".")
#meta <- read_csv(file = commandArgs(trailingOnly = TRUE)[3])
#meta <- read.table("/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/results/metadata_v2.csv", header=1, sep=",",  dec = ".")
#colnames(AF) <- samplenamesnew$V1

## make sure that order of AF data and Meta data match
#meta<-meta[match(colnames(AF),meta$sampleId),]

setwd(directory_path)

## make suuuuper simple function that fits a linear regression model and returns the p-values
lm.func <- function(x) {
  summary(lm(unlist(x)~y))[[4]][,4][2]
}

#print(meta)
pvalcsv <- data.frame(DATA$Chr, DATA$Pos)
facs <- c()
for ( i in colnames(meta)){
  y=meta[[i]]
  if (is.numeric(y) && length(table(y)) > 2){
    p.val<-apply(AF,1,lm.func)
    ID=paste0(i,".pval")
    DATA[[ID]]<-p.val
    pvalcsv[[ID]]<- p.val 
    facs <- append(facs,ID)
  } else {
    next
  }
}

#export DATA to directory
#pvalpath=paste(directoy_path, "/P_Values.csv")
write.csv(pvalcsv, "P_Values.csv")

## multiple testing problem??
y=runif(ncol(AF))
D=data.frame("X.1"=runif(nrow(AF)))
for (i in seq(2,ncol(AF),1)){
  D[[paste0("X.",i)]]<-runif(nrow(AF))
}

Bonf=0.05/nrow(DATA)
Test.p<-apply(D,1,lm.func)
Test.p[Test.p<Bonf]


pdf("Multtest_control3.pdf",
    width=15,
    height=5)
hist(Test.p,breaks=100)
#dev.off()

## Boferroni-correctd p-value threshold
Bonf=-log10(0.05/nrow(DATA))


#### plot with ggplot
PLOT<-facs %>%
  map(function(z){
    PLOT.df<-data.frame(Pos = DATA$Pos/1000000, P.val = -log10(DATA[[z]]))
    pl <- ggplot(PLOT.df, aes(x = Pos, y = P.val)) +
      geom_point(alpha=0.3) +
      xlab("Genomic Position on 3R [Mbp]") +
      ylab("-log10(P.value)")+
      geom_hline(yintercept=Bonf,
                 col="blue",
                 lty=2)+
      geom_hline(yintercept=-log10(0.05),
                 col="red",
                 lty=2)+
      ggtitle(z)+
      theme_bw()
    return(pl)
  }) %>%
  arrangeGrob(grobs = ., nrow=3) %>% #add ncol/nrow argument here
  grid.arrange()

ggsave("3R_Pvalues.png",
       PLOT,
       width=15,
       height=7)

