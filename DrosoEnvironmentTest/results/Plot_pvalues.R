
library(tidyverse)
library(gridExtra)

## set working directory
setwd("/media/inter/ssteindl/FC/DROSO_SANDBOX")

## read AlleleFrequency Matrix
DATA=read.table("MK/Subsampled_2L_10k.recode.af",
  header=T)
AF=DATA[,3:ncol(DATA)]

## read Metadata
meta=read.csv("MK/metadata.csv",
  header=T)

## make sure that order of AF data and Meta data match
meta<-meta[match(colnames(AF),meta$sampleId),]

## make suuuuper simple function that fits a linear regression model and returns the p-values
lm.func <- function(x) {
      summary(lm(unlist(x)~y))[[4]][,4][2]
}

## fit linear models for both biovariables to AF data in each row (i.e. SNP) and append corresponding p-value to new column in raw Data table
for ( i in c("bio1","bio12")){
  y=meta[[i]]
  p.val<-apply(AF,1,lm.func)
  ID=paste0(i,".pval")
  DATA[[ID]]<-p.val
}

Bonf=-log10(0.05/ncol(DATA))
## plot with ggplot
PLOT<-c("bio1.pval","bio12.pval") %>%
  map(function(z){
    PLOT.df<-data.frame(Pos = DATA$Pos/1000000, P.val = -log10(DATA[[z]]))
    pl <- ggplot(PLOT.df, aes(x = Pos, y = P.val)) +
      geom_point(alpha=0.3) +
      xlab("Genomic Position on 2L [Mbp]") +
      geom_hline(yintercept=Bonf,
        col="blue",
        lty=2)+
      ylab("-log10(P.value)")+
      ggtitle(z)+
      theme_bw()
    return(pl)
  }) %>%
  arrangeGrob(grobs = ., nrow=2) %>% #add ncol/nrow argument here
grid.arrange()

ggsave("MK/Pvalues.png",
  PLOT,
  width=15,
  height=5)


