args <- commandArgs(TRUE)
setwd(args[1])
library(readr)
library(ggplot2)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args){
  cat("Example: -working directory ")}

##zscorfile and/or p-value getting
#project = load.lfmmProject("genotypes_gradients.lfmmProject")
#folder <- paste(fol,"genotypes_gradients.lfmm")
#folder <- "/media/inter/ssteindl/FC/test/genotypes_gradients.lfmm/"

DATA<-read_table(args[4])

##from manual, 1:5 needs to be adjustet to number of repetitions parsed by user 
z.table = NULL
j=args[2]
var=args[5]

markerpos<- data.frame(DATA$Chr, DATA$Pos)

for (i in 1:args[3]){
  #fol<-paste("run",i,"/")
  #folder <- paste(fol,"genotypes_gradients.lfmm")
  file.name <- paste(var,"_run",i,"/genotypes_gradients.lfmm/K",j, "/run1/genotypes_r1_s1.3.zscore", sep="")
  z.table = cbind(z.table, read.table(file.name)[,1])
  }
z.score = apply(z.table, MARGIN = 1, median) #combines z-scores
lambda = median(z.score^2)/0.456
adjusted.p.values = pchisq(z.score^2/lambda, df = 1, lower = F) #re-adjust p-values
markerpos <- cbind(markerpos,adjusted.p.values)

png(paste(var,"_adjp.png"))
hist(adjusted.p.values, col = 3)
dev.off()

getwd()

#outfile <- paste(var,"_adjusted_pvals.csv")
write.csv(markerpos, paste(var,"_adjusted_pvals.csv"))

##plot the p values as mp

Bonf=-log10(0.05/ncol(DATA))

PLOT.df<-data.frame(Pos = DATA$Pos/1000000, P.val = -log10(adjusted.p.values))

pl <- ggplot(PLOT.df, aes(x = Pos, y = P.val)) +
geom_point(alpha=0.3) +
  xlab(paste("Genomic Position on 3R")) +
  geom_hline(yintercept=Bonf,
             col="blue",
             lty=2)+
  ylab(paste("-log10(P.value)"))+
  ggtitle(paste("p values for",var))+
  theme_bw()

ggsave(paste(var,"_ManhattanLEA",".png", sep=""),pl, width=15, height=5)