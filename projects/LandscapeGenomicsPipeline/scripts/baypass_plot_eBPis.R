args <- commandArgs(TRUE)
library(ggplot2)
library(readr)
file <- read.table(args[1],h=T)
arm=args[2]
#file <- read.table(\"/media/inter/ssteindl/FC/LandscapeGenomicsPipeline/results/BAYPASS/BayPass_summary_betai_reg.out\",h=T)
cov_info <- read.table(args[3])

for (cov in 1:max(file$COVARIABLE)){
  print(cov)
  file[file$COVARIABLE==cov,]$eBPis -> CovP
  plotdf <- data.frame(x=c(1:10000), y=CovP)
  pod.thresh=quantile(CovP,probs=0.95)
  PLOT <- ggplot(plotdf, aes(x=x, y=CovP)) + geom_point(alpha=0.3) + 
    geom_hline(yintercept=pod.thresh,col="blue",lty=2)+
    ylab("eBPis (-log10P)") + theme_bw() + ggtitle(paste0(cov_info[,cov])) +xlab(paste0("Genomic Position on ",arm)) 
  ggsave(paste0("results/",arm,"/BAYPASS/",cov_info[,cov],"eBPis.png"),
         PLOT,
         width=15,
         height=5)
}