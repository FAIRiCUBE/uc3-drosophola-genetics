args <- commandArgs(TRUE)
library(ggplot2)
library(readr)
#XtX=read.table(\"/media/inter/ssteindl/FC/LandscapeGenomicsPipeline/results/BAYPASS/BayPass_summary_pi_xtx.out\",h=T)$M_XtX
XtX=read.table(args[1],h=T)
arm=args[2]
snps <- length(XtX$log10.1.pval.)
pod.thresh=quantile(XtX$log10.1.pval.,probs=0.95)
plotdf <- data.frame(x=c(1:snps), y=XtX$log10.1.pval.)
PLOT <- ggplot(plotdf, aes(x=x, y=y)) + geom_point(alpha=0.3) + 
  geom_hline(yintercept=pod.thresh,col="blue",lty=2)+
  ylab("-log10_pValue") + theme_bw() + xlab(paste("Genomic Position", arm, sep=""))
ggsave(paste0("results/",arm,"/BAYPASS/XtXstats.png"),
       PLOT,
       width=15,
       height=5)



       
       