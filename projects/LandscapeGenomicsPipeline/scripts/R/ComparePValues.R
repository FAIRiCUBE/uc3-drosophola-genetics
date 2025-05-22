##load data
library(ggplot2)
args <- commandArgs(TRUE)

#var.names=args[1]
###from gm
##from baypass
#baypassfolder=args[2]

##from lea
#leafolder=read.table(args[3])

### 1. Plot p-value distribution of baypass for each covariable
##create paths to locate files 

#bayfile <- read.table("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_11/results/fullgenome/BAYPASS2/BayPass3_summary_betai_reg.out", header=TRUE)
#bayfile <- read.table(args[1], header=TRUE)
affile <- read.table(args[1],h=T)
#affile <- read.table("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_11/polymorphPOS.af",h=TRUE)
corr <- length(affile$Pos)
gmfolderpath <- args[2]
gm_file_all <- read.table(paste(gmfodlerpath, "/P_Values.csv", sep=""),header=TRUE, sep=",")
gm_file_all_arcsin <- read.table(paste(gmfodlerpath, "/P_Values.arcsin.csv", sep=""),header=TRUE, sep=",")

var1 <- paste(args[3])
var2 <- paste(args[4])

##histogram of p values for linear model 

gm_file_all <- read.table("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k_1/results/fullgenome/GM/P_Values.csv",h=T,sep=",")
gm_file_all_arcsin <- read.table("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k_1/results/fullgenome/GM/P_Values.arcsin.csv",h=T,sep=",")

gm_bio1 <- gm_file_all_arcsin$bio1.pval.arcsin
gm_nsat <- gm_file_all_arcsin$near_surface_air_temperature.pval.arcsin

pdf("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k_1/comparison/GM_pvals_raw.pdf",
    width=7,
    height=5)
hist(gm_bio1,breaks=100, main = paste("P-value distribution of 50.000 SNPs for Europe"), sub="Worldclim: Bio1", xlab=paste("p value"))
hist(gm_nsat,breaks=100, main = paste("P-value distribution of 50.000 SNPs for Europe"), sub="Copernicus: Near Surface Air Temperature", xlab=paste("p value"))
dev.off() 


gm_bio1_arcsin <- gm_file_all_arcsin$bio1.pval.arcsin
gm_nsat_arcsin <- gm_file_all_arcsin$near_surface_air_temperature.pval.arcsin
data <- as.data.frame(cbind(gm_bio1_arcsin, gm_nsat_arcsin))

pdf("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k_1/comparison/GM_pvals_arcsin.pdf",
    width=7,
    height=5)
hist(gm_bio1_arcsin,breaks=100, main = paste("P-value distribution of 50.000 SNPs for Europe"), sub="Worldclim: Bio1", xlab=paste("p value"))
hist(gm_nsat_arcsin,breaks=100, main = paste("P-value distribution of 50.000 SNPs for Europe"), sub="Copernicus: Near Surface Air Temperature", xlab=paste("p value"))
ggplot(data, aes(x = gm_bio1_arcsin, y = gm_nsat_arcsin)) +
  geom_point(alpha=0.5, size=0.5)  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle("Scatterplot NSAT vs. Bio1")
dev.off() 





####lea

lea_bio1 <- read.table("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k_1/results/fullgenome/LEA/bio1_adjusted_pvals.csv",h=T, sep=",")
lea_nsat <- read.table("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k_1/results/fullgenome/LEA/near_surface_air_temperature_adjusted_pvals.csv",h=T, sep=",")
lea_bio1_raw <- lea_bio1$adjusted.p.values
lea_nsat_raw <- lea_nsat$adjusted.p.values
data2 <- as.data.frame(cbind(lea_nsat_raw, lea_bio1_raw))

pdf("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k_1/comparison/LEA_pvals_raw.pdf",
    width=7,
    height=5)
hist(lea_bio1_raw,breaks=100, main = paste("P-value distribution (LFMM) "), sub="Worldclim: Bio1", xlab=paste("p value"))
hist(lea_nsat_raw,breaks=100, main = paste("P-value distribution (LFMM)"), sub="Copernicus: Near Surface Air Temperature", xlab=paste("p value"))
ggplot(data2, aes(x = lea_bio1_raw, y = lea_nsat_raw)) +
  geom_point(alpha=0.5, size=0.5)  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle("Scatterplot NSAT vs. Bio1")
dev.off()

data3 <- as.data.frame(cbind(data, data2))


pdf("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k_1/comparison/comparepvalsmethods.pdf",
    width=7,
    height=5)
ggplot(data3, aes(x = gm_bio1_arcsin, y = lea_bio1_raw)) +
  geom_point(alpha=0.5, size=0.5)  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle("Bio1 GM vs. LEA")
ggplot(data3, aes(x = gm_nsat_arcsin, y = lea_nsat_raw)) +
  geom_point(alpha=0.5, size=0.5)  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle("NSAT GM vs. LEA")
dev.off()


###trying to rank the values

rank1 <- rank(data$gm_bio1_arcsin)
rank2 <- rank(data$gm_nsat_arcsin)

# Create a rank plot
plot(rank1, gm_bio1_arcsin, type = "o", pch = 16, col = "blue",
     xlab = "Rank", ylab = "P-values",
     main = "Rank Plot of Two Sets of P-values")
points(rank2, gm_nsat_arcsin, type = "o", pch = 16, col = "red")
legend("topright", legend = c("Bio1", "NSAT"),
       col = c("blue", "red"), pch = 16, bty = "n")

##rank plot with ggplot

databio <- data.frame(Rank = rank1, Pvalue = gm_bio1_arcsin, Distribution = "Bio1")
datansat <- data.frame(Rank = rank2, Pvalue = gm_nsat_arcsin, Distribution = "NSAT")
combined_data <- rbind(databio, datansat)

# Create rank plot using ggplot
ggplot(combined_data, aes(x = Rank, y = Pvalue, color = Distribution)) +
  geom_line() +
  geom_point() +
  labs(x = "Rank", y = "P-values", title = "Rank Plot of Two Sets of P-values") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal()

##manhattan plots
library(ggplot2)
library(dplyr)
library(ggrepel)

gg.manhattan <- function(df, threshold, hlight, col, ylims, title, chr_spacing = 5) {
  df.tmp <- df %>%
    group_by(CHR) %>%
    summarise(chr_len = max(BP)) %>%
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    left_join(df, ., by = c("CHR" = "CHR")) %>%
    arrange(CHR, BP) %>%
    mutate(BPcum = BP + tot) %>%
    mutate(is_highlight = ifelse(SNP %in% hlight, "yes", "no")) %>%
    mutate(is_annotate = ifelse(P < threshold, "yes", "no"))
  
  axisdf <- df.tmp %>%
    group_by(CHR) %>%
    summarise(center = (max(BPcum) + min(BPcum)) / 2)
  
  ggplot(df.tmp, aes(x = BPcum, y = -log10(P))) +
    geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 2) +
    scale_color_manual(values = rep(col, 22)) +
    scale_x_continuous(
      breaks = axisdf$center,
      labels = axisdf$CHR,
      expand = c(0, chr_spacing)  # Adjust chr_spacing to control the spacing
    ) +
    scale_y_continuous(expand = c(0, 0), limits = ylims) +
    ggtitle(paste0(title)) +
    labs(x = "Chromosome") +
    geom_hline(yintercept = -log10(0.05/corr), color="red") +
    geom_hline(yintercept = -log10(0.01/corr), linetype = "dashed", color="blue") +
    #  geom_label_repel(data = df.tmp[df.tmp$is_annotate == "yes", ], aes(label = as.factor(SNP), alpha = 0.7), size = 5, force = 1.3) +
    theme_bw(base_size = 22) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}




###facet_grid()

neg_log_p_1<- as.numeric(var1)
neg_log_p_2<- as.numeric(var2)
p1 <- 10^(-neg_log_p_1)
p2 <- 10^(-neg_log_p_2)
#p1 <- neg_log_p_1[neg_log_p_1 < 1]
#p2 <- neg_log_p_2[neg_log_p_2 < 1]
folder_path=args[6]

pdf(paste0(folder_path, args[2], ".vs.",args[4], ".pdf"),
    width=7,
    height=5)
hist(p1,breaks=100, main = paste("Histogram of p-values for", args[4], "(BAYPASS)"), xlab=paste("p-values"))
hist(p2,breaks=100, main = paste("Histogram of p-values for", args[2], "(BAYPASS)"), xlab=paste("p-values"))

# Print a message indicating where the plot is saved
#pvals_bio1=10^(-as.numeric(bayfile$V10[bayfile$V1=="4"]))
#pvals_nsat=10^(-as.numeric(bayfile$V10[bayfile$V1=="5"]))

#qqplot(pvals_bio1, pvals_nsat, xlab = "bio1", ylab = "nsat", main = "BAYPASS p-values")
plot(p1, p2, main=paste("P-values BAYPASS"), xlab=paste(args[2]), ylab=paste(args[4]), pch=19, rgb=0,0,0,1)
abline(0,1)
#ggplot <- plot(p1,p2)
plotname=paste0(var1, "_",var2, "_qqplot.png")

lea_dir=args[7]
leafile1=read.table(paste(lea_dir, args[2],"_adjusted_pvals.csv", sep=""), sep=",")
leafile2=read.table(paste(lea_dir, args[4],"_adjusted_pvals.csv", sep=""), sep=",")
LEA_pval_1=as.numeric(leafile1$V4)
LEA_pval_2=as.numeric(leafile2$V4)

hist(LEA_pval_1,breaks=100, main = paste("Histogram of p-values for", args[4], "(LEA)"), xlab=paste("p-values"))
hist(LEA_pval_2,breaks=100, main = paste("Histogram of p-values for", args[2], "(LEA)"), xlab=paste("p-values"))

plot(LEA_pval_1, LEA_pval_2, main=paste("P-values LEA"), xlab=paste(args[2]), ylab=paste(args[4]))
abline(0,1)

## LINEAR 

##gm=read.table("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_9/results/fullgenome/GM/P_Values.csv", sep=",", header=TRUE)
###gm=read.table(args[?], sep=",")
##
##pv1 <- as.numeric(gm$near_surface_air_temperature.pval)
##pv2 <- as.numeric(gm$bio1.pval)
##
###v1 <- as.numeric(gm[[paste(args[2], ".pval", sep = "")]])
###v2 <- as.numeric(gm[[paste(args[4], ".pval", sep = "")]])
##
##hist(pv1,breaks=100, main = paste("Histogram of p-values for", args[2], "(LINEAR MODEL)"), xlab=paste("p-values"))
##hist(pv2,breaks=100, main = paste("Histogram of p-values for", args[4], "(LINEAR MODEL)"), xlab=paste("p-values"))



## MANHATTAN PLOT LEA

#LEA_file <- read.table("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_10/results/fullgenome/LEA/bio1_adjusted_pvals.csv", sep=",")
#affile <- read.table("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_10/polymorphPOS.af",h=T)
library(ggplot2)
library(dplyr)
library(ggrepel)


# Example function with increased spacing on the x-axis
gg.manhattan <- function(df, threshold, hlight, col, ylims, title, chr_spacing = 7) {
  df.tmp <- df %>%
    group_by(CHR) %>%
    summarise(chr_len = max(BP)) %>%
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    left_join(df, ., by = c("CHR" = "CHR")) %>%
    arrange(CHR, BP) %>%
    mutate(BPcum = BP + tot) %>%
    mutate(is_highlight = ifelse(SNP %in% hlight, "yes", "no")) %>%
    mutate(is_annotate = ifelse(P < threshold, "yes", "no"))
  
  axisdf <- df.tmp %>%
    group_by(CHR) %>%
    summarise(center = (max(BPcum) + min(BPcum)) / 2)
  
  ggplot(df.tmp, aes(x = BPcum, y = -log10(P))) +
    geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 2) +
    scale_color_manual(values = rep(col, 22)) +
    scale_x_continuous(
      breaks = axisdf$center,
      labels = axisdf$CHR,
      expand = c(0, chr_spacing)  # Adjust chr_spacing to control the spacing
    ) +
    scale_y_continuous(expand = c(0, 0), limits = ylims) +
    ggtitle(paste0(title)) +
    labs(x = "Chromosome") +
    geom_hline(yintercept = -log10(0.05/corr), color="red") +
    geom_hline(yintercept = -log10(0.01/corr), linetype = "dashed", color="blue") +
    #  geom_label_repel(data = df.tmp[df.tmp$is_annotate == "yes", ], aes(label = as.factor(SNP), alpha = 0.7), size = 5, force = 1.3) +
    theme_bw(base_size = 22) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}





##LEA MANHATTAN BIO1
affile <- read.table("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k_newSNPs/Subsampled3.af",h=T,sep="")
pod.thresh=quantile(na.exclude(lea_bio1_raw,probs=0.95))
td <-as.data.frame(cbind(lea_bio1_raw,affile$Chr, affile$Pos))
td$SNP <- paste(affile$Chr, affile$Pos, sep = "")
colnames(td) <- c("P", "CHR", "BP", "SNP")
td$P <- as.numeric(td$P)
td$BP <- as.numeric(td$BP)
td_nsat <-as.data.frame(cbind(lea_nsat_raw,affile$Chr, affile$Pos))
td_nsat$SNP <- paste(affile$Chr, affile$Pos, sep = "")
colnames(td_nsat) <- c("P", "CHR", "BP", "SNP")
td_nsat$P <- as.numeric(td_nsat$P)
td_nsat$BP <- as.numeric(td_nsat$BP)

mypalette <- c("#3c7b2e","#0d4b21","#318056","#6c7831","#2f6316") 
mypalette <- c("#51b44d",
                "#bca948",
                "#57a672",
                "#99b840",
                "#6c7831")


###GM MANHATTAN
td_gm_bio1 <-as.data.frame(cbind(gm_bio1_arcsin,affile$Chr, affile$Pos))
td_gm_bio1$SNP <- paste(affile$Chr, affile$Pos, sep = "")
colnames(td_gm_bio1) <- c("P", "CHR", "BP", "SNP")
td_gm_bio1$P <- as.numeric(td_gm_bio1$P)
td_gm_bio1$BP <- as.numeric(td_gm_bio1$BP)

td_gm_nsat <-as.data.frame(cbind(gm_nsat_arcsin,affile$Chr, affile$Pos))
td_gm_nsat$SNP <- paste(affile$Chr, affile$Pos, sep = "")
colnames(td_gm_nsat) <- c("P", "CHR", "BP", "SNP")
td_gm_nsat$P <- as.numeric(td_gm_nsat$P)
td_gm_nsat$BP <- as.numeric(td_gm_nsat$BP)

gg.manhattan(td, threshold=corr, hlight=mysnps, col=mypalette, ylims=c(0,10), title=paste("Manhattan-Plot for Bio1 (LEA)"))
gg.manhattan(td_nsat, threshold=corr, hlight=mysnps, col=mypalette, ylims=c(0,10), title=paste("Manhattan-Plot for NSAT (LEA)"))

gg.manhattan(td_gm_bio1, threshold=corr, hlight=mysnps, col=mypalette, ylims=c(0,10), title=paste("Manhattan-Plot for Bio1 (GM)"))
gg.manhattan(td_gm_nsat, threshold=corr, hlight=mysnps, col=mypalette, ylims=c(0,10), title=paste("Manhattan-Plot for NSAT (GM)"))


gg.manhattan(df = your_data,
             threshold = your_threshold,
             hlight = your_highlight,
             col = your_color,
             ylims = your_ylims,
             title = your_title,
             chr_order = c("chr1", "chr2", ..., "chrX", "chrY"))
pdf("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k_1/comparison/ManhattanPlots.pdf",
    width=7,
    height=5)
corr_td <- min(head(sort(-log10(td$P), decreasing = TRUE), 100))
gg.manhattan(td, threshold=corr_td, hlight=mysnps, col=mypalette, ylims=c(0,10), title=paste("Manhattan-Plot for Bio1 (LEA)"), chr_order = c("X", "2L","2R", "3L", "3R"))

gg.manhattan(td_nsat, threshold=corr, hlight=mysnps, col=mypalette, ylims=c(0,10), title=paste("Manhattan-Plot for NSAT (LEA)"))
gg.manhattan(td_gm_bio1, threshold=corr, hlight=mysnps, col=mypalette, ylims=c(0,10), title=paste("Manhattan-Plot for Bio1 (GM)"))
gg.manhattan(td_gm_nsat, threshold=corr, hlight=mysnps, col=mypalette, ylims=c(0,10), title=paste("Manhattan-Plot for NSAT (GM)"))
dev.off() 


###report figure B

td <- td[td$CHR != 4,]
td <- td[td$CHR != "Y",]
corr_td <- min(head(sort(-log10(td$P), decreasing = TRUE), 100))
td_nsat <- td_nsat[td_nsat$CHR != 4,]
td_nsat <- td_nsat[td_nsat$CHR != "Y",]
corr_td_nsat <- min(head(sort(-log10(td_nsat$P), decreasing = TRUE), 100))

pdf("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k_1/comparison/ReportFigureB.pdf",
    width=7,
    height=5)
plotX
plotY <- gg.manhattan(td, threshold=corr_td, hlight=mysnps, col=mypalette, ylims=c(0,10), title=paste("Manhattan-Plot for Bio1 (LEA)"), chr_order = c("X", "2L","2R", "3L", "3R"))
plotZ <- gg.manhattan(td_nsat, threshold=corr_td_nsat, hlight=mysnps, col=mypalette, ylims=c(0,10), title=paste("Manhattan-Plot for NSAT (LEA)"), chr_order = c("X", "2L","2R", "3L", "3R"))
dev.off() 

# Create the compound plot

# Assign positions to the plots
compound_plot <- compound_plot + plot_layout(
  ncol = 2, 
  heights = c(2, 1),  # The first plot spans two rows
  widths = c(1, 1)    # The left and right plots have the same width
)

# Print the compound plot
print(compound_plot)

####


pdf(paste0(folder_path, args[2], ".vs.",args[4], "_MANHATTAN.pdf"),
    width=15,
    height=7)
pod.thresh=quantile(na.exclude(LEA_pval_2,probs=0.95))
td2 <-as.data.frame(cbind(-log10(LEA_pval_2[0:-1]),affile$Chr, affile$Pos))

td2$SNP <- paste(affile$Chr, affile$Pos, sep = "")

colnames(td2) <- c("P", "CHR", "BP", "SNP")
td2$P <- as.numeric(td2$P)
td2$BP <- as.numeric(td2$BP)

#intersectLEA <- merge(td, td2, by = c("SNP"))
#mysnps <- intersectLEA$SNP
mysnps <- c("a","b")

gg.manhattan(td, threshold=corr, hlight=mysnps, col=mypalette, ylims=c(0,10), title=paste("Manhattan-Plot for ",args[2], "(LEA)"))
gg.manhattan(td2, threshold=corr, hlight=mysnps, col=mypalette, ylims=c(0,10), title=paste("Manhattan-Plot for ",args[4], "(LEA)"))


#######
#m <- ggplot(td, aes(x=V3, y=V1, color= V2)) +
#  geom_point(size=0.5) +
#  scale_colour_manual(values=c("#59c382",
#                               "#4aa5d4",
#                               "#41c45b",
#                               "#3abec8",
#                               "#155126",
#                               "#3abec8",
#                               "#155126")) +
#  geom_hline(yintercept=pod.thresh,col="blue",lty=2) +
#  labs(title="manhattan",
#       x="pos",
#       y="pval-log10",
#       color="Chromosme") + 
#  theme_minimal() + 
#  theme(legend.position = "none", 
#        axis.text.x = element_blank(),  
#        axis.text.y = element_blank(),  
#        axis.ticks.x = element_blank(), 
#        axis.ticks.y = element_blank()) 
#
#manh_facet <- m + facet_grid(.~V2, scales ="free_x", space="free_x")
#manh_facet


#BAYPASS MANHATTAN

B1 <- as.data.frame(cbind(-log(p1), affile$Chr, affile$Pos, c(seq(1:length(affile$Chr)))))
B2 <- as.data.frame(cbind(-log(p2), affile$Chr, affile$Pos, c(seq(1:length(affile$Chr)))))

colnames(B1) <- c("P", "CHR", "BP", "SNP")
colnames(B2) <- c("P", "CHR", "BP", "SNP")

B2$P <- as.numeric(B2$P)
B2$BP <- as.numeric(B2$BP)
B1$P <- as.numeric(B1$P)
B1$BP <- as.numeric(B1$BP)

#B1 <- B1[B1$P < 15.3,]
#B2 <- B2[B2$P < 15.3,]


###delete
plot(a,b,col=rgb(0,0,0,0.05),pch=16)
plot(c,d)
plot(a,c,col=rgb(0,0,0,0.01),pch=16)

###111
#1
plot(-log10(lea_bio1_raw), -log10(gm_bio1_arcsin),col=rgb(0,0,0,0.01),pch=16, xlim=c(0,20))
model <- lm(gm_bio1_arcsin ~ lea_bio1_raw)
rsquared <- summary(model)$r.squared

#1.1 DONE FOR NOW (R SQUARED MISSING)
ggplot(data, aes(x = -log10(lea_bio1_raw), y = -log10(gm_bio1_arcsin))) +
  geom_point(color = rgb(0, 0, 0, 0.1), shape = 16) +
  xlim(0, 10) + ggtitle(" Comparison of -log10 pvalues for Worldclim Bio1", ) +
  labs(x="Latent Factor Mixed Model", y="Linear Model") + geom_smooth(method = "lm", se = FALSE, linetype="dashed")+
  theme_bw() + # Apply black and white theme
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())
####122222
##1
plot(-log10(lea_nsat_raw), -log10(gm_nsat_arcsin),col=rgb(0,0,0,0.01),pch=16)
###1.1 NEEDS R SQUARED
ggplot(data, aes(x = -log10(lea_nsat_raw), y = -log10(gm_nsat_arcsin))) +
  geom_point(color = rgb(0, 0, 0, 0.1), shape = 16) +
  xlim(0, 10) + ggtitle(" Comparison of -log10 pvalues for Copernicus NSAT ", ) +
  labs(x="Latent Factor Mixed Model", y="Linear Model") + geom_smooth(method = "lm", se = FALSE, linetype="dashed")+
  theme_bw() + # Apply black and white theme
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())

model <- lm((-log10(lea_nsat_raw)) ~ (-log10(gm_nsat_arcsin)) )
rsquared <- summary(model)$r.squared
rsquared
####

#1
plot(-log10(lea_nsat_raw), -log10(lea_bio1_raw),col=rgb(0,0,0,0.1),pch=16, xlim=c(0,50))
abline(lm( -log10(lea_bio1_raw) ~ -log10(lea_nsat_raw)  ))
#1.1
ggplot(data, aes(x = -log10(lea_nsat_raw), y = -log10(lea_bio1_raw))) +
  geom_point(color = rgb(0, 0, 0, 0.1), shape = 16) +
  xlim(0, 10) + ggtitle(" Comparison of -log10 pvalues for Latent Factor Mixed Model") +
  labs(x="Copernicus NSAT", y="Worldclim Bio1") + geom_abline()

plotX <- ggplot(data, aes(x = -log10(lea_nsat_raw), y = -log10(lea_bio1_raw))) +
  geom_point(color = rgb(0, 0, 0, 0.1), shape = 16) +
  xlim(0, 10) + ggtitle(" Comparison of -log10 pvalues for Latent Factor Mixed Model", ) +
  labs(x="NSAT", y="Bio1") + geom_smooth(method = "lm", se = FALSE, linetype="dashed")+
  theme_bw() + # Apply black and white theme
  theme(
    plot.title = element_text(hjust = 0.5, size=18),
    axis.title.x = element_text(size = 18),  # Adjust x-axis label size
    axis.title.y = element_text(size = 18),  # Adjust y-axis label size
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )


##2xmh


####

ggplot(data, aes(x = -log10(lea_bio1_raw), y = -log10(gm_bio1_arcsin))) +
  geom_point(color = rgb(0, 0, 0, 0.01), shape = 16) +
  xlim(0, 20) + ggtitle("-log10 pvalues comparison for Bio1") + labs(x="Latent Factor Mixed Model ", y="Linear Model") +
  theme_bw() + geom_text(x = 10, y = 2, label = paste("R-squared =", round(rsquared, 3)), hjust = 0, size = 4) 
  

plot(gm_nsat_arcsin, gm_bio1_arcsin)

plot(c(1:length(e)),e,type="l",col="red")
lines(c(1:length(f)),f,col="green")


LEA_pvals <- read.csv("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k_newSNPs/results/fullgenome/LEA/bio1.rep.run1/bio1_LEA_pvals.csv")
