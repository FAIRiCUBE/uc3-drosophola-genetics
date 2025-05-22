###Compare Results from Linear Model and LEA
library(ggplot2)
library(gridExtra)

#1.EUROPE

##GM arcsin transformed
linear_pvalues <- read.csv("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k_1/results/fullgenome/GM/P_Values.arcsin.csv")
##GM  Bio1
gm_bio1 <- linear_pvalues$bio1.pval.arcsin
##GM NSAT 
gm_nsat <- linear_pvalues$near_surface_air_temperature.pval.arcsin

###LEA_BIO1

lea_bio1_full <- read.csv("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k_newSNPs/results/fullgenome/LEA/bio1.rep.run1/bio1_LEA_pvals.csv")
lea_bio1 <- lea_bio1_full$pv.pvalues

###LEA_NSAT
lea_nsat_full <- read.csv("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k_newSNPs/results/fullgenome/LEA/near_surface_air_temperature.rep.run1/near_surface_air_temperature_LEA_pvals.csv")
lea_nsat <- lea_nsat_full$pv.pvalues


data <- as.data.frame(cbind(gm_bio1, gm_nsat, lea_bio1, lea_nsat))

### BIO1 LEA VS GM
plot1 <- ggplot(data, aes(x = -log10(lea_bio1), y = -log10(gm_bio1))) +
  geom_point(color = rgb(0, 0, 0, 0.1), shape = 16) +
  ggtitle(" Comparison of -log 10 pvalues for Worldclim Bio1", ) +
  labs(x="Latent Factor Mixed Model", y="Linear Model") + geom_smooth(method = "lm", se = FALSE, linetype="dashed")+
  theme_bw() + # Apply black and white theme
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())

## BIO1 P VAL DISTRIBUTION LEA
plot2 <- ggplot(data, aes(x = lea_bio1)) +
  geom_histogram(breaks = seq(0, 1, by = 0.02), fill = "grey", color = "black") +
  labs(x = "p values", y = "Frequency", title = "P Value Distribution for Bio1 (LEA)") +
  theme_bw() 

## BIO1 P VAL DISTRIBUTION GM
plot3 <- ggplot(data, aes(x = gm_bio1)) +
  geom_histogram(breaks = seq(0, 1, by = 0.02), fill = "grey", color = "black") +
  labs(x = "p values", y = "Frequency", title = "P Value Distribution for Bio1 (Linear Model)") +
  theme_bw() 


### GM 
plot4 <- ggplot(data, aes(x = -log10(lea_nsat), y = -log10(gm_nsat))) +
  geom_point(color = rgb(0, 0, 0, 0.1), shape = 16) +
  ggtitle(" Comparison of -log10 pvalues for Copernicus NSAT", ) +
  labs(x="Latent Factor Mixed Model", y="Linear Model") + geom_smooth(method = "lm", se = FALSE, linetype="dashed")+
  theme_bw() + # Apply black and white theme
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())

plot5 <- ggplot(data, aes(x = lea_nsat)) +
  geom_histogram(breaks = seq(0, 1, by = 0.02), fill = "grey", color = "black") +
  labs(x = "p values", y = "Frequency", title = "P Value Distribution for NSAT (LEA)") +
  theme_bw() 

plot6 <- ggplot(data, aes(x = gm_nsat)) +
  geom_histogram(breaks = seq(0, 1, by = 0.02), fill = "grey", color = "black") +
  labs(x = "p values", y = "Frequency", title = "P Value Distribution for NSAT (Linear model)") +
  theme_bw() 


compound_plot <- grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, ncol = 3)


pdf("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k_1/comparison/Report_Europe_PlotsAlog.pdf",
    width=15,
    height=8)
compound_plot <- grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, ncol = 3)
dev.off()


###Manhattan PLots

# Example function with increased spacing on the x-axis
gg.manhattan <- function(df, threshold, hlight, col, ylims, title, chr_order, chr_spacing = 7) {
  # Convert CHR column to factor with custom levels
  df$CHR <- factor(df$CHR, levels = chr_order)
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
  
  # Calculate y-axis limits with a buffer
  ylim_buffer <- 3  # Adjust this buffer as needed
  max_y_value <- max(-log10(df.tmp$P)) + ylim_buffer
  ylims <- c(0, max_y_value)
  
  ggplot(df.tmp, aes(x = BPcum, y = -log10(P))) +
    geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 2) +
    scale_color_manual(values = rep(col, length(chr_order))) +
    scale_x_continuous(
      breaks = axisdf$center,
      labels = chr_order,
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


-log10()
plotX <- ggplot(data, aes(x = -log10(lea_nsat), y = -log10(lea_bio1))) +
  geom_point(color = rgb(0, 0, 0, 0.1), shape = 16) +
  ggtitle(" Comparison of -log10 pvalues for Latent Factor Mixed Model", ) +
  labs(x="NSAT", y="Bio1") + geom_smooth(method = "lm", se = FALSE, linetype="dashed")+
  theme_bw() + # Apply black and white theme
  theme(
    plot.title = element_text(hjust = 0.5, size=18),
    axis.title.x = element_text(size = 18),  # Adjust x-axis label size
    axis.title.y = element_text(size = 18),  # Adjust y-axis label size
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )






man_gm_bio1 <- as.data.frame(cbind(linear_pvalues$DATA.Chr, as.numeric(linear_pvalues$DATA.Pos), as.numeric(linear_pvalues$bio1.pval.arcsin)))
man_gm_bio1$SNP <- paste(man_gm_bio1$V1, man_gm_bio1$V2, sep = "")
colnames(man_gm_bio1) <- c("CHR", "BP", "P", "SNP")
man_gm_bio1 <- man_gm_bio1[!(man_gm_bio1$CHR %in% c('4', 'Y')), ]
as.numeric(man_gm_bio1$BP) -> man_gm_bio1$BP
as.numeric(man_gm_bio1$P) -> man_gm_bio1$P


man_gm_nsat <- as.data.frame(cbind(linear_pvalues$DATA.Chr, linear_pvalues$DATA.Pos, linear_pvalues$near_surface_air_temperature.pval.arcsin))
man_gm_nsat$SNP <- paste(man_gm_nsat$V1, man_gm_nsat$V2, sep = "")
colnames(man_gm_nsat) <- c("CHR", "BP", "P", "SNP")
man_gm_nsat <- man_gm_nsat[!(man_gm_nsat$CHR %in% c('4', 'Y')), ]
as.numeric(man_gm_nsat$BP) -> man_gm_nsat$BP
as.numeric(man_gm_nsat$P) -> man_gm_nsat$P

man_lea_bio1 <- as.data.frame(cbind(lea_bio1_full$V1, lea_bio1_full$V2, lea_bio1_full$pv.pvalues))
man_lea_bio1 $SNP <- paste(man_lea_bio1$V1, man_lea_bio1$V2, sep = "")
colnames(man_lea_bio1) <- c("CHR", "BP", "P", "SNP")
man_lea_bio1 <- man_lea_bio1[!(man_lea_bio1$CHR %in% c('4', 'Y')), ]
as.numeric(man_lea_bio1$BP) -> man_lea_bio1$BP
as.numeric(man_lea_bio1$P) -> man_lea_bio1$P

man_lea_nsat  <- as.data.frame(cbind(lea_nsat_full$V1, lea_nsat_full$V2, lea_nsat_full$pv.pvalues))
man_lea_nsat$SNP <- paste(lea_nsat_full$V1, lea_nsat_full$V2, sep = "")
colnames(man_lea_nsat) <- c("CHR", "BP", "P", "SNP")
man_lea_nsat <- man_lea_nsat[!(man_lea_nsat$CHR %in% c('4', 'Y')), ]
as.numeric(man_lea_nsat$BP) -> man_lea_nsat$BP
as.numeric(man_lea_nsat$P) -> man_lea_nsat$P

corr = nrow(man_gm_bio1)
pdf("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k_1/comparison/Report_Europe_FigureB.pdf",
    width=12,
    height=5)
plotX
plotY <- gg.manhattan(man_lea_bio1, threshold=corr, hlight=mysnps, col=mypalette, ylims=c(0,10), title=paste("Manhattan-Plot for Bio1 (LEA)"), chr_order = c("X", "2L","2R", "3L", "3R"))
plotZ <- gg.manhattan(man_lea_nsat, threshold=corr, hlight=mysnps, col=mypalette, ylims=c(0,10), title=paste("Manhattan-Plot for NSAT (LEA)"), chr_order = c("X", "2L","2R", "3L", "3R"))
plotY
plotZ
plotZ1 <- gg.manhattan(man_gm_bio1, threshold=corr, hlight=mysnps, col=mypalette, ylims=c(0,10), title=paste("Manhattan-Plot for Bio1 (Linear Model)"), chr_order = c("X", "2L","2R", "3L", "3R"))
plotZ2 <- gg.manhattan(man_gm_nsat, threshold=corr, hlight=mysnps, col=mypalette, ylims=c(0,10), title=paste("Manhattan-Plot for NSAT (Linear Model)"), chr_order = c("X", "2L","2R", "3L", "3R"))
plotZ1
plotZ2
dev.off() 

compound_plot <- ggarrange(
  plotX,
  ggarrange(plotY, plotZ, ncol = 1),
  nrow = 1
)
ggsave("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k_1/comparison/Report_Europe_FigureB1_1.pdf", compound_plot, width=15)


### 2. North AMERICA
na_linear <- read.csv("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_NorthAmerica50k/results/fullgenome/GM/P_Values.arcsin.csv")
na_gm_bio1 <- na_linear$bio1.pval.arcsin
na_gm_nsat <- na_linear$near_surface_air_temperature.pval.arcsin

na_lea_bio1_full <- read.csv("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_NorthAmerica50k/results/fullgenome/LEA/bio1.rep.run1/bio1_LEA_pvals.csv")
na_lea_bio1 <- na_lea_bio1_full$pv.pvalues

###LEA_NSAT
na_lea_nsat_full <- read.csv("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_NorthAmerica50k/results/fullgenome/LEA/near_surface_air_temperature.rep.run1/near_surface_air_temperature_LEA_pvals.csv")
na_lea_nsat <- na_lea_nsat_full$pv.pvalues


plot1.1 <- ggplot(data, aes(x = na_lea_bio1, y = na_gm_bio1)) +
  geom_point(color = rgb(0, 0, 0, 0.1), shape = 16) +
  ggtitle(" Comparison of pvalues for Worldclim Bio1", ) +
  labs(x="Latent Factor Mixed Model", y="Linear Model") + geom_smooth(method = "lm", se = FALSE, linetype="dashed")+
  theme_bw() + # Apply black and white theme
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())

## BIO1 P VAL DISTRIBUTION LEA
plot2.1 <- ggplot(data, aes(x = na_lea_bio1)) +
  geom_histogram(breaks = seq(0, 1, by = 0.02), fill = "grey", color = "black") +
  labs(x = "p values", y = "Frequency", title = "P Value Distribution for Bio1 (LEA)") +
  theme_bw() 

## BIO1 P VAL DISTRIBUTION GM
plot3.1 <- ggplot(data, aes(x = na_gm_bio1)) +
  geom_histogram(breaks = seq(0, 1, by = 0.02), fill = "grey", color = "black") +
  labs(x = "p values", y = "Frequency", title = "P Value Distribution for Bio1 (Linear Model)") +
  theme_bw() 


### GM 
plot4.1 <- ggplot(data, aes(x = na_lea_nsat, y = na_gm_nsat)) +
  geom_point(color = rgb(0, 0, 0, 0.1), shape = 16) +
  ggtitle(" Comparison of pvalues for Copernicus NSAT", ) +
  labs(x="Latent Factor Mixed Model", y="Linear Model") + geom_smooth(method = "lm", se = FALSE, linetype="dashed")+
  theme_bw() + # Apply black and white theme
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())

plot5.1 <- ggplot(data, aes(x = na_lea_nsat)) +
  geom_histogram(breaks = seq(0, 1, by = 0.02), fill = "grey", color = "black") +
  labs(x = "p values", y = "Frequency", title = "P Value Distribution for NSAT (LEA)") +
  theme_bw() 

plot6.1 <- ggplot(data, aes(x = na_gm_nsat)) +
  geom_histogram(breaks = seq(0, 1, by = 0.02), fill = "grey", color = "black") +
  labs(x = "p values", y = "Frequency", title = "P Value Distribution for NSAT (Linear model)") +
  theme_bw() 


data <- as.data.frame(cbind(na_gm_bio1, na_gm_nsat, na_lea_bio1, na_lea_nsat))


compound_plot <- grid.arrange(plot1.1, plot2.1, plot3.1, plot4.1, plot5.1, plot6.1, ncol = 3)

pdf("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k_1/comparison/Report_NorthAmerica_PlotsA.pdf",
    width=15,
    height=8)
compound_plot <- grid.arrange(plot1.1, plot2.1, plot3.1, plot4.1, plot5.1, plot6.1, ncol = 3)
dev.off()


man_gm_bio1 <- as.data.frame(cbind(na_linear$DATA.Chr, as.numeric(na_linear$DATA.Pos), as.numeric(na_linear$bio1.pval.arcsin)))
man_gm_bio1$SNP <- paste(man_gm_bio1$V1, man_gm_bio1$V2, sep = "")
colnames(man_gm_bio1) <- c("CHR", "BP", "P", "SNP")
man_gm_bio1 <- man_gm_bio1[!(man_gm_bio1$CHR %in% c('4', 'Y')), ]
as.numeric(man_gm_bio1$BP) -> man_gm_bio1$BP
as.numeric(man_gm_bio1$P) -> man_gm_bio1$P


man_gm_nsat <- as.data.frame(cbind(na_linear$DATA.Chr, na_linear$DATA.Pos, na_linear$near_surface_air_temperature.pval.arcsin))
man_gm_nsat$SNP <- paste(man_gm_nsat$V1, man_gm_nsat$V2, sep = "")
colnames(man_gm_nsat) <- c("CHR", "BP", "P", "SNP")
man_gm_nsat <- man_gm_nsat[!(man_gm_nsat$CHR %in% c('4', 'Y')), ]
as.numeric(man_gm_nsat$BP) -> man_gm_nsat$BP
as.numeric(man_gm_nsat$P) -> man_gm_nsat$P

man_lea_bio1 <- as.data.frame(cbind(na_lea_bio1_full$V1, na_lea_bio1_full$V2, na_lea_bio1_full$pv.pvalues))
man_lea_bio1 $SNP <- paste(man_lea_bio1$V1, man_lea_bio1$V2, sep = "")
colnames(man_lea_bio1) <- c("CHR", "BP", "P", "SNP")
man_lea_bio1 <- man_lea_bio1[!(man_lea_bio1$CHR %in% c('4', 'Y')), ]
as.numeric(man_lea_bio1$BP) -> man_lea_bio1$BP
as.numeric(man_lea_bio1$P) -> man_lea_bio1$P

man_lea_nsat  <- as.data.frame(cbind(na_lea_nsat_full$V1, na_lea_nsat_full$V2, na_lea_nsat_full$pv.pvalues))
man_lea_nsat$SNP <- paste(na_lea_nsat_full$V1, na_lea_nsat_full$V2, sep = "")
colnames(man_lea_nsat) <- c("CHR", "BP", "P", "SNP")
man_lea_nsat <- man_lea_nsat[!(man_lea_nsat$CHR %in% c('4', 'Y')), ]
as.numeric(man_lea_nsat$BP) -> man_lea_nsat$BP
as.numeric(man_lea_nsat$P) -> man_lea_nsat$P

corr = nrow(man_gm_bio1)

pdf("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k_1/comparison/Report_NorthAmerica_FigureB.pdf",
    width=12,
    height=5)
plotX
plotY <- gg.manhattan(man_lea_bio1, threshold=corr, hlight=mysnps, col=mypalette, ylims=c(0,10), title=paste("Manhattan-Plot for Bio1 (LEA)"), chr_order = c("X", "2L","2R", "3L", "3R"))
plotZ <- gg.manhattan(man_lea_nsat, threshold=corr, hlight=mysnps, col=mypalette, ylims=c(0,10), title=paste("Manhattan-Plot for NSAT (LEA)"), chr_order = c("X", "2L","2R", "3L", "3R"))
plotY
plotZ
plotZ1 <- gg.manhattan(man_gm_bio1, threshold=corr, hlight=mysnps, col=mypalette, ylims=c(0,10), title=paste("Manhattan-Plot for Bio1 (Linear Model)"), chr_order = c("X", "2L","2R", "3L", "3R"))
plotZ2 <- gg.manhattan(man_gm_nsat, threshold=corr, hlight=mysnps, col=mypalette, ylims=c(0,10), title=paste("Manhattan-Plot for NSAT (Linear Model)"), chr_order = c("X", "2L","2R", "3L", "3R"))
plotZ1
plotZ2
dev.off() 

compound_plot <- ggarrange(
  plotX,
  ggarrange(plotY, plotZ, ncol = 1),
  nrow = 1
)
ggsave("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k_1/comparison/Report_NorthAmerica_FigureB1.pdf", compound_plot, width=15)

