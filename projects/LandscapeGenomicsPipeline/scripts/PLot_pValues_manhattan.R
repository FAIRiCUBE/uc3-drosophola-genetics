# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Set the directory
setwd("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/results/fullgenome/GM")

# Specify the file name
file_name <- "/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/results/fullgenome/GM/P_Values.csv"

# Read the file
data <- read.csv(file_name)

# Create a directory to save plots if it doesn't exist
if (!dir.exists("plots")) {
  dir.create("plots")
}

plot_data <- data %>%
  select(DATA.Chr, DATA.Pos, ends_with(".pval")) %>%
  pivot_longer(cols = ends_with(".pval"), names_to = "pval_type", values_to = "pval_value")

Bonf=-log10(0.01/nrow(data))

# Loop through each unique p-value type
for (pval_col in unique(plot_data$pval_type)) {
  # Filter data for the current p-value type
  current_data <- data %>%
    select(DATA.Chr, DATA.Pos, !!sym(pval_col)) %>%
    rename(pval_value = !!sym(pval_col))
  
  # Generate a file name for the plot
  plot_file_name <- paste0("plots/", pval_col, "_manhattan_plot.png")
  
  # Create the plot
  p <- ggplot(current_data, aes(x = DATA.Pos, y = -log10(pval_value))) +
    geom_point(alpha = 0.5) +
    facet_wrap(~ DATA.Chr, scales = "free_x") +
    geom_hline(yintercept=Bonf,col="blue",lty=2) +
    geom_hline(yintercept=-log10(0.05),col="red",lty=2) +
    labs(title = paste("Manhattan Plot for", pval_col), x = "Position", y = "-log10(p-value)") +
    theme_bw()
  
  outliersBonf <- current_data %>% filter(pval_value < 10^(-Bonf))
  write.csv(outliersBonf, paste0("plots/LinearOutliers0.01_", pval_col, ".csv"), row.names = FALSE)

  # Save the plot
  ggsave(plot_file_name, plot = p)
}


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
