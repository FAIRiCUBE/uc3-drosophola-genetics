#install.packages('tidyverse')
options(repos = "https://cran.r-project.org/")
#install.packages('gridExtra')
library(tidyverse)
library(gridExtra)
library(ggplot2)
library(readr)
library(maps)
args <- commandArgs(TRUE)

## set working directory
#print(args[1])
#print(args[4])
setwd(args[1])
wd=args[1]
arm=args[4]
outpath=args[5]
meta <- read.table(args[3], header=1, sep=",",  dec = ".")
#

merged_df <- data.frame(Chr = character(0), Pos = numeric(0), Response = character(0))
csvfiles <- list()
plot_list <- list()



for (i in colnames(meta)) {
    print(colnames(meta))
    directory_path <- file.path(wd, "results", arm, "LEA", paste0(i, "rep1"))
    csv_file <- list.files(path = directory_path, pattern = "pvals.csv$", full.names = TRUE)
    csvfiles[[i]] <- csv_file
}

print(csvfiles)

# Loop through the list of file paths
for (i in names(csvfiles)) {
  csv_file <- csvfiles[[i]]
  if (!is.null(csv_file) && length(csv_file) > 0 && file.exists(csv_file)) {
    print(csv_file)
    pvalcsv_sub <- read.csv(csv_file)
    colnames(pvalcsv_sub) <- c("Response", "Chr", "Pos", i)
    pvalcsv_sub$Chr <- as.factor(pvalcsv_sub$Chr)
    print(colnames(merged_df))
    print(colnames(pvalcsv_sub))
    #print(pvalcsv_sub$Pvalue)
    #merged_df[[i]] <- pvalcsv_sub[]
    #merged_df <- cbind(pvalcsv_sub$Chr, pvalcsv_sub$Pos, pvalcsv_sub$i)
    merged_df <- merge(merged_df, pvalcsv_sub, by = c("Response","Chr", "Pos"), all = TRUE)
    # Process the dataframe as needed
    pl <- ggplot(merged_df, aes(x = Pos, y = -log10(.data[[i]]))) +
        geom_point(col = rgb(0, 0, 0, 0.1), pch = 16) +
        facet_grid(. ~ Chr, scales = "free_x", space = "free") +
        geom_hline(yintercept = -log10(0.05 / (nrow(merged_df) + 1)), colour = "blue") +
        xlab("Genomic Position [Mbp]") +
        ylab(paste0("-log10(P.value): ", i)) +
        theme_bw()
        ggtitle(paste0(i)) +
        theme(legend.position = "none")
    ggsave(paste0("LEA_Pvalues_",i,"_Manhattan.png"), pl, width = 12, height = 7)
    plot_list[[i]] <- pl
    
  } else {
    warning(sprintf("File for column '%s' does not exist or is invalid.", colnames(meta)[i]))
  }
}

#print(merged_df)

print(plot_list)

# Arrange the plots in a grid
PLOT_grid <- do.call(grid.arrange, c(plot_list, ncol = 1)) # Adjust nrow and ncol as needed

output=paste0(wd, "/results/", arm, "/LEA/LEA_Pvalues.png")

merge=paste0(wd, "/results/", arm, "/LEA/Merged_Pvalues.csv")
print(merged_df)

write.csv(merged_df, file = merge, row.names = FALSE)

# Save the combined plot to a file
ggsave(output, PLOT_grid, width = 15, height = 7)
ggsave(paste0(outpath,"/LEA_Pvalues.png"), PLOT_grid, width = 15, height = 7)



### read AlleleFrequency Matrix
#DATA=read.table(args[2],
#  header=3)
#AF=DATA[,3:ncol(DATA)]
#meta <- read.table(args[3], header=1, sep=",",  dec = ".")
#pvalcsv <- data.frame(DATA$Chr, DATA$Pos)
#pvalcsv2 <- data.frame(DATA$Chr, DATA$Pos)
### read Metadata
##meta=read.csv("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k/data/metadata.csv",header=T)
##DATA=read.table("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k/data/Subsampled_Europe_50k.recode.vcf.gz.af",header=3)
#cleaned_row_names <- gsub("\\.", "-", colnames(AF))
#cleaned_row_names_DATA <- gsub("\\.", "-", colnames(DATA))
#
#colnames(AF) <- cleaned_row_names
#colnames(DATA) <- cleaned_row_names_DATA
#
#
#actsamp<-intersect(cleaned_row_names, meta$sampleId)
#actsamp2<-intersect(cleaned_row_names_DATA, meta$sampleId)
#actsamp2 <- c("Chr", "Pos", actsamp2)
#                    
#
##meta <- read_csv(file = commandArgs(trailingOnly = TRUE)[3])
##meta <- read.table("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/testNSAT/data/NSAT.csv", header=1, sep=",",  dec = ".")
#
#
#
#AF <- AF[, actsamp, drop = FALSE]
#DATA <- DATA[, actsamp2, drop = FALSE]
##AF2 <- AF[meta$sampleId %in% unlist(colnames(AF)),]
#meta <- meta[meta$sampleId %in% colnames(AF),]
#
#setwd(directory_path)
#
### make suuuuper simple function that fits a linear regression model and returns the p-values
#lm.func <- function(x) {
#  summary(lm(unlist(x)~y))[[4]][,4][2]
#}
#
#print(colnames(meta))
##pvalcsv <- data.frame(DATA$Chr, DATA$Pos)
#facs <- c()
#facs2 <- c()
#for ( i in colnames(meta)){
#  y=meta[[i]]
#  print(y)
#  if (is.numeric(y) && length(table(y)) > 2){
#    p.val<-apply(AF,1,lm.func)
#    p.val.arcsin<-apply(asin(sqrt(AF**2)),1,lm.func)
#    ID=paste0(i,".pval")
#    ID2=paste0(i,".pval.arcsin")
#    DATA[[ID]]<-p.val
#    pvalcsv[[ID]]<- p.val 
#    facs <- append(facs,ID)
#    DATA[[ID2]]<-p.val.arcsin
#    pvalcsv2[[ID2]]<- p.val.arcsin
#    facs2 <- append(facs2,ID2)
#    write.csv(pvalcsv, paste(i,"_pvalues.csv", sep=""))
#    write.csv(pvalcsv2, paste(i,"_pvalues.arcsin.csv", sep=""))
#    pdf("Histograms_P_Values.pdf",
#    width=15,
#    height=5)
#    hist(as.numeric(p.val),breaks=100, main=paste0("Histogram of p-values", i), xlab="p-value")
#    hist(as.numeric(p.val.arcsin),breaks=100, main=paste0("Histogram of p-values", i), xlab="p-value")
#  } else {
#    next
#  }
#}
#
#
## Function to create Manhattan plot
#create_manhattan_plot <- function(data, pval_column, bonf_threshold, plot_filename, facet_by_chr = FALSE) {
#  data <- data %>% mutate(log_pval = -log10(!!sym(pval_column)))
#  p <- ggplot(data, aes(x = Pos / 1000000, y = log_pval)) +
#    geom_point(alpha = 0.3) +
#    xlab("Genomic Position [Mbp]") +
#    ylab("-log10(P.value)") +
#    geom_hline(yintercept = bonf_threshold, col = "blue", lty = 2) +
#    geom_hline(yintercept = -log10(0.05), col = "red", lty = 2) +
#    ggtitle(pval_column) +
#    theme_bw()
#  
#  if (facet_by_chr) {
#    p <- p + facet_wrap(~ DATA.Chr)
#  }
#  
#  ggsave(plot_filename, plot = p, width = 15, height = 7)
#}
#
## Bonferroni-corrected p-value threshold
#bonf_threshold <- -log10(0.05 / nrow(DATA))
#
## Identify p-value columns
#pval_columns <- colnames(pvalcsv)[grepl("\\.pval$", colnames(DATA))]
#
## Create and save Manhattan plots for each p-value column
#for (pval_column in pval_columns) {
#  plot_filename <- paste0(pval_column, "_manhattan.png")
#  create_manhattan_plot(DATA, pval_column, bonf_threshold, plot_filename, facet_by_chr = TRUE)
#}
