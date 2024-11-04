#install.packages('tidyverse')
library(tidyverse)
library(gridExtra)
library(ggplot2)
library(readr)
library(maps)
args <- commandArgs(TRUE)

#print(args[2])
#print(args[3])
#print(args[4])


## set working directory
setwd(args[1]) 
arm=args[4]
directory_path=paste("results/",arm, "/LinearRegressions", sep="")
dir.create(directory_path, recursive = TRUE)
## read AlleleFrequency Matrix
DATA=read.table(args[2],header=3)
#print(DATA)
AF=DATA[,3:ncol(DATA)]

rownames(AF) <- paste(DATA$Chr, DATA$Pos, sep=".")
envfile <- read.table(args[3], header=1, sep=",",  dec = ".")

samplesFreqs <- colnames(AF)
null_values <- read.csv("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/NUllValuesDetail.csv", header=FALSE)
##null_values in args

columns_with_only_na_string <- sapply(envfile, function(x) all(x == "na"))
na_string_columns <- names(envfile)[columns_with_only_na_string]
Rclean <- envfile[, colSums(envfile == "na") != nrow(envfile)]
Rclean <- as.data.frame(Rclean[2:ncol(Rclean)])
df_numeric <- as.data.frame(lapply(Rclean, as.numeric))
rownames(df_numeric) <- envfile$sample

# Nur Spaltennamen auswählen, die sowohl in null_values$V1 als auch in df_numeric vorhanden sind
valid_columns <- intersect(null_values$V1, colnames(df_numeric))

# Diese validen Spalten aus df_numeric auswählen
#new <- df_numeric[, valid_columns]
#sel <- intersect(rownames(df_numeric), samplesFreqs)

filtered_null_values <- null_values[null_values$V1 %in% colnames(df_numeric),]
named_list <- as.list(setNames(filtered_null_values$V2, filtered_null_values$V1))

null_values_list <- named_list

#df_new <- as.data.frame(lapply(seq_along(df_numeric), function(i) {
#  replace(new[[i]], new[[i]] == null_values_list[[names(new)[i]]], NA)
#}))

is_monomorphic <- function(column) {
  length(unique(column)) == 1}

# Apply the function to each column and get the names of monomorphic columns
monomorphic_columns <- sapply(df_numeric, is_monomorphic)
df_new <- df_numeric[, !monomorphic_columns]


colnames(df_new) <- colnames(df_numeric)
rownames(df_new) <- rownames(df_numeric)

#print(length(envfile$sample))
#print(length(samplesFreqs))
#rownames(df_new) <- envfile$sample.
#colnames(df_new) <- colnames(new)

core_europe <- df_new[!grepl("^UA|RU|BY|DK", rownames(df_new)), ]
na_count_col <- sapply(core_europe, function(x) sum(is.na(x)))
#absolutely no NAs
#cols_to_keep <- na_count_col <= 0
cols_to_keep <- na_count_col <= 10
remains <- core_europe[, cols_to_keep]
samplesFreqs <- gsub("\\.", "-", samplesFreqs)
ncol(remains)
nrow(remains)

sel <- intersect(rownames(remains), samplesFreqs)
#print(sel)
Env <- remains[sel,]

#print(length(rownames(remains)))

envfile <- Env
rownames(envfile) <- gsub("\\-", ".", rownames(envfile))


## read Metadata
#meta=read.csv("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k/data/metadata.csv",header=T)

actsamp<-intersect(rownames(envfile), colnames(AF))
AF <- AF[, actsamp, drop = FALSE]

Loci <- paste(AF$Chr, AF$Pos, sep=".")

annotations <- read.table("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/results/annotationdata/annotations.txt", header = FALSE)
annotLoci <- paste(annotations$V1, annotations$V2, sep=".")
colnames(annotations) <- c("Chr", "Pos", "Gene")
annotations <- unique(annotations)


#setwd(directory_path)
#print(colnames(AF))
#print(rownames(AF))
### make suuuuper simple function that fits a linear regression model and returns the p-values
lm.func <- function(x) {
  summary(lm(unlist(x)~y))[[4]][,4][2]
}


pvalcsv2 <- annotations
thres_env <- 0.05/length(DATA$Pos)
#print(colnames(envfile))
#pvalcsv <- data.frame(DATA$Chr, DATA$Pos)

facs2 <- c()
for ( i in colnames(envfile)){
  y=envfile[[i]]
  #print(y)
  facs <- c()
  pvalcsv <- annotations
  colnames(pvalcsv) <- c("Chr", "Pos", "Gene")
  if (is.numeric(y) && length(table(y)) > 2){
    #p.val<-apply(AF,1,lm.func)
    p.val.arcsin<-apply(asin(sqrt(AF**2)),1,lm.func)
    #ID=paste0(i,".pval")
    ID2=paste0(i,".pval.arcsin")
    #DATA[[ID]]<-p.val
    pvalcsv[[ID2]]<- p.val.arcsin 
    facs <- append(facs,ID2)
    DATA[[ID2]]<-p.val.arcsin
    pvalcsv2[[ID2]]<- p.val.arcsin
    facs2 <- append(facs2,ID2)
    #write.csv(pvalcsv, paste(directory_path,"/",i,"_pvalues_single.csv", sep=""), row.names= FALSE)
    write.csv(pvalcsv, paste("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullDataRun/results/AllSNPs/LinearRegressions/LinearRegr2/",i,"_pvalues.csv", sep=""), row.names=FALSE)
    #outliers <- na.exclude(pvalcsv[pvalcsv[[ID2]] < thres_env,])
    outliers <- na.omit(pvalcsv[pvalcsv[[ID2]] < thres_env,])
    #write.csv(outliers, paste(directory_path,"/",i,"_pvalues_outliers.csv", sep=""), row.names= FALSE)
    write.csv(outliers, paste("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullDataRun/results/AllSNPs/LinearRegressions/LinearRegr2/",i,"_pvalue_outliers.csv", sep=""), row.names=FALSE)
    #pdf(paste0("Histograms_P_Values_",i,".pdf"),
    #width=15,
    #height=5)
    #hist(as.numeric(p.val),breaks=100, main=paste0("Histogram of p-values", i), xlab="p-value")
    #hist(as.numeric(p.val.arcsin),breaks=100, main=paste0("Histogram of p-values", i), xlab="p-value")
    DATA <- as.data.frame(cbind(annotations, -log10(p.val.arcsin)))
    colnames(DATA) <- c("Chr", "Pos", "Gene", "pvalslog")
    Bonf <- -log10(thres_env)
    pl <- ggplot(DATA, aes(x = Pos, y = pvalslog)) +
        geom_point(aes(color = pvalslog > Bonf), alpha = 0.5) +  # Color points based on pvalslog > 10
        facet_wrap(~ Chr, scales = "free_x") +
        geom_hline(yintercept = Bonf, col = "blue", lty = 2) +
        geom_hline(yintercept = -log10(0.05), col = "red", lty = 2) +
        labs(title = paste("Manhattan Plot for", "pval"), x = "Position", y = "-log10(p-value)") +
        theme_bw() +
        geom_text(aes(label = ifelse(pvalslog > Bonf, Gene, '')), vjust = -0.5, check_overlap = TRUE) +  # Only label points where pvalslog > 10
        scale_color_manual(values = c("black", "red")) +
        theme(legend.position = "none")
    # Export as PNG
    png(filename = paste0(directory_path,"/",i, "_Pvalues.png"), width = 800, height = 600)
    #png(filename = paste0("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullDataRun/results/AllSNPs/LinearRegressions/LinearRegr2/",i, "_Pvalues.png"), width = 800, height = 600)
    print(pl)
    dev.off()
  } else {
    next
  }
}

write.csv(pvalcsv2, paste(directory_path,"/Merged_pvalues.arcsin.csv", sep=""), row.names= FALSE)

##export DATA to directory
##pvalpath=paste(directoy_path, "/P_Values.csv")
#write.csv(pvalcsv, "P_Values.csv")
#write.csv(pvalcsv2, "P_Values.arcsin.csv")
#
### multiple testing problem??
#y=runif(ncol(AF))
#D=data.frame("X.1"=runif(nrow(AF)))
#for (i in seq(2,ncol(AF),1)){
#  D[[paste0("X.",i)]]<-runif(nrow(AF))
#}
#
#Bonf=0.05/nrow(DATA)
#Test.p<-apply(D,1,lm.func)
#Test.p[Test.p<Bonf]
#
#
#pdf("Multtest_control.pdf",
#    width=15,
#    height=5)
#hist(Test.p,breaks=100, main="Histogram of p-values overall", xlab="p-value")
##dev.off()
#
#outliersBonf <- pvalcsv2[which(pvalcsv2$lat.pval.arcsin < Bonf),]
#OlBonf <- as.data.frame(paste0(outliersBonf$DATA.Chr,".", outliersBonf$DATA.Pos))
#write.csv(OlBonf, "LinearOutliers_ArcSin.csv")
#
### Boferroni-correctd p-value threshold
##Bonf=-log10(0.05/nrow(DATA))
#
###### plot with ggplot
##PLOT<-facs %>%
##  map(function(z){
##    PLOT.df<-data.frame(Pos = DATA$Pos/1000000, P.val = -log10(DATA[[z]]))
##    pl <- ggplot(PLOT.df, aes(x = Pos, y = P.val)) +
##      geom_point(alpha=0.3) +
##      xlab("Genomic Position on 3R [Mbp]") +
##      ylab("-log10(P.value)")+
##      geom_hline(yintercept=Bonf,
##                 col="blue",
##                 lty=2)+
##      geom_hline(yintercept=-log10(0.05),
##                 col="red",
##                 lty=2)+
##      ggtitle(z)+
##      theme_bw()
##    return(pl)
##  }) %>%
##  arrangeGrob(grobs = ., nrow=3) %>% #add ncol/nrow argument here
##  grid.arrange()
##
##ggsave("3R_Pvalues.png",
##       PLOT,
##       width=15,
##       height=7)
#
### Create a ggplot function for mapping
##ggplot_func <- function(PLOT.df, Bonf, title) {
##  ggplot(PLOT.df, aes(x = Pos, y = P.val)) +
##    geom_point(alpha = 0.3) +
##    xlab("Genomic Position [Mbp]") +
##    ylab("-log10(P.value)") +
##    geom_hline(yintercept = Bonf, col = "blue", lty = 2) +
##    geom_hline(yintercept = -log10(0.05), col = "red", lty = 2) +
##    ggtitle(title) +
##    theme_bw()
##}
#
### Apply ggplot function to each factor
##PLOT_list <- lapply(facs, function(z) {
##  PLOT.df <- data.frame(Pos = DATA$Pos / 1000000, P.val = -log10(DATA[[z]]))
##  ggplot_func(PLOT.df, Bonf, z)
##})
##
### Arrange and save the plots
##pdf("3R_Pvalues.pdf", width = 12, height = 7)
##grid.arrange(grobs = PLOT_list, nrow = 3)
##dev.off()
#
##for (i in seq_along(facs)) {
##  z <- facs[i]
##  PLOT.df <- data.frame(Pos = DATA$Pos / 100000, P.val = -log10(DATA[[z]]))
##  pl <- ggplot_func(PLOT.df, Bonf, z)
##  
##  # Export as PNG
##  png(filename = paste0(z, "_Pvalues.png"), width = 800, height = 600)
##  print(pl)
##  dev.off()
##}
#






