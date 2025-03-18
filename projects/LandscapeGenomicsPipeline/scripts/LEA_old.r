##This scripts purpose will be applying LEA LFMM to a VCF and parallelize workflow via manual z-score pulling

args <- commandArgs(TRUE)
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args){
  cat("Example: -working directory ")}

library(LEA)
library(readr)
library(gtools)
library(tidyverse)

rep=args[5]

repdir=paste0(args[1],"/",args[4])
dir.create(repdir)
setwd(repdir)

metadata <- read.csv(args[3])

is_monomorphic <- function(row) {
  return(length(unique(row)) == 1)
}

# Function to identify and remove monomorphic rows
remove_monomorphic_rows <- function(genotype_data) {
  monomorphic_rows <- apply(genotype_data, 1, is_monomorphic)
  polymorphic_data <- genotype_data[!monomorphic_rows, ]
  return(polymorphic_data)
}

DATA_full <- read.table(args[2], h=T, sep="\t")
rownames(DATA_full) <- paste0(DATA_full$Chr, DATA_full$Pos)
SUB <- DATA_full[3:length(DATA_full)]
###extra filter 
freq_mean <- rowMeans(SUB)
SUB <- SUB[-which(freq_mean>=0.95 | freq_mean<=0.05),]
####
samples <- gsub("\\.", "-", colnames(DATA_full))
colnames(SUB) <- gsub("\\.", "-", colnames(SUB))
ss_data <- metadata[metadata$sample %in% samples,]
s_data <- na.exclude(ss_data)
SUB2 <- SUB[colnames(SUB) %in% s_data$sample]
SUB2 <- asin(sqrt(SUB2))
SUB2 <- remove_monomorphic_rows(SUB2)


s_data <- s_data[s_data$sample %in% colnames(SUB2), ]
DATA_full <- DATA_full[rownames(DATA_full) %in% rownames(SUB2), ]



varname=paste(args[4])
#print(colnames(s_data))
#varname <- gsub("\\+", ".", varname_or)

print(varname)

if (!(varname %in% colnames(s_data))) {
    stop(paste("Column", varname, "not found in s_data"))
}



print(paste0("After removing monomorphic rows,", nrow(SUB2), "SNPs  are being analysed"))


geno_object <- paste0(repdir,"/genotypes.lfmm")

if (file.exists(geno_object)) {
  # If the file exists, read it
  tr <- read.lfmm(geno_object)
  print("File exists and has been opened.")
} else {
  tr <- t(as.matrix(SUB2))
  #tr_asin <- asin(sqrt(tr))
  tr <- round(tr,2)
  print("DIMENSIONS OF GENO-LFMM")
  print(nrow(tr))
  print(ncol(tr))
  # Write the matrix to an LFMM file
  write.lfmm(tr,geno_object)
  print("A new geno.lffm.file has been created.")
}

print(dim(tr))
first_line <- readLines(geno_object, n = 1)
print(length(strsplit(first_line, " ")[[1]]))

var_index <- match(varname, colnames(s_data))
print(varname)
f <- as.matrix(s_data[,varname])  
write.env(f, paste0(repdir, "/gradients.env"))

run_lfmm <- function(data, env, K_range = 1:10) {
    project <- snmf(data, K = K_range, entropy = TRUE, repetitions = 2, project = "new")
    CE <- sapply(K_range, function(K) mean(cross.entropy(project, K = K)))
    K <- which.min(CE)
    mod.lfmm2 <- lfmm2(data, env, K = K)
    return(mod.lfmm2)
}

print('CALCULATING LFMM2')
mod.lfmm2 <- run_lfmm(geno_object, f)

print('TESTING THE LFMM2 OBJECT')           
pv <- lfmm2.test(object = mod.lfmm2, input = tr, env = f, linear = TRUE)


print('EXPORTING P-VALUES TO CSV FILE')     
markerpos <- data.frame(
  Chr = DATA_full$Chr,
  Pos = DATA_full$Pos,
  Pval = as.numeric(pv$pvalues)  # Make sure the p-values are numeric
)
                 
#colnames(markerpos) <- c("Chr","Pos", "Pval")
                 
write.csv(markerpos, paste(repdir,"/",varname,"_LEA_pvals.csv", sep=""), row.names=FALSE)

Bonf=0.05/(nrow(markerpos)+1)
#markerpos <- as.data.frame(markerpos)
outliers <- markerpos[markerpos$Pval < Bonf,]
write.csv(outliers, paste(repdir,"/",varname,"_LEA_pvals_outliers.csv", sep=""), row.names=FALSE)


print("LOADING ANNOTATIONS")

annotations <- read.table("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/results/annotationdata/annotations.txt", header = FALSE)
annotLoci <- paste(annotations$V1, annotations$V2, sep=".")
colnames(annotations) <- c("Chr", "Pos", "Gene")
u_annotations <- unique(annotations)
DATA <- merge(u_annotations, markerpos, by= c("Chr","Pos"))
                 
print("PLOT MANHATTAN")

pl <- ggplot(DATA, aes(x = Pos, y = -log10(Pval))) +
    geom_point(aes(color = -log10(Pval) > -log10(Bonf)), alpha = 0.5) +  # Color points based on pvalslog > 10
    facet_wrap(~ Chr, scales = "free_x") +
    geom_hline(yintercept = -log10(Bonf), col = "blue", lty = 2) +
    geom_hline(yintercept = -log10(0.05), col = "red", lty = 2) +
    labs(title = paste("Manhattan Plot for", "pval"), x = "Position", y = "-log10(p-value)") +
    theme_bw() +
    geom_text(aes(label = ifelse(-log10(Pval) > -log10(Bonf), Gene, '')), vjust = -0.5, check_overlap = TRUE) +  # Only label points where pvalslog > 10
    scale_color_manual(values = c("black", "red")) +
    theme(legend.position = "none")
    
png(filename = paste0(repdir,"/",varname, "_Pvalues.png"), width = 800, height = 600)
print(pl)
dev.off()