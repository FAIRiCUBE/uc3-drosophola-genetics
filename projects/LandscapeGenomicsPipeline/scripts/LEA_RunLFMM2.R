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
print(args[1])
print(args[2])
print(args[3])
print(args[4])
print(args[5])

rep=args[5]

repdir=paste0(args[1],"/",args[4])
dir.create(repdir)
print("READING METADATA")
metadata <- read.csv(args[3])

is_monomorphic <- function(row) {
  return(length(unique(row)) == 1)
}

## Function to identify and remove monomorphic rows
remove_monomorphic_rows <- function(genotype_data) {
  monomorphic_rows <- apply(genotype_data, 1, is_monomorphic)
  polymorphic_data <- genotype_data[!monomorphic_rows, ]
  return(polymorphic_data)
}


geno_object <- paste0(args[1],"/genotypes.lfmm")
Loci_list <-  paste0(args[1],"/Loci.csv", sep="")

print("READING ALLELE FREQUENCIES")
DATA_full <- read.table(args[2], h=T, sep="\t")
rownames(DATA_full) <- paste(DATA_full$Chr, DATA_full$Pos, sep=".")
#SUB <- DATA_full[3:length(DATA_full)]
SUB <- DATA_full
samples <- gsub("\\.", "-", colnames(DATA_full))
colnames(SUB) <- gsub("\\.", "-", colnames(SUB))
ss_data <- metadata[metadata$sample %in% samples,]
s_data <- na.exclude(ss_data)
SUB2 <- SUB[colnames(SUB) %in% s_data$sample]
SUB2 <- remove_monomorphic_rows(SUB2)
SUB3 <- SUB2[,3:ncol(SUB2)]
s_data <- s_data[s_data$sample %in% colnames(SUB3),]

split_parts <- strsplit(rownames(SUB3), "\\.")
Chr <- sapply(split_parts, function(x) x[1])  
Pos <- sapply(split_parts, function(x) x[2])  
print(nrow(SUB3))
print(ncol(SUB3))

print("Loci")
Loci <- as.data.frame(cbind(Chr, Pos))

print(nrow(Loci))
#colnames(Loci) <- c("Chr", "Pos")

if (file.exists(geno_object) && file.info(geno_object)$size > 0){
  # If the file exists, read it
  print("FILES EXIST, PLEASE CHECK FOR CORRECTNESS")
  tr <- read.lfmm(geno_object)
  Loci <- read.csv(Loci_list)
  print("File exists and has been opened.")
} else {
  write.csv(Loci, Loci_list,row.names = FALSE)
  #write.csv(SUB2, paste0(args[1],"/CleanedEnv.csv", sep=","),row.names = FALSE)
  tr <- t(as.matrix(SUB3))
  tr <- round(tr,4)
  print("DIMENSIONS OF GENO")
  print(nrow(tr))
  print(ncol(tr))
  # Write the matrix to an LFMM file
  write.lfmm(tr,geno_object)
  print("A new geno.lffm.file has been created.")
}
#
varname_old=paste(args[4])
varname=gsub("\\+0", "\\.0", varname_old)
#var=grep(args[4], colnames(s_data))

#print(colnames(s_data))
print(varname)
#print(var)
var_index <- match(varname, colnames(s_data))
#print(colnames(s_data))
print(varname)
#print(s_data[,varname])
#if (is.na(var_index)) {
#    stop("The specified column does not exist in s_data.")
#} else {
f <- as.matrix(s_data[, var_index])  # Use the index to extract the column
#}


print(f)
write.env(f, paste0(repdir, "/gradients.env"))
#
run_lfmm <- function(data, env, K_range = 1:2) {
    project <- snmf(data, K = K_range, entropy = TRUE, repetitions = 1, project = "new")
    CE <- sapply(K_range, function(K) mean(cross.entropy(project, K = K)))
    K <- which.min(CE)
    print(paste("Number of K: ", K))
    mod.lfmm2 <- lfmm2(data, env, K = K)
    return(mod.lfmm2)
}

print(ncol(f))
print(nrow(f))

print('CALCULATING LFMM2')
mod.lfmm2 <- run_lfmm(geno_object, f)
print('SUCCESFULLY CALCULATED MOD.LFMM2')

print('TESTING THE LFMM2 OBJECT')           
pv <- lfmm2.test(object = mod.lfmm2, input = tr, env = f, linear = TRUE)


print('EXPORTING P-VALUES TO CSV FILE')
print(pv$pvalues)
markerpos <- data.frame(
  Chr = SUB2$Chr,
  Pos = SUB2$Pos,
  Pval = as.numeric(pv$pvalues)  # Make sure the p-values are numeric
)

write.csv(markerpos, paste(repdir,"/",varname,"_LEA_pvals_all.csv", sep=""),row.names = FALSE)

Bonf=0.05/(nrow(markerpos)+1)
#markerpos <- as.data.frame(markerpos)
#colnames(markerpos) <- c("Pos", "Pval")
outliers <- markerpos[markerpos$Pval < Bonf,]
write.csv(outliers, paste(repdir,"/",varname,"_LEA_pvals_outliers.csv", sep=""),row.names = FALSE)


print("LOADING ANNOTATIONS")

annotations <- read.table("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/results/annotationdata/annotations.txt", header = FALSE)
annotLoci <- paste(annotations$V1, annotations$V2, sep=".")
colnames(annotations) <- c("Chr", "Pos", "Gene")
u_annotations <- unique(annotations)
#annotations <- merge(u_annotations, markerpos, by = c("Chr", "Pos"), all.x = TRUE)


#plogv <- -log10(markerpos$Pval)

DATA <- merge(u_annotations, markerpos, by= c("Chr","Pos"))
#colnames(DATA) <- c("Chr", "Pos", "Gene", "pvalslog")

print("PLOT MANHATTAN")
#plot_list <- list()
#
#pl <- ggplot(markerpos, aes(x = Pos, y = -log10(Pval))) +
#  geom_point(col = rgb(0, 0, 0, 0.1), pch = 16) +
#  facet_grid(. ~ DATA.Chr, scales = "free_x", space = "free") +
#  geom_hline(yintercept = -log10(0.05 / (nrow(markerpos) + 1)), colour = "blue") +
#  xlab("Genomic Position [Mbp]") +
#  ylab(paste0("-log10(P.value): ", varname)) +
#  theme_bw()
#  ggtitle(paste0("TEST")) +
#  theme(legend.position = "none")
#plot_list.append(pl)
# Arrange the plots in a grid
#PLOT_grid <- do.call(grid.arrange, c(plot_list, ncol = 1)) # Adjust nrow and ncol as needed

# Save the combined plot to a file
#ggsave(paste0(repdir,"/P_values.png"), PLOT_grid, width = 15, height = 7)      

#pl <- ggplot(DATA, aes(x = Pos, y = Pval)) +
#        geom_point(aes(color = Pval > Bonf), alpha = 0.5) +  # Color points based on pvalslog > 10
#        facet_wrap(~ Chr, scales = "free_x") +
#        geom_hline(yintercept = Bonf, col = "blue", lty = 2) +
#        geom_hline(yintercept = -log10(0.05), col = "red", lty = 2) +
#        labs(title = paste("Manhattan Plot for", "pval"), x = "Position", y = "-log10(p-value)") +
#        theme_bw() +
#        geom_text(aes(label = ifelse(Pval > Bonf, Gene, '')), vjust = -0.5, check_overlap = TRUE) +  # Only label points where pvalslog > 10
#        scale_color_manual(values = c("black", "red")) +
#        theme(legend.position = "none")
#

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