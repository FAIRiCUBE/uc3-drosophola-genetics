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
samples <- gsub("\\.", "-", colnames(DATA_full))
colnames(SUB) <- gsub("\\.", "-", colnames(SUB))
ss_data <- metadata[metadata$sampleId %in% samples,]
s_data <- na.exclude(ss_data)
SUB2 <- SUB[colnames(SUB) %in% s_data$sampleId]
SUB2 <- remove_monomorphic_rows(SUB2)

geno_object <- paste0(args[1],"/genotypes.lfmm")

if (file.exists(geno_object)) {
  # If the file exists, read it
  tr <- read.lfmm(geno_object)
  print("File exists and has been opened.")
} else {
  tr <- t(as.matrix(SUB2))
  tr <- round(tr,4)
  print("DIMENSIONS OF GENO")
  print(nrow(tr))
  print(ncol(tr))
  # Write the matrix to an LFMM file
  write.lfmm(tr,geno_object)
  print("A new geno.lffm.file has been created.")
}

varname=paste(args[4])

var=grep(args[4], colnames(s_data))
f <- as.matrix(s_data[,var])
write.env(f, paste0(repdir, "/gradients.env"))

run_lfmm <- function(data, env, K_range = 1:1) {
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


print('EXPORTIN P-VALUES TO CSV FILE')     
markerpos <- cbind(rownames(SUB2),c(pv$pvalues))
write.csv(markerpos, paste(repdir,"/",varname,"_LEA_pvals_noname.csv", sep=""))

Bonf=0.05/(nrow(markerpos)+1)
markerpos <- as.data.frame(markerpos)
colnames(markerpos) <- c("Pos", "Pval")
outliers <- markerpos[markerpos$Pval < Bonf,]
write.csv(outliers, paste(repdir,"/",varname,"_LEA_pvals_outliers.csv", sep=""))

                 
print("PLOT MANHATTAN")
plot_list <- list()

pl <- ggplot(markerpos, aes(x = Pos, y = -log10(Pval))) +
  geom_point(col = rgb(0, 0, 0, 0.1), pch = 16) +
  facet_grid(. ~ DATA.Chr, scales = "free_x", space = "free") +
  geom_hline(yintercept = -log10(0.05 / (nrow(markerpos) + 1)), colour = "blue") +
  xlab("Genomic Position [Mbp]") +
  ylab(paste0("-log10(P.value): ", z)) +
  theme_bw()
  ggtitle(paste0("TEST")) +
  theme(legend.position = "none")
plot_list.append(pl)

# Arrange the plots in a grid
PLOT_grid <- do.call(grid.arrange, c(plot_list, ncol = 1)) # Adjust nrow and ncol as needed

# Save the combined plot to a file
ggsave(paste0(repdir,"/P_values.png"), PLOT_grid, width = 15, height = 7)      
