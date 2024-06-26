##This scripts purpose will be applying LEA LFMM to a VCF and parallelize workflow via manual z-score pulling

args <- commandArgs(TRUE)
## Default setting when no arguments passed
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
rep=args[5]

repdir=paste0(args[1],"/",args[4])
dir.create(repdir)
#setwd(repdir)
fol <- paste(args[1])
#print(SUB[,1])
metadata <- read.csv(args[3])
#DATA_full <- read.table(args[2], h=T)

DATA_full <- read.table(args[2], h=T, sep="\t")
SUB <- DATA_full[3:length(DATA_full)]
samples <- colnames(SUB)

#print(SUB)
samples <- gsub("\\.", "-", samples)
colnames(SUB) <- gsub("\\.", "-", colnames(SUB))

#metadata2 <- read.csv("/media/inter/ssteindl/FC/DEST_freeze1/populationInfo/worldClim/dest.worldclim.csv")
s_data <- metadata[metadata$sampleId %in% samples,]
s_data <- s_data[s_data$V5!="",]
SUB2 <- SUB[colnames(SUB) %in% s_data$sampleId]
####intersection of samples and metadata needs to be checked!

#print(s_data$sampleId)
print("TEST")
#print(colnames(SUB2))
#write lfmm genotypes 
#tr <- t(as.matrix(SUB2))
#tr <- round(tr,4)
#print("DIMENSIONS OF GENO")
#print(nrow(tr))
#print(ncol(tr))
#
### one #

geno_object <- paste0(args[1],"/genotypes.lfmm")
if (file.exists(geno_object)) {
  # If the file exists, read it
  lfmm_geno <- read.lfmm(geno_object)
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




#
##### NEWLY ADDED ####
#pc = pca("genotypes.lfmm")
#tw = tracy.widom(pc)
#a=stars.pval(tw$pvalues)
#plot(tw$percentage, pch = 19, col = "darkblue", cex = .8)
#nK<-quick.elbow(tw$eigenvalues, low=0.02, max.pc = 0.9)
#print("NK")
#print(nK)
##genotype=lfmm2geno("genotypes.lfmm")

print(args[4])
#write environmental lfmm with K number of latent factos and r number of repetitions
var=grep(args[4], colnames(s_data))
varname=paste(args[4])
##K=args[5]
r=args[5]
####this factor needs to be possible as parameter
##factor=s_data[,3]
f <- s_data[,var]
#print(factor)

env_path <- paste0(repdir,"/gradients.env")
write.env(f, env_path)

run_lfmm <- function(data, env, K_range = 1:10) {
    project <- snmf(data, K = K_range, entropy = TRUE, repetitions = 2, project = "new")
    CE <- sapply(K_range, function(K) mean(cross.entropy(project, K = K)))
    K <- which.min(CE)
    mod.lfmm2 <- lfmm2(data, env, K = K)
    return(mod.lfmm2)
}


print('CALCULATING LFMM2')
#mod2 <- lfmm2("genotypes.lfmm", "gradients.env", K = nK)
mod.lfmm2 <- run_lfmm(geno_object, env_path)
               
pv <- lfmm2.test(object = mod.lfmm2, input = geno_object, env = env_path, linear = TRUE)
print(pv$pvalues)
#print(tr)
#mod2 <- lfmm2(input = tr, env = f, K = 5)
#print(mod2)
#pv <- lfmm2.test(object = mod2, input = tr, env = f, linear = TRUE)

plot(-log10(pv$pvalues), col = "grey", cex = .4, pch = 19)
#points(target, -log10(pv$pvalues[target]), col = "red")

markerpos <- as.data.frame(cbind(DATA_full$Chr, DATA_full$Pos))
markerpos <- cbind(markerpos,pv$pvalues)

write.csv(markerpos, paste(varname,"_LEA_pvals.csv", sep=""))
                 
print("PLOT MANHATTAN")
print(facs)
plot_list <- list()

markerpos$DATA.Chr <- as.factor(markerpos$V1)

#print(pvalcsv$DATA.Chr)
#PLOT.df <- data.frame(Chr = DATA$Chr, Pos = DATA$Pos / 1000000, P.val = -log10(DATA[[z]]))
pl <- ggplot(pvalcsv, aes(x = DATA.Pos, y = -log10(.data[[z]]))) +
  geom_point(col = rgb(0, 0, 0, 0.1), pch = 16) +
  facet_grid(. ~ DATA.Chr, scales = "free_x", space = "free") +
  geom_hline(yintercept = -log10(0.05 / (nrow(pvalcsv) + 1)), colour = "blue") +
  xlab("Genomic Position [Mbp]") +
  ylab(paste0("-log10(P.value): ", z)) +
  theme_bw()
  ggtitle(paste0("tttt")) +
  theme(legend.position = "none")
plot_list.append(pl)

print(plot_list)

# Arrange the plots in a grid
PLOT_grid <- do.call(grid.arrange, c(plot_list, ncol = 1)) # Adjust nrow and ncol as needed

# Save the combined plot to a file
ggsave("P_values.png", PLOT_grid, width = 15, height = 7)      
