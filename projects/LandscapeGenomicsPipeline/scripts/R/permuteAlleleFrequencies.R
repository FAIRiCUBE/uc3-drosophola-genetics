##PermuteAF

# Load Allele Frequencies
af_file="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/results/fullgenome/Subsampled_fullgenome.final_DP15.af"
af <- read.table(af_file, header=TRUE, na.strings = NaN)
af2 <- as.data.frame(t(af %>% dplyr::select(3:ncol(af))))
## Create an array with numbers from 1 to 300
#arrow_array <- 1:ncol(af)
#
#permutationMatrix <- c()
## Generate a random resampling of the array
#for (i in 1:100) {
#  resampled_arrows <- sample(arrow_array, size = length(arrow_array), replace = TRUE)
#  permutationMatrix[[i]] <- resampled_arrows
#}

arrow_array <- rownames(af2)
#
permutationPops <- c()
# Generate a random resampling of the array
for (i in 1:100) {
  resampled_arrows <- sample(arrow_array, size = length(arrow_array), replace = FALSE)
  permutationPops[[i]] <- resampled_arrows
}

write.table(permutationPops, "/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/PermutationMatrix.rds", col.names = FALSE, row.names = FALSE)

##### Version 2 

arrow_array <- rownames(af2)
permutationPops <- list()

# Function to generate a valid permutation efficiently
generate_valid_permutation <- function(arr) {
  shuffled <- sample(arr)  # Initial shuffle
  idx <- 1:length(arr)

  # Check and swap conflicts
  for (i in seq_along(arr)) {
    if (substr(arr[i], 1, 2) == substr(shuffled[i], 1, 2)) {
      swap_idx <- which(substr(arr, 1, 2) != substr(shuffled[i], 1, 2))
      swap_idx <- swap_idx[swap_idx != i]  # Avoid swapping with itself
      if (length(swap_idx) > 0) {
        j <- sample(swap_idx, 1)  # Pick a random valid swap
        shuffled[c(i, j)] <- shuffled[c(j, i)]
      }
    }
  }
  return(shuffled)
}

# Generate permutations
for (i in 1:100) {
  permutationPops[[i]] <- generate_valid_permutation(arrow_array)
}

write.table(permutationPops, "/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/PermutationMatrix2.rds", col.names = FALSE, row.names = FALSE)
