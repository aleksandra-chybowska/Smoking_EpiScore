library(tidyverse)

# Read in weights of predictors

setwd("<cluster_home_dir>/Smoking/Replication/")

predictors = read.csv("Predictors_long.csv")
mCigarette = subset(predictors, Predictor == "mCigarette")
Daniel = subset(predictors, Predictor == "McCartney et al.")
BayesR = subset(predictors, Predictor == "Trejo et al.")
BayesR = subset(BayesR, Beta != 0)
EpiSmokEr = read.csv("elliot.tsv")


# Calculate overlaps with mCigarette
mCig_Daniel = intersect(Daniel$CpG, mCigarette$CpG)  # Overlap with McCartney (19)
mCig_Bayes = intersect(BayesR$CpG, mCigarette$CpG)    # Overlap with Trejo
mCig_EpiSmokEr = intersect(EpiSmokEr$CpG, mCigarette$CpG)

# Print results
cat("Overlap between mCigarette and Daniel: ", length(mCig_Daniel), "\n")
cat("Overlap between mCigarette and BayesR: ", length(mCig_Bayes), "\n")
cat("Overlap between mCigarette and EpiSmokEr: ", length(mCig_EpiSmokEr), "\n")

# Identify unique CpG sites in mCigarette
unique_mCigarette_Daniel = setdiff(mCigarette$CpG, Daniel$CpG)
unique_mCigarette_Bayes = setdiff(mCigarette$CpG, BayesR$CpG)
unique_mCigarette_EpiSmokEr = setdiff(mCigarette$CpG, EpiSmokEr$CpG)

# Unique CpG sites in mCigarette compared to all predictors
all_other_CpGs = unique(c(Daniel$CpG, BayesR$CpG, EpiSmokEr$CpG))
unique_mCigarette_all = setdiff(mCigarette$CpG, all_other_CpGs)

# Print unique sites for each comparison
cat("Number of unique mCigarette CpGs (not in Daniel): ", length(unique_mCigarette_Daniel), "\n")
cat("Unique mCigarette CpGs (not in Daniel):\n", unique_mCigarette_Daniel, "\n")

cat("Number of unique mCigarette CpGs (not in BayesR): ", length(unique_mCigarette_Bayes), "\n")
cat("Unique mCigarette CpGs (not in BayesR):\n", unique_mCigarette_Bayes, "\n")

cat("Number of unique mCigarette CpGs (not in EpiSmokEr): ", length(unique_mCigarette_EpiSmokEr), "\n")
cat("Unique mCigarette CpGs (not in EpiSmokEr):\n", unique_mCigarette_EpiSmokEr, "\n")

results = data.frame(CpG = mCigarette$CpG, Comment = "")

determine_overlap <- function(site, mCig_Daniel, mCig_Bayes, mCig_EpiSmokEr) {
  overlap_with <- c()  # Initialize an empty vector to store overlaps

  # Check for overlaps
  if (site %in% mCig_Daniel) {
    overlap_with <- c(overlap_with, "McCartney et al.")
  }
  if (site %in% mCig_Bayes) {
    overlap_with <- c(overlap_with, "BayesR")
  }
  if (site %in% mCig_EpiSmokEr) {
    overlap_with <- c(overlap_with, "EpiSmokEr")
  }

  # Determine if the site is unique or has overlaps
  if (length(overlap_with) == 0) {
    return("Unique to mCigarette")
  } else {
    return(paste0("Overlap with ", paste(overlap_with, collapse = ", ")))
  }
}

results$Comment <- sapply(results$CpG, function(site) {
  determine_overlap(site, mCig_Daniel, mCig_Bayes, mCig_EpiSmokEr)
})
intercept = data.frame(CpG = "(Intercept)", Beta = -8.93e-17, Predictor = "mCigarette")
mCigarette = rbind(intercept, mCigarette)
mCigarette = mCigarette[c("CpG", "Beta")]
df = merge(mCigarette, results, by.x = "CpG", by.y = "CpG", all.x = TRUE)
write.csv(df, "mCigarette_overlap_with_other_predictors.csv", row.names = FALSE)
