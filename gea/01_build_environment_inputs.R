args <- commandArgs(trailingOnly = TRUE)
env_file <- if (length(args) >= 1) args[1] else "Sheep_21_Env_Variables_1km.csv"
fam_file <- if (length(args) >= 2) args[2] else "sheep_220_final.fam"

env <- read.csv(env_file, stringsAsFactors = FALSE, check.names = FALSE)
names(env)[1] <- "Pop"
fam <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
names(fam)[1:2] <- c("FID", "IID")

# Exact variables that reproduce the LFMM environmental PC1 used in the study.
lfmm_variables <- c(
  "Elevation", "BIO1", "BIO2", "BIO3", "BIO5",
  "BIO9", "BIO12", "BIO14", "BIO15", "BIO17"
)
if (!all(lfmm_variables %in% names(env))) stop("Required environment columns are missing")

pca <- prcomp(env[, lfmm_variables], center = TRUE, scale. = TRUE)
pc1 <- pca$x[, 1]

# Resolve the arbitrary PCA sign so that larger PC1 corresponds to elevation.
if (cor(pc1, env$Elevation) < 0) {
  pc1 <- -pc1
  pca$rotation[, 1] <- -pca$rotation[, 1]
}

pop_pc1 <- data.frame(Pop = env$Pop, PC1 = pc1)
write.table(pop_pc1$PC1, "New_Gradients_PC1.env", row.names = FALSE, col.names = FALSE)
write.table(
  data.frame(variable = rownames(pca$rotation), loading_PC1 = pca$rotation[, 1]),
  "environment_PC1_loadings.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)

fam$Pop <- gsub("[0-9]+", "", fam$IID)
fam$Pop[fam$Pop == "GB"] <- "GBW"
sample_pc1 <- pop_pc1$PC1[match(fam$Pop, pop_pc1$Pop)]
if (anyNA(sample_pc1)) stop("At least one sample population was not matched to environment data")
if (length(sample_pc1) != 220) stop("Expected 220 sample-level PC1 values")

write.table(sample_pc1, "New_Gradients_PC1_220.env", row.names = FALSE, col.names = FALSE)
write.table(
  data.frame(IID = fam$IID, Pop = fam$Pop, Environmental_PC1 = sample_pc1),
  "New_Gradients_PC1_220_with_ID.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)

cat("LFMM PC1 variables:", paste(lfmm_variables, collapse = ", "), "\n")
cat("PC1 variance explained:", summary(pca)$importance[2, 1], "\n")

