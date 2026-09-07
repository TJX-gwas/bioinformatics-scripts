suppressPackageStartupMessages(library(rMVP))

geno_prefix <- "mvp.data"
pheno_file <- "phenotype_final_merged.txt"
pca_file <- "My_PCA_Covariates.txt"
ind_file <- paste0(geno_prefix, ".geno.ind")

geno_ids <- read.table(ind_file, stringsAsFactors = FALSE)[[1]]
phenotype <- read.table(pheno_file, header = TRUE, sep = "\t", check.names = FALSE)
pca <- read.table(pca_file, header = FALSE, stringsAsFactors = FALSE)

names(phenotype)[1] <- "Taxa"
pca <- pca[, c(2, 3:7)]
names(pca) <- c("Taxa", paste0("PC", 1:5))

if (anyDuplicated(phenotype$Taxa) || anyDuplicated(pca$Taxa)) {
  stop("Duplicate sample IDs detected")
}
if (!setequal(geno_ids, phenotype$Taxa) ||
    !setequal(geno_ids, pca$Taxa)) {
  stop("Genotype, phenotype, and PCA sample sets do not match")
}

phenotype <- phenotype[match(geno_ids, phenotype$Taxa), ]
pca <- pca[match(geno_ids, pca$Taxa), ]
stopifnot(
  identical(geno_ids, phenotype$Taxa),
  identical(geno_ids, pca$Taxa)
)

covariate <- as.matrix(cbind(pca[, -1, drop = FALSE]))
storage.mode(covariate) <- "numeric"

genotype <- bigmemory::attach.big.matrix(paste0(geno_prefix, ".geno.desc"))
map <- read.table(paste0(geno_prefix, ".geno.map"), header = TRUE)
trait_names <- names(phenotype)[-1]

cat("rMVP version:", as.character(packageVersion("rMVP")), "\n")
cat("Samples:", length(geno_ids), "\n")
cat("Markers:", nrow(map), "\n")
cat("Covariates:", paste(colnames(covariate), collapse = ", "), "\n")

for (trait in trait_names) {
  out_dir <- paste0("MVP_Result_", trait)
  check_file <- file.path(out_dir, paste0(trait, ".FarmCPU.csv"))
  if (file.exists(check_file)) {
    cat("Skipping completed trait:", trait, "\n")
    next
  }
  dir.create(out_dir, showWarnings = FALSE)
  phe <- phenotype[, c("Taxa", trait), drop = FALSE]

  MVP(
    phe = phe,
    geno = genotype,
    map = map,
    CV.GLM = covariate,
    CV.MLM = covariate,
    CV.FarmCPU = covariate,
    method = c("GLM", "MLM", "FarmCPU"),
    ncpus = 8,
    outpath = out_dir,
    file.output = TRUE,
    p.threshold = 0.05
  )
  gc()
}
