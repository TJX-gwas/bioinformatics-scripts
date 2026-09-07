args <- commandArgs(trailingOnly = TRUE)
geno_ind_file <- if (length(args) >= 1) args[1] else "mvp.data.geno.ind"
fam_file <- if (length(args) >= 2) args[2] else "sheep_220_final.fam"
pheno_file <- if (length(args) >= 3) args[3] else "phenotype_final_merged.txt"
pca_file <- if (length(args) >= 4) args[4] else "My_PCA_Covariates.txt"

geno_ids <- read.table(geno_ind_file, stringsAsFactors = FALSE)[[1]]
fam_ids <- read.table(fam_file, stringsAsFactors = FALSE)[[2]]
pheno <- read.table(pheno_file, header = TRUE, sep = "\t", check.names = FALSE)
pca <- read.table(pca_file, header = FALSE, stringsAsFactors = FALSE)

pheno_ids <- pheno[[1]]
pca_ids <- pca[[2]]

sources <- list(
  mvp_ind = geno_ids,
  fam = fam_ids,
  phenotype = pheno_ids,
  pca = pca_ids
)

for (nm in names(sources)) {
  ids <- sources[[nm]]
  if (length(ids) != 220) stop(nm, " has ", length(ids), " samples, expected 220")
  if (anyDuplicated(ids)) stop(nm, " contains duplicate sample IDs")
  if (!setequal(geno_ids, ids)) {
    stop(nm, " sample set differs from mvp.data.geno.ind")
  }
  if (!identical(geno_ids, ids)) {
    stop(nm, " sample order differs from mvp.data.geno.ind")
  }
}

cat("All four sample sources contain the same 220 IDs in the same order.\n")

