suppressPackageStartupMessages(library(LEA))

genotype_file <- "sheep_220_final.lfmm"
environment_file <- "New_Gradients_PC1_220.env"
bim_file <- "sheep_220_final.bim"

env <- scan(environment_file, quiet = TRUE)
if (length(env) != 220) stop("LFMM environment file must contain 220 values")

set.seed(20260729)
model <- lfmm2(input = genotype_file, env = environment_file, K = 6)
test <- lfmm2.test(
  object = model,
  input = genotype_file,
  env = environment_file,
  linear = TRUE
)

p_raw <- as.numeric(test$pvalues)
p_for_chisq <- pmin(pmax(p_raw, .Machine$double.xmin), 1 - .Machine$double.eps)
chisq <- qchisq(p_for_chisq, df = 1, lower.tail = FALSE)
lambda_gc <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
p_gc <- pchisq(chisq / lambda_gc, df = 1, lower.tail = FALSE)
fdr_bh <- p.adjust(p_gc, method = "BH")

bim <- read.table(bim_file, header = FALSE, stringsAsFactors = FALSE)
names(bim) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")
if (length(p_raw) != nrow(bim)) {
  stop("LFMM P-value count differs from the final BIM marker count")
}

results <- data.frame(
  CHR = bim$CHR,
  SNP = bim$SNP,
  BP = bim$BP,
  A1 = bim$A1,
  A2 = bim$A2,
  P_raw = p_raw,
  P_GC = p_gc,
  FDR_BH = fdr_bh
)
write.csv(results, "LFMM_PC1_K6_all_SNPs.csv", row.names = FALSE, quote = FALSE)

significant <- results[!is.na(results$FDR_BH) & results$FDR_BH < 0.05, ]
write.csv(
  significant,
  "LFMM_PC1_K6_GIF_BH_FDR05_Significant_SNPs.csv",
  row.names = FALSE,
  quote = FALSE
)

cat("R version:", R.version.string, "\n")
cat("LEA version:", as.character(packageVersion("LEA")), "\n")
cat("K: 6\n")
cat("Lambda GC:", lambda_gc, "\n")
cat("FDR < 0.05 SNPs:", nrow(significant), "\n")

