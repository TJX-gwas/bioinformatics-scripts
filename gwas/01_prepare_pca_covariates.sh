#!/usr/bin/env bash
set -euo pipefail

# This reconstructs the commands recorded in the successful PLINK logs.
PLINK=${PLINK:-plink}
INPUT_PREFIX=${INPUT_PREFIX:-Sheep_QC_Ultra}

"$PLINK" --bfile "$INPUT_PREFIX" \
  --chr-set 26 \
  --allow-extra-chr \
  --recode vcf-iid \
  --out Sheep_For_LD

"$PLINK" --vcf Sheep_For_LD.vcf \
  --chr-set 26 \
  --allow-extra-chr \
  --chr 1-26 \
  --geno 0.05 \
  --maf 0.05 \
  --indep-pairwise 50 10 0.5 \
  --out moderate_pruning

"$PLINK" --vcf Sheep_For_LD.vcf \
  --chr-set 26 \
  --allow-extra-chr \
  --extract moderate_pruning.prune.in \
  --pca 5 \
  --out Sheep_PCA_5PCs

test "$(wc -l < Sheep_PCA_5PCs.eigenvec)" -eq 220
test "$(wc -l < moderate_pruning.prune.in)" -eq 5495
