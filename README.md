## Analysis Workflow and Parameters

GWAS was performed using rMVP, with the first five genetic principal components included as covariates. FarmCPU results were used for the primary analyses reported in the manuscript.

For GEA, the first principal component derived from environmental variables was used as the explanatory variable. Each population-level environmental score was assigned to individuals from that population, and the resulting vector was aligned with the sample order in the genotype file. LFMM association results were corrected for genomic inflation, followed by Benjamini–Hochberg correction for multiple testing. Candidate loci were identified at a false discovery rate (FDR)<0.05.

## Usage Notes

This repository provides GWAS and GEA scripts associated with the study. It does not include the complete data-preprocessing workflow or the study input datasets. Before running the scripts, users should verify the input formats, sample identifiers, and sample order, and configure file paths for their computing environment.
