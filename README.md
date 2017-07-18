# CBSD_Trancriptomics
This repository contains scripts and additional information for the manuscript:

Working title:
  ## Genomic selection for Cassava Response to CBSV. Advances on incorporating transcriptomics data to the prediction model
  
  List of scripts used in this manuscript:
  
  ### Describing_Population.Rmd (R notebook)
  1) Mean LDscore plot (Figure 1a)
  2) LD decay plot (Figure 1b)
  3) PCA plot (Using the SNPrelate R package) (Figure 1c)
  4) Reference Allele frequencies by population (Figure 1d)

  ### Impute2_genotypes.Rmd (R notebook)
  1) Concordance values for imputation (Cross validationshowing the accuracy of imputation per chromsome segment)
  2) Alternate allele frequency 
  3) info/AR2 distribution
  4) info/AR2 by Allele Frequency
  5) Distribution of the markers, comparing the expected natural variation as measured by Whole Genome Sequencing (Ramu et al.) with Imputation to Whole genome sequence using Beagle and Impute2. (Figure 2a)
  
  ### Impact_imputation_level.Rmd (R notebook)
  1) Comparison Of Genomic Prediction Accuracies when using GBS, BEAGLE and IMPUTE2 markers (Figure 3a).
  2) Comparison of Genomic Prediction Accuracies when using GBS, IMPUTE2 and subsets of IMPUTE2 markers (Figure 3b).
  
  ### Accounting_for_QTLs (R notebook)
  1) Prediction Accuracies using two kernel GBLUP for three traits. Comparing prediction accuracies when using chromosomes known to include QTLs and random chromosomes. (Figure S4a)
  2) Prediction Accuracies using a three kernel model. Chr4 and Chr11 (QTL + chromosomes) plus a random genetic effect representing the rest of the genome. (Figure S4b)
  3) Figure 4
