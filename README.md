## Analysis scripts

All scripts are located in:

r_analysis/

The analysis is structured as a pipeline:

### Compartment assignment

- 01_assign_compartments_to_bins.R  
  Assigns A/B compartments from Hi-C eigenvectors (including optional 30% thresholds)

- 02_assign_compartments_to_genes.R  
  Maps compartment assignments to genes

---

### Core compartment results

- 03_compartment_distributions.R  
  Generates the main figure showing A/B distribution across tissues and chromosomes

- 05_upset_compartments.R  
  Produces the UpSet plot showing compartment combinations for duplicated genes

---

### Ohnolog analysis

- 06_ohnolog_preprocessing.R  
  Builds mapping between LOC and Ensembl IDs and constructs ohnolog gene sets

- 07_ohnolog_analysis.R  
  Performs:
  - Ohnolog vs singleton classification  
  - Pair-level compartment analysis  
  - Regulatory divergence (theta-shift) analysis  + figure
  - Statistical testing (Fisher and Wilcoxon)  
  - GO enrichment analysis
