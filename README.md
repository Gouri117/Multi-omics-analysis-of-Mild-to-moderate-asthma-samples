# T2–ICS-300: Transcriptomic Analysis of Type 2 Asthma and Inhaled Corticosteroid Response
## Overview
This project investigates Type 2 (T2) asthma heterogeneity and inhaled corticosteroid (ICS) response using bulk RNA-seq data from ~300 induced sputum samples.
The goal is to identify gene expression patterns and co-expression modules associated with:
* T2-high vs T2-low asthma
* ICS use vs no ICS
* Interaction effects between T2 status and ICS treatment

The analysis combines differential expression, network-based methods (WGCNA), and functional enrichment to move beyond single-gene signals toward biologically interpretable modules.

## Dataset
* Sample type: Induced sputum
* Sample size: ~300 subjects
* Data: Bulk RNA-seq (raw counts)

* Metadata includes:
- T2 status (T2-high / T2-low)
- ICS use (Yes / No)
* Clinical covariates (e.g., asthma status, demographics)

Raw data files are not publicly included due to privacy restrictions.

## Analysis Workflow
1. Preprocessing
* Quality control and filtering of low-count genes
* Variance-stabilizing transformation (VST)
* Sample and gene QC checks

2. Differential Expression
* T2-high vs T2-low
* ICS vs No ICS
* Stratified contrasts (e.g., T2-high + ICS vs T2-high without ICS)
* Visualization using volcano plots and heatmaps

3. Co-expression Network Analysis (WGCNA)
i) Construction of scale-free gene co-expression networks
ii) Identification of gene modules
iii) Calculation of:
- Module eigengenes (ME)
- Module membership (MM)
- Gene significance (GS)
iv) Correlation of modules with:
- T2 status
- ICS use
- Combined phenotypes
v) Functional Enrichment
* GO and pathway enrichment for key modules
* Biological interpretation of module-level behavior
* Identification of hub genes within clinically relevant modules

## Key Questions Addressed

1. How does T2 inflammation shape transcriptomic profiles in asthma?

2. Which gene networks are associated with ICS use?

3. Are there modules that distinguish steroid-responsive vs non-responsive patients?

Can network-level signals explain heterogeneity beyond DE genes?
