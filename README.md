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
* Preprocessing
  1. Quality control and filtering of low-count genes
  2. Variance-stabilizing transformation (VST)
  3. Sample and gene QC checks

* Differential Expression
  1. T2-high vs T2-low
  2. ICS vs No ICS
  3. Stratified contrasts (e.g., T2-high + ICS vs T2-high without ICS)
  4. Visualization using volcano plots and heatmaps

* Co-expression Network Analysis (WGCNA)
  1. Construction of scale-free gene co-expression networks
  2. Identification of gene modules
  3. Calculation of:
     * Module eigengenes (ME)
     * Module membership (MM)
     * Gene significance (GS)
* Correlation of modules with:
 1. T2 status
 2. ICS use
 3. Combined phenotypes
    
* Functional Enrichment
  1. GO and pathway enrichment for key modules
  2. Biological interpretation of module-level behavior
  3. Identification of hub genes within clinically relevant modules

## Key Questions Addressed

1. How does T2 inflammation shape transcriptomic profiles in asthma?

2. Which gene networks are associated with ICS use?

3. Are there modules that distinguish steroid-responsive vs non-responsive patients?

Can network-level signals explain heterogeneity beyond DE genes?
