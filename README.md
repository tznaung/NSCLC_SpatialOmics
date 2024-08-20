# NSCLC Immunotherapy Response Prediction using Multi-Omics Spatial Analysis
This repository contains the code and analysis pipeline developed for a research project to understand the tumor immune microenvironment (TIME) in non-small cell lung cancer (NSCLC) using a multi-omics approach. Repository contents are 1) LASSO modelling, 2) data processing pipelines, and visualizations and plots. We aim to develop a spatially resolved predictive model to determine immunotherapy outcomes in NSCLC, utilizing advanced spatial proteomic and whole transcriptome profiling techniques.
# Background
NSCLC has an increasing number of targeted and systemic therapies, with subsets of patients benefiting from long-term durable responses. A comprehensive molecular characterization of the TIME is essential to understanding and predicting these therapeutic outcomes. We hypothesize that an integrated multi-omics approach will uncover novel interactions within the NSCLC TIME and identify biomarkers predictive of immunotherapy outcomes, thus advancing precision oncology.
# Methods
1. Spatial Proteomics: We used Phenocycler Fusion (PCF) for single-cell resolution spatial mapping of protein expression in the TIME.
2. Spatial Whole Transcriptome Profiling: We employed Digital Spatial Profiling (DSP-GeoMx-WTA) for multi-cellular readout of gene expression at cellular compartment resolution.
# Data Analysis
We applied a Least Absolute Shrinkage and Selection Operator (LASSO) model to derive gene signatures from cell type signatures and predict treatment outcomes.
# Cohorts
We analyzed tissue samples from two independent cohorts of advanced NSCLC patients treated with PD-1-based immunotherapies from Yale and the University of Queensland.
# Key Finding
1. Resistance associated cell types
2. Response associated cell types
3. Gene signature models
# Repository Contents
1. scripts/: Contains scripts for processing and analyzing spatial omics data.
2. data/: Input datasets (normalized and processed)
3. README.md: This file.
# Installation
1. R 4.0+
2. Required R packages: ggplot2, glmnet, survival, survminer, tidyverse, edgeR, pROC, anndata, umap, ggpubr, RColorBrewer
# License
This project is licensed under the MIT License - see the LICENSE file for details.
# Citation
Please cite our related research paper if you use or build upon this work.
Contact
For questions or issues, please contact Thazin Nwe Aung at thazin.aung@yale.edu OR tznaung@gmail.com.

