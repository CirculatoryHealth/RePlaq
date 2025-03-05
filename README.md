<!--  ## AE_TEMPLATE -->
<!--   -->
<!--  This is a "Lookup request" template.  -->
<!--  The naming of the repository should follow the convention we use in > Trello, _etc._, _e.g._  "AE_20190910_008_JHILLEBRANDS_SDEJAGER_TEMS_TIE2". -->
<!--   -->
<!--  This template includes some (standard) codes/scripts for: -->
<!--   -->
<!--  - baseline tables, sample selections -->
<!--  - SNP lookups, GWAS, or gene-based lookups -->
<!--  - bulk RNAseq analyses -->
<!--  - scRNAseq projections and lookups -->

<!--  Provide a title.         -->
[AE_TEMPLATE](https://github.com/CirculatoryHealth/AE_TEMPLATE)<img align="right" height="200" src=images/AE_Genomics_2010.png>
============

<!--  Provide details on the people involved in the project.         -->
*Collaborators*

[First name] [Initials] [Last name]
[First name] [Initials] [Last name]
[First name] [Initials] [Last name]

*Athero-Express Team*

Sander W. van der Laan, 
Michal Mokry, 
Ernest Diez Benavente, 
Hester den Ruijter, 
Dominique de Kleijn,
Gert Jan de Borst, 
Gerard Pasterkamp.

<!--  Provide the project ID.         -->
**Project ID** [`AE_[YYYYMMDD]_[PROJECTNUMBER]_[LEADCOLLABORATOR]_[PROJECTNAME]`]

# Background
<!--  Provide some background, study design, results, etc.         -->
Collaboration to study common variants and (single-cell) gene expression in relation to atherosclerotic plaques characteristics. 


## Study design

We will test the hypothesis that common variants in and expression of (a) gene(s)-of-interest are associated with plaque characteristics. We will use data from the **Athero-Express Biobank Study**.

These are the questions we will address: 

- Are any of the _variants_ associated to plaque characteristics?
- Is the _gene expression_ correlated to characteristics of plaques?
- In which _cell types_ are the target genes expressed? 


## Athero-Express Biobank Study

We have bulk RNAseq (n ± 1,100 samples) and single-cell RNAseq data (n ± 46), genome-wide methylation (Illumina 450K) in n ± 600, as well as overlapping genetic data for ±2,000 individuals with extensive histological plaque characterisation. 


### Genetic analyses

For the genetic analyses we will perform regression analyses adjusted for age, sex (where applicable) and principal components. So, we will apply the following model:

We will perform regression analyses adjusted for age, sex (where applicable) and principal components. 

- model 1: `phenotype ~ age + sex + chip-used + PC1 + PC2 + year-of-surgery`

Optionally, we can also run sex-stratified analyses:

- model 1m: `phenotype ~ age + chip-used + PC1 + PC2 + year-of-surgery` males only
- model 1f: `phenotype ~ age + chip-used + PC1 + PC2 + year-of-surgery` females only

The phenotypes are:

- `calcification`, coded `Calc.bin` no/minor vs. moderate/heavy staining
- `collagen`, coded `Collagen.bin` no/minor vs. moderate/heavy staining
- `fat10`, coded `Fat.bin_10` no/<10% fat vs. >10% fat
- `fat40`, coded `Fat.bin_40` no/<40% fat vs. >40% fat
- `intraplaque hemorrhage`, coded `IPH.bin` no vs. yes
- `macrophages (CD68)`, coded `macmean0` mean of computer-assisted calculation CD68<sup>+</sup> region of interest
- `smooth muscle cells (alpha-actin)`, coded `smcmean0` mean of computer-assisted calculation SMA<sup>+</sup> region of interest
- `intraplaque vessel density (CD34)`, coded `vessel_density` manually counted CD34<sup>+</sup> cells per 3-4 hotspots
- `mast cells`, coded `Mast_cells_plaque` manually counted mast cell tryptase<sup>+</sup> cells (https://academic.oup.com/eurheartj/article/34/48/3699/484981) [Note: low sample size]
- `neutrophils (CD66b)`, coded `neutrophils` manually counted CD66b<sup>+</sup> cells (https://pubmed.ncbi.nlm.nih.gov/20595650/) [Note: low sample size]
- `plaque vulnerability index`, scaled from 0 to 4, where 0 is most stable, and 4 is least stable plaque phenotype.

Continuous variables were inverse-rank normal transformated, indicated by `_rankNorm`. 

`slideToolKit` derived histological plaque phenotypes are not yet available. 

**Figure 1: Genotyped individuals in the Athero-Express Biobank Study**
![Genotyped individuals in the Athero-Express Biobank Study](PLOTS/20240611.overlap.AEDB_AEGS123.UpSetR.png)


### Whole-plaque RNAseq

For the expression analysis we used carotid plaque-derived bulk RNAseq data and queried it for the gene list. Below a graph showing the overall expression of the genes (not all are in the data) compared to the mean expression of 1,000 randomly picked genes. 

**Figure 2A: Overall expression of target genes in carotid plaques from the Athero-Express Biobank Study**
![Overall expression of target genes in carotid plaques from the Athero-Express Biobank Study](PLOTS/20240611.TargetExpression_vs_1000genes.png)

**Figure 2B: Overall expression of target genes by experiment in carotid plaques from the Athero-Express Biobank Study**
![Overall expression of target genes by experiment in carotid plaques from the Athero-Express Biobank Study](PLOTS/20240611.TargetExpression_vs_1000genes_by_AERNAStudy.png)

We assessed the correlation with plaque characteristics (mentioned above) and secondary major adverse cardiovascular events (MACE [major]) at 30 days and 3 years after CEA. 


### Single cell RNAseq

We projected target genes to the single-cell RNAseq data derived from 46 carotid plaque samples. We identified cell communities (Figure 2), mapped and projected target gene expression to the cell communities (Figure 3). 

**Figure 3: Cell communities identified in carotid plaques from the Athero-Express Biobank Study**
![Cell communities identified in carotid plaques from the Athero-Express Biobank Study](PLOTS/20240611.UMAP.png)

**Figure 4: Dotplot showing expression of target genes per cell type in carotid plaques from the Athero-Express Biobank Study**
![Dotplot showing expression of target genes per cell type in carotid plaques from the Athero-Express Biobank Study](PLOTS/20240611.DotPlot.Targets.png)


# Where do I start?

You can load this project in RStudio by opening the file called 'RePlaq.Rproj'.

## Project structure

<!--  You can add rows to this table, using "|" to separate columns.         -->
File                      | Description                      | Usage         
------------------------- | -------------------------------- | --------------
README.md                 | Description of project           | Human editable
RePlaq.Rproj              | Project file                     | Loads project 
LICENSE                   | User permissions                 | Read only     
.worcs                    | WORCS metadata YAML              | Read only     
prepare_data.R            | Script to process raw data       | Human editable
manuscript/manuscript.Rmd | Source code for paper            | Human editable
manuscript/references.bib | BibTex references for manuscript | Human editable
renv.lock                 | Reproducible R environment       | Read only     
File                                    | Description                          | Usage         
--------------------------------------- | ------------------------------------ | --------------
README.md                               | Description of project               | Human editable
AE_TEMPLATE.Rproj                       | Project file                         | Loads project
LICENSE                                 | User permissions                     | Read only
.worcs                                  | WORCS metadata YAML                  | Read only
renv.lock                               | Reproducible R environment           | Read only
images                                  | image directory for project          | Human editable
BASELINE                                | Baseline characteristics directory   | Human editable
OUTPUT                                  | Output directory                     | Human editable
PLOTS                                   | Some results                         | Human editable
SNP                                     | SNP analysis directory               | Human editable
scripts                                 | Scripts directory                    | Human editable
targets                                 | Directory containing list of targets | Human editable
manuscript                              | Source code for paper                | Human editable
manuscript/references.bib               | Manuscript BibTex references         | Human editable
packages.bib                            | BibTex references for packages used  | Human editable
references.bib                          | BibTex references                    | Human editable
preregistration                         | Preregistered hypotheses             | Human editable
preregistration/preregistration.rmd     | Preregistered document.              | Human editable
1_AEDB.CEA.baseline.Rmd                | Preparing data, baseline table       | Human editable
2_SNP_analyses.Rmd                     | Preparing SNP analyses, b38               | Human editable
2_SNP_analyses_b37.Rmd                     | Preparing SNP analyses, b37               | Human editable
3_1_bulkRNAseq.preparation.Rmd          | Preparing bulk RNAseq analyses       | Human editable
3_2_bulkRNAseq.exploration.Rmd          | Exploration RNAseq.                  | Human editable
3_3_bulkRNAseq.main_analysis.Rmd        | Main RNAseq analyses                 | Human editable
3_4_bulkRNAseq.additional_figures.Rmd   | Additional RNAseq figures            | Human editable
4_scRNAseq.Rmd                          | Single-cell RNAseq analyses          | Human editable
5_DNAm_AEMS450K_preparation.Rmd         | Prepare dataset for DNAm analyses          | Human editable

<!--  You can consider adding the following to this file:                    -->
<!--  * A citation reference for your project                                -->
<!--  * Contact information for questions/comments                           -->
<!--  * How people can offer to contribute to the project                    -->
<!--  * A contributor code of conduct, https://www.contributor-covenant.org/ -->

## Create a target list

You should use targets.xlsx to create a list of target genes and variants. This list will be used in the analyses. 

### Variants

The genetic variants should be b37 (legacy HRC r1.1 merge with 1000G p3v5 imputed data) or b38 (TOPMed release 3, freeze 10) formatted. For b37 you only need the rsID (optional), chromosome and position. For b38 you need the chromosome and position, and reference and alternate alleles. The per-variant lookup is done using the `VariantID` which is a combination of chromosome and position (b37) or chromosome, position, reference allele, and alternate allele (b38).

For example, for b37 you could have: rs1234567 1 1234567 where `VariantID` will automatically be coded as `1:1234567`. For b38 you could have: rs1234567 1 1234567 A T where `VariantID` will automatically be coded as `chr1:1234567:A:T`.

To check the reference allele use Bravo: https://bravo.sph.umich.edu. 

### Genes

The bulk RNAseq and scRNAseq are mapped against b38 and Ensembl 86, so genes should match these versions. In 9 out of 10 of the cases all genes are found. When a gene is not found this is likely due to a different gene name or a gene that is not in the data. In the former case you can use the `biomaRt` package to find the correct gene name, or simply use [www.genecards.org](https://www.genecards.org) to check out synonyms.

### 1. AEDB.CEA.baseline.Rmd

This notebook is used to create a `RDS`-formatted datasets with the selection of samples.

### 2. SNP_analyses.Rmd

This notebook is used to prepare the SNP analyses using the output of the `1_AEDB.CEA.baseline.Rmd` notebook.

A few notes on the SNP analyses. These make use of `GWASToolKit` which is a collection of scripts to make analyses with `SNPTEST` (currently `v2.5.4-beta3`) easier and streamline a lot of default and bottersome steps in the process. 

We are using `vcf.gz` files; these are bgzipped and tabix indexed. Since we are using `vcf`-files we have to indicate the genotype-field to have `SNPTEST` handle genotype probabilities instead of hardcoded genotypes. At the moment this is hardcoded in `GWASToolKit`. A future release of `GWASToolKit`, dubbed `GWASToolKit2`, will enable using `SNPTEST 2.5.6`, `plink2` and `regenie` (as each have their practical or theoretical statistical pro's and cons), and make the genotype-field and many new options a parameter while providing an easier commandline interface.

### 3.1 bulkRNAseq.preparation.Rmd

This notebook is used to prepare the bulk RNAseq analyses using the output of the `3_1_bulkRNAseq.preparation.Rmd` notebook.

### 3.2 bulkRNAseq.exploration.Rmd

This notebook used the output from `3_1_bulkRNAseq.preparation.Rmd` to explore the expression of the target gene(s) compared to some variables and the expression of 1,000 randomly picked genes.

### 3.3 bulkRNAseq.main_analysis.Rmd

This is the main-analysis notebook uses the output from `3_1_bulkRNAseq.preparation.Rmd` and `3_ _2_bulkRNAseq.exploration.Rmd`. Here linear and logistic regression analyses are executed using two models, `model 1` is a 'simple' model correcting for age, sex, year-of-surgery and RNAseq experiment (AERNAS1 and AERNAS2), `model 2` is the same as `model 1` but extended with a selection of traditional risk factors. 

### 3.4 bulkRNAseq.additional_figures.Rmd

Here the output from `3_3_bulkRNAseq.main_analysis.Rmd` is used to create additional informative figures that could be used for publication. 

### 4. scRNAseq.Rmd

Here the expression of the target gene(s) is inspected and visualised in the single-cell RNAseq dataset using output from `1_AEDB.CEA.baseline.Rmd`. 

### Additional utility scripts

A few key utility `python3` scripts are available. Please refer to `python3 scripts/[given_script].py --help` for more information of the respective functions.

#### parse_molqtl.py

The target variant(s) are inspected and looked up in the molQTL results (version 1/preliminary, n±600) using `parse_molqtl.py`. The following molQTL analyses results are queried.

- _cis_-acting expression QTL (eQTL) nominal results
- _trans_-acting expression QTL (eQTL) permuted (100-10,000x) results
- _cis_-acting methylation QTL (mQTL) permuted (100-10,000x) results
- _trans_-acting mQTL permuted (100-10,000x) results

> Note that this is currently limited to b37 based results. In the near future this will be updated to b38 with deeper sequenced RNAseq data. Sex- and smoking-interaction and stratified analyses will also be added, in addition to other stratified analyses for type 2 diabetes, kidney function, BMI, and medication use.

#### liftover.py

This script converts genomic coordinates between genome builds using the `pyliftover` package.

#### alt_models.py

This script runs logistic (`--logistic`) or linear regression (`--linear`) in PLINK v1.9 for a given SNP (`--snp`) and chromosome (`--chr`) and a given list of phenotypes (`--phenotypes`) and covariates (`--covariates`). Includes option for interaction analyses (`--interaction`). Can be run locally or on the Utrecht high-performance compute cluster.

#### alt_models_summarize.py

This script will concatenate the results of the logistic or linear regression models (`--glm`) for a given SNP (`--snp`) and a list of phenotypes (`--phenotypes`). Should be run directly on the data - there is no job-submission.


### manuscript & preregistration

In the spirit of Open Science you could opt to preregister your work, for instance through [OSF](https://osf.io). Alternatively, you can simply preregister your work by making this repository public and add the `preregistration.Rmd` filled in to the best of your knowledge and abilities. It is a document to register your intend. Research is at times intractable, and so don't worry too much about changes in the future. This is exactly the point: the research process is iterative and you may (have to) adapt your methods and/or your approach. Through this document you can document these thoughts and changes. 

Likewise, you can opt to write your manuscript - at least the draft - in `markdown` or `LaTeX` using `manuscript.Rmd`. 

# Reproducibility

This project uses the Workflow for Open Reproducible Code in Science (WORCS) to ensure transparency and reproducibility according to the principles of [FAIR](https://www.go-fair.org/fair-principles/). The workflow is designed to meet the principles of Open Science throughout a research project. 

To learn how WORCS helps researchers meet the TOP-guidelines and FAIR principles, read the preprint at https://osf.io/zcvbs/.

This repository works with `renv`. This is a package that creates a reproducible environment for your project. It is a good practice to use `renv` to ensure that your project is reproducible. However, it is not mandatory and it is up to the user to decide whether to use it or not.

## WORCS: Advice for authors

* To get started with `worcs`, see the [setup vignette](https://cjvanlissa.github.io/worcs/articles/setup.html)
* For detailed information about the steps of the WORCS workflow, see the [workflow vignette](https://cjvanlissa.github.io/worcs/articles/workflow.html)

## WORCS: Advice for readers

Please refer to the vignette on [reproducing a WORCS project]() for step by step advice.
<!-- If your project deviates from the steps outlined in the vignette on     -->
<!-- reproducing a WORCS project, please provide your own advice for         -->
<!-- readers here.                                                           -->

# Acknowledgements
Dr. Sander W. van der Laan is funded through EU H2020 TO_AITION (grant number: 848146), EU HORIZON NextGen (grant number: 101136962), EU HORIZON MIRACLE (grant number: 101115381), and Health~Holland PPP Allowance ‘Getting the Perfect Image’.

We are thankful for the support of the Leducq Fondation ‘PlaqOmics’ and ‘AtheroGen’, and the Chan Zuckerberg Initiative ‘MetaPlaq’. The research for this contribution was made possible by the AI for Health working group of the [EWUU alliance](https://aiforhealth.ewuu.nl/). The collaborative project ‘Getting the Perfect Image’ was co-financed through use of PPP Allowance awarded by Health~Holland, Top Sector Life Sciences & Health, to stimulate public-private partnerships.

Plaque samples are derived from endarterectomies as part of the [Athero-Express Biobank Study](https://doi.org/10.1007/s10564-004-2304-6) which is an ongoing study in the UMC Utrecht. We would like to thank all the (former) employees involved in the Athero-Express Biobank Study of the Departments of Surgery of the St. Antonius Hospital Nieuwegein and University Medical Center Utrecht for their continuing work. Lastly, we would like to thank all participants of the Athero-Express Biobank Study; without you these kinds of studies would not be possible.

The framework was based on the [`WORCS` package](https://osf.io/zcvbs/).

## Disclosures
Dr. Sander W. van der Laan has received Roche funding for unrelated work.

<a href='https://uefconnect.uef.fi/en/group/miracle/'><img src='images/UEF_Miracle_Logo-07.png' align="center" height="75" /></a> <a href='https://www.to-aition.eu'><img src='images/to_aition.png' align="center" height="75" /></a> <a href='https://www.health-holland.com'><img src='images/logo_NL_HealthHollland_Wit-Oranje_RGB.png' align="center" height="35" /></a> <a href='https://www.nextgentools.eu'><img src='images/NextGen_1_Red.png' align="center" height="35" /></a> <a href='https://www.era-cvd.eu'><img src='images/ERA_CVD_Logo_CMYK.png' align="center" height="75" /></a> <a href=''><img src='images/leducq-logo-large.png' align="center" height="75" /></a> <a href='https://www.fondationleducq.org'><img src='images/leducq-logo-small.png' align="center" height="75" /></a> <a href='https://osf.io/zcvbs/'><img src='images/worcs_icon.png' align="center" height="75" /></a> <a href='https://doi.org/10.1007/s10564-004-2304-6'><img src='images/AE_Genomics_2010.png' align="center" height="100" /></a>

#### Changes log
    
    _Version:_      v1.4.5</br>
    _Last update:_  2024-10-04</br>
    _Written by:_   Sander W. van der Laan (s.w.vanderlaan [at] gmail [dot] com).
    
    **MoSCoW To-Do List**
    The things we Must, Should, Could, and Would have given the time we have.
    _M_
    
    - [] merge datasets for variants on chromosome Y and MT in vcf.gz format, bgzipped and tabix indexed
    - [x] add methylation data preparation workflow
    - [] add methylation data analysis workflow for target gene(s)
    - [] add OLINK (and other protein measurement techniques) workflow for target gene(s)
    - [x] add proper molQTL lookup script to automatically add molQTL lookup
    - [] streamlining standard steps and functions to reduce lines of code
    
    _S_

    - [] update the workflow to churn out reports
    - [x] edit molQTL lookup script to include gene-to-variant and variant-to-gene lookups
    - [x] edit molQTL lookup script to annotate target-list
    - [x] edit molQTL lookup script include variant-to-gene-to-cpg and gene-to-cpg-to-variants lookups
    
    _C_

    _W_

    **Changes log**
    * v1.4.6 Updates to molQTL lookup script; accordingly the `targets.xlsx` file changed. Dropped all output to reduce size. Updates to the readme.me. Added wiki-pages with further explanations.
    * v1.4.5 Fixed an issue with the selection variable in KeyTable to fix the SNP_analyses.Rmd downstream script-steps.
    * v1.4.4 Added correct reference to KeyTable for b38 data. Added 2_SNP_analyses_b37.Rmd for b37 analyses. Added DNAM preparation analysis notebook.
    * v1.4.3 Fixed an issue with the KeyTable where the samples selected in 2018 were no not selected. These samples, _i.e._, a given surgery and corresponding studynumber, are duplicate patients but _different_ plaques, hence the phenotyping could be different. 
    * v1.4.2 Added DNAm data preparation. Fixed a few minor issues.
    * v1.4.1 Fixed an issue where the UPID used at the time of genotyping and QC was not properly included. Fixed an issue where the `STUDY_TYPE` and the `STUDY_ARTERY` were not properly included. Added an export to `PLINK`-format. Added a proper selection on informed consent.
    * v1.4.0 Major update. Added corrected AEDB. Added b38 genetic data. Added updated parsing of molQTL results. Additional per-notebooke fixes.
    * v1.3.4 Textual fixes to the RNAseq preparation notebook. Fixed baseline table writing. 
    * v1.3.3 Fixes to the renv(). Small textual fixes to the notebooks. Added more information to this readme. Reorganize a few things.
    * v1.3.2 Textual fixes. Fixes to the bulk RNAseq workflow.
    * v1.3.1 Fixed baseline table writing. Added additional saving options (raw, normalized, and log-transformed data) for the bulk RNAseq.
    * v1.3.0 Some script changes. Update to AEDB. Update to RNAseq (deeper sequencing). 
    * v1.2.2 Some organisational updates. 
    * v1.2.1 Fixed some references in the README.md. 
    * v1.2.0 There were some critical fixes in the way each notebook starts. Updates to the way the local system is set up and packages are loaded. Updates to functions. Updates to the organisation of the various notebooks. 
    * v1.1.0 Major update to WORCS system. 
    * v1.0.6 Small bug fixes. 
    * v1.0.5 Added png for overlap-figure.
    * v1.0.5 Removed obsolete references to objects.
    * v1.0.4 Fixed a mistake in the chr X sample-file creation. Now the order matches the chr X data.
    * v1.0.3 Fixed weight of files (limit of 10Mb per file for templates). Renamed entire repo.
    * v1.0.2 Added sex-specific .sample-files. Added GWASToolKit input-files.
    * v1.0.0 Initial version. Add 'plaque vulnerability index', Fixed baseline table, added codes, and results. Created sample-files. 

--------------

#### Creative Commons BY-NC-ND 4.0
##### Copyright (c) 1979-2024 Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com.
