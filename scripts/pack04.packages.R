################################################################################
#                                PACKAGES TO LOAD                              #
################################################################################

cat("\n* General packages...\n")
install.packages.auto("readr")
install.packages.auto("optparse")
install.packages.auto("tools")
install.packages.auto("dplyr")
install.packages.auto("tidyr")
install.packages.auto("tidylog")
library("tidylog", warn.conflicts = FALSE)
install.packages.auto("naniar")

# To get 'data.table' with 'fwrite' to be able to directly write gzipped-files
# Ref: https://stackoverflow.com/questions/42788401/is-possible-to-use-fwrite-from-data-table-with-gzfile
# install.packages("data.table", repos = "https://Rdatatable.gitlab.io/data.table")
library(data.table)

install.packages.auto("tidyverse")
install.packages.auto("knitr")
install.packages.auto("DT")

install.packages.auto("haven")
install.packages.auto("tableone")

# Install the devtools package from Hadley Wickham
install.packages.auto('devtools')

cat("\n* Genomic packages...\n")
install.packages.auto("org.Hs.eg.db")
install.packages.auto("mygene")
install.packages.auto("EnhancedVolcano")

install.packages.auto("GenomicFeatures")
install.packages.auto("bumphunter")
install.packages.auto("minfi")
install.packages.auto("SummarizedExperiment")
install.packages.auto("IlluminaHumanMethylation450kmanifest")
install.packages.auto("IlluminaHumanMethylation450kanno.ilmn12.hg19")
install.packages.auto("FDb.InfiniumMethylation.hg19")
install.packages.auto("TxDb.Hsapiens.UCSC.hg19.knownGene")
install.packages.auto("org.Hs.eg.db")
install.packages.auto("AnnotationDbi")
install.packages.auto("EnsDb.Hsapiens.v86")

# for plotting
install.packages.auto("pheatmap")
install.packages.auto("qqman")
install.packages.auto("forestplot")

# for meta-analysis
install.packages.auto("meta")
install.packages.auto("bacon")

# The actual DNAmArray package
cat("\n* DNAmArray package...\n")
# Also refer to: 
# - https://molepi.github.io/DNAmArray_workflow/index.html
# - https://github.com/molepi/DNAmArray
# - https://github.com/bbmri-nl/BBMRIomics
# For DNAmArray 'pls' is required
# install.packages.auto("pls")
library(devtools)
# install_github("molepi/DNAmArray", force = FALSE)
library(DNAmArray)
# install_github("molepi/omicsPrint", ref = "R3.4", force = FALSE)
# install_github("molepi/omicsPrint", force = FALSE)
library(omicsPrint)
# BBMRIomics requires 'VariantAnnotation'
# install.packages.auto("VariantAnnotation")
library(VariantAnnotation)
install.packages.auto("RPostgreSQL")
# install_github("bbmri-nl/BBMRIomics", subdir = "BBMRIomics", force = FALSE)
library(BBMRIomics)

# Load MultiAssayExperiment
cat("\n* MultiAssayExperiment package...\n")
# if (!require("BiocManager"))
#     install.packages("BiocManager")
# BiocManager::install("MultiAssayExperiment")

library(MultiAssayExperiment)
library(GenomicRanges)
library(SummarizedExperiment)
library(RaggedExperiment)
