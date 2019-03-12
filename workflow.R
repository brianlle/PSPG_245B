#R3.5.2
#Adaped from Dr. Bin Chen's drug repositioning pipeline
#This is a simplified version of a pipeline used to predict drug hits for a given disease
#Depending on your context, you can also compare different phenotypes, responders vs. non-responders, mutation vs. wild type, etc

setwd("~/PSPG_245_Test/")

#install packages needed to run the code
source("http://www.bioconductor.org/biocLite.R")
biocLite("impute")
biocLite("siggenes")
biocLite("RankProd")
biocLite("preprocessCore")
biocLite("GEOquery")
biocLite("qvalue")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db", version = "3.8")

install.packages(c("pheatmap", "gplots", "ggplot2", "RColorBrewer", "plyr"))

###############################
#parameters
disease <- "breast_cancer" #no space

#method to compute disease signatures. by default, we use SAM.
method_id <- 3 #2: rankprod, 3: siggenes sam, 

q_thresh <- 0.05 #fdr; if there are no differentially expressed genes, can try loosening the threshold
fold_change <- 1 #fold change cutoff

#disease signature file
dz_sig.file <- paste(disease, "/dz_signature_cmap.txt", sep="")
###############################

#create a folder for the disease
if (!dir.exists(disease)) dir.create(disease)

#create disease signatures
#need to identify a dataset used to create disease signatures
#in the real case, you need to find a more robust to validate disease signatures before making drug predictions.

#update location of the disease dataset
disease_data_filepath <- "data/HiSeqV2"

#This code snippet takes in the disease gene expression data and generates a differential gene expression signature
#Also outputs a heatmap visualization of the separation of cases vs. controls (stored in disease folder)
source("code/create_dz_signature_from_TCGA.R")


#Predict drugs using connectivity map data with generated disease signature
cmd <- paste("Rscript code/predict_drugs_cmap.R", disease, "cmap", paste(disease, "/dz_signature_cmap.txt", sep=""))
system(cmd)

#Analyze predictions: apply filtering to drug hits, generate visualization of best hits and worst hits
source("code/drug_repurpose.R")

#all the results will be stored in the disease folder
