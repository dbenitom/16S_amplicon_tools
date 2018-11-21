############################################################
### DADA2 - 16 S analysis from paired-end illumina reads ###
############################################################

# Script by: David Benito-Merino (dbenito@mpi-bremen.de).
# Based on the DADA2 tutorial: https://benjjneb.github.io/dada2/tutorial.html
# 
# Script for BACTERIAL 16S V3-V4 amplicons (long overlap of R1 and R2).
# 
# DADA2 works on paired-end Illumina reads in which non-biological sequences (primers, adapters, etc.) have been removed.
# Use the aphros pipeline ampliconNGS until step 1 included (primer clipping) and copy the clipped files in the working directory.

require(dada2)
require(Rcpp)

rm(list=ls())
dir <- "D:/16S_projects/20181109_GuayEnrich_16Stags_Nov2017/Bacteria_dada2/" # path to the directory containing the ptrim files.
setwd(dir)
getwd()

list.files(dir)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
R1_list <- sort(list.files(dir, pattern="_R1.fastq", full.names = TRUE))
R2_list <- sort(list.files(dir, pattern="_R2.fastq", full.names = TRUE))


# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq.
# In our case the samples are named "1", "2", etc. Conserve the fileID_R1 and fileID_R2 to refer to original sample names.
sample.names <- sapply(strsplit(basename(R1_list), "_"), `[`, 1)


# Visualisation of the reads' quality:
quality_plot_R1 <- plotQualityProfile(R1_list[1:length(R1_list)])
quality_plot_R2 <- plotQualityProfile(R2_list[1:length(R2_list)])

# Place length filtered & quality trimmed files in the filtered/ subdirectory
R1_filtered <- file.path(dir, "filtered", paste0(sample.names, "_R1_filt.fastq.gz"))
R2_filtered <- file.path(dir, "filtered", paste0(sample.names, "_R2_filt.fastq.gz"))

# Some filtering parameters (default values shown):
# CHECK HELP FOR MORE PARAMETERS.
# If no R2 files are provided, the files are processed as single-reads.
# compress=T - compress to .fastq.gz
# truncQ=2 - truncate reads that have quality score =/<truncQ.
# truncLen=0 - remove reads shorter than the value. This may be modified based on read qualities and the overlap of R1 and R2. Default is 0.
# trimLeft=0 / trimRight=0 - number of nt to trim from the start/end of each read.
# maxLen=Inf - remove reads longer than defined. Enforced BEFORE trimming & truncation. Recommended for short amplicons.
# minLen=20 - remove reads shorter than value. Enforced AFTER trimming & truncation.
# maxN=0 - after truncation, remove reads with more Ns (generic base) than value.
# minQ=0 - after truncation, remove reads with quality score <minQ.
# maxEE=c(Inf,Inf) - after truncation, remove reads with higher expected errors than value. High values mean more relaxed threshold.
out <- filterAndTrim(R1_list, R1_filtered, R2_list, R2_filtered,
                     compress=TRUE,
                     truncQ=5,
                     truncLen=0, # Why would we truncate the reads? Bad quality?
                     maxLen=500,
                     minLen=100,
                     maxN=0,
                     maxEE=c(2,2),
                     multithread=FALSE) # On Windows set multithread=FALSE

head(out)

# Calclate and plot error rates for each nucleotide substitution (along 100 bp sequence):
R1_error <- learnErrors(R1_filtered, multithread=FALSE)
R2_error <- learnErrors(R2_filtered , multithread=FALSE)

error_plot_R1 <- plotErrors(R1_error, nominalQ=TRUE)
error_plot_R2 <- plotErrors(R2_error, nominalQ=TRUE)

# Dereplication:
R1_dereplicated <- derepFastq(R1_filtered, verbose=TRUE)
R2_dereplicated <- derepFastq(R2_filtered, verbose=TRUE)
# Use the same sample names in the dereplicated data:
names(R1_dereplicated) <- sample.names
names(R2_dereplicated) <- sample.names
# For very large datasets see https://benjjneb.github.io/dada2/bigdata.html

# Sample inference algorithm:
# pool=FALSE - by default samples are processed individually
# pool=TRUE - pooling to increase sensitivity (https://benjjneb.github.io/dada2/pool.html#pooling-for-sample-inference)
# pool="pseudo" - pseudo-pooled where samples are still processed individually (https://benjjneb.github.io/dada2/pseudo.html#Pseudo-pooling)
R1_dada <- dada(R1_dereplicated, err=R1_error, multithread=FALSE, pool=TRUE)
R2_dada <- dada(R2_dereplicated, err=R2_error, multithread=FALSE, pool=TRUE)

R1_dada[[1]]
R2_dada[[1]]

# Merging
## Non-overlapping reads supported with the parameter: mergePairs(..., justConcatenate=TRUE
R1R2_merged <- mergePairs(R1_dada, R1_dereplicated, R2_dada, R2_dereplicated, verbose=TRUE)
# Inspect the merged data.frame from the first sample
head(R1R2_merged[[1]])
# ATTENTION: not many reads were merged, check upstream parameters.
## Solved when filtering parameter was changed from truncLen=c(100,100) to truncLen=0.

# Construct the amplicon sequence variant table (ASV)
# This is a high resolution version of the OTU table.
ASV_table <- makeSequenceTable(R1R2_merged)
dim(ASV_table)
# Inspect distribution of sequence lengths
table(nchar(getSequences(ASV_table)))
## ATTENTION: not so many sequences!!! Check upstream parameters.

# Remove chimaeras
ASV_nochim <- removeBimeraDenovo(ASV_table, method="consensus", multithread=FALSE, verbose=TRUE)
dim(ASV_nochim)
# Proportion of non-chimaeric sequences:
sum(ASV_nochim)/sum(ASV_table)

# Read counts throughout the pipeline:
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(R1_dada, getN), sapply(R2_dada, getN), sapply(R1R2_merged, getN), rowSums(ASV_nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
# ATTENTION: This table looks strange from filtering on!

# Taxonomic classification:
# ATTENTION: download DB from dada2 github
ASV_taxonomy <- assignTaxonomy(ASV_nochim, "D:/Scripts/dada2/silva_nr_v132_train_set.fa.gz", multithread=FALSE)
# Add species:
ASV_taxonomy <- addSpecies(ASV_taxonomy, "D:/Scripts/dada2/silva_species_assignment_v132.fa.gz")

##############################################################################################

# Analyse and plot results with phyloseq.

# Import results to phyloseq:
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

# Extract the sample and ASV names:
samples.out <- rownames(ASV_nochim)
ASVs <- colnames(ASV_nochim)
# ASVs ID table:
ASVs_ID <- cbind(ASVs, paste("asv", c(1:ncol(ASV_nochim)), sep=""))

# rename the ASV to asv#:
colnames(ASV_nochim) <- paste("asv", c(1:ncol(ASV_nochim)), sep="")
rownames(ASV_taxonomy) <- paste("asv", c(1:nrow(ASV_taxonomy)), sep="")
ASV_taxonomy[is.na(ASV_taxonomy[,1])] <- "Unclassified" # Replace empty taxons (domain/kingdom level) with "Unclassified".

# Add sample names:
head (samples.out)
samples.out <- cbind(samples.out, c("Control_70", "Sed_6-10", "Methane_37", "Control_37", "Hexadecane_37",
                                    "Hexadecane_70", "Methane_70", "Hexadecane_50", "Sed_3-5", "Sed_0-2"))
colnames(samples.out) <- c("ID", "Sample")
rownames(samples.out) <- samples.out[,1] # Row names are samples IDs.
samples.out <- as.data.frame(samples.out)

OTU_phyloseq <- otu_table(ASV_nochim, taxa_are_rows = FALSE)
SAMPLE_phyloseq <- sample_data(samples.out)
TAX_phyloseq <- tax_table(ASV_taxonomy)

dada2_phyloseq <- phyloseq(OTU_phyloseq, TAX_phyloseq, SAMPLE_phyloseq) # Create a phyloseq object.

dada2_phyloseq_prune <- prune_samples(sample_names(dada2_phyloseq)!="10", dada2_phyloseq)

plot_bar(dada2_phyloseq_prune, fill="Phylum")

##############################################################################################

# Create a list object with 3 elements:
# Object "X" with 3 elements:
# Element "OTU": matrix otu# vs samples --> contingency table (without total).
# Element "TAX": matrix otu# vs taxonomic path (domain, phylum, class, order, family, genus) (dim = otus x 6)
# Element "ALIGN": matrix otu# vs amplicons seeds. Columns: accnos (sequence/cluster identifier)
#                                                           align (quality of the grouping)
#                                                           path (taxonomic path in ";" separated format)
rownames(ASV_table_3) <- paste("otu", c(1:nrow(ASV_table_3)), sep="")
ASV_table_3 <- ASV_table_3[,3:(ncol(ASV_table_3)-1)]

ASV_taxonomy_3 <- as.data.frame(ASV_taxonomy_3)
ASV_tax_3 <- separate(ASV_taxonomy_3, V3, into=c("domain", "class", "order", "family", "genus"), sep=";")
  sread.csv(textConnection(ASV_taxonomy_3[,3]), sep=";")

X <- list(OTU=ASV_table_3, TAX=ASV_tax_3, ALIGN=ASV_taxonomy_3)

##############################################################################################

# Prepare data for plotting.

# IMPORTANT NOTE:
# ASV_table is a dataframe in which the columns are the OTU sequences and the rows are abundances in each sample.
# ASV_taxonomy is a dataframe in which the rows are the OTU sequences and the columns are the taxonomic paths.
# ASV_table needs to be transposed to be plotted with the PlotAbund.R script by Christiane Hassenrueck. 
# After transposing, add a first column with OTU ID and a second column with amplicon ID (OTU sequence).
require(readr)
require(tidyverse)
ASV_table <- as.data.frame(ASV_table)
OTU_ID <- colnames(ASV_table)
ASV_table_2 <- t(ASV_table) # Transpose the matrix.
total <- rowSums(ASV_table_2)
ASV_table_3 <- cbind(c(1:nrow(ASV_table_2)), OTU_ID, ASV_table_2, total)
rownames(ASV_table_3) <- c()
colnames(ASV_table_3) <- c("OTU", "amplicon", c(1:(ncol(ASV_table_3)-3)), "total")
ASV_table_3 <- as.data.frame(ASV_table_3)
write_tsv(ASV_table_3, "ASV_table_3.tsv")
X <- read.table("ASV_table_3.tsv", h=T, sep="\t")
# ATTENTION: do we need a first column with OTU IDs??
ASV_taxonomy_2 <- ASV_taxonomy[,1:6] # Remove the 2x species columns
colnames(ASV_taxonomy_2) <- c()
rownames(ASV_taxonomy_2) <- c()
ASV_taxonomy_2 <- as.data.frame(ASV_taxonomy_2)
ASV_taxonomy_2 <- unite(ASV_taxonomy_2, sep=";", col="V1", remove=FALSE)
ASV_taxonomy_2 <- ASV_taxonomy_2[,1]
ASV_taxonomy_2 <- gsub("NA;", "", ASV_taxonomy_2)
ASV_taxonomy_2 <- gsub("NA", "", ASV_taxonomy_2)
ASV_taxonomy_3 <- cbind(OTU_ID, rep(100, 135), ASV_taxonomy_2)
ASV_taxonomy_3[ASV_taxonomy_3=="NA"] <- "Unclassified;" # Replace empty taxons with "Unclassified".
colnames(ASV_taxonomy_3) <- c()

save.image(file="20181120_dada2_bacteria.RData")