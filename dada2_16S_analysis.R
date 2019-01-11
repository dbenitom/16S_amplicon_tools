############################################################
### DADA2 - 16 S analysis from paired-end illumina reads ###
############################################################

# Script by: David Benito-Merino (dbenito@mpi-bremen.de).
# Based on the DADA2 tutorial: https://benjjneb.github.io/dada2/tutorial.html
# 
# Script for BACTERIAL 16S V3-V4 amplicons (long overlap of R1 and R2).
# V3-V4 amplicon length: 344 bp.
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
quality_plot_R1 <- plotQualityProfile(R1_list[1:length(R1_list)]); quality_plot_R1
quality_plot_R2 <- plotQualityProfile(R2_list[1:length(R2_list)]); quality_plot_R2

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
      # trimLeft=c(x,y) can be used to remove primers of length x (R1) and length y (R2).
# maxLen=Inf - remove reads longer than defined. Enforced BEFORE trimming & truncation. Recommended for short amplicons.
# minLen=20 - remove reads shorter than value. Enforced AFTER trimming & truncation.
# maxN=0 - after truncation, remove reads with more Ns (generic base) than value.
# minQ=0 - after truncation, remove reads with quality score <minQ.
# maxEE=c(Inf,Inf) - after truncation, remove reads with higher expected errors than value. High values mean more relaxed threshold.
out <- filterAndTrim(R1_list, R1_filtered, R2_list, R2_filtered,
                     compress=TRUE,
                     truncQ=5,
                     maxLen=500,
                     minLen=100,
                     maxN=0,
                     maxEE=c(2,2),
                     truncLen=c(250,200), # We truncate the length based on the quality graphs.
                     verbose=TRUE,
                     multithread=FALSE) # On Windows set multithread=FALSE

head(out)

save.image(file="20181120_dada2_bacteria.RData")

# Calculate and plot error rates for each nucleotide substitution (along 100 bp sequence):
R1_error <- learnErrors(R1_filtered, multithread=FALSE)
R2_error <- learnErrors(R2_filtered , multithread=FALSE)

error_plot_R1 <- plotErrors(R1_error, nominalQ=TRUE)
error_plot_R2 <- plotErrors(R2_error, nominalQ=TRUE)

save.image(file="20181120_dada2_bacteria.RData")

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

R1_dada_rep <- dada(R1_dereplicated, err=R1_error, multithread=FALSE) # Without pooling.
R2_dada_rep <- dada(R2_dereplicated, err=R2_error, multithread=FALSE)

R1_dada_rep[[1]]
R2_dada_rep[[1]]

save.image(file="20181120_dada2_bacteria.RData")

# Merging
## Non-overlapping reads supported with the parameter: mergePairs(..., justConcatenate=TRUE
R1R2_merged <- mergePairs(R1_dada, R1_dereplicated, R2_dada, R2_dereplicated, verbose=TRUE)
# Inspect the merged data.frame from the first sample
head(R1R2_merged[[1]])
# Also for non-pooled dada inference:
R1R2_merged_rep <- mergePairs(R1_dada_rep, R1_dereplicated, R2_dada_rep, R2_dereplicated, verbose=TRUE)
# Inspect the merged data.frame from the first sample
head(R1R2_merged_rep[[1]])

# Construct the amplicon sequence variant table (ASV)
# This is a high resolution version of the OTU table.
ASV_table <- makeSequenceTable(R1R2_merged)
dim(ASV_table)
# Inspect distribution of sequence lengths
table(nchar(getSequences(ASV_table)))

# Also for non-pooled data:
ASV_table_rep <- makeSequenceTable(R1R2_merged_rep)
dim(ASV_table_rep)
# Inspect distribution of sequence lengths
table(nchar(getSequences(ASV_table_rep)))

# Remove chimaeras
ASV_nochim <- removeBimeraDenovo(ASV_table, method="consensus", multithread=FALSE, verbose=TRUE)
dim(ASV_nochim)
# Proportion of non-chimaeric sequences:
sum(ASV_nochim)/sum(ASV_table)

# For non-pooled data.
ASV_nochim_rep <- removeBimeraDenovo(ASV_table_rep, method="consensus", multithread=FALSE, verbose=TRUE)
dim(ASV_nochim_rep)
# Proportion of non-chimaeric sequences:
sum(ASV_nochim_rep)/sum(ASV_table_rep)

# Read counts throughout the pipeline:
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(R1_dada, getN), sapply(R2_dada, getN), sapply(R1R2_merged, getN), rowSums(ASV_nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

track_rep <- cbind(out, sapply(R1_dada_rep, getN), sapply(R2_dada_rep, getN), sapply(R1R2_merged_rep, getN), rowSums(ASV_nochim_rep))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track_rep) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_rep) <- sample.names
head(track_rep)

# Taxonomic classification:
# The database in the correct format can be found in the dada2 website.
ASV_taxonomy <- assignTaxonomy(ASV_nochim, "D:/Scripts/dada2/silva_nr_v132_train_set.fa.gz", multithread=FALSE)
# Add species:
ASV_taxonomy <- addSpecies(ASV_taxonomy, "D:/Scripts/dada2/silva_species_assignment_v132.fa.gz")

# Also for non-pooled data:
ASV_taxonomy_rep <- assignTaxonomy(ASV_nochim_rep, "D:/Scripts/dada2/silva_nr_v132_train_set.fa.gz", multithread=FALSE)
# Add species:
ASV_taxonomy_rep <- addSpecies(ASV_taxonomy_rep, "D:/Scripts/dada2/silva_species_assignment_v132.fa.gz")


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

plot_bar(dada2_phyloseq, fill="Phylum", x="sample_Sample")

# Remove archaeal reads: 
dada2_phyloseq_prune <- prune_taxa((tax_table(dada2_phyloseq))=="Bacteria", dada2_phyloseq)

# Obtain only the 10 most abundant taxa:
f1 <- filterfun_sample(topk(10))
abundant_list <- genefilter_sample(dada2_phyloseq, f1)
dada2_phyloseq_10 <- prune_taxa(abundant_list, dada2_phyloseq)
dim(tax_table(dada2_phyloseq_10))

# Remove Archaea, Eukarya and NA from tax_table
dada2_phyloseq <- subset_taxa(dada2_phyloseq, Kingdom=="Bacteria")

# Obtain the 10 most abundant taxa per sample:
RelAbun_ASV <- otu_table(dada2_phyloseq) / rowSums(otu_table(dada2_phyloseq)) # relative abundance of each ASV per sample.
colnames(RelAbun_ASV)
rownames(RelAbun_ASV)
RelAbun_ASV <- as.data.frame(RelAbun_ASV)


k <- 5 # number of taxa to obtain
mx <- t(apply(RelAbun_ASV, 1, 
              function(x) names(Rel_abun_ASV)[sort(head(order(x, decreasing=TRUE), k))]
              ))

list_mx <- unique(as.character(mx))
all_ASV <- colnames(otu_table(dada2_phyloseq))
my_ASV <- all_ASV[!(all_ASV %in% list_mx)] # What is this?
dada2_phyloseq_prune <- prune_taxa(my_ASV, dada2_phyloseq)

phylumGlommed <-  tax_glom(dada2_phyloseq_trial, "Phylum")
genusGlommed <-  tax_glom(dada2_phyloseq_trial, "Genus")
orderGlommed <- tax_glom(dada2_phyloseq_trial, "Order")

plot_bar(orderGlommed, fill="Order", x="sample_Sample") +
  theme_classic() +
  theme(legend.position = "none")+ 
  geom_bar(position="fill")

plot_bar(dada2_phyloseq_trial, fill="Phylum", x="sample_Sample")
plot_bar(phylumGlommed, fill="Phylum", x="sample_Sample")

# To classify the NAs as a higher known taxon:
# Alternatively
NameTax <- function(x, ind){
  if(is.na(x[ind])){
    x[ind] <- x[ind]
  } else {
    if(ind==1){x[ind] <- paste("d", x[ind], sep="_")} else{                  # Domain
      if(ind==2){x[ind] <- paste("p", x[ind], sep="_")} else{                # Phylum
        if(ind==3){x[ind] <- paste("c", x[ind], sep="_")} else{              # Class
          if(ind==4){x[ind] <- paste("o", x[ind], sep="_")} else{            # Order
            if(ind==5){x[ind] <- paste("f", x[ind], sep="_")} else{          # Family
              if(ind==6){x[ind] <- paste("g", x[ind], sep="_")} else{        # Genus
                if(ind==7){x[ind] <- paste("s", x[ind-1], x[ind], sep="_")}  # Species
              }
            }
          }
        }
      }
    }
  }
}

tax.tab <- data.frame(tax_table(dada2_phyloseq_trial))

for (i in 1:7) {
  tax_table(dada2_phyloseq_trial)[,i] <- apply(tax.tab, 1, NameTax, ind=i)
}

ModifyTax <- function(x,ind){
  #   xth row in the dataframe
  #   ind taxonomy level to change
  if(is.na(x[ind])){
    nonNa <- which(!is.na(x[-ind])) # which taxa are not NA excepting the one we're interested in.
    maxNonNa <- max(nonNa)
    x[ind] <- x[maxNonNa]
  }else{x[ind] <- x[ind]}
}

for (i in 1:7) {
  tax_table(dada2_phyloseq_trial)[,i] <- apply(tax.tab,1,ModifyTax,ind=i)
}

# Transform to relative abundance and filter out ASVs with low variances:
phylumGlommed <-  tax_glom(dada2_phyloseq_trial, "Phylum")
physeq_RA <- transform_sample_counts(phylumGlommed, function(x) x/sum(x))
physeq_F <- filter_taxa(physeq_RA, function(x) var(x) > 1e-05, TRUE) # Keep OTUs with variance greater than 0.0001 

plot_bar(physeq_RA, fill="Phylum", x="sample_Sample")


###################
#Transform into relative abundances + filter to a mean threshold
physeq2 = filter_taxa(dada2_phyloseq_trial, function(x) mean(x) > 0.1, TRUE)
physeq3 = transform_sample_counts(physeq2, function(x) x / sum(x) )

# create dataframe from phyloseq object for stomach samples... to graph samples without condensed phyla
stom <- psmelt(physeq3)

#Condense low abundance taxa into an "Other" category for the barplot
# Turn all OTUs into phylum counts
glom <- tax_glom(physeq3, taxrank ="Order")
glom # should list # taxa as # phyla
data <- psmelt(glom) # create dataframe from phyloseq object
data$Phylum <- as.character(data$Phylum) #convert to character

#simple way to rename phyla with < 1% abundance
data$Phylum[data$Abundance < 0.01] <- "< 1% abund."
# I also tried using this method, based on median abundance threshold... either way, I still end up with sample that do not reach 100% on the bar graph
# group dataframe by Phylum, calculate median rel. abundance

medians <- ddply(data, ~Phylum, function(x) c(median=median(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%

remainder <- medians[medians$median <= 0.01,]$Phylum
# list of low abundance phyla

remainder
# change their name to "Phyla < 1% abund."

data[data$Phylum %in% remainder,]$Phylum <- "Phyla < 1% abund."
#rename phyla with < 1% relative abundance
data$Phylum[data$Abundance < 0.01] <- "Phyla < 1% abund."

#plot with condensed phyla into "< 1% abund" category
p <- ggplot(data=data, aes(x=Sample, y=Abundance, fill=Phylum))
p + geom_bar(aes(), stat="identity", position="stack", color="black")

g <- transform_sample_counts(dada2_phyloseq_trial, function(x) sort(x, decreasing=T))
otu_table(g)

x1 <- tax_glom(dada2_phyloseq_trial, "Genus" )#glom to genus-level
TopNOTUs <- names(sort(taxa_sums(x1), TRUE)[1:20]) #generate list of top 20 OTUs
x3 = transform_sample_counts(x1, function(x) x/sum(x)) #generate rel abundance by group
x4 <- prune_taxa(TopNOTUs, x3) #prune less abundant taxa

plot_bar(x4, fill="Genus", x="sample_Sample") +
  theme(legend.position = "bottom", legend.text = element_text(size=5))

plot_bar(genusGlommed, fill="Genus", x="sample_Sample") +
  theme(legend.position = "none", legend.text = element_text(size=5))

save.image(file="20181120_dada2_bacteria.RData")