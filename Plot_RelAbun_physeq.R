# This script works separating the least abundant phyla in a separate category

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