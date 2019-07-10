
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)

# microbial analyses
# Read in genome-size normalised table, apply abundance filter
otu_normalized <- read.csv("Rout/Rout_post_genome_normalisation_190704.csv", header = TRUE, row.names = 1)
colnames(otu_normalized) <- gsub("X","", colnames(otu_normalized))

# Sum rows to obtain sum of relative abundance
otu_sumRelAbundance <- data.frame("Relative_Abundance_Sum"=rowSums(otu_normalized))
# Divide relative abundances by total abundance in sample to aquire per taxa proportional abundance
otu_normalized_propAbundance <- otu_normalized/otu_sumRelAbundance$Relative_Abundance_Sum

# Define Filtering function and cutoff value
prop_cutoff <- 0.0003
filter_taxa <- function(x) {
  if (x > prop_cutoff) {
    
    return(x)
  } else {
    return(0)
  }
}

# Filter out taxa with less than 0.03% relative abundance in sample.
otu_abundance_filtered <- apply(otu_normalized_propAbundance,c(1,2),filter_taxa)
# Transform back from relative abundance
otu_abundance_filtered <- otu_abundance_filtered*otu_sumRelAbundance$Relative_Abundance_Sum

# Add Pseudocount of (1 read)/(mean genome size) to avoid log(0)
genome_sizes = read.table(file = "Rout/average_genome_sizes.txt")
mgs <- mean(genome_sizes[,2]) #mean genome size
mgs
otu_abundance_filt_pcount <- apply(otu_abundance_filtered,c(1,2),function(x){x+1/mgs})
# Appy CLR
otu_CLR <- t(apply(otu_abundance_filt_pcount, 1, function(x){log(x)-mean(log(x))}))

head(otu_abundance_filt_pcount)

write.csv(otu_abundance_filt_pcount, "Rout/Jena_Calc_CLR_abundance_filt_190708")




