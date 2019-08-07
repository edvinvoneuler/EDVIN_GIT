library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)

# Read in genome-size normalised table, apply abundance filter
otu_normalized <- read.delim("kraken2_otu_table_merged-individuals_bracken_190807.txt", header = TRUE, row.names = 1, skip = 1) 
row.names(otu_normalized) <- gsub("X", "", row.names(otu_normalized))
colnames(otu_normalized) <- gsub("kraken2_report_bracken", "", colnames(otu_normalized))

# Sum columns to obtain sum of relative abundance per sample
otu_sumRelAbundance <- data.frame("Relative_Abundance_Sum"=colSums(otu_normalized))
# Divide abundances by total abundance in sample to aquire per taxa proportional abundance
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

#filter out human tax ID and rows/columns that sum to 0
dim(otu_abundance_filtered)
otu_abundance_filtered_drop <- otu_abundance_filtered[,colnames(otu_abundance_filtered) != "9606"]
dim(otu_abundance_filtered_drop)

otu_abundance_filtered_drop <- otu_abundance_filtered_drop[rowSums(otu_abundance_filtered_drop)>0,]
dim(otu_abundance_filtered_drop)

otu_abundance_filtered_drop <- otu_abundance_filtered_drop[,colSums(otu_abundance_filtered_drop)>0]
dim(otu_abundance_filtered_drop)

length(setdiff(row.names(otu_abundance_filtered),row.names(otu_abundance_filtered_drop))) # 3854
setdiff(colnames(otu_abundance_filtered),colnames(otu_abundance_filtered_drop)) #0

write.table(otu_abundance_filtered_drop, "Rout_filter/Jena_Calculus_bracken_abundance_filt_190807.txt")
