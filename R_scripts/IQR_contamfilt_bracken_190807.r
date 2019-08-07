library(dplyr)
library(tidyverse)
library(ggplot2)

#Read files, remove subset-specific parts of samplenames
kraken_collapsed <- read.delim(file = "kraken2_otu_table_merged-individuals_bracken_190807.txt", skip=1)
colnames(kraken_collapsed) <- gsub(pattern = "kraken2_report_bracken",replacement = "",x = colnames(kraken_collapsed))

kraken_uncollapsed <- read.delim(file = "kraken2_otu_table_uncollapsed_bracken_190807.txt", skip=1)
colnames(kraken_uncollapsed) <- gsub(pattern = "_unmapped_uncollapsed_kraken2_report_bracken", replacement = "",x = colnames(kraken_uncollapsed))

# Subset into two dataframes of shared taxas and their abundances 
shared_taxa_collapsed <- subset(kraken_collapsed, kraken_collapsed$X.OTU.ID %in% kraken_uncollapsed$X.OTU.ID)
shared_taxa_uncollapsed <- subset(kraken_uncollapsed, kraken_uncollapsed$X.OTU.ID %in% kraken_collapsed$X.OTU.ID)

# Sort dataframes to facilitate subtraction, add OTU ID as rownames.
shared_taxa_collapsed <- shared_taxa_collapsed[order(shared_taxa_collapsed$X.OTU.ID),]
row.names(shared_taxa_collapsed) <- shared_taxa_collapsed$X.OTU.ID
shared_taxa_uncollapsed <- shared_taxa_uncollapsed[order(shared_taxa_uncollapsed$X.OTU.ID),]
row.names(shared_taxa_uncollapsed) <- shared_taxa_uncollapsed$X.OTU.ID

# Subtract abundances to obtain distance between collapsed and uncollapsed composition.
vert_dist <- shared_taxa_uncollapsed - shared_taxa_collapsed

# Sanity checks
print("Colnames of both original dataframes in same order:")
all(colnames(shared_taxa_collapsed)==colnames(shared_taxa_uncollapsed))=="TRUE"
print("OTU IDs sum to 0:")
all(vert_dist$X.OTU.ID==0)

# Replace OTU ID column as it was now 0.
vert_dist[,1] <- row.names(vert_dist)

# flatten Abundances, IQR takes a vector
vert_dist_flat <- gather(vert_dist, Sample, Abundance, 2:22)

# Outlier-cutoff with all samples
outlier_cutoff <- quantile(vert_dist_flat$Abundance, 3/4) + 1.5 * IQR(vert_dist_flat$Abundance)

print("All taxa from uncollapsed and collapsed overlap")
print("Q3:")
quantile(vert_dist_flat$Abundance, 3/4)
print("IQR:")
IQR(vert_dist_flat$Abundance)
print("Cutoff:")
outlier_cutoff

contams <- filter(vert_dist_flat, vert_dist_flat$Abundance > outlier_cutoff)
contams_unique <- unique(contams$X.OTU.ID)
print("Number of contaminants")
length(contams_unique)

# Subset orignal table with taxa not in list of contams.
kraken_collapsed_IQRfilt <- subset(kraken_collapsed, !(kraken_collapsed$X.OTU.ID %in% contams_unique))

# Sanity check to make sure a bunch of taxa didn't magically dissapear.
dim(kraken_collapsed)[1]-dim(kraken_collapsed_IQRfilt)[1]==length(contams_unique)

write((paste0("#IQR-outlier filtered kraken table")), "Rout_filter/kraken_bracken_IQR_filtered.tsv")
write.table(kraken_collapsed_IQRfilt, "Rout_filter/kraken_bracken_IQR_filtered.tsv", quote=FALSE ,row.names=FALSE, append=TRUE)
