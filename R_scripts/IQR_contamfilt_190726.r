
library(dplyr)
library(tidyverse)
library(ggplot2)

#Read files, remove subset-specific parts of samplenames
kraken_collapsed <- read.delim(file = "/crex/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_JENA_GORILLA_181009/P7_kraken2_merged-indiviuals_190717/kraken2_otu_table_merged_individuals_1908718.txt", skip=1)
colnames(kraken_collapsed) <- gsub(pattern = "kraken2_report",replacement = "",x = colnames(kraken_collapsed))

kraken_uncollapsed <- read.delim(file = "/crex/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_JENA_GORILLA_181009/P7_kraken2_uncollapsed_190722/kraken2_otu_table_uncollapsed_190722.txt", skip=1)
colnames(kraken_uncollapsed) <- gsub(pattern = "_unmapped_uncollapsed_kraken2_report", replacement = "",x = colnames(kraken_uncollapsed))

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
vert_dist_flat <- gather(vert_dist, Taxa, Abundance, 2:22)
IQR_cutoff <- IQR(vert_dist_flat$Abundance)

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

# Outlier-cutoff without zero abundance samples
vert_dist_flat_nozeros <- subset(vert_dist_flat, vert_dist_flat$Abundance>0)
outlier_cutoff_nozeros <- quantile(vert_dist_flat_nozeros$Abundance, 3/4) + 1.5 * IQR(vert_dist_flat_nozeros$Abundance)


print("No zero-distance taxa")
print("Q3:")
quantile(vert_dist_flat_nozeros$Abundance, 3/4)
print("IQR:")
IQR(vert_dist_flat_nozeros$Abundance)
print("Cutoff:")
outlier_cutoff_nozeros

contams_nozeros <- filter(vert_dist_flat_nozeros, vert_dist_flat_nozeros$Abundance > outlier_cutoff_nozeros)
contams_unique_nozeros <- unique(contams_nozeros$X.OTU.ID)
print("Number of contaminants:")
length(contams_unique_nozeros)

# Subset orignal table with taxa not in list of contams.
kraken_collapsed_IQRfilt <- subset(kraken_collapsed, !(kraken_collapsed$X.OTU.ID %in% contams_unique_nozeros))

# Sanity check to make sure a bunch of taxa didn't magically dissapear.
dim(kraken_collapsed)[1]-dim(kraken_collapsed_IQRfilt)[1]==length(contams_unique_nozeros)

write((paste0("#IQR-outlier filtered kraken table")), "/crex/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_JENA_GORILLA_181009/SCRIPTS/EDVIN_GIT/R_scripts/Rout_filter/kraken_IQR_filtered.tsv")
write.table(kraken_collapsed_IQRfilt, "/crex/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_JENA_GORILLA_181009/SCRIPTS/EDVIN_GIT/R_scripts/Rout_filter/kraken_IQR_filtered.tsv", quote=FALSE ,row.names=FALSE, append=TRUE)
