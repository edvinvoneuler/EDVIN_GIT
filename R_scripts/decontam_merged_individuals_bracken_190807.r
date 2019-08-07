
library(decontam)
library(ggplot2)
library(phyloseq)


############## Functions for extracting Taxonomic IDs #################

taxonomy <- read.csv("taxa_tables/TaxonID_TaxonName_Kraken2_DB_20190301.csv")
kk.taxonomy <- read.csv("taxa_tables/TaxonName_TaxonPath_Kraken2_DB_20190301.csv")

# current plan: use taxonomy to get scientific name from taxon ID, then use kk.taxonomy to get full name by grep of name
# if genus, probably multiple hits, but first hit should be correct one. Except where there are stupid names
taxID2tax_f <- function(id) {
  taxname <- as.character(taxonomy[which(taxonomy$Taxon_ID == id), "Taxon_name"])
  if (length(taxname) > 0) {
    taxpath <- grep(paste0("__",taxname), kk.taxonomy$V1, value = TRUE)
    taxpathN <- length(taxpath)
    taxpath <- taxpath[1] #first hit
    lowest.level <- regmatches(as.character(taxpath),regexpr("\\|[^\\|]*$",as.character(taxpath)))
    level <- substr(lowest.level, 2, 2)
  } else {
    taxname <- "NA"
    taxpath <- 0
    taxpathN <- 0
    level <- "NA"
  }
  return(c(id, taxname, taxpath, taxpathN, level))
}

taxID2tax_LIST_f <- function(idlist) {
  tax.df <- data.frame(do.call(rbind.data.frame, sapply(idlist, taxID2tax_f)))
  colnames(tax.df) <- c("Taxon_ID", "Taxon_Name","Taxon_Path","No.Path.hits", "Taxon_Level")
  tax.df$Taxon_ID <- as.integer(levels(tax.df$Taxon_ID))[tax.df$Taxon_ID]
  tax.df$Taxon_Name <- as.character(tax.df$Taxon_Name)
  tax.df$Taxon_Path <- as.character(tax.df$Taxon_Path)
  tax.df$No.Path.hits <- as.integer(levels(tax.df$No.Path.hits))[tax.df$No.Path.hits]
  tax.df$Taxon_Level <- as.character(tax.df$Taxon_Level)
  return(tax.df)
}

kraken_otu <- read.delim("kraken2_otu_table_merged-individuals_bracken_190807.txt", skip = 1)
row.names(kraken_otu) <- kraken_otu$X.OTUID
colnames(kraken_otu) <- gsub(pattern = "kraken2_report_bracken", replacement = "", x = colnames(kraken_otu))
kraken_otu_decontam <- subset(kraken_otu, select = -1)

# Transform OTU table into Sample-by-Taxa integer matrix 
kraken_otu_decontam <- t(kraken_otu_decontam)
otu_matrix <- as.matrix(kraken_otu_decontam)

# Transform raw read counts into integer Sample-by-Count integer vector 
read_count <-as.data.frame(c(20843383, 4675633, 11582643, 17702434, 23121319, 15570486, 2195920, 30143901, 14600405, 2902850, 5926113, 4485752, 11232681, 18575472, 3082577, 8868720, 6377647,15657295,481956,6837100,13406134))
row.names(read_count) <- colnames(kraken_otu)[2:ncol(kraken_otu)]
colnames(read_count) <- "Reads"

# Estimate contaminants based on distribution of taxa in samples based on read count compared to 
# the conc. of total DNA  Note that I use total # of reads as an estimate of DNA concentration.

contam_df <- isContaminant(seqtab = otu_matrix, conc = read_count$Reads, method = "frequency")
nrow(contam_df[contam_df$contaminant == TRUE,]) #29
nrow(contam_df[contam_df$contaminant == FALSE,]) #5662

## TODO: Maybe redo the old analysis with read counts to see
## if the old results with qPCR data shows a large discrepancy with using reads.

# Create list linking taxonomix IDs to Taxonomic Paths.
contam_taxlist <- taxID2tax_LIST_f(subset(kraken_otu$X.OTU.ID, contam_df$contaminant))

png("figures/Rplot_allcontams_bracken_decontam.png", units = "in", width = 11.7, height = 8.3, res = 600)
plot_frequency(seqtab = otu_matrix, taxa = row.names(contam_df[contam_df$contaminant == TRUE,]), conc = read_count$Reads) + 
  xlab("Total read count in sample")
dev.off()

#Remove contaminants from original table
otu_filtered <- kraken_otu[!contam_df$contaminant,]
length(setdiff(kraken_otu$X.OTU.ID,otu_filtered$X.OTU.ID)) #Sanity Check (=29)
write.table(otu_filtered, quote = FALSE, file = "Rout_filter/Jena_Calculus_bracken_decontam_filt_190807.tsv")
