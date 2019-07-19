
library(decontam)
library(ggplot2)
library(phyloseq)


############## Functions for extracting Taxonomic IDs #################

taxonomy <- read.csv("TaxonID_TaxonName_Kraken2_DB_20190301.csv")
kk.taxonomy <- read.csv("TaxonName_TaxonPath_Kraken2_DB_20190301.csv")

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

kraken_otu <- read.delim("kraken2_otu_table_merged_190703.txt")
row.names(kraken_otu) <- kraken_otu$X.OTUID
colnames(kraken_otu) <- gsub(pattern = "_kraken2_report", replacement = "", x = colnames(kraken_otu))
kraken_otu_decontam <- subset(kraken_otu, select = -1)

# Transform OTU table into Sample-by-Taxa integer matrix 
kraken_otu_decontam <- t(kraken_otu_decontam)
otu_matrix <- as.matrix(kraken_otu_decontam)

# Transform raw read counts into integer Sample-by-Count integer vector 
read_count <- read.delim("raw_read_count.tsv", row.names = 1)
row.names(read_count) <- gsub("_passed_unmapped.fastq","", row.names(read_count))
read_vector <- as.matrix(read_count)

# Rename row to compensate for very weird naming anomaly in the samples
row.names(read_vector)[row.names(read_vector)=='MTM009.A0101'] <- 'MTM009A1'

### This section makes up for the fact that some samples were single end sequences and they
### were merged before kraken2 so in order for the read counts to make sense they must be summed.

# Make vector of single ended samples
single_ends <- subset(read_vector, !(row.names(read_vector) %in% row.names(otu_matrix)))
# Subset the sample without the single ended samples
read_vector <- subset(read_vector, (row.names(read_vector) %in% row.names(otu_matrix)))

# Add single ended read counts into their PE counterparts 
for(i in 1:nrow(single_ends))
{
    j <- 1

    while(TRUE)
    {
        if(grepl(row.names(single_ends)[i], row.names(read_vector)[j]))
        {
            print("MATCH")
            print(paste0(row.names(single_ends)[i],", read count:", single_ends[i]))
            print(paste0(row.names(read_vector)[j],", read count PE: ", read_vector[j]))            
            read_vector[j] = read_vector[j]+single_ends[i]
            print(paste0("Reads after merging:", read_vector[j]))
            cat("\n\n")
            break
        print(j)
        }
    j <- j+1
    }
    
}

dim(read_vector)
dim(otu_matrix)

# Estimate contaminants based on distribution of taxa in samples based on read count compared to 
# the conc. of total DNA  Note that I use total # of reads as an estimate of DNA concentration.

contam_df <- isContaminant(seqtab = otu_matrix, conc = read_vector, method = "frequency")
nrow(contam_df[contam_df$contaminant == TRUE,]) #60
nrow(contam_df[contam_df$contaminant == FALSE,]) # 6264

## TODO: Maybe redo the old analysis with read counts to see
## if the old results with qPCR data shows a large discrepancy with using reads.

# Create list linking taxonomix IDs to Taxonomic Paths.
contam_taxlist <- taxID2tax_LIST_f(row.names(subset(contam_df, contaminant)))

contam_taxlist

read_df <- as.data.frame(read_vector)
plot_frequency(seqtab = otu_matrix, taxa = c("37692"), conc = read_df[,1]) + 
  xlab("Total read count in sample")

