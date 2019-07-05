
library(ggplot2)
library(reshape2)

# Read normalized abundance table
otu_kraken <- read.csv(file = "Rout_post_genome_normalisation_190704.csv", row.names = 1)
colnames(otu_kraken) <- gsub("X","",colnames(otu_kraken))
otu_kraken <- t(otu_kraken)


# Create DF mapping calculus samples to "sink" for SourceTracker to recognize as data of interest.
# Create Env column for description of sample.
sourcetracker_map <- data.frame(colnames(otu_kraken))
colnames(sourcetracker_map)[1] <- "Sample"
sourcetracker_map$SourceSink <- "sink"
sourcetracker_map$Env <- "Dental Calculus Gorilla"
sourcetracker_map

# Create corresponding DF of contaminants
##TODO: use uncollapsed files from Jena to make my own set of potential contaminants.
sourcetracker_sources <- read.csv(file = "contams/kraken2_otu_table_sources_190301.txt",sep = "\t", skip = 1)
colnames(sourcetracker_sources) <- gsub("kraken2_report","", colnames(sourcetracker_sources))
colnames(sourcetracker_sources) <- gsub("_rmhuman.fastq","", colnames(sourcetracker_sources))
colnames(sourcetracker_sources)[1] <- "OTU_ID"

sourcetracker_source_map <- data.frame(colnames(sourcetracker_sources))
sourcetracker_source_map <- data.frame(sourcetracker_source_map[-c(1),])
colnames(sourcetracker_source_map)[1] <- "Sample"
sourcetracker_source_map$SourceSink <- "source"


# Extract Env description from sample ID.
envmap_f <- function(sid) {
  sv <- strsplit(sid, "\\_")[[1]]
  env <- paste(sv[2:length(sv)], collapse = "_")
  return(env)
}
sourcetracker_source_map$Env <- sapply(as.character(sourcetracker_source_map[,1]), envmap_f)

sourcetracker_merged_map <- rbind(sourcetracker_map,sourcetracker_source_map)

# Make OTU_IDs to a column instead of rownames
otu_kraken_norownames <- cbind(rownames(otu_kraken), data.frame(otu_kraken, row.names=NULL))
colnames(otu_kraken_norownames)[1] <- "OTU_ID"


dim(otu_kraken_norownames)
dim(sourcetracker_sources)
# merge kraken OTU table with contaminant OTU table, join on OTU_ID

sourcetracker_sources$OTU_ID <- as.factor(sourcetracker_sources$OTU_ID)

otu_merged <- merge(otu_kraken_norownames,sourcetracker_sources, by="OTU_ID", all = TRUE)
otu_merged[is.na(otu_merged)] <- 0
dim(otu_merged)


