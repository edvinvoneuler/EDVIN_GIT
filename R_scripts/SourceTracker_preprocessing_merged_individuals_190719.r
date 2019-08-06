# Read normalized abundance table
otu_kraken <- read.delim(file = "/crex/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_JENA_GORILLA_181009/P7_kraken2_merged-indiviuals_190717/kraken2_otu_table_merged_individuals_1908718.txt", sep="\t", skip=1, row.names=1)
colnames(otu_kraken) <- gsub("kraken2_report","", colnames(otu_kraken))
colnames(otu_kraken) <- gsub("X.","",colnames(otu_kraken))
head(otu_kraken)

# Create DF mapping calculus samples to "sink" for SourceTracker to recognize as data of interest.
# Create Env column for description of sample.
sourcetracker_map <- data.frame(colnames(otu_kraken))
colnames(sourcetracker_map)[1] <- "SampleID"
sourcetracker_map$SourceSink <- "sink"
sourcetracker_map$Env <- "Dental_Calculus_Gorilla"
#sourcetracker_map

# Create corresponding DF of contaminants
##TODO: use uncollapsed files from Jena to make my own set of potential contaminants(?)
sourcetracker_sources <- read.csv(file = "contams/kraken2_otu_table_sources_190301.txt",sep = "\t", skip = 1)
colnames(sourcetracker_sources) <- gsub("kraken2_report","", colnames(sourcetracker_sources))
colnames(sourcetracker_sources) <- gsub("_rmhuman.fastq","", colnames(sourcetracker_sources))
colnames(sourcetracker_sources)[1] <- "OTUID"

sourcetracker_source_map <- data.frame(colnames(sourcetracker_sources))
sourcetracker_source_map <- data.frame(sourcetracker_source_map[-c(1),])
colnames(sourcetracker_source_map)[1] <- "SampleID"
sourcetracker_source_map$SourceSink <- "source"


# Extract Env description from sample ID.
envmap_f <- function(sid) {
  sv <- strsplit(sid, "\\_")[[1]]
  env <- paste(sv[2:length(sv)], collapse = "_")
  return(env)
}
sourcetracker_source_map$Env <- sapply(as.character(sourcetracker_source_map[,1]), envmap_f)

sourcetracker_merged_map <- rbind(sourcetracker_source_map, sourcetracker_map)

# Make OTU_IDs to a column instead of rownames
otu_kraken_norownames <- cbind(rownames(otu_kraken), data.frame(otu_kraken, row.names=NULL))
colnames(otu_kraken_norownames)[1] <- "OTUID"

# Sanity checks
dim(otu_kraken_norownames)
dim(sourcetracker_sources)

# merge kraken OTU table with contaminant OTU table, join on OTUID
sourcetracker_sources$OTUID <- as.factor(sourcetracker_sources$OTUID)

otu_merged <- merge(sourcetracker_sources,otu_kraken_norownames, by="OTUID", all = TRUE)
colnames(otu_merged)[1] <- "#OTU ID"
otu_merged[is.na(otu_merged)] <- 0

#filter out human tax ID and rows/columns that sum to 0
otu_merged_filter <- otu_merged[otu_merged$`#OTU ID` != "9606",]
otu_merged_filter <- otu_merged[,c("#OTU ID",names(which(colSums(otu_merged_filter[,2:ncol(otu_merged_filter)])>0)))]
otu_merged_filter <- otu_merged_filter[rowSums(otu_merged_filter[,c(2:ncol(otu_merged_filter))])>0,]

dim(otu_merged)
setdiff(names(otu_merged),names(otu_merged_filter)) # 0
setdiff(row.names(otu_merged),row.names(otu_merged_filter)) #0
identical(otu_merged_filter,otu_merged) # True
# Nothing removed

cat(paste0(c("# Kraken biom table \n")), file = "Rout_merged_individuals/Jena_Calculus_SourceTracker_otu_allcontams_190719.txt")
write.table(otu_merged, "Rout_merged_individuals/Jena_Calculus_SourceTracker_otu_allcontams_190719.txt", sep = "\t",row.names = FALSE, quote = FALSE, append = TRUE)
write.table(sourcetracker_merged_map, "Rout_merged_individuals/Jena_Calculus_SourceTracker_map_allcontams_190719.txt", sep = "\t",row.names = FALSE, quote = FALSE)




