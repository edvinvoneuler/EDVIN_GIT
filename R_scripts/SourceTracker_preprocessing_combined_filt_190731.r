IQR_filt <- read.table("Rout_filter/kraken_IQR_filtered.tsv", header = TRUE, row.names = 1)
abundance_filt <- read.table("Rout_filter/Jena_Calculus_abundance_filt_1900730.tsv", header = TRUE)
abundance_filt <- as.data.frame(t(abundance_filt))
row.names(abundance_filt) <- gsub("X","",row.names(abundance_filt))
decontam_filt <- read.table("Rout_filter/Jena_Calculus_decontam_filt_1900731.tsv", header = TRUE, row.names = 2)
decontam_filt <- decontam_filt[-1]

tax_ids_filt_1 <- row.names(abundance_filt)
# Not negating the IQR filtering only left me with 100 taxa, all assigned to "unknown" by sourcetracker. 
tax_ids_filt_2 <- tax_ids_filt_1[!(tax_ids_filt_1 %in% row.names(IQR_filt))] 
## TODO: Check IQR filtering did I mess up and accidentally output contaminants?
tax_ids_filt_3 <- tax_ids_filt_2[tax_ids_filt_2 %in% row.names(decontam_filt)]

length(tax_ids_filt_3) # 1251 left

# Read normalized abundance table
otu_kraken <- read.delim(file = "/crex/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_JENA_GORILLA_181009/P7_kraken2_merged-indiviuals_190717/kraken2_otu_table_merged_individuals_1908718.txt",sep = "\t", skip = 1,header = TRUE)
colnames(otu_kraken) <- gsub("kraken2_report","", colnames(otu_kraken))
colnames(otu_kraken)[1] <- "OTUID"
head(otu_kraken)

# Subset the kraken table with the tax ids left after filtering
otu_kraken <- otu_kraken[otu_kraken$OTUID%in%tax_ids_filt_3,]

# Create DF mapping calculus samples to "sink" for SourceTracker to recognize as data of interest.
# Create Env column for description of sample.
sourcetracker_map <- data.frame(colnames(otu_kraken))
# Remove OTUID from rows
sourcetracker_map <- as.data.frame(sourcetracker_map[2:nrow(sourcetracker_map),])
colnames(sourcetracker_map)[1] <- "SampleID"
sourcetracker_map$SourceSink <- "sink"
sourcetracker_map$Env <- "Dental_Calculus_Gorilla"
#sourcetracker_map

# Create corresponding DF of contaminants
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

# Sanity checks
dim(otu_kraken)
dim(sourcetracker_sources)

# merge kraken OTU table with contaminant OTU table, join on OTUID
#sourcetracker_sources$OTUID <- as.factor(sourcetracker_sources$OTUID)

otu_merged <- merge(sourcetracker_sources,otu_kraken, by="OTUID", all = TRUE)
colnames(otu_merged)[1] <- "#OTU ID"
otu_merged[is.na(otu_merged)] <- 0

cat(paste0(c("# Kraken biom table IQR, abundance and decontam filtered \n")), file = "Rout_filter/Jena_Calculus_combined_filt_SourceTracker_otu_allcontams_190731.txt")
write.table(otu_merged, "Rout_filter/Jena_Calculus_combined_filt_SourceTracker_otu_allcontams_190731.txt", sep = "\t",row.names = FALSE, quote = FALSE, append = TRUE)
write.table(sourcetracker_merged_map, "Rout_filter/Jena_Calculus_combined_filt_SourceTracker_map_allcontams_190731.txt", sep = "\t",row.names = FALSE, quote = FALSE)
