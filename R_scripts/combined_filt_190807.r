IQR_filt <- read.table("Rout_filter/kraken_bracken_IQR_filtered.tsv", header = TRUE, row.names = 1)
abundance_filt <- read.table("Rout_filter/Jena_Calculus_bracken_abundance_filt_190807.txt", header = TRUE)
row.names(abundance_filt) <- gsub("X","",row.names(abundance_filt))
decontam_filt <- read.table("Rout_filter/Jena_Calculus_bracken_decontam_filt_190807.tsv", header = TRUE, row.names = 2)
decontam_filt <- decontam_filt[-1]

tax_ids_filt_1 <- row.names(abundance_filt)
tax_ids_filt_2 <- tax_ids_filt_1[tax_ids_filt_1 %in% row.names(IQR_filt)] 
tax_ids_filt_3 <- tax_ids_filt_2[tax_ids_filt_2 %in% row.names(decontam_filt)]

length(tax_ids_filt_3) # 1813 left

# Read normalized abundance table
otu_kraken <- read.delim(file = "kraken2_otu_table_merged-individuals_bracken_190807.txt",sep = "\t", skip = 1,header = TRUE)
#otu_kraken <- read.delim(file = "/crex/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_JENA_GORILLA_181009/P7_kraken2_merged-indiviuals_190717/kraken2_otu_table_merged_individuals_1908718.txt",sep = "\t", skip = 1,header = TRUE)
colnames(otu_kraken) <- gsub("kraken2_report","", colnames(otu_kraken))
colnames(otu_kraken)[1] <- "OTUID"
head(otu_kraken)

# Subset the kraken table with the tax ids left after filtering
otu_kraken_filtererd <- otu_kraken[otu_kraken$OTUID%in%tax_ids_filt_3,]

write.table(otu_kraken_filtererd, "Rout_filter/kraken_otu_all_filtered_merged-individuals.tsv", quote = FALSE, row.names = FALSE)