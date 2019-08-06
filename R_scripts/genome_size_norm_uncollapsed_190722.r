library(dplyr)
library(reshape2)
library(tidyr)

# Functions

taxonomy <- read.csv("/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_FIRSTSCREEN_181220/R_analyses/TaxonID_TaxonName_Kraken2_DB_20190301.csv")
kk.taxonomy <- read.csv("/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_FIRSTSCREEN_181220/R_analyses/TaxonName_TaxonPath_Kraken2_DB_20190301.csv")

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


kraken2_otu <- read.delim(file="/home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/P7_kraken2_uncollapsed_190722/kraken2_otu_table_uncollapsed_190722.txt", sep="\t", skip=1)
colnames(kraken2_otu) <- append("Taxon_ID",gsub("\\_.*","",colnames(kraken2_otu[,2:ncol(kraken2_otu)])))


##### Genome size normalisation #####
#prokaryotes <- read.delim(file = "/proj/sllstore2017021/nobackup/JAELLE/REFERENCES/genome_reports_mod190215/prokaryotes.txt")
prokaryotes <- read.delim(file = "taxa_tables/prokaryotes.txt")
unique(prokaryotes$Status)
colnames(prokaryotes)
prok.sizes <- dplyr::select(prokaryotes, X.Organism.Name, TaxID, Status, Size..Mb.)
prok.sizes <- subset(prok.sizes, Status == "Complete Genome" | Status == "Complete")

kraken2_otu <- gather(kraken2_otu,"Sample",value="Raw_Abundance",-1)

average.sizes <- aggregate(prok.sizes[,4], list(prok.sizes$TaxID), mean)
colnames(average.sizes) <- c("Taxon_ID", "Average_Genome_Size")

for (i in 1:length(kraken2_otu$Taxon_ID)){
  if (kraken2_otu$Taxon_ID[i] %in% average.sizes$Taxon_ID){
    kraken2_otu$Average_Genome_Size[i] <- average.sizes[which(average.sizes$Taxon_ID == kraken2_otu$Taxon_ID[i]),"Average_Genome_Size"]}
  else{kraken2_otu$Average_Genome_Size[i] <- NA}
}


kraken2.nodata <- subset(kraken2_otu, is.na(Average_Genome_Size))
kraken2.nodata <- as.data.frame(unique(kraken2.nodata$Taxon_ID))
write.csv(kraken2.nodata, "tax_reports_uncollapsed/190722_no_data.csv", quote = FALSE, row.names = FALSE)
kraken2.nodata <- read.delim("tax_reports_uncollapsed/tax_report_prok.txt")
kraken2.nodata <- dplyr::select(kraken2.nodata, taxid, taxname)
kraken2.nodata <- kraken2.nodata[1:(length(kraken2.nodata$taxid)-2),]
length(kraken2.nodata$taxid)

##average sizes based on genus name##

prok.sizes.name <- prok.sizes
prok.sizes.name$X.Organism.Name <- gsub(" .*$", "", prok.sizes.name$X.Organism.Name)
prok.sizes.name$X.Organism.Name <- gsub("[\\[,]", "", prok.sizes.name$X.Organism.Name)
prok.sizes.name$X.Organism.Name <- gsub("[\\],]", "", prok.sizes.name$X.Organism.Name)
prok.sizes.name$X.Organism.Name <- gsub("[\\',]", "", prok.sizes.name$X.Organism.Name)

average.sizes.name <- aggregate(prok.sizes[, 4],list(prok.sizes.name$X.Organism.Name), mean)
colnames(average.sizes.name) <- c("Genus", "Average_Genome_Size")

kraken2.nodata$taxname <- gsub(" .*$", "", kraken2.nodata$taxname)
kraken2.nodata$taxname <- gsub("[\\[,]", "", kraken2.nodata$taxname)
kraken2.nodata$taxname <- gsub("[\\],]", "", kraken2.nodata$taxname)
kraken2.nodata$taxname <- gsub("[\\',]", "", kraken2.nodata$taxname)


for (i in 1:length(kraken2.nodata$taxname)){
  if (kraken2.nodata$taxname[i] %in% average.sizes.name$Genus){
    kraken2.nodata$Average_Genome_Size[i] <- average.sizes.name[which(average.sizes.name$Genus == kraken2.nodata$taxname[i]),"Average_Genome_Size"]}
  else{kraken2.nodata$Average_Genome_Size[i] <- NA}
}

sum(is.na(kraken2.nodata$Average_Genome_Size)) #159

for (i in 1:length(kraken2_otu$Taxon_ID)){
  if (kraken2_otu$Taxon_ID[i] %in% kraken2.nodata$taxid){
    kraken2_otu$Average_Genome_Size[i] <- kraken2.nodata[which(kraken2.nodata$taxid == kraken2_otu$Taxon_ID[i]),"Average_Genome_Size"]}
  else{kraken2_otu$Average_Genome_Size[i] <- kraken2_otu$Average_Genome_Size[i]}
}

kraken2_otu$Ratio <- kraken2_otu$Raw_Abundance/kraken2_otu$Average_Genome_Size


kraken2.nodata2 <- subset(kraken2.nodata, is.na(Average_Genome_Size))
length(kraken2.nodata2$taxid) #159

viruses <- read.delim(file = "taxa_tables/viruses.txt")

unique(viruses$Status)
colnames(viruses)
virus.sizes <- dplyr::select(viruses, X.Organism.Name, TaxID, Status, Size..Kb.)
virus.sizes <- subset(virus.sizes, Status == "Complete Genome" | Status == "Complete")
virus.sizes$Size..Kb. <- virus.sizes$Size..Kb.*0.001
colnames(virus.sizes)[colnames(virus.sizes)=="Size..Kb."] <- "Size..Mb."

average.virus.sizes <- aggregate(virus.sizes[,4], list(virus.sizes$TaxID), mean)
colnames(average.virus.sizes) <- c("Taxon_ID", "Average_Genome_Size")

for (i in 1:length(kraken2_otu$Taxon_ID)){
  if (kraken2_otu$Taxon_ID[i] %in% average.virus.sizes$Taxon_ID){
    kraken2_otu$Average_Genome_Size[i] <- average.virus.sizes[which(average.virus.sizes$Taxon_ID == kraken2_otu$Taxon_ID[i]),"Average_Genome_Size"]}
  else{kraken2_otu$Average_Genome_Size[i] <- kraken2_otu$Average_Genome_Size[i]}
}

kraken2.nodata3 <- subset(kraken2_otu, is.na(Average_Genome_Size))
kraken2.nodata3 <- dplyr::select(kraken2.nodata3, Taxon_ID, Average_Genome_Size)
kraken2.nodata3 <- unique( kraken2.nodata3[ , 1:2 ] )
length(kraken2.nodata3$Taxon_ID) #133

write.csv(kraken2.nodata3$Taxon_ID, "tax_reports_uncollapsed/190722_no_data_viruses.csv", quote = FALSE, row.names = FALSE)
kraken2.nodata3 <- read.delim("tax_reports_uncollapsed/taxreports_no_data_viruses.txt")
kraken2.nodata3 <- dplyr::select(kraken2.nodata3, taxid, taxname)
kraken2.nodata3 <- kraken2.nodata3[1:(length(kraken2.nodata3$taxid)-2),]
length(kraken2.nodata3$taxid) #133

for (i in 1:length(kraken2.nodata3$taxname)){
  if (kraken2.nodata3$taxname[i] %in% average.sizes.name$Genus){
    kraken2.nodata3$Average_Genome_Size[i] <- average.sizes.name[which(average.sizes.name$Genus == kraken2.nodata3$taxname[i]),"Average_Genome_Size"]}
  else{kraken2.nodata3$Average_Genome_Size[i] <- NA}
}

for (i in 1:length(kraken2_otu$Taxon_ID)){
  if (kraken2_otu$Taxon_ID[i] %in% kraken2.nodata3$taxid){
    kraken2_otu$Average_Genome_Size[i] <- kraken2.nodata3[which(kraken2.nodata3$taxid == kraken2_otu$Taxon_ID[i]),"Average_Genome_Size"]}
  else{kraken2_otu$Average_Genome_Size[i] <- kraken2_otu$Average_Genome_Size[i]}
}

kraken2.NA <- kraken2_otu[which(is.na(kraken2_otu$Average_Genome_Size)),]
length(unique(kraken2.NA$Taxon_ID)) #133 taxa for which I could not obtain average genome size.
kraken2.NA.taxa <- as.data.frame(unique(kraken2.NA$Taxon_ID))

##As it stands, I've removed all prokaryotes and viruses for which complete genomes exist.
##In the case of genus level, or strains without sequenced complete genomes, I assign them the average of the genome sizes of that genus.
##For the 133 taxa that remain, no complete genome of a taxon in their genus exists.



###REST###

prok.sizes.rest <- dplyr::select(prokaryotes, X.Organism.Name, TaxID, Status, Size..Mb.)
prok.sizes.rest <- subset(prok.sizes.rest, !(Status == "Complete Genome" | Status == "Complete"))

average.sizes.rest <- aggregate(prok.sizes.rest[,4], list(prok.sizes.rest$TaxID), mean)
colnames(average.sizes.rest) <- c("Taxon_ID", "Average_Genome_Size")

for (i in 1:length(kraken2_otu$Taxon_ID)){
  if (kraken2_otu$Taxon_ID[i] %in% average.sizes.rest$Taxon_ID){
    kraken2_otu$Average_Genome_Size[i] <- average.sizes.rest[which(average.sizes.rest$Taxon_ID == kraken2_otu$Taxon_ID[i]),"Average_Genome_Size"]}
  else{kraken2_otu$Average_Genome_Size[i] <- kraken2_otu$Average_Genome_Size[i]}
}

kraken2.NA.2 <- kraken2_otu[which(is.na(kraken2_otu$Average_Genome_Size)),]
length(unique(kraken2.NA.2$Taxon_ID)) #73

kraken2.nodata4 <- subset(kraken2_otu, is.na(Average_Genome_Size))
kraken2.nodata4 <- as.data.frame(unique(kraken2.nodata4$Taxon_ID))
write.csv(kraken2.nodata4, "tax_reports_uncollapsed/190722_no_data_rest.csv", quote = FALSE, row.names = FALSE)
kraken2.nodata4 <- read.delim("tax_reports_uncollapsed/tax_report_rest.txt")
kraken2.nodata4 <- dplyr::select(kraken2.nodata4, taxid, taxname)
kraken2.nodata4 <- kraken2.nodata4[1:(length(kraken2.nodata4$taxid)-2),]
length(kraken2.nodata4$taxid) #73

##average sizes based on genus name##

prok.sizes.name.rest <- prok.sizes.rest
prok.sizes.name.rest$X.Organism.Name <- gsub(" .*$", "", prok.sizes.name.rest$X.Organism.Name)
prok.sizes.name.rest$X.Organism.Name <- gsub("[\\[,]", "", prok.sizes.name.rest$X.Organism.Name)
prok.sizes.name.rest$X.Organism.Name <- gsub("[\\],]", "", prok.sizes.name.rest$X.Organism.Name)
prok.sizes.name.rest$X.Organism.Name <- gsub("[\\',]", "", prok.sizes.name.rest$X.Organism.Name)

average.sizes.name.rest <- aggregate(prok.sizes.rest[,4], list(prok.sizes.name.rest$X.Organism.Name), mean)
colnames(average.sizes.name.rest) <- c("Genus", "Average_Genome_Size")

kraken2.nodata4$taxname <- gsub(" .*$", "", kraken2.nodata4$taxname)
kraken2.nodata4$taxname <- gsub("[\\[,]", "", kraken2.nodata4$taxname)
kraken2.nodata4$taxname <- gsub("[\\],]", "", kraken2.nodata4$taxname)
kraken2.nodata4$taxname <- gsub("[\\',]", "", kraken2.nodata4$taxname)


for (i in 1:length(kraken2.nodata4$taxname)){
  if (kraken2.nodata4$taxname[i] %in% average.sizes.name.rest$Genus){
    kraken2.nodata4$Average_Genome_Size[i] <- average.sizes.name.rest[which(average.sizes.name.rest$Genus == kraken2.nodata4$taxname[i]),"Average_Genome_Size"]}
  else{kraken2.nodata4$Average_Genome_Size[i] <- NA}
}

sum(is.na(kraken2.nodata4$Average_Genome_Size)) #45

for (i in 1:length(kraken2_otu$Taxon_ID)){
  if (kraken2_otu$Taxon_ID[i] %in% kraken2.nodata4$taxid){
    kraken2_otu$Average_Genome_Size[i] <- kraken2.nodata4[which(kraken2.nodata4$taxid == kraken2_otu$Taxon_ID[i]),"Average_Genome_Size"]}
  else{kraken2_otu$Average_Genome_Size[i] <- kraken2_otu$Average_Genome_Size[i]}
}

kraken2_otu$Ratio <- kraken2_otu$Raw_Abundance / kraken2_otu$Average_Genome_Size

kraken2.nodata5 <- subset(kraken2.nodata4, is.na(Average_Genome_Size))
length(kraken2.nodata5$taxid) #45

virus.sizes.rest <- dplyr::select(viruses, X.Organism.Name, TaxID, Status, Size..Kb.)
virus.sizes.rest <- subset(virus.sizes.rest, !(Status == "Complete Genome" | Status == "Complete"))
virus.sizes.rest$Size..Kb. <- virus.sizes.rest$Size..Kb.*0.001
colnames(virus.sizes.rest)[colnames(virus.sizes.rest)=="Size..Kb."] <- "Size..Mb."

average.virus.sizes.rest <- aggregate(virus.sizes.rest[,4], list(virus.sizes.rest$TaxID), mean)
colnames(average.virus.sizes.rest) <- c("Taxon_ID", "Average_Genome_Size")

for (i in 1:length(kraken2_otu$Taxon_ID)){
  if (kraken2_otu$Taxon_ID[i] %in% average.virus.sizes.rest$Taxon_ID){
    kraken2_otu$Average_Genome_Size[i] <- average.virus.sizes.rest[which(average.virus.sizes.rest$Taxon_ID == kraken2_otu$Taxon_ID[i]),"Average_Genome_Size"]}
  else{kraken2_otu$Average_Genome_Size[i] <- kraken2_otu$Average_Genome_Size[i]}
}

kraken2.nodata6 <- subset(kraken2_otu, is.na(Average_Genome_Size))
kraken2.nodata6 <- dplyr::select(kraken2.nodata6, Taxon_ID, Average_Genome_Size)
kraken2.nodata6 <- unique( kraken2.nodata6[ , 1:2 ] )
length(kraken2.nodata6$Taxon_ID) #41

write.csv(kraken2.nodata6$Taxon_ID, "tax_reports_uncollapsed/190722_no_data_viruses_rest.csv", quote = FALSE, row.names = FALSE)
kraken2.nodata6 <- read.delim("tax_reports_uncollapsed/tax_report_virus_rest.txt")
kraken2.nodata6 <- dplyr::select(kraken2.nodata6, taxid, taxname)
kraken2.nodata6 <- kraken2.nodata6[1:(length(kraken2.nodata6$taxid)-2),]

for (i in 1:length(kraken2.nodata6$taxname)){
  if (kraken2.nodata6$taxname[i] %in% average.virus.sizes.rest$Genus){
    kraken2.nodata6$Average_Genome_Size[i] <- average.virus.sizes.rest[which(average.virus.sizes.rest$Genus == kraken2.nodata6$taxname[i]),"Average_Genome_Size"]}
  else{kraken2.nodata6$Average_Genome_Size[i] <- NA}
}

virus.sizes.name.rest <- virus.sizes.rest
virus.sizes.name.rest$X.Organism.Name <- gsub(" .*$", "", virus.sizes.name.rest$X.Organism.Name)
virus.sizes.name.rest$X.Organism.Name <- gsub("[\\[,]", "", virus.sizes.name.rest$X.Organism.Name)
virus.sizes.name.rest$X.Organism.Name <- gsub("[\\],]", "", virus.sizes.name.rest$X.Organism.Name)
virus.sizes.name.rest$X.Organism.Name <- gsub("[\\',]", "", virus.sizes.name.rest$X.Organism.Name)

kraken2.nodata6$taxname <- gsub(" .*$", "", kraken2.nodata6$taxname)
kraken2.nodata6$taxname <- gsub("[\\[,]", "", kraken2.nodata6$taxname)
kraken2.nodata6$taxname <- gsub("[\\],]", "", kraken2.nodata6$taxname)
kraken2.nodata6$taxname <- gsub("[\\',]", "", kraken2.nodata6$taxname)

virus.sizes.name <- virus.sizes
virus.sizes.name$X.Organism.Name <- gsub(" .*$", "", virus.sizes.name$X.Organism.Name)
virus.sizes.name$X.Organism.Name <- gsub("[\\[,]", "", virus.sizes.name$X.Organism.Name)
virus.sizes.name$X.Organism.Name <- gsub("[\\],]", "", virus.sizes.name$X.Organism.Name)
virus.sizes.name$X.Organism.Name <- gsub("[\\',]", "", virus.sizes.name$X.Organism.Name)

average.virus.sizes.name <- aggregate(virus.sizes[,4], list(virus.sizes.name$X.Organism.Name), mean)
colnames(average.virus.sizes.name) <- c("Genus", "Average_Genome_Size")


average.virus.sizes.name.rest <- aggregate(virus.sizes.rest[,4], list(virus.sizes.name.rest$X.Organism.Name), mean)
colnames(average.virus.sizes.name.rest) <- c("Genus", "Average_Genome_Size")

for (i in 1:length(kraken2.nodata6$taxname)){
  if (kraken2.nodata6$taxname[i] %in% average.virus.sizes.name$Genus){
    kraken2.nodata6$Average_Genome_Size[i] <- average.virus.sizes.name[which(average.virus.sizes.name$Genus == kraken2.nodata6$taxname[i]),"Average_Genome_Size"]}
  else{kraken2.nodata6$Average_Genome_Size[i] <- NA}
}

for (i in 1:length(kraken2.nodata6$taxname)){
  if (kraken2.nodata6$taxname[i] %in% average.virus.sizes.name.rest$Genus){
    kraken2.nodata6$Average_Genome_Size[i] <- average.virus.sizes.name.rest[which(average.virus.sizes.name.rest$Genus == kraken2.nodata6$taxname[i]),"Average_Genome_Size"]}
  else{kraken2.nodata6$Average_Genome_Size[i] <- kraken2.nodata6$Average_Genome_Size[i]}
}

for (i in 1:length(kraken2_otu$Taxon_ID)){
  if (kraken2_otu$Taxon_ID[i] %in% kraken2.nodata6$taxid){
    kraken2_otu$Average_Genome_Size[i] <- kraken2.nodata6[which(kraken2.nodata6$taxid == kraken2_otu$Taxon_ID[i]),"Average_Genome_Size"]}
  else{kraken2_otu$Average_Genome_Size[i] <- kraken2_otu$Average_Genome_Size[i]}
}


kraken2.NA <- kraken2_otu[which(is.na(kraken2_otu$Average_Genome_Size)),]
length(unique(kraken2.NA$Taxon_ID)) #34
kraken2.NA.taxa <- as.data.frame(unique(kraken2.NA$Taxon_ID))
write.csv(kraken2.NA.taxa, "Rout_uncollapsed/Taxa_without_average_genome_size_190722.csv")


kraken2_otu$Ratio <- kraken2_otu$Raw_Abundance / kraken2_otu$Average_Genome_Size
sum(is.na(kraken2_otu$Average_Genome_Size)) #714
kraken2_otu <- na.omit(kraken2_otu)


#make table with just average genome size and taxon
mean.av.genome.size <- unique(kraken2_otu[,c(1,4)])
write.table(mean.av.genome.size, "Rout_uncollapsed/average_genome_sizes.txt")


kraken2_otu_out <- dplyr::select(kraken2_otu, Taxon_ID, Sample, Ratio)
kraken2_otu_out <- acast(kraken2_otu_out, Sample~Taxon_ID, value.var="Ratio")
write.csv(kraken2_otu_out, "Rout_uncollapsed/Rout_post_genome_normalisation_uncollapsed_190722.csv", quote = FALSE, row.names = TRUE)

