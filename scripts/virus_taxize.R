# Grab virus taxonomy from taxize

require(dplyr);packageVersion("dplyr")#1.1.4
require(tidyr); packageVersion("tidyr")#1.3.0
require(taxize); packageVersion("taxize")#0.9.99

viruses <- read.csv("../raw_data/VirusNames_translation_Feb23_2024.csv")

df <- data.frame(Virus_ICTV=sort(unique(viruses$Virus_ICTV)))

if (!file.exists("../raw_data/viruses_taxize.rds")){

	tax <- classification(df$Virus_ICTV, db="ncbi")
	saveRDS(tax,"../raw_data/viruses_taxize.rds")

} else { tax <- readRDS("../raw_data/viruses_taxize.rds")}

taxdf <- dplyr::tbl_df(cbind(tax))
taxdf <- taxdf[, -grep("id$", colnames(taxdf))]
names(taxdf)[names(taxdf) == "species"] <- "Virus_NCBI"
names(taxdf)[names(taxdf) == "clade"] <- "realm"

taxdf$Virus_ICTV <- df$Virus_ICTV
taxdf <- taxdf %>% select(-c(subphylum, subfamily, no.rank, suborder, subgenus, query))
taxdf <- select(taxdf, Virus_ICTV, superkingdom:Virus_NCBI)

write.csv(taxdf, "../clean_data/Virus_taxonomy.csv", row.names=FALSE)
