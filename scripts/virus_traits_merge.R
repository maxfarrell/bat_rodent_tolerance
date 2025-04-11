# virus_traits_merge.R
# merge experimental level data with mollentze 2020 reservoir data and guth 2022 CFR data

require(dplyr);packageVersion("dplyr")#1.1.4

primary <- read.csv("../raw_data/primary_screen.csv")
speciesonly <- read.csv("../raw_data/species_only_data.csv")
dat <- read.csv("../raw_data/individual_data.csv")
vtax <- read.csv("../clean_data/Virus_taxonomy.csv")
viruses <- read.csv("../raw_data/VirusNames_translation_Feb23_2024.csv")
viruses <- left_join(viruses, vtax)

primary <- left_join(primary, viruses)
speciesonly <- left_join(speciesonly, viruses)
dat <- left_join(dat, viruses)

# remove NA viruses
primary <- primary[!is.na(primary$Virus_ICTV),]
speciesonly <- speciesonly[!is.na(speciesonly$Virus_ICTV),]
dat <- dat[!is.na(dat$Virus_ICTV),]

# subset and join names
cols <- intersect(names(dat), names(primary))
dat <- dat[,cols]
primary <- primary[,cols]
dat <- full_join(dat,primary)
dat <- unique(dat)

cols <- intersect(names(dat), names(speciesonly))
speciesonly <- speciesonly[,cols]
dat <- full_join(dat,speciesonly)
dat <- unique(dat)


# load reservoir data
mollentze <- read.csv("../raw_data/mollentze_VirusReservoirs.csv")
guth <- read.csv("../raw_data/guth_stringent_data.csv")

guth$SppName_ICTV_MSL2018b[grep("Middle East", guth$SppName_ICTV_MSL2018b)] <- "Middle East respiratory syndrome-related coronavirus"

reservoir_dat <- full_join(mollentze, guth)


# But perhaps Guth has information for Mollentze's "ORPHAN" viruses
reservoirs <- reservoir_dat %>% select(SppName_ICTV_MSL2018b, Reservoir, hOrder) %>%
								unique() %>% arrange(SppName_ICTV_MSL2018b)

# View(reservoirs[reservoirs$Reservoir%in%"Orphan" & !is.na(reservoirs$hOrder) ,])
# Guth adds reservoir data for Banzi, Ilheus, MERS, Zika
# using Mollentze as source

# leaving the option below to FILL THESE IN MANUALLY - however, I'm not sure of the confidence in these (e.g. MERS)

reservoir_dat <- mollentze

reservoir_dat$Virus_ICTV <- reservoir_dat$SppName_ICTV_MSL2018b
intersect(dat$Virus_ICTV, reservoir_dat$Virus_ICTV)# 59
setdiff(dat$Virus_ICTV, reservoir_dat$Virus_ICTV)# 16

reservoir_dat$Virus_ICTV[grep("West Caucasian",reservoir_dat$Virus_ICTV)] <- "Lyssavirus caucasicus"
reservoir_dat$Virus_ICTV[grep("Lagos bat lyssavirus",reservoir_dat$Virus_ICTV)] <- "Lyssavirus lagos"
reservoir_dat$Virus_ICTV[grep("Khujand",reservoir_dat$Virus_ICTV)] <- "Lyssavirus khujand"
reservoir_dat$Virus_ICTV[grep("Irkut",reservoir_dat$Virus_ICTV)] <- "Lyssavirus irkut"
reservoir_dat$Virus_ICTV[grep("Aravan",reservoir_dat$Virus_ICTV)] <- "Lyssavirus aravan"
reservoir_dat$Virus_ICTV[grep("Aravan",reservoir_dat$Virus_ICTV)] <- "Lyssavirus aravan"
reservoir_dat$Virus_ICTV[grep("Rabies",reservoir_dat$Virus_ICTV)] <- "Lyssavirus rabies"
reservoir_dat$Virus_ICTV[grep("Australian bat",reservoir_dat$Virus_ICTV)] <- "Lyssavirus australis"
reservoir_dat$Virus_ICTV[grep("Colorado tick fever",reservoir_dat$Virus_ICTV)] <- "Colorado tick fever coltivirus"

setdiff(dat$Virus_ICTV, reservoir_dat$Virus_ICTV)# 8

# Mollentze 2020  and Guth 2022 are missing:

# Tahyna orthobunyavirus
# La Crosse orthobunyavirus
# Trivittatus orthobunyavirus
# Potosi orthobunyavirus
# Rio grande phlebovirus
# Spondweni virus (Mollentze has UniversalName, but not ICTV2018 name)
reservoir_dat$Virus_ICTV[reservoir_dat$UniversalName=="Spondweni virus"] <- "Spondweni virus"
# Snowshoe hare orthobunyavirus
# Pacui virus

# reservoir_dat$Reservoir[reservoir_dat$Virus_ICTV%in%"Banzi virus"] <- "Rodent"
# reservoir_dat$Reservoir[reservoir_dat$Virus_ICTV%in%"Ilheus virus"] <- "Avian"
# reservoir_dat$Reservoir[reservoir_dat$Virus_ICTV%in%"Middle East respiratory syndrome-related coronavirus"] <- "Chiroptera"
# reservoir_dat$Reservoir[reservoir_dat$Virus_ICTV%in%"Zika virus"] <- "NonHumanPrimate"

vtraits <- left_join(ungroup(dat), reservoir_dat)
vtraits <- vtraits %>% select(Virus_ICTV, SppName_ICTV_MSL2018b, Reservoir, Vector.borne, Reference_Reservoir, Reference_Vector) %>% unique()
vtraits <- vtraits[!is.na(vtraits$Virus_ICTV),]

vtraits$Virus_ICTV[vtraits$Reservoir%in%"MULTIPLE"]
# Usutu virus

length(unique(vtraits$Virus_ICTV[vtraits$Reservoir=="Rodent"]))#16
length(unique(vtraits$Virus_ICTV[vtraits$Reservoir=="Pterobat"]))#8
length(unique(vtraits$Virus_ICTV[vtraits$Reservoir=="Vespbat"]))#11

vtraits <- arrange(vtraits, Virus_ICTV) 

# adding CFR and other virus traits from guth
vtraits <- left_join(vtraits, select(guth, c(SppName_ICTV_MSL2018b, CFR_avg, human.trans, IsVectorBorne, GenomeType, ExclusiveCytoplasmicReplication))) %>% unique()
# looks like influenza A has two CFRs
# lower CFR is for H1N1; higher is for H5N1 and H7N9 
# none of these are included in our exerimental infection data
# because of large range, set as CFR as NA
vtraits$CFR_avg[vtraits$Virus_ICTV=="Influenza A virus"] <- NA
vtraits <- unique(vtraits)

write.csv(vtraits, "../clean_data/Virus_traits.csv",row.names=FALSE)
