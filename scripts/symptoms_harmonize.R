# symptoms_harmonize.R

require(dplyr);packageVersion("dplyr")#1.1.4
require(tidyr); packageVersion("tidyr")#1.3.0
require(corrplot); packageVersion("corrplot")#0.91
require(ggplot2);packageVersion("ggplot2")#3.4.4
require(tm);packageVersion("tm")#0.7.8
require(wordcloud);packageVersion("wordcloud")#2.6


primary <- read.csv("../raw_data/primary_screen.csv")
indiv <- read.csv("../raw_data/individual_data.csv")
speciesonly <- read.csv("../raw_data/species_only_data.csv")
vtax <- read.csv("../clean_data/Virus_taxonomy.csv")
viruses <- read.csv("../raw_data/VirusNames_translation_Feb23_2024.csv")
viruses <- dplyr::left_join(viruses, vtax)
hosts <- read.csv("../raw_data/HostNames_NCBI_Upham.csv")
hosts$Host_Upham <- gsub(" ", "_", hosts$Host_Upham)
#subset to only host reported and upham supertree
hosts <- hosts %>% select(Host_reported, Host_Upham) %>% unique()

# construct column for individualID
indiv <- indiv %>%                              # Create numbering variable
      group_by(PaperID) %>%
        mutate(IndividualID = paste(row_number(),unique(PaperID), sep="_"))

# add names
primary <- dplyr::left_join(primary, viruses)
primary <- dplyr::left_join(primary, hosts)
speciesonly <- dplyr::left_join(speciesonly, viruses)
speciesonly <- dplyr::left_join(speciesonly, hosts)
indiv <- dplyr::left_join(indiv, viruses)
indiv <- dplyr::left_join(indiv, hosts)

# remove genus only hosts
primary <- primary[!is.na(primary$Host_Upham),]
speciesonly <- speciesonly[!is.na(speciesonly$Host_Upham),]
indiv <- indiv[!is.na(indiv$Host_Upham),]


# remove NA viruses
indiv <- indiv[!is.na(indiv$Virus_ICTV),]
primary <- primary[!is.na(primary$Virus_ICTV),]
speciesonly <- speciesonly[!is.na(speciesonly$Virus_ICTV),]

# Harmonizing Organ_tissue_damaged
organs_indiv <- select(indiv, IndividualID, Organ_tissue_damaged)
organs_speconly <- select(speciesonly, PaperID, Organ_tissue_damaged)
organs <- full_join(organs_indiv, organs_speconly)

# making minor adjustments
organs$Organ_tissue_damaged[grep("arteries and veins", organs$Organ_tissue_damaged)] <- 
		"stomach; jejunum; mesentery; kidney; spleen; endothelial cells; tunica media"

organs_long <- organs %>% mutate(organ = strsplit(tolower(as.character(Organ_tissue_damaged)),"; *")) %>%
						unnest(organ)

# View(sort(unique(organs_long$organ)))

# cleaning and collapsing
organs_long$organ[organs_long$organ=="alveoli"] <- "lung"
organs_long$organ[organs_long$organ=="bronchioles"] <- "lung"
organs_long$organ[organs_long$organ=="caudal nose"] <- "nose"
organs_long$organ[organs_long$organ=="cilia"] <- "respiratory cilia"
organs_long$organ[organs_long$organ=="cornea"] <- "eye"
organs_long$organ[organs_long$organ=="duodenum"] <- "intestine"
organs_long$organ[grep("enterocyte", organs_long$organ)] <- "intestine"
organs_long$organ[organs_long$organ=="hipocampus"] <- "hippocampus"
organs_long$organ[organs_long$organ=="jejunum"] <- "intestine"
organs_long$organ[organs_long$organ=="kidneys"] <- "kidney"
organs_long$organ[organs_long$organ=="large intestine"] <- "intestine"
organs_long$organ[organs_long$organ=="lip"] <- "lips"
organs_long$organ[organs_long$organ=="lip"] <- "lips"
organs_long$organ[organs_long$organ=="lungs"] <- "lung"
organs_long$organ[organs_long$organ=="lungs (multifocal interstitial pneumonia)"] <- "lung"
organs_long$organ[organs_long$organ=="lymph node"] <- "lymph nodes"
organs_long$organ[organs_long$organ=="lymph node"] <- "lymph nodes"
organs_long$organ[grep("mild degenerative",organs_long$organ)] <- "intestine"
organs_long$organ[organs_long$organ=="nasal mucosa"] <- "nasal cavity"
organs_long$organ[organs_long$organ=="respiratory epithelium"] <- "respiratory tract"
organs_long$organ[organs_long$organ=="right cranial lungs lobe"] <- "lung"
organs_long$organ[organs_long$organ=="rostral nose"] <- "nasal cavity"
organs_long$organ[organs_long$organ=="small intestine"] <- "intestine"
organs_long$organ[organs_long$organ=="testicle"] <- "testes"
organs_long$organ[organs_long$organ=="upper lip"] <- "lip"

# NOTE: collapsing has not been done to a set level (e.g. organ)
# therefore, if investigating number of unique organs, this should not be summed without further collapsing!

# View(sort(unique(organs_long$organ)))
 
# internal / external distinction
organs_long$external_damage <- 0

external_orgs <- c("digit", "extremities", "eye", "groin", "gums", "lip,", "lips",
					"mouth", "nasa cavity", "nose", "oral cavity", "skin",
					"tongue")

organs_long$external_damage[organs_long$organ %in% external_orgs] <- 1

# internal is everything NOT external except NA and "none"
organs_long$internal_damage <- 1 - organs_long$external_damage
organs_long$internal_damage[organs_long$organ%in%"none"] <- 0


# neural tropism
organs_long$neural_damage <- 0
neural_orgs <- c("brain", "brainstem", "cerebellum", "hippocampus", "spinal cord")
organs_long$neural_damage[organs_long$organ %in% neural_orgs] <- 1


# respiratory tropism
organs_long$respiratory_damage <- 0
respiratory_orgs <- c("lung", "respiratory cilia", "respiratory tract")
organs_long$respiratory_damage[organs_long$organ %in% respiratory_orgs] <- 1


# reproductive tropism
organs_long$reproductive_damage <- 0
reproduction_orgs <- c("endometrial gland (distended)", "ovaries",
						"prostate", "testes", "uterus")
organs_long$reproductive_damage[organs_long$organ %in% reproduction_orgs] <- 1


# digestive tropism
organs_long$digestive_damage <- 0
digestive_orgs <- c("colon", "intestines", "rectum", "stomach")
organs_long$digestive_damage[organs_long$organ %in% digestive_orgs] <- 1


# Put NAs back then merge with original data
organs_long$external_damage [is.na(organs_long$Organ_tissue_damaged)] <- NA
organs_long$internal_damage [is.na(organs_long$Organ_tissue_damaged)] <- NA
organs_long$neural_damage [is.na(organs_long$Organ_tissue_damaged)] <- NA
organs_long$respiratory_damage [is.na(organs_long$Organ_tissue_damaged)] <- NA
organs_long$reproductive_damage [is.na(organs_long$Organ_tissue_damaged)] <- NA
organs_long$digestive_damage [is.na(organs_long$Organ_tissue_damaged)] <- NA

organs_long <- select(organs_long, -Organ_tissue_damaged) %>% unique()
organs_long <- select(organs_long, -organ) %>% unique()


# converting back to one row per individual
organs_wide <- as.data.frame(organs_long) %>% group_by(PaperID, IndividualID) %>%
								summarise(external_damage 		= (sum(external_damage, na.rm=TRUE)>0)*1,
										  internal_damage 		= (sum(internal_damage, na.rm=TRUE)>0)*1,
										  neural_damage   		= (sum(neural_damage, na.rm=TRUE)>0)*1,
										  respiratory_damage 	= (sum(respiratory_damage, na.rm=TRUE)>0)*1,
										  reproductive_damage 	= (sum(reproductive_damage, na.rm=TRUE)>0)*1,
										  digestive_damage 		= (sum(digestive_damage, na.rm=TRUE)>0)*1,
										  )

# merge with original data
speciesonly <- left_join(speciesonly, organs_wide)
indiv <- left_join(indiv, organs_wide)


# merging individual and species only data
dat <- full_join(indiv, speciesonly)


# manually collapse disease manifestations (Behaviour + Physical_symptoms) 
# based on keywords (grep)
# creating new columns

# first view word clouds

behaviours <- Corpus(VectorSource(unique(dat$Behaviour)))
behaviours <- tm_map(behaviours, stripWhitespace) 
behaviours <- tm_map(behaviours, content_transformer(tolower))
behaviours <- tm_map(behaviours, removeWords, stopwords("english"))
# behaviours <- tm_map(behaviours, stemDocument) # optional in this case
behaviours <- tm_map(behaviours, removeNumbers) # optional in this case
behaviours <- tm_map(behaviours, removePunctuation)
# tuning
words_to_remove <- c("none","note","see","cage","woodwork")
behaviours <- tm_map(behaviours, removeWords, words_to_remove)

# png(file = "../plots_tables/wordcloud_behaviours.png", width = 1500, height = 1500, units = "px", res = 300, bg = "white")
# wordcloud(behaviours, scale=c(3,0.5), max.words=100, random.order=FALSE, rot.per=0.35, use.r.layout=FALSE, colors=brewer.pal(8, "Dark2"))
# dev.off()


symptoms <- Corpus(VectorSource(unique(dat$Physical_symptoms)))
symptoms <- tm_map(symptoms, stripWhitespace) 
symptoms <- tm_map(symptoms, content_transformer(tolower))
symptoms <- tm_map(symptoms, removeWords, stopwords("english"))
# symptoms <- tm_map(symptoms, stemDocument) # optional in this case
symptoms <- tm_map(symptoms, removeNumbers) # optional in this case
symptoms <- tm_map(symptoms, removePunctuation)
# tuning
words_to_remove <- c("none","note","see","cage","woodwork")
symptoms <- tm_map(symptoms, removeWords, words_to_remove)

# png(file = "../plots_tables/wordcloud_symptoms.png", width = 1500, height = 1500, units = "px", res = 300, bg = "white")
# wordcloud(symptoms, scale=c(3,0.5), max.words=100, random.order=FALSE, rot.per=0.35, use.r.layout=FALSE, colors=brewer.pal(8, "Dark2"))
# dev.off()


# Behaviour collapse
# core aspects: aggression/biting, lethargy/depression, vocalizations, incoordination,
# View(sort(c(unique(dat$Behaviour))))
# sensitivity, appetite loss

# "reclusive bevahiour"? 
# "salivation"
# "unresponsive to touch"


# key terms
aggression_terms <- c("aggression", "aggressive", "agitation", "biting", "vicious", "irritable", "furious", "fury", "agitation")

lethargy_terms <- c("depression", "prostration", "lethargy", "listless", "sleepy", "reluctance to move")

# odd or increased vocalizations
vocalization_terms <- c("atypical vocalizations","odd vocalization", "vocalization", "dysphasia", "vocal", "hissing sounds")

incoordination_terms <- c("hunching", "staggering", "head tilt", "incoordination")
# including with movement terms below (!)

# increase sensitivity to stimulti or over-excitement - maybe not necessary?
sensitivity_terms <- c("excitement", "fly violently", "flying violently", "irritability to light")

# loss of appetite ./ refusal to eat or drink
appetite_terms <- c("inappetance", "inapetence", "ate very little","reduced food intake", "refused food", "dehydration")



# Physical symptoms collapse
# View(sort(c(unique(dat$Physical_symptoms))))
# some behaviour terms are included: "aggressive","agitation","lethargy",

# core aspects: 
# lesions (skin?), weight loss, movement inhibition, pneumonia,
# nasal/ocular discharge

# "splenitis"
# "hyperplasia"
# "altered reflexes"
# "bruising" / "hematomas"
# hyperaesthesia
# "nystagmus" - uncontrolled movement of the eyes (maybe movement inhibition)
# incontinence
# cellular infiltrates in heart / other cardiac stuff? (tiger heart)
# "crusty noses"
# "skin pustules"
# "laboured breathing"
# "vesicular rash"
# "pustular rash"
# "skin pustules"
# "enlarged lymph nodes"
# "green coloured diarrhoea"
# "enlarged discoloured liver"
# "tachypnea"

skin_issue_terms <- c("skin lesions", "cutaneous lesion", "hypopigmentation", "tongue lesion",
					"small lesion on lip", "lesions near both eyes", "lesions in the skin",
					"corneal lesion", "epithelial hyperplasia", "multiple dermal foci",
					"lesion on gum", "lesions on tongue", "skin pustules", "petechial rash",
					"vesicle on tongue", "vesicular rash", "pustular rash", "lesions on the head") 

weight_loss_terms <- c("weight loss", "anorexia")

movement_terms <- c("paralysis","tremor","spasms", "ataxia","weakness","paresis", 
					"nystagmus", "seizure","muscle weakness", "trembling", "coma", 
					"unable to fly", "unable to move", "inability to stand", 
					"altered reflexes")

# include incoordination terms
movement_terms <- c(movement_terms, incoordination_terms)


eye_nose_infection_terms <- c("rhinitis", "mild conjunctivitis", "nasal discharge", "ocular discharge", "red puffy nose", 
		"crusty nose", "nasal crust", "bloody nose")

pneumonia_terms <- c("pneumonia", "respiratory disease")


internal_organ_damage_terms <- c("splenitis", "cardiomyocyte necrosis", "hemmorhage in lungs",
							"lesions on brain", "lesions on testes", "lesions on salivary gland", 
							"hepatocellular degeneration", 
							"ulceration of villus tips", "multifocally fused villi", "degeneration of villus cells")

# specifially note observations of damage to nervous, respiratory, redproductive, or digestive organs
# merge with data on organ damage from above

nervous_terms <- c("lesions on brain")
respiratory_terms <- c("hemmorhage in lungs")
reproductive_terms <- c("lesions on testes")
digestive_terms <- c("ulceration of villus tips", "multifocally fused villi", "degeneration of villus cells")


# Constructing binary variables

# constructing single "Symptoms" column that collapses behaviour and physical columns
# this is because some behaviours like aggression are included in physical symptoms 


# individual level data
symptoms <- paste0(dat$Behaviour, ";", dat$Physical_symptoms)

dat$aggression <- unlist(sapply(symptoms, FUN=function(x) any(grep(paste0(aggression_terms,collapse="|"), x))*1 ))
dat$lethargy <- unlist(sapply(symptoms, FUN=function(x) any(grep(paste0(lethargy_terms,collapse="|"), x))*1 ))
dat$vocalization <- unlist(sapply(symptoms, FUN=function(x) any(grep(paste0(vocalization_terms,collapse="|"), x))*1 ))
dat$sensitivity <- unlist(sapply(symptoms, FUN=function(x) any(grep(paste0(sensitivity_terms,collapse="|"), x))*1 ))
dat$appetite <- unlist(sapply(symptoms, FUN=function(x) any(grep(paste0(appetite_terms,collapse="|"), x))*1 ))
dat$skin_issue <- unlist(sapply(symptoms, FUN=function(x) any(grep(paste0(skin_issue_terms,collapse="|"), x))*1 ))
dat$weight_loss <- unlist(sapply(symptoms, FUN=function(x) any(grep(paste0(weight_loss_terms,collapse="|"), x))*1 ))
dat$movement <- unlist(sapply(symptoms, FUN=function(x) any(grep(paste0(movement_terms,collapse="|"), x))*1 ))
dat$eye_nose <- unlist(sapply(symptoms, FUN=function(x) any(grep(paste0(eye_nose_infection_terms,collapse="|"), x))*1 ))
dat$pneumonia <- unlist(sapply(symptoms, FUN=function(x) any(grep(paste0(pneumonia_terms,collapse="|"), x))*1 ))

# terms to integrate with observed organ damage column
dat$obs_nervous <- unlist(sapply(symptoms, FUN=function(x) any(grep(paste0(nervous_terms,collapse="|"), x))*1 ))
dat$obs_respiratory <- unlist(sapply(symptoms, FUN=function(x) any(grep(paste0(respiratory_terms,collapse="|"), x))*1 ))
dat$obs_reproductive <- unlist(sapply(symptoms, FUN=function(x) any(grep(paste0(reproductive_terms,collapse="|"), x))*1 ))
dat$obs_digestive <- unlist(sapply(symptoms, FUN=function(x) any(grep(paste0(digestive_terms,collapse="|"), x))*1 ))
dat$obs_internal_organ <- unlist(sapply(symptoms, FUN=function(x) any(grep(paste0(internal_organ_damage_terms,collapse="|"), x))*1 ))


# merge with observed organ damage
dat$internal_damage[dat$obs_internal_organ==1] <- 1
dat$neural_damage[dat$obs_nervous==1] <- 1
dat$respiratory_damage[dat$obs_respiratory==1] <- 1
dat$reproductive_damage[dat$obs_reproductive==1] <- 1
dat$digestive_damage[dat$obs_digestive==1] <- 1



symptoms_dat <- dat %>% ungroup() %>% select(Virus_ICTV, Host_Upham,
								aggression, lethargy, vocalization, 
								sensitivity, appetite, skin_issue,
								weight_loss, movement, eye_nose, pneumonia,
								internal_damage) %>% unique()

symptoms_dat <- symptoms_dat %>% select(-Virus_ICTV, -Host_Upham)

# replace with zeros for corellation plot
symptoms_dat <- mutate_all(symptoms_dat, ~coalesce(.,0))

res <- cor(symptoms_dat)

clean_names <- c("Aggression", "Lethargy", "Vocalizations", "Increased sensitivity",
					"Loss of appetite", "Skin rashes / lesions", "Weight loss", "Impaired movement",
					"Eye / nose infections", "Pneumonia", "Internal organ damage")

colnames(res) <- clean_names
rownames(res) <- clean_names

# pdf("../plots_tables/symptoms_hostXvirus.pdf", width = 10, height=8)
# corrplot(res, type="upper", tl.col="black", tl.srt=45)
# dev.off()

# corrplot(res, type="upper", tl.col="black", tl.srt=45, diag=FALSE)
# corrplot(res, type="upper", tl.col="black", tl.srt=45, order = 'AOE')
# corrplot(res, type="upper", tl.col="black", tl.srt=45, method = 'shade', order = 'AOE', diag = FALSE)

pdf("../plots_tables/symptoms_hostXvirus.pdf", width = 10, height=8)
corrplot(res, type="upper", tl.col="black", tl.srt=45, method = 'color', order = 'AOE', diag = FALSE)
dev.off()

max(res[upper.tri(res)])# max corellation is 0.458


symptoms_dat <- dat %>% ungroup() %>% select(Virus_ICTV, Host_Upham, ClinicalDisease,
								aggression, lethargy, vocalization, 
								sensitivity, appetite, skin_issue,
								weight_loss, movement, internal_damage, 
								eye_nose, pneumonia) %>% unique()

symptoms_dat <- symptoms_dat %>% select(-Virus_ICTV, -Host_Upham) %>% unique()

# View(symptoms_dat[symptoms_dat$ClinicalDisease=="N",])
# some of the individuals reported as no Clinical Disease did indeed display symptoms... 
# odd - this questions whether we have a consistent definition of clinical disease
# alternative is to assume any individual with symptoms has signs of "clinical disease"...


# make smaller symptoms data frame for use elsewhere
dat_small <- select(dat, PaperID, IndividualID, Virus_ICTV, Virus_strain_isolate, Host_Upham, Host_order, 
					aggression, lethargy, vocalization, 
								sensitivity, appetite, skin_issue,
								weight_loss, movement, 
								eye_nose, pneumonia,
								neural_damage, respiratory_damage, 
								reproductive_damage, digestive_damage,
								external_damage, internal_damage) %>% unique()


# re-doing symptoms plot across individuals 
symptoms_dat <- dat %>% filter(!is.na(IndividualID)) %>% select(PaperID, IndividualID, Virus_ICTV, Host_Upham, 
								aggression, lethargy, vocalization, 
								sensitivity, appetite, skin_issue,
								weight_loss, movement, 
								eye_nose, pneumonia, internal_damage) %>% unique() %>% 
				ungroup() %>% select(-PaperID, -IndividualID, -Virus_ICTV, -Host_Upham)

symptoms_dat <- mutate_all(symptoms_dat, ~coalesce(.,0))

res <- cor(symptoms_dat)
colnames(res) <- clean_names
rownames(res) <- clean_names
corrplot(res, type="upper", tl.col="black", tl.srt=45, method = 'color', order = 'AOE', diag = FALSE)

pdf("../plots_tables/symptoms_individuals.pdf", width = 10, height=8)
corrplot(res, type="upper", tl.col="black", tl.srt=45, method = 'color', order = 'AOE', diag = FALSE)
dev.off()

max(res[upper.tri(res)])# max corellation is 0.597



##################################
## constructing a severity rank ##
##################################

# exploring author-reported severity
# View(unique(sort(dat$severity_reported)))
sum(!is.na(dat$severity_reported))/nrow(dat)# only 8.65% of entries have data

# Histopathology information
# sort(unique(dat$Histopathology))
# medium, mild, moderate, none, severe

# View(dat[dat$Histopathology%in%c("mild"),])
# View(dat[dat$Histopathology%in%c("moderate"),])


# for the sake of this section we assume if something is unreported it is an NA
dat$internal_damage[is.na(dat$internal_damage)] <- 0
dat$external_damage[is.na(dat$external_damage)] <- 0
dat$neural_damage[is.na(dat$neural_damage)] <- 0
dat$respiratory_damage[is.na(dat$respiratory_damage)] <- 0
dat$reproductive_damage[is.na(dat$reproductive_damage)] <- 0
dat$digestive_damage[is.na(dat$digestive_damage)] <- 0


dat$severity_rank <- NA


# 1 - No disease
# authors say "no clinical disease" or do not mention clinical disease, and no other symptoms reported

dat$severity_rank[	(dat$ClinicalDisease%in%"N" | is.na(dat$ClinicalDisease))  
					 & dat$aggression%in%0
					 & dat$lethargy%in%0 
					 & dat$vocalization%in%0 
					 & dat$sensitivity%in%0
					 & dat$appetite%in%0 
					 & dat$skin_issue%in%0
					 & dat$weight_loss%in%0
					 & dat$movement%in%0
					 & dat$eye_nose%in%0
					 & dat$pneumonia%in%0
					 & dat$internal_damage%in%0 
					 & dat$external_damage%in%0 
					 & (dat$Mortality_YN%in%"N" | is.na(dat$Mortality_YN))
					 ] <- 1


# 2 - Subclinical disease			
# authors say no clinical disease, but there is mention of some symptoms (assuming these are mild)
# OR  authors say clinical disease, but there is no mention of other symptons (assuming mild clinical disease)
dat$severity_rank[	(dat$ClinicalDisease%in%"N" 
					 & (dat$aggression%in%1
					 | dat$lethargy%in%1
					 | dat$vocalization%in%1 
					 | dat$sensitivity%in%1
					 | dat$appetite%in%1 
					 | dat$skin_issue%in%1
					 | dat$weight_loss%in%1
					 | dat$movement%in%1
					 | dat$eye_nose%in%1
					 | dat$pneumonia%in%1
					 | dat$external_damage%in%1
					 | dat$internal_damage%in%1
					 ))
					 ] <- 2


# 3 - Mild to moderate disease 
# authors say clinical disease, but there is no mention of other symptons (assuming mild clinical disease)					 
# OR authors do not say no clinical disease (either reported clinical disease or unreported)					 
# OR authors report mild minimal or mild to moderate disease 
# OR changes in eye_nose, skin_issue, lethargy, loss of appetite, weight loss, aggression, vocalization, 
					 # sensitivity, pneumonia
					 # no movement impairment

dat$severity_rank[ ((dat$ClinicalDisease%in%"Y")  
					 & dat$aggression%in%0
					 & dat$lethargy%in%0 
					 & dat$vocalization%in%0 
					 & dat$sensitivity%in%0
					 & dat$appetite%in%0 
					 & dat$skin_issue%in%0
					 & dat$weight_loss%in%0
					 & dat$movement%in%0
					 & dat$eye_nose%in%0
					 & dat$pneumonia%in%0
					 & dat$internal_damage%in%0 
					 & dat$external_damage%in%0 
					 & (dat$Mortality_YN%in%"N" | is.na(dat$Mortality_YN))
					 )
					 
					 |

					 (dat$ClinicalDisease%in%"Y" | is.na(dat$ClinicalDisease))
					 & (
					   1:nrow(dat)%in%grep(c("mild | minimal"), dat$severity_reported, ignore.case=TRUE)
					 | 1:nrow(dat)%in%grep(c("mild | moderate | medium"), dat$Histopathology, ignore.case=TRUE)
					 | dat$eye_nose%in%1
					 | dat$lethargy%in%1
					 | dat$skin_issue%in%1
					 | dat$appetite%in%1
					 | dat$external_damage%in%1
					 | dat$weight_loss%in%1
					 | dat$aggression%in%1
					 | dat$vocalization%in%1 
					 | dat$sensitivity%in%1
					 | dat$pneumonia%in%1
					 )
					 & (dat$movement%in%0
					 # & dat$internal_damage%in%0
					 & (dat$Mortality_YN%in%"N" | is.na(dat$Mortality_YN))
					 )
					 ] <- 3


# 4 - Severe disease
# authors report severe disease
# OR signs of internal organ damage OR imparied modement
# OR the disease manifestation was too severe that the authors culled the animal for ethical reasons 	


dat$severity_rank[   (dat$ClinicalDisease%in%"Y" | is.na(dat$ClinicalDisease))
					 &
					 (1:nrow(dat)%in%grep(c("severe"), dat$severity_reported, ignore.case=TRUE)
					 | 1:nrow(dat)%in%grep(c("severe"), dat$Histopathology, ignore.case=TRUE)
					 | dat$internal_damage%in%1 | dat$movement%in%1 | dat$Mortality_reason=="Culled_ethics")
					 & (dat$Mortality_YN%in%"N" | is.na(dat$Mortality_YN))

					 ] <- 4


# 5 - death
dat$severity_rank[  dat$Mortality_YN%in%"Y" ] <- 5

table(dat$severity_rank)
# are all entries are ranked?
# nrow(dat[is.na(dat$severity_rank),])
# yes


# set severity rank to NA for paper-host-virus triples where pathology is not reported, 
# and individuals are not susceptible
dat <- left_join(dat, unique(primary[,c("PaperID","Host_order","PathologyReported_YN")]))

# what does severeity look like for studies we say pathology is not reported?
unique(dat$severity_rank[dat$Susceptible_YN%in%"Y" & dat$PathologyReported_YN%in%"N"])
# 1
# these should be set to NA then removed 


# joining severity rank info
# to avoid duplicate entries, set IndividualID for speciesonly data from NA to the paperID
# construct column for individualID
dat[is.na(dat$IndividualID),] <- dat[is.na(dat$IndividualID),] %>%                              # Create numbering variable
      group_by(PaperID) %>%
        mutate(IndividualID = paste(row_number(),unique(PaperID), sep="_"))

dat_small[is.na(dat_small$IndividualID),] <- dat_small[is.na(dat_small$IndividualID),] %>%                              # Create numbering variable
      group_by(PaperID) %>%
        mutate(IndividualID = paste(row_number(),unique(PaperID), sep="_"))

dat_small <- left_join(dat_small, unique(select(dat, PaperID, Virus_strain_isolate, Host_order, IndividualID, severity_rank, Susceptible_YN, PathologyReported_YN)))


# setting disease severity for entries with non-susceptible individuals (all zeros)
dat_small$severity_rank[!dat_small$Susceptible_YN%in%"Y"] <- NA

# setting disease severity for entries with no pathology reported
dat_small$severity_rank[dat_small$PathologyReported_YN=="N"] <- NA
dat_small <- select(dat_small, -PathologyReported_YN)

# removing entries with severity_rank is NA
dat_small <- dat_small[!is.na(dat_small$severity_rank),]

table(dat_small$severity_rank[dat_small$Susceptible_YN%in%"Y"])
sum(table(dat_small$severity_rank[dat_small$Susceptible_YN%in%"Y"]))#2030

# setting Individual_ID back to NA for species only data
dat_small$IndividualID[dat_small$PaperID%in%speciesonly$PaperID] <- NA

# writing csv
write.csv(dat_small,"../clean_data/symptoms_severity.csv", row.names=FALSE)


# plot distribution

# 5 color viridis 'inferno'
colors <- (c("#fcffa4","#f98e09","#bc3754","#57106e","#000004"))

# distribution across individuals
p <- ggplot(dat_small[dat_small$Susceptible_YN%in%"Y" & !is.na(dat_small$IndividualID),], aes(x=as.factor(severity_rank))) + geom_bar(fill=colors) + ylab("Number of individuals") + 
				xlab("Severity rank") + 
				scale_x_discrete(labels=as.character(sort(unique(dat_small$severity_rank))), breaks=sort(unique(dat$severity_rank))) + 
				theme_minimal() + scale_y_continuous(limits=c(0,1400),breaks=seq(0,1500,200)) +  guides(fill = guide_legend(reverse=TRUE))

p

ggsave("../plots_tables/severity_distribution_individuals.pdf", width=4, height=3)
ggsave("../plots_tables/severity_distribution_individuals.png", width=4, height=3)


# stacked barplot split by host order
# p + facet_wrap(~Host_order)

nrow(dat_small[!is.na(dat_small$severity_rank),])

# 5 color viridis 'inferno'
colors <- rev(c("#fcffa4","#f98e09","#bc3754","#57106e","#000004"))

p2 <- ggplot(dat_small[!is.na(dat_small$severity_rank) & !is.na(dat_small$IndividualID),], aes(fill=forcats::fct_rev(as.factor(severity_rank)), y=Host_order)) + 
		geom_bar(position="fill", stat="count") +
		scale_fill_manual(values=colors, name="Severity rank   ") + scale_color_manual(values=colors) +
		theme_minimal() + theme(legend.position = "bottom", legend.direction="horizontal", 
								legend.justification="right", legend.box.margin=margin(-15,-10,-10,-10),
								legend.spacing.x = unit(-0.1, 'cm'), legend.title=element_text(size=10),) + 
        guides(fill = guide_legend(reverse=TRUE, title.vjust=0.2, title.hjust=0.5, title.position="left",
        							label.position="top", label.vjust=-2)) +    
		ylab("") + xlab("Proportion of individuals")  

p2

ggsave("../plots_tables/severity_proportion_hostOrder.pdf", width=6, height=3)
ggsave("../plots_tables/severity_proportion_hostOrder.png", width=6, height=3)
