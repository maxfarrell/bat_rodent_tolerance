# virion_merge.R

require(dplyr);packageVersion("dplyr")#1.1.4
require(vroom);packageVersion("vroom")#1.5.5
require(ape);packageVersion("ape")#5.5
require(PhyloMeasures);packageVersion("PhyloMeasures")#2.1


# load experiment data
primary <- read.csv("../raw_data/primary_screen.csv")
vtax <- read.csv("../clean_data/Virus_taxonomy.csv")
viruses <- read.csv("../raw_data/VirusNames_translation_Feb23_2024.csv")
viruses <- dplyr::left_join(viruses, vtax)
hosts <- read.csv("../raw_data/HostNames_NCBI_Upham.csv")
hosts$Host_Upham <- gsub(" ", "_", hosts$Host_Upham)
#subset to only host reported and upham
hosts <- hosts %>% select(Host_reported, Host_Upham) %>% unique()

# add names
primary <- dplyr::left_join(primary, viruses)
primary <- dplyr::left_join(primary, hosts)

# remove NA viruses
primary <- primary[!is.na(primary$Virus_ICTV),]

dat <- primary

# load virion
virion <- vroom("../raw_data/Virion.csv.gz", delim="\t")

missing <- setdiff(dat$Virus_NCBI_TaxID, virion$VirusTaxID)
missing
# only two missing..
unique(dat$Virus_ICTV[dat$Virus_NCBI_TaxID%in%missing])

grep("Trivi", virion$Virus)
grep("Potosi", virion$Virus)

# investigating viruses which have fewer reported bat/rat hosts in VIRION
# compared to our experimental data
# note that some of our primary screen data include experiments where individuals
# were found to be not susceptible 
# unique(virion$Host[grep("banzi", virion$Virus)])
# unique(virion$Host[grep("bwamba", virion$Virus)])
# unique(virion$Host[grep("pongola", virion$Virus)])
# unique(virion$Host[grep("tai forest", virion$Virus)])

# virion to upham harmonization
virion <- virion[virion$VirusTaxID %in% dat$Virus_NCBI_TaxID,]

# Upham mammal supertree
tree <- read.nexus("../raw_data/upham_tree_666.nex")
tax <- read.csv("../raw_data/taxonomy_mamPhy_5911species.csv")

tree$tip.label <- tolower(gsub("_"," ", tree$tip.label))
tax$Species_Name <- tolower(gsub("_"," ", tax$Species_Name))
virion_mamm <- virion[virion$HostClass=="mammalia",]

setdiff(virion_mamm$Host, tree$tip.label)# 66 hosts missing :)
missing <- setdiff(virion_mamm$Host, tree$tip.label)

# harmonizing host names to Upham (using clover translation table)
trans_t <- read.csv("../raw_data/CLOVER/mammal_phylo_translations.csv")

intersect(setdiff(missing, tree$tip.label), trans_t$Host)# get 44 with trans table
lookup <- setNames(tolower(trans_t$Host_Upham), trans_t$Host)

# remove mysterious NA hosts
virion_mamm <- virion_mamm[!is.na(virion_mamm$Host),]

virion_mamm$Host[!virion_mamm$Host%in%tree$tip.label]
virion_mamm$Host_Upham <- lookup[virion_mamm$Host]

# any NAs created re-insert original host name
virion_mamm$Host_Upham[is.na(virion_mamm$Host_Upham)] <- virion_mamm$Host[is.na(virion_mamm$Host_Upham)]

setdiff(virion_mamm$Host_Upham, tree$tip.label)
# only 21

# giraffa giraffa is a subspecies of g. camelopardalis
virion_mamm$Host_Upham[virion_mamm$Host_Upham=="giraffa giraffa"] <- "giraffa camelopardalis"

# cricetulus griseus is a subspecies of c. barabensis
virion_mamm$Host_Upham[virion_mamm$Host_Upham=="cricetulus griseus"] <- "cricetulus barabensis"

# molossus ater synonym of m. rufus (https://www.gbif.org/species/5218728)
virion_mamm$Host_Upham[virion_mamm$Host_Upham=="molossus ater"] <- "molossus rufus"

# loxodonta cyclotis was split from l africana but Upham does not recognize this
virion_mamm$Host_Upham[virion_mamm$Host_Upham=="loxodonta cyclotis"] <- "loxodonta africana"

# hypsugo pulveratus is in MDD, but not upham tree...
# moved from Pipistrellus to Hypsugo
# upham tree still has pipistrellus pulveratus
virion_mamm$Host_Upham[virion_mamm$Host_Upham=="hypsugo pulveratus"] <- "pipistrellus pulveratus"

# microtus obscurus
# was split from arvalis
virion_mamm$Host_Upham[virion_mamm$Host_Upham=="microtus obscurus"] <- "microtus arvalis"

# rhinolophus monoceros
# has been considered a synonym of R. pusillus, although recent publications have retained the species as distinct
virion_mamm$Host_Upham[virion_mamm$Host_Upham=="rhinolophus monoceros"] <- "rhinolophus pusillus"

# laephotis capensis
# moved from Neoromicia to Laephotis; includes the recently described melckorum
virion_mamm$Host_Upham[virion_mamm$Host_Upham=="laephotis capensis"] <- "neoromicia capensis"

# eothenomys eleusis
# split from E. melanogaster
virion_mamm$Host_Upham[virion_mamm$Host_Upham=="eothenomys eleusis"] <- "eothenomys melanogaster"

# dobsonia magna
# includes magna for the time being, but D. moluccensis likely represents a species complex
virion_mamm$Host_Upham[virion_mamm$Host_Upham=="dobsonia magna"] <- "dobsonia moluccensis"

# rhabdomys dilectus
# previously included R. chakae
# but r chakae is not in upham
# only rhabdomys un upham is Rhabdomys pumilio

# virion_mamm$Host[grep("rhabdomys", virion_mamm$Host)]
# virion_mamm[grep("rhabdomys", virion_mamm$Host),]
# # only monkeypox
# sort(unique(virion_mamm$Host[grep("monkeypox", virion_mamm$Virus)]))#107 hosts
# sort(unique(virion_mamm$Host[virion_mamm$Virus %in% "monkeypox virus" & virion_mamm$HostOrder%in%"rodentia"]))# 39 rodents

# in this case, I will collapse these two species - 
# very similar biogeography and since this collapse only impacts MPX
# which has a large number of rodent hosts, this should be a negligible influence 
virion_mamm$Host_Upham[virion_mamm$Host_Upham=="rhabdomys dilectus"] <- "rhabdomys pumilio"

# macaca brunnescens
# nominal name of Macaca ochreata
virion_mamm$Host_Upham[virion_mamm$Host_Upham=="macaca brunnescens"] <- "macaca ochreata"

# cercopithecus doggetti
# synonym of Cercopithecus_mitis
virion_mamm$Host_Upham[virion_mamm$Host_Upham=="cercopithecus doggetti"] <- "cercopithecus mitis"

# cercopithecus roloway
# MSW reports this is sometimes considered a subspecies of C. diana
virion_mamm$Host_Upham[virion_mamm$Host_Upham=="cercopithecus roloway"] <- "cercopithecus diana"

# cricetomys ansorgei
# upham says it was previously a synonym of c. emini OR c. gambianus, both are in the tree
# both are also in virion_mamm
# virion_mamm$Host[grep("cricetomys",virion_mamm$Host)] %>% sort()
# virion_mamm$Virus[grep("cricetomys",virion_mamm$Host)] %>% sort()
# again all are for mokeypox virus, so perhaps doesn't matter so much if we collapse one
# only gambiatus is in our data
virion_mamm$Host_Upham[virion_mamm$Host_Upham=="cricetomys ansorgei"] <- "cricetomys gambianus"

# allochrocebus preussi
# moved from Cercopithecus to Allochrocebus
virion_mamm$Host_Upham[virion_mamm$Host_Upham=="allochrocebus preussi"] <- "cercopithecus preussi"

# macaca speciosa
# known as stump tailed macaque (Macaca arctoides, synonym m. speciosus)
virion_mamm$Host_Upham[virion_mamm$Host_Upham=="macaca speciosa"] <- "macaca arctoides"

# cercopithecus kandti
# synonym of c. mitis
virion_mamm$Host_Upham[virion_mamm$Host_Upham=="cercopithecus kandti"] <- "cercopithecus mitis"

setdiff(virion_mamm$Host_Upham, tree$tip.label)

# only hybrids and "marmosets"

virion_mamm <- virion_mamm[virion_mamm$Host_Upham%in%tree$tip.label,]

#replace virion mammals with virion_mamm
virion <- virion[!virion$HostClass%in%"mammalia",]
virion <- plyr::rbind.fill(virion, virion_mamm)

# collapse mammal names only (make Host as Host_Upham when not NA)
virion$Host[!is.na(virion$Host_Upham)] <- virion$Host_Upham[!is.na(virion$Host_Upham)]

length(unique(virion$Virus))#72
length(unique(virion$Host))#1488

virion <- virion %>% rename(Virus_NCBI_TaxID = VirusTaxID)

# Filter spurious non-mammalian and non-ave host-virus associations
# only for viruses in our experimental data ***

to_remove <- virion %>% filter(Virus == "australian bat lyssavirus",
							HostClass == "actinopteri") %>% pull(NCBIAccession)

virion <- virion %>% filter(!NCBIAccession %in% to_remove)

# before merging we have to clean dat to have upham names
# and taxonomy
dat$Host_Upham <- tolower(gsub("_"," ", dat$Host_Upham))
setdiff(dat$Host_Upham, virion$Host_Upham)# 10 mammals (all in upham tree)

dat$Host <- dat$Host_Upham

dat <- dat %>% rename(HostOrder = Host_order)
dat$HostOrder <- tolower(dat$HostOrder)
dat$HostClass <- "mammalia"

# we want to remove the virion Virus names because they are inconsistent
# and cause merge error (e.g. lyssavirus naming)
virion <- virion %>% select(-Virus)

# keep track of which data are experimental infections prior to merge
dat$Experiment <- 1

dat <- full_join(virion, dat)

dat$Susceptible_primary_YN[is.na(dat$Susceptible_primary_YN)] <- "Y"

# experimental infections with host to genus only
to_remove <- which(dat$HostClass=="mammalia" & is.na(dat$Host_Upham))
dat <- dat[-to_remove,]

# calculate host range metrics, whether they infect humans / bats / rodents

hr <- dat %>% filter(Susceptible_primary_YN=="Y" & !is.na(Host)) %>% 
					group_by(Virus_NCBI_TaxID) %>%
					mutate( n_classes = n_distinct(HostClass, na.rm = T),
							n_mammals = n_distinct(Host_Upham[HostClass=="mammalia"], na.rm = T),
							n_orders = n_distinct(HostOrder, na.rm = T),
							hostRichness = n_distinct(Host, na.rm = T),
							batHosts 	= n_distinct(Host_Upham[HostOrder=="chiroptera"], na.rm = T),
							rodentHosts = n_distinct(Host_Upham[HostOrder=="rodentia"], na.rm = T),
							zoonotic	= any(Host%in%"homo sapiens")*1) %>%							
					select(Virus_NCBI_TaxID, hostRichness, n_mammals, n_orders, n_classes, 
							batHosts, rodentHosts, zoonotic) %>%
					unique()

# View(hr)

# add back Virus_ICTV names
viruses <- read.csv("../raw_data/VirusNames_translation_Feb23_2024.csv")
viruses <- select(viruses, Virus_ICTV, Virus_NCBI_TaxID) %>% unique()

hr <- left_join(viruses, hr) %>% arrange(Virus_ICTV)

hr <- hr[!is.na(hr$Virus_ICTV),]

write.csv(hr, "../clean_data/host_range_data.csv", row.names=FALSE)


# export full h-v list for data plot

dat <- left_join(dat, viruses)

edge_list <- dat %>% select(Host, Virus_NCBI_TaxID) %>% 
						unique() %>% arrange(Virus_NCBI_TaxID, Host)

edge_list <- left_join(edge_list, viruses)

# write.csv(edge_list, "../clean_data/edge_list.csv", row.names=FALSE) 


edge_list_sm <- dat %>% filter(HostOrder%in%c("chiroptera", "rodentia")) %>%
						select(Host_Upham, Virus_NCBI_TaxID) %>% 
						unique() %>% arrange(Virus_NCBI_TaxID, Host_Upham)

edge_list_sm <- left_join(edge_list_sm, viruses)

edge_list_sm <- edge_list_sm %>% filter(!is.na(Host_Upham))

# write.csv(edge_list_sm, "../clean_data/edge_list_bats_rodents.csv", row.names=FALSE) 





# calculate phylo metrics / evolIso / distance to reservoir

# calculating diversity metrics

edge_list_mamm <- dat %>% filter(HostClass=="mammalia" & Susceptible_primary_YN=="Y") %>%
						select(Host, Virus_NCBI_TaxID) %>% 
						unique() %>% arrange(Virus_NCBI_TaxID, Host)

names(edge_list_mamm)

com <- table(edge_list_mamm$Host, edge_list_mamm$Virus_NCBI_TaxID)
com[com>1] <- 1
dim(com)# 834.. lost one because it only appears as a non-susceptible experimental host
rownames(com)

setdiff(rownames(com), tree$tip.label)

# subset tree to only hosts in com
tree_sub <- drop.tip(tree, setdiff(tree$tip.label,rownames(com)))


# re-order com_mamm to match tree
com <- com[tree_sub$tip.label,]

all(rownames(com)==tree_sub$tip.label)#TRUE

# transpose to use mpd.query
com_t <- t(com)

pd <- pd.query(tree_sub, com_t, standardize = FALSE, null.model="uniform", 
				reps=1000, seed=1989)

ses_pd <- pd.query(tree_sub, com_t, standardize = TRUE, null.model="uniform", 
				reps=1000, seed=1989)

mpd <- mpd.query(tree_sub, com_t, standardize = FALSE, null.model="uniform", 
				reps=1000, seed=1989)

ses_mpd <- mpd.query(tree_sub, com_t, standardize = TRUE, null.model="uniform", 
				reps=1000, seed=1989)

spec_vals <- data.frame(Virus_NCBI_TaxID=as.numeric(rownames(com_t)), pd=pd, ses_pd=ses_pd, mpd=mpd, ses_mpd=ses_mpd)


# calculate evolutionary isolation

evoiso <- dat %>% filter(Experiment==1) %>% 
						select(Host_Upham, Virus_NCBI_TaxID) %>%
						filter(!is.na(Host_Upham)) %>% unique()

dim(evoiso)#242

names(evoiso)

for(i in 1:nrow(evoiso)){

  focal_host <- evoiso$Host_Upham[i]
  para <- evoiso[i,2]
  hosts <- unique(c(focal_host, rownames(com)[which(com[,colnames(com)%in%para]==1)]))
  
  if(length(hosts)==1) {
                evoiso$EvoIso[i] <-0
  } else {
      pruned_tree <- drop.tip(tree,tree$tip.label[-match(hosts, tree$tip.label)])
      d_mat_mamm  <- cophenetic(pruned_tree)
      d_mat_mamm  <- d_mat_mamm[rownames(d_mat_mamm)==focal_host, colnames(d_mat_mamm)!=focal_host]
      evoiso$EvoIso[i] <- mean(d_mat_mamm)
  }
}

evoiso[evoiso$EvoIso>200,]

spec_vals <- left_join(spec_vals, evoiso)

hr_2 <- left_join(hr, spec_vals)


# now add information on reservoirs....
# average distance to the known reservoir clade...
# to do this, pull in list of reservoirs, identify single representative taxa
# calculate distance from this taxa to the infected animal
# if the infected animal is in the same taxonomic group, set this to 0
# (this is necessary because within a clade this will be impacted by choice of representative species)

reservoir <- read.csv("../clean_data/Virus_traits.csv")

reservoir <- select(reservoir, Virus_ICTV, Reservoir) %>% unique()

# get representative species
reservoir$Reservoir[reservoir$Reservoir==NA] <- "Orphan"
reservoir$Reservoir[reservoir$Reservoir%in%
				c("Anseriformes","Charadriiformes","Columbiformes",
					"Galliformes","Passeriformes","Pelecaniformes",
					"Suliformes")] <- "AVES"

reservoir$Reservoir[reservoir$Reservoir%in%"Vespbat"] <- "VESPERTILIONIDAE"
reservoir$Reservoir[reservoir$Reservoir%in%"Pterobat"] <- "PTEROPODIDAE"
reservoir$Reservoir[reservoir$Reservoir%in%"Artiodactyl"] <- "CETARTIODACTYLA"
reservoir$Reservoir[reservoir$Reservoir%in%"Perissodactyla"] <- "PERISSODACTYLA"
reservoir$Reservoir[reservoir$Reservoir%in%"Carnivore"] <- "CARNIVORA"
reservoir$Reservoir[reservoir$Reservoir%in%"Human"] <- "homo sapiens"
reservoir$Reservoir[reservoir$Reservoir%in%"NonHumanPrimate"] <- "PRIMATES"
reservoir$Reservoir[reservoir$Reservoir%in%"Rodent"] <- "RODENTIA"
reservoir$Reservoir[reservoir$Reservoir%in%"Diprotodontia"] <- "DIPROTODONTIA"

res_reps <- data.frame(Reservoir=sort(unique(reservoir$Reservoir)))
res_reps$rep_species <- NA
res_reps$rep_species[res_reps$Reservoir%in%"homo sapiens"] <- "homo sapiens"

for (i in seq_along(res_reps$Reservoir)){

	if (res_reps$Reservoir[i]%in%tax$ord) res_reps$rep_species[i] <- sample(tax$Species_Name[tax$ord%in%res_reps$Reservoir[i]], 1)

	if (res_reps$Reservoir[i]%in%tax$fam) res_reps$rep_species[i] <- sample(tax$Species_Name[tax$fam%in%res_reps$Reservoir[i]], 1)

}

res_reps
reservoir <- left_join(reservoir, res_reps)

hr_3 <- left_join(hr_2, reservoir)

hr_3$res_dist <- NA

hr_3[is.na(hr_3$Host_Upham),]
# One host is NA because it wasn't reported to species level
hr_3 <- hr_3[!is.na(hr_3$Host_Upham),]


for(i in 1:nrow(hr_3)){

	if(!is.na(hr_3$rep_species[i])){

		focal_host <- hr_3$Host_Upham[i]
		para <- hr_3$Virus_NCBI_TaxID[i]
		hosts <- unique(c(focal_host, hr_3$rep_species[i]))
		  
		  if(length(hosts)==1) {
		                hr_3$res_dist[i] <-1
		  } else {

		      pruned_tree <- drop.tip(tree, tree$tip.label[-match(hosts, tree$tip.label)])
		      d_mat_mamm  <- cophenetic(pruned_tree)
		      d_mat_mamm  <- d_mat_mamm[rownames(d_mat_mamm)==focal_host, colnames(d_mat_mamm)!=focal_host]
		      hr_3$res_dist[i] <- mean(d_mat_mamm)
		  }
	}
}

# now set distances to zero when experimental host is from the same order as reservoir
tax2 <- tax
tax2 <- select(tax2, Species_Name, ord, fam)
names(tax2)[1] <- "Host_Upham"
hr_3 <- left_join(hr_3, tax2)

hr_3$res_dist[!is.na(hr_3$res_dist) & hr_3$Reservoir==hr_3$ord] <- 0
hr_3$res_dist[!is.na(hr_3$res_dist) & hr_3$Reservoir==hr_3$fam] <- 0

# manually filling in others
# mrca for mammals/ birds was ~ 310 MYA
hr_3$res_dist[hr_3$Reservoir%in%"AVES"] <- 620

# average reservoir distances for viruses with multiple reported reservoirs
mean_res_dist <- hr_3 %>% filter(!is.na(res_dist)) %>% 
							select(Virus_ICTV, Reservoir, res_dist) %>% 
							unique() %>% group_by(Virus_ICTV) %>% 
							mutate(mean_reservoir_dist = mean(res_dist))


hr_3 <- left_join(hr_3, mean_res_dist)

hr_3 <- select(hr_3, -c(Reservoir, rep_species, res_dist, ord, fam))
hr_3 <- unique(hr_3)

# re-format Host_Upham

hr_3$Host_Upham <- Hmisc::capitalize(hr_3$Host_Upham)
hr_3$Host_Upham <- gsub(" ", "_", hr_3$Host_Upham)

# overwrite original host range data
write.csv(hr_3, "../clean_data/host_range_data.csv", row.names=FALSE)



