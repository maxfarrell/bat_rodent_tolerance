# dose_harmonize.R

require(dplyr);packageVersion("dplyr")#1.1.4
require(ggplot2);packageVersion("ggplot2")#3.4.4
require(cowplot);packageVersion("cowplot")#1.1.1

require(brms);packageVersion("brms")#2.21.6
require(bayesplot);packageVersion("bayesplot")#1.8.1
require(ape);packageVersion("ape")#5.5


# data
primary <- read.csv("../raw_data/primary_screen.csv")
dat <- read.csv("../raw_data/individual_data.csv")
vtax <- read.csv("../clean_data/Virus_taxonomy.csv")
viruses <- read.csv("../raw_data/VirusNames_translation_Feb23_2024.csv")
viruses <- dplyr::left_join(viruses, vtax)
hosts <- read.csv("../raw_data/HostNames_NCBI_Upham.csv")
hosts$Host_Upham <- gsub(" ", "_", hosts$Host_Upham)
#subset to only host reported and upham
hosts <- hosts %>% select(Host_reported, Host_Upham) %>% unique()

# construct column for individualID
dat <- dat %>%                              # Create numbering variable
		  group_by(PaperID) %>%
  		  mutate(IndividualID = paste(row_number(),unique(PaperID), sep="_"))


# add names
primary <- dplyr::left_join(primary, viruses)
primary <- dplyr::left_join(primary, hosts)
dat <- dplyr::left_join(dat, hosts)
dat <- dplyr::left_join(dat, viruses)

# remove NA viruses
primary <- primary[!is.na(primary$Virus_ICTV),]
dat <- dat[!is.na(dat$Virus_ICTV),]

# adding mean body size
# via COMBINE https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecy.3344

combine <- read.csv("../raw_data/COMBINE/trait_data_imputed.csv")
combine$phylacine_binomial <- gsub(" ","_",combine$phylacine_binomial)

# intersect(dat$Host_Upham, combine$phylacine_binomial)#82
# setdiff(dat$Host_Upham, combine$phylacine_binomial)# Tamias_amoenus, Liomys_salvini, NA
combine$phylacine_binomial[grep("Neotamias_amoenus",combine$phylacine_binomial)] <- "Tamias_amoenus"
combine$phylacine_binomial[grep("Heteromys_salvini",combine$phylacine_binomial)] <- "Liomys_salvini"

combine$Host_Upham <- combine$phylacine_binomial
bodymass <- select(combine, Host_Upham, adult_mass_g) %>% unique()

# for duplicated host species take average of mass
bodymass <- bodymass %>% group_by(Host_Upham) %>% mutate(mean_mass=mean(adult_mass_g))
bodymass$adult_mass_g <- bodymass$mean_mass
bodymass <- select(bodymass, -mean_mass)
bodymass <- unique(bodymass)

dat <- left_join(dat, bodymass)


# cleaning dose amount and converting to numeric
dat$Dose_amount[dat$Dose_amount=="Unknown"] <- NA

dat$Dose_amount <- sapply(dat$Dose_amount, function(txt) eval(parse(text=txt)))

# harmonizing dose unit
sort(unique(dat$Dose_unit))
dat$Dose_unit[dat$Dose_unit%in%
				c("Adult mouse LD50","Median MICLD50","Median_MICLD50", "MICLD50",
				"MICLD50_newborns","MICLD50_suckling","MIPLD50","Mouse LD50")] <- "Mouse LD50"

dat$Dose_unit[dat$Dose_unit%in%
				c("Mean_TCID50","Median_TCID50","TCID50")] <- "TCID50"


sort(table(dat$Dose_unit), decreasing=TRUE)
# 1262 Mouse LD50
# 567 PFU
# 540 TCID50
# 116 FFU


# Harmonize Inoculation Route
dat <- dat %>% tidyr::separate(col=Inoculation_route, 
									into=c("Route_type","Route_location"), 
									sep=" [(]", extra="merge")

sort(unique(dat$Route_type))
sort(unique(dat$Route_location))

dat$Route_type[grep("Intranasal", dat$Route_location)] <- "Subcutaneous; Intranasal"
dat$Route_type[grep("Contact", dat$Route_location)] <- "Bite; Contact"

sort(unique(dat$Route_type))



# less stringent collapsing of route types

dat_rt_lrg <- dat

# lower / delayed severity (localized & peripheral)
dat_rt_lrg$Route_type[dat_rt_lrg$Route_type%in%
				c("Subcutaneous", "Intradermal","Subdermal",
					"Footpad")] <- "Subcutaneous / Subdermal / Intradermal"

# Intraperitoneal
dat_rt_lrg$Route_type[dat_rt_lrg$Route_type%in%
				c("Intraperitoneal")] <- "Intraperitoneal"

# Intramuscular
dat_rt_lrg$Route_type[dat_rt_lrg$Route_type%in%
				c("Intramuscular")] <- "Intramuscular"

# Intraperitoneal
dat_rt_lrg$Route_type[dat_rt_lrg$Route_type%in%
				c("Intraperitoneal")] <- "Intraperitoneal"

# Intravenous
dat_rt_lrg$Route_type[dat_rt_lrg$Route_type%in%
				c("Intravenous")] <- "Intravenous"

# Intracardial
dat_rt_lrg$Route_type[dat_rt_lrg$Route_type%in%
				c("Intracardial")] <- "Intracardial"

# Injection into brain:
dat_rt_lrg$Route_type[dat_rt_lrg$Route_type%in%
				c("Intracerebral","Intracranial")] <- "Injection (brain)"


# Oral / nasal inoculation:
dat_rt_lrg$Route_type[dat_rt_lrg$Route_type%in%
				c("Oronasal","Oral","Intranasal", "Nasal")] <- "Oronasal"

# Simluated "Natural" infections
dat_rt_lrg$Route_type[dat_rt_lrg$Route_type%in%
				c("Transplacental","Contact","Skin","Vector",
					"Bite","Bite; Contact", "Aerosol","Natural")] <- "Natural"

# setting others (multiple routes or vague "injected") to NA
dat_rt_lrg$Route_type[dat_rt_lrg$Route_type%in%
				c("Injected - no further details","Subcutaneous; Intranasal",
					"Intracardial; Intraperitoneal", "Subcutaneous; Intraperitoneal",
					"Intraperitoneal; Intravenous")] <- NA

sort(table(dat_rt_lrg$Route_type), decreasing=TRUE)

sum(is.na(dat_rt_lrg$Route_type))#122



colors <- c("#648FFF","#785EF0","#DC267F","#FE6100","#FFB000", "#949698")


dat_rt_lrg_plot <- dat_rt_lrg[!is.na(dat_rt_lrg$Route_type),]


rt_stack <- ggplot(dat_rt_lrg_plot, aes(fill=Route_type,  x=Host_order)) + 
		geom_bar(position="fill", stat="count") +
		# geom_point(position = position_jitter(seed = 1, width = 0.2), 
		# 	alpha=0.4) + 
		# scale_y_continuous(trans='log10') +
		scale_fill_manual(values=colors) + scale_color_manual(values=colors) +
		# theme(legend.position = "none",legend.title=element_blank()) + 
		xlab("") + ylab("Proportion of individuals")  

rt_stack

# View(dat_rt_lrg_plot)


# merging with severity data to plot severity by inoculation route

dat_symptoms <- read.csv("../clean_data/symptoms_severity.csv")

dat_rt_lrg_plot <- left_join(dat_symptoms, dat_rt_lrg_plot)
dat_rt_lrg_plot <- dat_rt_lrg_plot[!is.na(dat_rt_lrg_plot$Dose_unit),]

colours_BR <- c("#1B85BF", "#AB1808")

p5 <- ggplot(dat_rt_lrg_plot, aes(y=severity_rank, x=Route_type, fill=Host_order, color=Host_order, shape=Host_order)) + 
		# geom_violin(alpha=0.4) +
		geom_point(position = position_jitter(seed = 1, height = 0.15), 
			alpha=0.4) + 
		# scale_x_continuous(trans='log10') +
		scale_fill_manual(values=colours_BR) + scale_color_manual(values=colours_BR) +
		# theme(legend.position = "none") + 
	    labs(color="  Host order", shape="  Host order", fill="  Host order") + ylab("Severity") + xlab("Inoculation Route")
		# facet_grid(cols=vars(Route_type)) + theme(legend.key=element_blank()) 
				# stat_summary(fun = "mean",
    #            geom = "crossbar", 
    #            width = 0.7,
    #            colour = "black")
p5


# barplot version 

# 1. Pre-calculate the sample sizes per bar (Host + Route combination)
bar_labels <- dat_rt_lrg_plot %>%
  group_by(Host_order, Route_type) %>%
  summarize(
    n_individuals = n(),
    n_studies = n_distinct(PaperID),
    .groups = "drop"
  ) %>%
  # Create a clean string label to print on the plot
  mutate(label_text = paste0("Individuals = ", n_individuals, "\nStudies = ", n_studies))

# 2. Plot with text overlays
p5_new <- ggplot(dat_rt_lrg_plot, aes(x = Route_type)) +
  # Use fill mapping inside geom_bar so it only applies to the bars, not the text labels
  geom_bar(aes(fill = factor(severity_rank)), position = "fill", width = 0.7) + 
  
  # Add the text layer using our pre-calculated label data frame
  geom_text(
    data = bar_labels,
    aes(x = Route_type, y = 1.02, label = label_text),
    inherit.aes = FALSE,   # Prevents it from looking for severity_rank in bar_labels
    vjust = 0,             # Aligns text so it sits on top of the y-coordinate
    size = 3,              # Adjust text size here
    lineheight = 0.85      # Tightens the vertical space between the N and S lines
  ) +
  
  facet_wrap(~ Host_order) +
  
  # Expand the upper limit of the Y axis slightly (to 1.1) so the text labels don't get cut off
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.1), breaks = seq(0, 1, 0.25)) +
  
  scale_fill_brewer(palette = "YlOrRd") + 
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(1.5, "lines")
  ) +
  labs(
    y = "Proportion of Individuals",
    x = "Inoculation Route",
    fill = "Severity Rank"
  )

# Render the plot
p5_new

ggsave("../plots_tables/severity_by_inoculationroute.pdf", p5_new, width=15, height=7)

# # saving as dat
# dat <- dat_rt_lrg_plot



# setting new harmonized inoculation routes as main dat
dat <- dat_rt_lrg


# adjusting for body mass
# dat_sm$Dose_mass <- log10(dat_sm$Dose_amount) / log10(dat_sm$adult_mass_g)
# dat_sm$Dose_mass <- log10(dat_sm$Dose_amount) / (dat_sm$adult_mass_g)
dat$Dose_mass <- (dat$Dose_amount) / (dat$adult_mass_g)
hist(dat$Dose_mass)


# Graphs with only major dose units
dat_sm <- filter(dat, Dose_unit%in%c("Mouse LD50","PFU","TCID50","FFU"))

# keeping only sucessful infections  
dat_sm <- filter(dat_sm, Susceptible_YN%in%c("Y"))

dat_sm$Dose_unit <- factor(dat_sm$Dose_unit, levels = c("Mouse LD50", "PFU", "TCID50","FFU"))



colors <- c("#648FFF","#785EF0","#DC267F","#FE6100","#FFB000")

p <- ggplot(dat_sm, aes(Host_order, Dose_amount, fill=factor(Dose_unit), color=factor(Dose_unit))) + 
		geom_violin(alpha=0.4) +
		geom_point(position = position_jitter(seed = 1, width = 0.2), 
			alpha=0.4) + 
		scale_y_continuous(trans='log10') +
		scale_fill_manual(values=colors) + scale_color_manual(values=colors) +
		theme(legend.position = "none") + xlab("Host order") + ylab("Dose") + 
		facet_grid(cols=vars(Dose_unit)) + 
				stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.7,
               colour = "black")
p

ggsave("../plots_tables/Doses_by_type.pdf", p, width=8, height=5)
ggsave("../plots_tables/Doses_by_type.png", p, width=8, height=5)


# per body size
p2 <- ggplot(dat_sm, aes(Host_order, Dose_mass, fill=Dose_unit, color=Dose_unit)) + 
		geom_violin(alpha=0.4) +
		geom_point(position = position_jitter(seed = 1, width = 0.2), 
			alpha=0.4) + 
		scale_y_continuous(trans='log10') +
		scale_fill_manual(values=colors) + scale_color_manual(values=colors) +
		# theme_bw() +
		theme(legend.position = "none") + xlab("Host order") + ylab("Dose/g") + 
		facet_grid(cols=vars(Dose_unit)) + 
				stat_summary(fun = "mean",   
            geom = "crossbar", 
               width = 0.7,
               colour = "black")
p2

ggsave("../plots_tables/Dose_per_g_by_type.pdf", p2, width=8, height=5)
ggsave("../plots_tables/Dose_per_g_by_type.png", p2, width=8, height=5)


# Dose per viral family

p3 <- ggplot(dat_sm, aes(Host_order, Dose_mass, fill=factor(Dose_unit), color=factor(Dose_unit))) + 
		geom_violin(alpha=0.4) +
		# geom_point(position = position_jitter(seed = 1, width = 0.1), 
		# 	alpha=0.4) + 
		scale_y_continuous(trans='log10') +
		scale_fill_manual(values=colors) + scale_color_manual(values=colors) +
		theme(legend.position = "bottom", legend.title=element_blank()) + ylab("Dose/g") + xlab("") + 
		facet_grid(cols=vars(family)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
		# + 
		# 		stat_summary(fun = "mean",
  #              geom = "crossbar", 
  #              inherit.aes = TRUE)
p3

ggsave("../plots_tables/Dose_by_viral_fam.pdf", p3, height=4, width=16)
ggsave("../plots_tables/Dose_by_viral_fam.png", p3, height=4, width=16)


p4 <- ggplot(dat_sm, aes(Host_order, Dose_mass, fill=factor(Dose_unit), color=factor(Dose_unit))) + 
		geom_violin(alpha=0.4) +
		# geom_point(position = position_jitter(seed = 1, width = 0.1), 
		# 	alpha=0.4) + 
		scale_y_continuous(trans='log10') +
		scale_fill_manual(values=colors) + scale_color_manual(values=colors) +
		theme(legend.position = "bottom", legend.title=element_blank()) + ylab("Dose/g") + xlab("") + 
		facet_grid(cols=vars(order)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
		# + 
		# 		stat_summary(fun = "mean",
  #              geom = "crossbar", 
  #              inherit.aes = TRUE)
p4

ggsave("../plots_tables/Dose_by_viral_order.pdf", p4, height=4, width=12)
ggsave("../plots_tables/Dose_by_viral_order.png", p4, height=4, width=12)

dose_dat <- select(dat, PaperID, IndividualID, Virus_ICTV, Host_Upham,
					adult_mass_g, Dose_amount, Dose_unit, Dose_mass, Route_type) %>% unique()

dim(dose_dat)#2885

# updating to reflect inoculation data is in this
write.csv(dose_dat, "../clean_data/dose_inoculation_data.csv", row.names=FALSE)


# modelling dose differences between bats and rodents

bayesplot_theme_set(theme_bw())

dat_dose <- dose_dat

dat_dose$Host_name <- dat_dose$Host_Upham
dat_dose$Virus_name <- dat_dose$Virus_ICTV


# host data
tree <- read.nexus("../raw_data/upham_tree_666.nex")
hosts <- read.csv("../raw_data/HostNames_NCBI_Upham.csv")
hosts$Host_Upham <- gsub(" ", "_", hosts$Host_Upham)

# virus names and taxonomy
vtax <- read.csv("../clean_data/Virus_taxonomy.csv")
viruses <- read.csv("../raw_data/VirusNames_translation_Feb23_2024.csv")
viruses <- left_join(viruses, vtax)

#subset to only host reported and upham
hosts <- hosts %>% select(Host_reported, Host_Upham) %>% unique()

# make subset tree to sampled hosts
host_tree <- drop.tip(tree, setdiff(tree$tip.label, dat$Host_Upham))

# make virus tree
vtax <- as.data.frame(unclass(vtax), stringsAsFactors=TRUE)
frm <- ~superkingdom/realm/kingdom/phylum/class/order/family/genus/Virus_ICTV
vtree <- as.phylo(frm, data = vtax, collapse=FALSE)
vtree$edge.length <- rep(1, nrow(vtree$edge))

# include only viruses in subset data (e.g. after removing mole)
vtree <- drop.tip(vtree, setdiff(vtree$tip.label, dat$Virus_ICTV))
# plot(vtree)

# phylogenetic correlation structure
host_cov <- vcv(host_tree, corr=TRUE)
virus_cov <- vcv(vtree, corr=TRUE)


# adding host_order

dat_dose <- left_join(dat_dose, select(dat, c(Host_Upham, Host_order), ))



# merging with severity data to plot severity by dose

dat_symptoms <- read.csv("../clean_data/symptoms_severity.csv")

dat <- left_join(dat_symptoms, dat_sm)
dat <- dat[!is.na(dat$Dose_unit),]

colours_BR <- c("#1B85BF", "#AB1808")

p5 <- ggplot(dat, aes(y=severity_rank, x=Dose_mass, fill=Host_order, color=Host_order, shape=Host_order)) + 
		# geom_violin(alpha=0.4) +
		geom_point(position = position_jitter(seed = 1, height = 0.15), 
			alpha=0.4) + 
		scale_x_continuous(trans='log10') +
		scale_fill_manual(values=colours_BR) + scale_color_manual(values=colours_BR) +
		# theme(legend.position = "none") + 
	    labs(color="  Host order", shape="  Host order", fill="  Host order") + ylab("Severity") + xlab("Dose/g") + 
		facet_grid(cols=vars(Dose_unit)) + theme(legend.key=element_blank()) 
				# stat_summary(fun = "mean",
    #            geom = "crossbar", 
    #            width = 0.7,
    #            colour = "black")
p5

ggsave("../plots_tables/severity_by_dose.pdf", p5, width=10, height=4)
ggsave("../plots_tables/severity_by_dose.png", p5, width=10, height=4)
