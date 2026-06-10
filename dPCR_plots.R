
library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(janitor)

sample_mapping <- read_delim("~/Desktop/Spartina/Spartina2025/Rennes_Sampling.tsv") %>% 
  dplyr::select(`Sample number`, "Species") %>%
  set_colnames(c("User sample ID", "species"))

root_weights <- readxl::read_excel("JMF-2508-06_rootweights.xlsx") %>%
  mutate(`User sample ID` = str_split_i(`User sample ID`, " ", 1)) %>% 
  filter(`Sample description` == "roots") %>%
  dplyr::select(`User sample ID`, `Biomass in g`)

# read in ddpcr results and change relevant column to numeric
ddpcr <- read.csv("~/Desktop/Spartina/Spartina2025/ddPCR/DOME-SD-SPARTINA16S-RUN1_analysis_31_03_2026_08_33_15_UTC+02_00.csv", skip = 2)
ddpcr$Conc...cp.µL...undiluted.sample. %<>% as.numeric()

# read in PCR results and JMF samplesheet to cross-reference sample names
pcr <- readxl::read_xlsx("~/Desktop/Spartina/Spartina2025/ddPCR/KatiedPCR.xlsx")
jmf <- readxl::read_xlsx("~/Desktop/Spartina/Spartina2025/ddPCR/JMF-2508-06-1.xlsx")
pcr$`Sample.NTC.Control` <- pcr$`JMF sample ID` %>% str_sub(-2,-1) %>% as.integer() %>% as.character()

concentrations <- inner_join(ddpcr, pcr, by="Sample.NTC.Control") %>% 
  dplyr::select("User sample ID", `Conc...cp.µL...undiluted.sample.`) %>% 
  mutate(`User sample ID` = str_remove(`User sample ID`, " r"))


results <- inner_join(concentrations, sample_mapping, by="User sample ID") %>%
  mutate(species=case_when(species == "Spartina anglica" ~ "Sporobolus anglicus",
                           species == "Spartina alternifllora" ~ "Sporobolus alterniflorus",
                           species == "Spartina maritima" ~ "Sporobolus maritimus")) %>%
  rename("Conc...cp.µL...undiluted.sample." = "copies_per_ul") %>%
  left_join(root_weights)

# calculate copies/ul for the reaction volume 
# using these forumlas from Astrid at DOME
# value from software * (reaction volume/template volume) * dilution factor = cp/µl in original sample tube
# For example:  2038 cp/µl * (12 µl/2 µl) * 100 = 1 222 800 cp/µl in original sample tube
dilution_factor <- 100000
results %<>% mutate(copies_per_ul_converted=(copies_per_ul * (12/2) * dilution_factor)/`Biomass in g`)

png("FigureS2_16S_copies_barplot.png", height=500, width=500)
ggplot(results, aes(x=species, y=(copies_per_ul_converted), fill=species)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("brown2", "palegreen3", "dodgerblue2")) +
  theme(axis.title = element_text(size=20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=20),
        strip.background =element_rect(fill="skyblue2"),
        legend.title = element_blank(),
        legend.text = element_text(size=15)) +
  ylab("16S Copies per Microlitre")
dev.off()



##################################################################
#.                  checking dPCR with Astrids help              #
##################################################################

# Read in RDS and sample data
phylo_rennes <- readRDS("/Users/katieemelianova/Desktop/Spartina/JMF_results/JMF-2508-06_16S_raw_phyloseq.rds")
sample_info <- read_tsv("/Users/katieemelianova/Desktop/Spartina/Spartina2025//Rennes_Sampling.tsv") %>% 
  dplyr::select(Locality, Species, `Sample number`)


  


# abundance and prevalence filter
phylo_rennes <- tax_filter(
  phylo_rennes,
  min_prevalence = 10,
  prev_detection_threshold = 5,
  min_total_abundance = 5,
  min_sample_abundance = 5)

# sample names currently in format "1058 r", split this and set sample ID as the number 
# and compartment as the letter translated to the compartment name
phylo_rennes@sam_data$User_sample_ID_number <- phylo_rennes@sam_data %>% data.frame() %>% pull(User_sample_ID) %>% str_split_i(" ", 1)
phylo_rennes@sam_data$compartment <- phylo_rennes@sam_data %>% data.frame() %>% pull(User_sample_ID) %>% str_split_i(" ", 2) %>% 
  str_replace("root", "r") %>%
  str_replace("rhizome", "z") %>%
  str_replace("soil", "s") %>%
  str_replace("r", "Root") %>%
  str_replace("z", "Rhizome") %>%
  str_replace("s", "Rhizosphere") %>% 
  replace_na('Unknown')
phylo_rennes@sam_data$Locality <- phylo_rennes@sam_data %>% data.frame() %>% left_join(sample_info, by = c("User_sample_ID_number" = "Sample number")) %>% pull(Locality) %>% replace_na('Unknown')
phylo_rennes@sam_data$Species <- phylo_rennes@sam_data %>% data.frame() %>% left_join(sample_info, by = c("User_sample_ID_number" = "Sample number")) %>% pull(Species) %>% replace_na('Unknown')
phylo_rennes@sam_data$Species %<>% 
  str_replace("Spartina alternifllora", "Sporobolus alterniflorus") %>%
  str_replace("Spartina maritima", "Sporobolus maritimus") %>%
  str_replace("Spartina anglica", "Sporobolus anglicus")

# get root only and select oiut only mitochondrial seqs
phylo_rennes <- subset_samples(phylo_rennes, compartment == "Root")
phylo_rennes_plastid <- subset_taxa(phylo_rennes, Family %in% c("Mitochondria", "Chloroplast") | Order %in% c("Mitochondria", "Chloroplast"))
phylo_rennes_bacteria <- subset_taxa(phylo_rennes, !(Family %in% c("Mitochondria", "Chloroplast") | Order %in% c("Mitochondria", "Chloroplast")))


# get relative abundance
phylo_rennes_prop_plastid <- transform_sample_counts(phylo_rennes_plastid, function(otu) otu/sum(otu))
phylo_rennes_prop_bacteria <- transform_sample_counts(phylo_rennes_bacteria, function(otu) otu/sum(otu))


png("spartina_absolute_plastid_abundance.png", height=500, width=500)
phylo_rennes_prop_plastid %>% ps_melt() %>% 
  dplyr::select(OTU, Abundance, User_sample_ID_number, sample_Species) %>%
  rename("User_sample_ID_number" = "User sample ID") %>%
  left_join(results, by = "User sample ID", relationship = "many-to-many") %>%
  mutate(absolute_abundance = Abundance * copies_per_ul_converted) %>%
  drop_na() %>%
  ggplot(aes(x=sample_Species, y=log(absolute_abundance))) + 
  geom_boxplot(aes(fill = sample_Species)) +
  ylab("Absolute Abundance Plastid") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5),
        axis.text = element_text(size=15),
        axis.title.y = element_text(size = 15),
        legend.title=element_blank()) +
  scale_fill_manual(values=c("brown2", "palegreen3", "dodgerblue2"))
dev.off()

png("spartina_absolute_bacteria_abundance.png", height=500, width=500)
phylo_rennes_prop_bacteria %>% ps_melt() %>% 
  dplyr::select(OTU, Abundance, User_sample_ID_number, sample_Species) %>%
  rename("User_sample_ID_number" = "User sample ID") %>%
  left_join(results, by = "User sample ID", relationship = "many-to-many") %>%
  mutate(absolute_abundance = Abundance * copies_per_ul_converted) %>%
  drop_na() %>%
  ggplot(aes(x=sample_Species, y=log(absolute_abundance))) + 
  geom_boxplot(aes(fill = sample_Species)) +
  ylab("Absolute Abundance Bacteria") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5),
        axis.text = element_text(size=15),
        axis.title.y = element_text(size = 15),
        legend.title=element_blank()) +
  scale_fill_manual(values=c("brown2", "palegreen3", "dodgerblue2")) 
dev.off()














