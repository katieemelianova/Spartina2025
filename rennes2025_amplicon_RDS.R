library(microViz)
library(DESeq2)

#####################################################
#           Read in RDS and sample data             # 
#####################################################

phylo_rennes <- readRDS("/Users/katieemelianova/Desktop/Spartina/JMF_results/JMF-2508-06_16S_raw_phyloseq.rds")
sample_info <- read_tsv("/Users/katieemelianova/Desktop/Spartina/Spartina2025//Rennes_Sampling.tsv") %>% 
  dplyr::select(Locality, Species, `Sample number`)


######################################################################################
#    filter on prevalence, set metadata and remove chloroplast and mitochondria      # 
######################################################################################

# abundance and prevalence filter
phylo_rennes <- tax_filter(
  phylo_rennes,
  min_prevalence = 10,
  prev_detection_threshold = 5,
  min_total_abundance = 5,
  min_sample_abundance = 5)

# sample names currently in format "1058 r", split this and set sample ID as the number 
# and compartment as the letter translated to the compartment name
phylo_rennes@sam_data
phylo_rennes@sam_data$User_sample_ID_number <- phylo_rennes@sam_data %>% data.frame() %>% pull(User_sample_ID) %>% str_split_i(" ", 1)
phylo_rennes@sam_data$compartment <- phylo_rennes@sam_data %>% data.frame() %>% pull(User_sample_ID) %>% str_split_i(" ", 2) %>% 
  str_replace("root", "r") %>%
  str_replace("rhizome", "z") %>%
  str_replace("soil", "s") %>%
  str_replace("r", "root") %>%
  str_replace("z", "rhizome") %>%
  str_replace("s", "soil") %>% 
  replace_na('Unknown')

# set locality and species as metadata slots
phylo_rennes@sam_data$Locality <- phylo_rennes@sam_data %>% data.frame() %>% left_join(sample_info, by = c("User_sample_ID_number" = "Sample number")) %>% pull(Locality) %>% replace_na('Unknown')
phylo_rennes@sam_data$Species <- phylo_rennes@sam_data %>% data.frame() %>% left_join(sample_info, by = c("User_sample_ID_number" = "Sample number")) %>% pull(Species) %>% replace_na('Unknown')

# set one phyloseq object as chloroplast only to see if we can detect phylogenetic signal
phylo_rennes_chloroplast <- subset_samples(phylo_rennes, Species != "Unknown")
phylo_rennes_chloroplast <- subset_taxa(phylo_rennes_chloroplast, Order %in% c("Chloroplast"))

# set another phyloseq object where we do not have chloroplast or mitochondria in for main analysis
phylo_rennes <- subset_samples(phylo_rennes, Species != "Unknown")
phylo_rennes <- subset_taxa(phylo_rennes, !(Family %in% c("Mitochondria", "Chloroplast"))) %>% subset_taxa(!(Order %in% c("Mitochondria", "Chloroplast")))


######################################
#          Plot richness             # 
######################################

plot_richness(phylo_rennes, x="compartment", measures=c("Shannon"), color="Species") +
  geom_point(size = 3.5) +
  theme(strip.text.x = element_text(size=25),
        axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size=25),
        legend.text = element_text(size=20),
        axis.title.x = element_blank(),
        legend.title = element_blank()) 


######################################
#         plot ordination            # 
######################################

# transform sample to relative abundance
phylo_rennes_prop <- transform_sample_counts(phylo_rennes, function(otu) otu/sum(otu))
ord.nmds.bray_jr <- ordinate(phylo_rennes_prop, method="NMDS", distance="bray")

# plot ordination plot
plot_ordination(phylo_rennes_prop, ord.nmds.bray_jr, color="Species", title="Bray NMDS", shape="compartment") + 
  geom_point(size = 7) +
  theme(strip.text.x = element_text(size=25),
        axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size=25),
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        panel.background = element_blank()) +
  ggtitle("")


# now do a bray curtis plot but using onyly the chloroplast amplicons
# the purpose of this is to confirm the maternal parent of spartina anglica as spartina alterniflora
# also the plot can be used to confirm the expectation that no host dna is included in the 
# rhizosphere/soil samples, as here do not cluster by phylogenetic expectations
# the plot can be included as supplementary
phylo_rennes_chloroplast_prop <- transform_sample_counts(phylo_rennes_chloroplast, function(otu) otu/sum(otu))
ord.nmds.bray_chloroplast <- ordinate(phylo_rennes_chloroplast_prop, method="NMDS", distance="bray")

plot_ordination(phylo_rennes_chloroplast_prop, ord.nmds.bray_chloroplast, color="Species", title="Bray NMDS", shape="compartment") + 
  geom_point(size = 7) +
  theme(strip.text.x = element_text(size=25),
        axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size=25),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20)) +
  ggtitle("")




#####################################################
#           differential abundance analysis         # 
#####################################################

# geometric means function
# doing it this way because see here:
# https://github.com/joey711/phyloseq/issues/445
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


run_deseq <- function(phylo_object){
  phylo_deseq <- phylo_object %>% phyloseq_to_deseq2(~ 0 + Species)
  phylo_deseq_geomeans <- apply(counts(phylo_deseq), 1, gm_mean)
  phylo_deseq = estimateSizeFactors(phylo_deseq, geoMeans = phylo_deseq_geomeans)
  phylo_deseq = DESeq(phylo_deseq, fitType="local")
  return(phylo_deseq)
}

# convert to deseq object specifying variable to test DE on per compartment
# then run deseq, using geo means
# reason for geo means function: https://github.com/joey711/phyloseq/issues/445
phylo_rennes_deseq_root <- subset_samples(phylo_rennes, compartment == "root") %>% run_deseq()
phylo_rennes_deseq_rhizome <- subset_samples(phylo_rennes, compartment == "rhizome") %>% run_deseq()
phylo_rennes_deseq_soil <- subset_samples(phylo_rennes, compartment == "soil") %>% run_deseq()

alpha = 0.01
resultsNames(phylo_rennes_deseq)


results(phylo_rennes_deseq, list("SpeciesSpartina.alternifllora", "SpeciesSpartina.anglica"))

# get significant results sorted by adjusted pvalue
res_alterniflora_maritima = results(phylo_rennes_deseq, contrast=c("Species", "Spartina alternifllora","Spartina maritima")) %>% data.frame() %>% arrange(padj) %>% filter(padj < alpha) %>% rownames_to_column(var="amplicon") %>% left_join(tax_table(phylo_rennes) %>% data.frame() %>% rownames_to_column(var="amplicon"))
res_anglica_maritima = results(phylo_rennes_deseq, contrast=c("Species", "Spartina anglica","Spartina maritima")) %>% data.frame() %>% arrange(padj) %>% filter(padj < alpha) %>% rownames_to_column(var="amplicon") %>% left_join(tax_table(phylo_rennes) %>% data.frame() %>% rownames_to_column(var="amplicon"))
res_anglica_alterniflora = results(phylo_rennes_deseq, contrast=c("Species", "Spartina alternifllora","Spartina anglica")) %>% data.frame() %>% arrange(padj) %>% filter(padj < alpha) %>% rownames_to_column(var="amplicon") %>% left_join(tax_table(phylo_rennes) %>% data.frame() %>% rownames_to_column(var="amplicon"))


results(phylo_rennes_deseq, contrast=c("Species", "Spartina alternifllora","Spartina maritima")) 




res_alterniflora_maritima %>% pull(Family) %>% table() %>% sort(decreasing = TRUE) %>% head(5) %>% data.frame()
res_anglica_maritima %>% pull(Family) %>% table() %>% sort(decreasing = TRUE) %>% head(5) %>% data.frame()
res_anglica_alterniflora %>% pull(Family) %>% table() %>% sort(decreasing = TRUE) %>% head(5) %>% data.frame()

subset_taxa(phylo_rennes_prop, Family == "Desulfosarcinaceae") %>% 
  subset_samples(compartment == "root") %>%
  subset_taxa(!(is.na(Family))) %>% 
  tax_glom("Family") %>%
  plot_bar(fill="Family") + facet_wrap(~sample_Species, scales="free_x", ncol=3) + 
  theme(strip.text.x = element_text(size=20),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size=25),
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.ticks.x = element_blank()) +
  ylab("Relative Abundance")








