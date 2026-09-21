library(microViz)
library(DESeq2)
library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyverse)
library(DESeq2)
library(phyloseq)
library(patchwork)


######################################################################
#          functional annotation info for taxa from Rolando paper    #
######################################################################


functional <- readxl::read_xlsx("/Users/katieemelianova/Desktop/Spartina/JR_functionally_annotated_genera.xlsx", skip=2) %>% 
  dplyr::select("S/Sulfate reducing genera", "Sulfur oxidizing genera", "Iron oxidizing genera", "Nitrifying genera") %>%
  set_colnames(c("sulfate_reducers", "sulfur_oxidisers", "iron_oxidisers", "nitrifiers"))

sulfate_reducers <- functional %>% dplyr::select(sulfate_reducers) %>% drop_na() %>% pull() %>% str_replace_all("_", " ")
sulfur_oxidisers <- functional %>% dplyr::select(sulfur_oxidisers) %>% drop_na() %>% pull() %>% str_replace_all("_", " ")
iron_oxidisers <- functional %>% dplyr::select(iron_oxidisers) %>% drop_na() %>% pull() %>% str_replace_all("_", " ")
nitrifiers <- functional %>% dplyr::select(nitrifiers) %>% drop_na() %>% pull() %>% str_replace_all("_", " ")

# removing arcbacter because only one of them is a sulfur oxidiser
sulfur_oxidisers <- sulfur_oxidisers[sulfur_oxidisers != "Arcobacter"]

#####################################################
#           Read in RDS and sample data             # 
#####################################################

phylo_rennes <- readRDS("/Users/katieemelianova/Desktop/Spartina/JMF_results/JMF-2508-06_16S_raw_phyloseq.rds")
sample_info <- read_tsv("/Users/katieemelianova/Desktop/Spartina/Spartina2025//Rennes_Sampling.tsv") %>% 
  dplyr::select(Locality, Species, `Sample number`)

##############################################
#    filter on prevalence, set metadata      # 
##############################################

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
  str_replace("r", "Root") %>%
  str_replace("z", "Rhizome") %>%
  str_replace("s", "Rhizosphere") %>% 
  replace_na('Unknown')

# set locality and species as metadata slots
phylo_rennes@sam_data$Locality <- phylo_rennes@sam_data %>% data.frame() %>% left_join(sample_info, by = c("User_sample_ID_number" = "Sample number")) %>% pull(Locality) %>% replace_na('Unknown')
phylo_rennes@sam_data$Species <- phylo_rennes@sam_data %>% data.frame() %>% left_join(sample_info, by = c("User_sample_ID_number" = "Sample number")) %>% pull(Species) %>% replace_na('Unknown')
phylo_rennes@sam_data$Species %<>% 
  str_replace("Spartina alternifllora", "Sporobolus alterniflorus") %>%
  str_replace("Spartina maritima", "Sporobolus maritimus") %>%
  str_replace("Spartina anglica", "Sporobolus anglicus")



######################################
#         plot ordination            # 
######################################

# transform sample to relative abundance and remove chloroplast and mitochonria and unknowns
phylo_rennes <- subset_samples(phylo_rennes, compartment != "Unknown")
phylo_rennes <- subset_samples(phylo_rennes, Species != "Unknown")
phylo_rennes <- subset_taxa(phylo_rennes, !(Family %in% c("Mitochondria", "Chloroplast"))) %>% subset_taxa(!(Order %in% c("Mitochondria", "Chloroplast")))
phylo_rennes <- prune_samples(sample_sums(phylo_rennes) >= 2000, phylo_rennes)
set.seed(1)
phylo_rennes <- rarefy_even_depth(phylo_rennes, sample.size = min(sample_sums(phylo_rennes)), replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
phylo_rennes_prop <- transform_sample_counts(phylo_rennes, function(otu) otu/sum(otu))


ord.nmds.bray_jr <- ordinate(phylo_rennes_prop, method="NMDS", distance="bray")

## plot ordination plot
ordination_plot <- plot_ordination(phylo_rennes_prop, ord.nmds.bray_jr, color="Species", title="Bray NMDS", shape="compartment") + 
  geom_point(size = 8) +
  theme(strip.text.x = element_text(size=30),
        axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=25),
        axis.title = element_text(size=30),
        legend.title = element_blank(),
        legend.text = element_text(size=32),
        panel.background = element_blank()) +
  ggtitle("") + 
  scale_colour_manual(values=c("brown2", "palegreen3", "dodgerblue2")) +
  scale_colour_discrete(guide = guide_legend(label.theme = element_text(angle = 0, face = "italic")))






#################################
#           shannon             #
#################################



alpha_df <- estimate_richness(phylo_rennes, measures = c("Shannon", "Observed"))

# test if data is normally distributed
alpha_df$Shannon %>% shapiro.test

# no it is not, so use a wilcoxon rank sum test to test for difference between groups


# first get metadata and merge in data
meta_df <- data.frame(sample_data(phylo_rennes))
alpha_merged <- cbind(alpha_df, meta_df) %>% dplyr::select(Shannon, compartment, Locality, Species)

# do the test
pairwise.wilcox.test(alpha_merged$Shannon, alpha_merged$Species, p.adjust.method = "bonf")

alpha_merged %>% group_by(compartment, Species) %>% summarise(meanalpha=mean(Shannon))


alpha_diversity <- 
  phylo_rennes %>%
  {
    sample_data(.)$Species_facet <-
      paste0("italic('", sample_data(.)$Species, "')")
    .
  } %>%
  plot_richness(x = "compartment", measures = c("Shannon")) + 
  geom_boxplot(aes(fill = Species)) +
  facet_wrap(~Species_facet, labeller = label_parsed) +
  theme(axis.text.x = element_text(size=20, angle = 45, hjust=1, vjust=1),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=30),
        axis.title = element_text(size=30),
        strip.text.x = element_text(size=20),
        legend.position = "none") + 
  xlab("Compartment")
  scale_fill_manual(values=c("brown2", "palegreen3", "dodgerblue2"))

png("FigureS1_alpha_diversity.png", height=600, width=800)
alpha_diversity
dev.off()


###########################################################
# compare dispersion of rhizosphere to root and rhizome   #
###########################################################


# dispersion 
anglicus <- phylo_rennes %>% subset_samples(Species %in% c("Sporobolus anglicus"))
anglicus_metadata <- as(sample_data(anglicus), "data.frame")
anglicus_betadisper <- distance(anglicus, method = "bray") %>% betadisper(anglicus_metadata$compartment) %>% TukeyHSD()

maritimus <- phylo_rennes %>% subset_samples(Species %in% c("Sporobolus maritimus"))
maritimus_metadata <- as(sample_data(maritimus), "data.frame")
maritimus_betadisper <- distance(maritimus, method = "bray") %>% betadisper(maritimus_metadata$compartment) %>% TukeyHSD()

alterniflorus <- phylo_rennes %>% subset_samples(Species %in% c("Sporobolus alterniflorus"))
alterniflorus_metadata <- as(sample_data(alterniflorus), "data.frame")
alterniflorus_betadisper <- distance(alterniflorus, method = "bray") %>% betadisper(alterniflorus_metadata$compartment) %>% TukeyHSD()

# PUT RESULTS INTO TABLE SUPPLEMENTARY

# usign this to get results to copy and paste into supplementary tables
rbind((anglicus_betadisper$group %>% data.frame() %>% mutate(species ="Sporobolus anglicus")),
      (maritimus_betadisper$group %>% data.frame() %>% mutate(species ="Sporobolus maritimus")),
      (alterniflorus_betadisper$group %>% data.frame() %>% mutate(species ="Sporobolus alterniflorus"))) %>%
  rownames_to_column("comparison") %>%
  writexl::write_xlsx("betadispersiuon.xlsx")
  

########################################################################################
# test for difference between species and compartment host associated tissue clusters   #
########################################################################################

# first run adonis to test fr differences between compartments, species, and an interaction of both

#phylo_rennes_plant <- phylo_rennes %>% subset_samples(compartment %in% c("Root", "Rhizome"))
metadata <- as(sample_data(phylo_rennes), "data.frame")
perm_design <- how(nperm = 999, blocks = metadata$User_sample_ID_number)
dist_matrix <- as.matrix(distance(phylo_rennes, method = "bray"))
#make sure order is same
metadata <- metadata[rownames(dist_matrix), ]

permanova_result <- adonis2(dist_matrix ~ Species * compartment, data = metadata, permutations = perm_design, by="terms")
print(permanova_result)

# all three significant. Now we can see, per compartment, whether species are significantly different:

compartments <- unique(metadata$compartment)

pairwise_species_by_compartment <- lapply(compartments, function(comp) {
  sub_meta_comp <- metadata[metadata$compartment == comp, ]
  species_list <- unique(sub_meta_comp$Species)
  pairs <- combn(species_list, 2, simplify = FALSE)
  
  res_list <- lapply(pairs, function(pair) {
    sub_meta <- sub_meta_comp[sub_meta_comp$Species %in% pair, ]
    sub_dist <- as.dist(dist_matrix[rownames(sub_meta), rownames(sub_meta)])
    res <- adonis2(sub_dist ~ Species, data = sub_meta, permutations = 999)
    data.frame(
      compartment = comp,
      comparison = paste(pair, collapse = " vs "),
      R2 = res$R2[1],
      F  = res$F[1],
      p  = res$`Pr(>F)`[1]
    )
  })
  do.call(rbind, res_list)
})


pairwise_species_by_compartment_df <- do.call(rbind, pairwise_species_by_compartment)
pairwise_species_by_compartment_df$p_adj <- p.adjust(pairwise_species_by_compartment_df$p, method = "BH")
pairwise_species_by_compartment_df %>%
  writexl::write_xlsx("pairwise_permanova.xlsx")

















#####################################################
#           differential abundance analysis         # 
#####################################################

# geometric means function
# doing it this way because see here:
# https://github.com/joey711/phyloseq/issues/445
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

run_deseq <- function(phylo_object, design_term){
  formula_parsed<-paste("~", design_term)
  phylo_deseq <- phylo_object %>% phyloseq_to_deseq2(as.formula(formula_parsed))
  phylo_deseq_geomeans <- apply(counts(phylo_deseq), 1, gm_mean)
  phylo_deseq = estimateSizeFactors(phylo_deseq, geoMeans = phylo_deseq_geomeans)
  phylo_deseq = DESeq(phylo_deseq, fitType="local")
  return(phylo_deseq)
}

annotate_deseq_results <- function(deseq_result, phylo_object){
  annotated <- deseq_result %>% 
    data.frame() %>% 
    arrange(padj) %>% 
    filter(padj < alpha) %>% 
    rownames_to_column(var="amplicon") %>% 
    left_join(tax_table(phylo_object) %>% data.frame() %>% rownames_to_column(var="amplicon"))
}

# convert to deseq object specifying variable to test DE on per compartment
# then run deseq, using geo means
# reason for geo means function: https://github.com/joey711/phyloseq/issues/445
phylo_rennes_deseq_root <- subset_samples(phylo_rennes, compartment == "Root") %>% run_deseq("0 + Species")
phylo_rennes_deseq_rhizome <- subset_samples(phylo_rennes, compartment == "Rhizome") %>% run_deseq("0 + Species")
phylo_rennes_deseq_soil <- subset_samples(phylo_rennes, compartment == "Rhizosphere") %>% run_deseq("0 + Species")

phylo_rennes_deseq_alt_root_rhizosphere <- subset_samples(phylo_rennes, Species == "Sporobolus alterniflorus" & compartment %in% c("Root", "Rhizosphere")) %>% run_deseq("0 + compartment")
phylo_rennes_deseq_ang_root_rhizosphere <- subset_samples(phylo_rennes, Species == "Sporobolus anglicus" & compartment %in% c("Root", "Rhizosphere")) %>% run_deseq("0 + compartment")
phylo_rennes_deseq_mar_root_rhizosphere <- subset_samples(phylo_rennes, Species == "Sporobolus maritimus" & compartment %in% c("Root", "Rhizosphere")) %>% run_deseq("0 + compartment")



# set pvalue threshold
alpha = 0.01

# this shows the comparisons which is handy when specifying them later - the format is always a bit funny
resultsNames(phylo_rennes_deseq_root)



resultsNames(phylo_rennes_deseq_ang_root_rhizosphere)


# now we can get the per tissue results for contrasts between species

### TODO! might need to filter taxa by their abundance in N samples per group
# like in RNAseq
# whta I would do is get the OTU table and then just subset OTU by abundance in each species/Locality grouping, one by one
# then I would get all these OTUs and subset the phyloseq object using those

####################################
#      alterniflora vs anglica     #
####################################

contrast <- list("SpeciesSporobolus.alterniflorus", "SpeciesSporobolus.anglicus")

root_alt_ang <- results(phylo_rennes_deseq_root, contrast) %>% annotate_deseq_results(phylo_rennes)
rhiz_alt_ang <- results(phylo_rennes_deseq_rhizome, contrast) %>% annotate_deseq_results(phylo_rennes)
soil_alt_ang <- results(phylo_rennes_deseq_soil, contrast) %>% annotate_deseq_results(phylo_rennes)


####################################
#      alterniflora vs maritima    #
####################################
contrast <- list("SpeciesSporobolus.alterniflorus", "SpeciesSporobolus.maritimus")

root_alt_mar <- results(phylo_rennes_deseq_root, contrast) %>% annotate_deseq_results(phylo_rennes)
rhiz_alt_mar <- results(phylo_rennes_deseq_rhizome, contrast) %>% annotate_deseq_results(phylo_rennes)
soil_alt_mar <- results(phylo_rennes_deseq_soil, contrast) %>% annotate_deseq_results(phylo_rennes)


####################################
#      anglica vs maritima         #
####################################
contrast <- list("SpeciesSporobolus.anglicus", "SpeciesSporobolus.maritimus")

root_ang_mar <- results(phylo_rennes_deseq_root, contrast) %>% annotate_deseq_results(phylo_rennes)
rhiz_ang_mar <- results(phylo_rennes_deseq_rhizome, contrast) %>% annotate_deseq_results(phylo_rennes)
soil_ang_mar <- results(phylo_rennes_deseq_soil, contrast) %>% annotate_deseq_results(phylo_rennes)




#############################################
#      alterniflora root vs rhizosphere     #
#############################################

alt_root_rhizosphere <- phylo_rennes_deseq_alt_root_rhizosphere %>% results() %>% annotate_deseq_results(phylo_rennes)

#############################################
#       anglica root vs rhizosphere         #
#############################################

ang_root_rhizosphere <- phylo_rennes_deseq_ang_root_rhizosphere %>% results() %>% annotate_deseq_results(phylo_rennes)

#############################################
#       maritima root vs rhizosphere        #
#############################################

mar_root_rhizosphere <- phylo_rennes_deseq_mar_root_rhizosphere %>% results() %>% annotate_deseq_results(phylo_rennes)



# alterniflora has sedimenticolaceae enriched in root relative to rhjizosphere in USA
# it is also in Europe
# and now we can confirm it is also in  maritima and anglica
# next we can take the rtolando list of sulfur oxidisers and reducers and nitrifiers 
# and ask how many of each is enriched in root relative to rhizosphere per species
# and make a graph

annotate_functional_genus <- function(results_table, species){
  freq_table <- results_table$Genus %>% data.frame() %>% set_colnames(c("Genus"))
  freq_table %<>% mutate(annot_func = case_when(Genus %in% sulfate_reducers ~ "sulfate reducers",
                                                Genus %in% sulfur_oxidisers ~ "sulfur oxidisers",
                                                Genus %in% iron_oxidisers ~ "iron oxidisers",
                                                Genus %in% nitrifiers ~ "nitrifiers")) %>%
    mutate(Species=species) %>%
    drop_na()
  return(freq_table)
}

# this bit I get each resuots table of DA amplicons and keep only those which are annotated as one of the four functions (sox, sred, feox, nit)
ang_root_rhizosphere_functions <- ang_root_rhizosphere %>% annotate_functional_genus("Sporobolus anglicus")
mar_root_rhizosphere_functions <- mar_root_rhizosphere %>% annotate_functional_genus("Sporobolus maritimus")
alt_root_rhizosphere_functions <- alt_root_rhizosphere %>% annotate_functional_genus("Sporobolus alterniflorus")
















functional_da_asv_plot <- rbind(ang_root_rhizosphere_functions,
      mar_root_rhizosphere_functions,
      alt_root_rhizosphere_functions) %>%
  ggplot(aes(x=annot_func, fill=Species)) + 
  geom_histogram(stat="count") +
  facet_wrap(~Species) +
  xlab("Functional Grouping") +
  ylab("Observed richness in root associated ASVs") + 
  scale_fill_manual(values=c("brown2", "palegreen3", "dodgerblue2")) + 
  theme(strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    axis.title = element_text(size=30),
    axis.text = element_text(size=25),
    legend.title = element_blank(),
    legend.text = element_text(size=31, face="italic"),
    legend.key.size = unit(0.45,"cm"))






###########################################################################
#  plot differentially abundant Genus between root and soil per species.  #
###########################################################################

global_theme <- theme(strip.text.x = element_text(size = 30),
                     axis.title = element_text(size=30),
                     axis.text.x = element_blank(),
                     axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.y = element_text(size=25),
                     strip.text = element_text(size=30),
                     legend.title = element_blank(),
                     legend.text = element_text(size=30, face="italic"),
                     legend.justification = "left",
                     legend.key.size = unit(0.85,"cm"))


scale_colours <- scale_fill_manual(
  values = c(
    "Sedimenticola"        = "#FFAD0AFF",
    "Sulfurimonas"         = "#1BB6AFFF",
    "Sulfurovum"           = "#D72000FF",
    "Candidatus Thiodiazotropha" = "#132157FF",  # use EXACT string from tax_table(phylo_rennes_prop)[, "Genus"]
    "Thiolapillus"         = "pink"
  ),
  breaks = c("Sedimenticola", "Sulfurimonas", "Sulfurovum",
             "Candidatus Thiodiazotropha", "Thiolapillus"),
  labels = c("Sedimenticola", "Sulfurimonas", "Sulfurovum",
             "C. Thiodiazotropha", "Thiolapillus")
)


# a quick note - the > or < sign matters here for which way the comparisonis shown (up in root vs rhizosphere or vice versa)


# alterniflora
alterniflora_root_associated <- prune_taxa(alt_root_rhizosphere %>% filter(log2FoldChange > 0) %>% pull(amplicon), phylo_rennes_prop) %>% # get the amplicons which are DA and prune to include only those
  subset_samples(compartment %in% c("Root") & Species == "Sporobolus alterniflorus") %>% # then subset the object to include only sampes which are in root and rhizosphere
  subset_taxa(Genus %in% (alt_root_rhizosphere_functions %>% filter(annot_func %in% c("sulfur oxidisers")) %>% pull(Genus))) %>% # then of those DA amplicons in root or rhizosphere samples, select only those where the Genus matches one in sulfur oxidiser list
  tax_glom("Genus") %>%
  {
    sample_data(.)$Species_facet <-
      paste0("italic('", sample_data(.)$Species, "')")
    .
  } %>%
  plot_bar(fill="Genus") + 
  facet_wrap(~Species_facet, scales="free_x", ncol=3, labeller = label_parsed) +
  global_theme  +
  ylim(0, 0.2) +
  scale_colours +
  ylab("")
  



# maritima
maritima_root_associated <- prune_taxa(mar_root_rhizosphere %>% filter(log2FoldChange > 0) %>% pull(amplicon), phylo_rennes_prop) %>%
  subset_samples(compartment %in% c("Root") & Species == "Sporobolus maritimus") %>%
  subset_taxa(Genus %in% (mar_root_rhizosphere_functions %>% filter(annot_func %in% c("sulfur oxidisers")) %>% pull(Genus))) %>%
  tax_glom("Genus") %>%
  {
    sample_data(.)$Species_facet <-
      paste0("italic('", sample_data(.)$Species, "')")
    .
  } %>%
  plot_bar(fill="Genus") + 
  facet_wrap(~Species_facet, scales="free_x", ncol=3, labeller = label_parsed) +
  global_theme  +
  ylim(0, 0.2) +
  scale_colours + 
  ylab("Root Relative Abundance")




# anglica
anglica_root_associated <- prune_taxa(ang_root_rhizosphere %>% filter(log2FoldChange > 0) %>% pull(amplicon), phylo_rennes_prop) %>%
  subset_samples(compartment %in% c("Root") & Species == "Sporobolus anglicus") %>%
  subset_taxa(Genus %in% (ang_root_rhizosphere_functions %>% filter(annot_func %in% c("sulfur oxidisers")) %>% pull(Genus))) %>%
  tax_glom("Genus") %>%
  {
    sample_data(.)$Species_facet <-
      paste0("italic('", sample_data(.)$Species, "')")
    .
  } %>%
  plot_bar(fill="Genus") + 
  facet_wrap(~Species_facet, scales="free_x", ncol=3, labeller = label_parsed) +
  global_theme  +
  scale_colours + 
  theme(axis.title.x = element_text(size=30)) + 
  ylim(0, 0.2) +
  xlab("Sample") +
  ylab("")


rel_abundance_plot <- (alterniflora_root_associated / maritima_root_associated / anglica_root_associated) +  plot_layout(heights = c(1, 1, 1))

png("Figure1_panel.png", width=1900, height=1400)
(rel_abundance_plot | (functional_da_asv_plot / ordination_plot)) +
  plot_layout(widths = c(1.1, 1)) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 35))
dev.off() 



####################################
#.     END OF MANUSCRIPT FIGURES.  #
####################################














###########################################################################
#  plot differentially abundant Genus between root and soil per species.  #
###########################################################################

# an example of a command to get a significantly lower abundance ASVs in anglica compared to maritima
root_ang_mar %>% filter(log2FoldChange < 0) %>% pull(amplicon)

# take this same command and plug it into an abundance plot by species
prune_taxa(soil_alt_mar %>% filter(log2FoldChange < 0) %>% pull(amplicon), phylo_rennes) %>%
  tax_glom("Family") %>%
  plot_bar(fill="sample_Species") + facet_wrap(~sample_Species, scales="free_x", ncol=3)

# a nice sanity check to see that ASVs differentially abundant in barplot format indeed look so

# now we have differential abundance for ASVs
# we can ask which families are most represented in ASVs which are most abundant

# questions

# 1. how many ASVs are differentially abundant between each species pair?
# Answer = anglica and maritima have far fewer differentially abundant taxa than the other two comparisons

#pdf("differentially_abundant_ASV.pdf", height=7, width=12)
png("differentially_abundant_ASV.png", height=500, width=700)
data.frame(number_genus_diff_abundant = c(root_ang_mar$Family %>% length(), soil_ang_mar$Family %>% length(), rhiz_ang_mar$Family %>% length(),
                                          root_alt_mar$Family %>% length(), soil_alt_mar$Family %>% length(), rhiz_alt_mar$Family %>% length(),
                                          root_alt_ang$Family %>% length(), soil_alt_ang$Family %>% length(), rhiz_alt_ang$Family %>% length()),
           compartment=rep(c("root", "rhizosphere", "rhizome"), 3),
           comparison=c(rep("anglica vs maritima", 3), rep("anterniflora vs maritima", 3), rep("anterniflora vs anglica", 3))) %>%
  ggplot(aes(y=comparison, x=number_genus_diff_abundant, fill=compartment)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c("skyblue2", "salmon1", "khaki3"), labels=c("Rhizome", "Rhizosphere", "Root")) + 
                      xlab("Differentially Abundant ASVs") + 
                      ylab("Comparison") + 
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=16))
dev.off()


# 2. which families show differential abundance in each tissue between hybrid and parents but not between parents
# not present in alterniflora - maritima but present in anglica-mariitma and anglica-alterniflora
# i.e. transgressive hybrid bacterial taxa
# Answer=
# root: none
# rhizome: Cyclobacteriaceae, Desulfocapsaceae, Kiritimatiellaceae
# soil: Sulfurovaceae

# root
setdiff(intersect(root_ang_mar$Family, root_alt_ang$Family), root_alt_mar$Family)
# rhizome 
setdiff(intersect(rhiz_ang_mar$Family, rhiz_alt_ang$Family), rhiz_alt_mar$Family)
#soil
setdiff(intersect(soil_ang_mar$Family, soil_alt_ang$Family), soil_alt_mar$Family)



# 3. which families are present in all species across all localities?

# Answer = Desulfobacteraceae, Desulfocapsaceae, Flavobacteriaceae, Kiritimatiellaceae, Sedimenticolaceae, Thioalkalispiraceae"
# melt the relative abundance phyloseq object to easily access data
phylo_rennes_prop_melt <- phylo_rennes_prop %>% tax_glom("Family") %>% filter_taxa(function(x) mean(x) > 0.01, TRUE) %>% psmelt() %>% dplyr::select(c("OTU", "Sample", "Abundance" , "User_sample_ID_number", "compartment", "Locality", "sample_Species", "Domain", "Phylum", "Class", "Order", "Family"))

get_present_in_all <- function(melted_phyloseq, comp, mean_abundance_min){
  to_return <- melted_phyloseq %>% filter(compartment == comp & Abundance > mean_abundance_min) %>% 
    group_by(Family, sample_Species)  %>%
    summarise(length(unique(Locality))) %>% # how many localities is each family found in per species?
    filter(`length(unique(Locality))` == 2) %>% #get families which are present in both localities per species
    ungroup() %>%
    group_by(Family) %>%
    summarise(length(sample_Species)) %>% # how many species are both-locality families present in?
    filter(`length(sample_Species)` == 3) %>% #get families in which all three species are present in both localities
    drop_na() %>%
    pull(Family)
  return(to_return)
}


# get families present in both locvalities for all three species
present_in_all_root <- get_present_in_all(phylo_rennes_prop_melt, "root", 0.01)
present_in_all_rhizome <- get_present_in_all(phylo_rennes_prop_melt, "rhizome", 0.01)
present_in_all_soil <- get_present_in_all(phylo_rennes_prop_melt, "soil", 0.01)


# boxplot
#min_abnundance0.1_families_root <- phylo_rennes_prop_melt %>% 
#  filter(Family %in% present_in_all_root & compartment == "root" & Abundance > 0.01) %>%
#  dplyr::select(Abundance, Family, sample_Species, Locality) %>%
#  ggplot(aes(x=Family, y=log(Abundance), fill=sample_Species)) + 
#  geom_boxplot(width = 0.6) +
#  scale_fill_manual(values=c("brown2", "palegreen3", "dodgerblue2")) +
#  facet_wrap(~Family, scales="free_x", ncol=5) +
#  theme(strip.text.x = element_text(size = 13),
#        axis.title = element_text(size=14),
#        axis.text.x = element_blank(),
#        axis.title.x = element_blank(),
#        axis.text.y = element_text(size=15),
#        strip.background =element_rect(fill="khaki3"),
#        strip.text = element_text(size=16),
#        legend.position="none") +
#  ylab("Log(Root Abundance)")

#min_abnundance0.1_families_rhizome <- phylo_rennes_prop_melt %>% 
#  filter(Family %in% present_in_all_rhizome & compartment == "rhizome" & Abundance > 0.01) %>%
#  dplyr::select(Abundance, Family, sample_Species, Locality) %>%
#  ggplot(aes(x=Family, y=log(Abundance), fill=sample_Species)) + 
#  geom_boxplot(width = 0.6) +
#  scale_fill_manual(values=c("brown2", "palegreen3", "dodgerblue2")) +
#  facet_wrap(~Family, scales="free_x", ncol=5) +
#  theme(strip.text.x = element_text(size = 13),
#        axis.title = element_text(size=14),
#        axis.text.x = element_blank(),
#        axis.title.x = element_blank(),
#        axis.text.y = element_text(size=15),
#        strip.background =element_rect(fill="skyblue2"),
#        strip.text = element_text(size=16)) +
#  ylab("Log(Rhizome Abundance)")

#min_abnundance0.1_families_rhizosphere <- phylo_rennes_prop_melt %>% 
#  filter(Family %in% present_in_all_soil & compartment == "soil" & Abundance > 0.01) %>%
#  dplyr::select(Abundance, Family, sample_Species, Locality) %>%
#  ggplot(aes(x=Family, y=log(Abundance), fill=sample_Species)) + 
#  geom_boxplot(width = 0.6) +
#  scale_fill_manual(values=c("brown2", "palegreen3", "dodgerblue2")) +
#  facet_wrap(~Family, scales="free_x", ncol=5) +
#  theme(strip.text.x = element_text(size = 13),
#        axis.title = element_text(size=14),
#        axis.text.x = element_blank(),
#        axis.title.x = element_blank(),
#        axis.text.y = element_text(size=15),
#        strip.background =element_rect(fill="salmon1"),
#        legend.position="none",
#        strip.text = element_text(size=16)) +
#  ylab("Log(Rhizosphere Abundance)")
#
#pdf("min_abnundance0.1_families_compartments.pdf", height=11, width=12)
#png("min_abnundance0.1_families_compartments.png", height=100, width=100)
#(min_abnundance0.1_families_root / min_abnundance0.1_families_rhizome / min_abnundance0.1_families_rhizosphere) +  plot_layout(heights = c(2, 2, 1))
#dev.off()




# 4. which families arefound only in one species?

# having had a play around with this I'm not sure its wise to proceed as below
# present in all species at appreciable amounts is fine, but asking which species are found in one species but not in others requires filering for abundance and this seems to be getting into problems of variable dropout rates, taxa being reported because they are just one side of a threshold cutoff
# also it seems as though solutions I have played about with like sd of the means of species per family have strayed into differential abundance territory which I have already done
# so I think either I can leave this or consult with someione frim DOME


get_one_species_only <- function(melted_phyloseq, comp, mean_abundance_min){
  to_return <- melted_phyloseq %>% filter(compartment == comp) %>% 
    group_by(Family, sample_Species, Locality) %>%
    summarise(mean_abundance=mean(Abundance)) %>% # get mean abundance of each fanily oper species per locality
    filter(mean_abundance > mean_abundance_min) %>%
    summarise(length(unique(Locality))) %>% # how many localities is each family found in per species?
    filter(`length(unique(Locality))` == 2) %>% #get families which are present in both localities per species
    ungroup() %>%
    group_by(Family) %>%
    summarise(length(sample_Species)) %>% # how many species are both-locality families present in?
    filter(`length(sample_Species)` == 1) %>%
    pull(Family)
  return(to_return)
}

one_species_only_root <- get_one_species_only(phylo_rennes_prop_melt, "root", 0.01)
one_species_only_rhizome <- get_one_species_only(phylo_rennes_prop_melt, "rhizome", 0.01)
one_species_only_soil <- get_one_species_only(phylo_rennes_prop_melt, "soil", 0.01)


# barplot
subset_taxa(phylo_rennes_prop_family_abfilt, Family %in% rhiz_alt_ang$Family) %>% 
  subset_samples(compartment == "root") %>%
  subset_taxa(!(is.na(Family))) %>% 
  tax_glom("Family") %>%
  plot_bar(fill="Family") + facet_wrap(~sample_Species+Locality, scales="free_x", ncol=2) + 
  theme(strip.text.x = element_text(size=20),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size=25),
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.ticks.x = element_blank()) +
  ylab("Relative Abundance")





#############################################
#     USA rlando analysis comoarison
#############################################


phylo_rolando <- readRDS("/Users/katieemelianova/Desktop/Spartina/Spartina2025/JR_amplicons.rds")
phylo_rolando <- tax_filter(
  phylo_rolando,
  min_prevalence = 10,
  prev_detection_threshold = 5,
  min_total_abundance = 5,
  min_sample_abundance = 5)


phylo_rolando_deseq <- phylo_rolando %>% run_deseq("compartment")
test <- results(phylo_rolando_deseq) %>% annotate_deseq_results_jr(phylo_rolando)




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






#https://github.com/joey711/phyloseq/issues/763 
# make this into a funxction and do something like this for the oter species
alt_root_sampleids <- subset_samples(phylo_rennes, compartment == "root" & Species == "Spartina alternifllora") %>% sample_data() %>% rownames()
alt_root_samples <- phylo_rennes@otu_table[,alt_root_sampleids]
to_keep_otus <- alt_root_samples[rowSums(alt_root_samples > 5) >= 3, ] %>% rownames()
my_subset <- subset(otu_table(phylo_rennes), rownames(otu_table(phylo_rennes)) %in% to_keep_otus)
new_physeq <- merge_phyloseq(my_subset, tax_table(phylo_rennes), sample_data(phylo_rennes))



