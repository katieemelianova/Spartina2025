library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(phyloseq)
library(ggsignif)
library(microViz)
library(DESeq2)
library(png)
library(tidyverse)
library(glue)
library(ggtext)


soil <- readxl::read_xlsx("/Users/katieemelianova/Desktop/Spartina/Spartina2025/spartina_soil_chemistry_results.xlsx")
soil <- read_delim("/Users/katieemelianova/Desktop/Spartina/Spartina2025/spartina_photosynq_results.txt") %>% 
  dplyr::select("Which sample number?", "Which species?", "Which locality?") %>%
  dplyr::rename("sample" = `Which sample number?`, "species" = "Which species?", "locality" = "Which locality?") %>%
  mutate(species=case_when(species == "anglica" ~"S. anglicus",
                           species == "maritima" ~ "S. maritimus",
                           species == "alterniflora" ~ "S. alterniflorus")) %>%
  mutate(species = if_else(species %in% c("S. alterniflorus", "S. maritimus", "S. anglicus"), glue("<i>{species}</i>"), species)) %>%
  right_join(soil) %>%
  mutate(species=case_when(sample %in% c("C1", "C2", "C3") ~ "Control",
                           !(sample %in% c("C1", "C2", "C3")) ~ species)) %>% 
  set_colnames(c("sample", "species", "locality", "labnum", "ph", "EC", "dryweight", "nitrogen", "carbon")) %>%
  filter(!(is.na(species)))

common_theme <- theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1),
                      legend.position="none",
                      axis.title.x = element_blank(),
                      axis.title = element_text(size=30),
                      axis.text = element_text(size=25),
                      panel.background = element_blank(),
                      plot.margin = unit(c(1.5, 0, 0, 0),
                                         "inches")
                      #plot.margin = margin(1.5, 1.5, 1.5, 1.5, "cm")
                      )




nitrogen_plot <- soil %>% 
  mutate(species = fct_relevel(species, "Control", "<i>S. alterniflorus</i>", "<i>S. anglicus</i>", "<i>S. maritimus</i>")) %>%
  ggplot(aes(x=species, y=nitrogen, fill=species)) + 
  geom_boxplot(width = 0.6) +
  common_theme +
  scale_fill_manual(values=c("grey89", "brown2", "palegreen3", "dodgerblue2")) +
  ylab("Total Nitrogen (g/kg-1)") +
  geom_signif(
    test = "t.test",
    comparisons = list(c("S. anglicus", "S. alterniflorus"), c("S. alterniflorus", "S. maritimus"), c("S. anglicus", "S. maritimus")),
    map_signif_level = TRUE, textsize = 6) +
  theme(axis.text.x = ggtext::element_markdown())
  
carbon_plot <- soil %>% 
  mutate(species = fct_relevel(species, "Control", "<i>S. alterniflorus</i>", "<i>S. anglicus</i>", "<i>S. maritimus</i>")) %>%
  ggplot(aes(x=species, y=carbon, fill=species)) + 
  geom_boxplot(width = 0.6) +
  common_theme +
  scale_fill_manual(values=c("grey89", "brown2", "palegreen3", "dodgerblue2")) +
  ylab("Soil Organic Carbon (g/kg-1)") +
  geom_signif(
    test = "t.test",
    comparisons = list(c("S. anglicus", "S. alterniflorus"), c("S. alterniflorus", "S. maritimus"), c("S. anglicus", "S. maritimus")),
    map_signif_level = TRUE, textsize = 6) +
  theme(axis.text.x = ggtext::element_markdown())

cnratio_plot <- soil %>% 
  mutate(species = fct_relevel(species, "Control", "<i>S. alterniflorus</i>", "<i>S. anglicus</i>", "<i>S. maritimus</i>")) %>%
  ggplot(aes(x=species, y=carbon/nitrogen, fill=species)) + 
  geom_boxplot(width = 0.6) +
  common_theme +
  scale_fill_manual(values=c("grey89", "brown2", "palegreen3", "dodgerblue2")) +
  ylab("Carbon/Nitrogen Ratio") +
  geom_signif(
    test = "t.test",
    comparisons = list(c("S. anglicus", "S. alterniflorus"), c("S. alterniflorus", "S. maritimus"), c("S. anglicus", "S. maritimus")),
    map_signif_level = TRUE, textsize = 6) +
  theme(axis.text.x = ggtext::element_markdown())

ph_plot <- soil %>% 
  mutate(species = fct_relevel(species, "Control", "<i>S. alterniflorus</i>", "<i>S. anglicus</i>", "<i>S. maritimus</i>")) %>%
  ggplot(aes(x=species, y=ph, fill=species)) + 
  geom_boxplot(width = 0.6) +
  common_theme +
  scale_fill_manual(values=c("grey89", "brown2", "palegreen3", "dodgerblue2")) +
  ylab("pH") +
  geom_signif(
    test = "t.test",
    comparisons = list(c("S. anglicus", "S. alterniflorus"), c("S. alterniflorus", "S. maritimus"), c("S. anglicus", "S. maritimus")),
    map_signif_level = TRUE, textsize = 6) +
  theme(axis.text.x = ggtext::element_markdown())



############################################
#     composite with greenhouse photo      #  
############################################

#greenhouse2_cropped <- readPNG("greenhouse2_cropped.png")
test1 <- ggplot() + 
  theme(panel.background = element_blank(),
        plot.margin = margin(0,0,0,0, "cm")) + 
  cowplot::draw_image("greenhouse2.png")

test2 <- ggplot() + 
  theme(panel.background = element_blank(),
        plot.margin = margin(0,0,0,0, "cm")) + 
  cowplot::draw_image("Maritima_LocoalMendon.png")
  


png("Figure2_soil_chemistry_composite_greenhouse_pic.png", width=1600, height=1400)
chemistry <- (nitrogen_plot | carbon_plot | ph_plot)
photos <- (test1 | plot_spacer() | test2) + plot_layout(widths = c(5, -1.1, 5))
(chemistry / plot_spacer() /photos) + 
  plot_layout(heights = c(0.8, -0.01, 1)) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 35))
dev.off()






############################################
#              soil microbial              #
############################################



phylo_greenhouse <- readRDS("/Users/katieemelianova/Desktop/Spartina/JMF_results/Results_2026-03-03/JMF-2512-16_M7TD7_16S_rRNA_V4_CPR_samples_refined_gte_1000_phyloseq.rds") %>%
  tax_filter(min_prevalence = 10,
             prev_detection_threshold = 5,
             min_total_abundance = 5,
             min_sample_abundance = 5) %>%
  subset_taxa(!(Family %in% c("Mitochondria", "Chloroplast"))) %>% subset_taxa(!(Order %in% c("Mitochondria", "Chloroplast")))




phylo_greenhouse@sam_data$compartment <- ifelse(endsWith(phylo_greenhouse@sam_data$Sample_description, "R"), "Rhizosphere", "Bulk Soil")
phylo_greenhouse@sam_data$user_sample_id <- substr(phylo_greenhouse@sam_data$Sample_description, 1,4)



sample_mapping <- readxl::read_xlsx("/Users/katieemelianova/Desktop/Spartina/Spartina2025/rennes_sampling_spreadsheet.xlsx") %>%
  dplyr::select("Sample number", "Species", "Locality") %>%
  set_colnames(c("user_sample_id", "Species", "Locality"))
sample_mapping$user_sample_id %<>% as.integer() %>% as.character()




phylo_greenhouse@sam_data$Species <- left_join(data.frame(phylo_greenhouse@sam_data), sample_mapping, by = "user_sample_id") %>% pull(Species)
phylo_greenhouse@sam_data$Species <- ifelse(is.na(phylo_greenhouse@sam_data$Species), "Control", phylo_greenhouse@sam_data$Species)
phylo_greenhouse@sam_data$Locality <- left_join(data.frame(phylo_greenhouse@sam_data), sample_mapping, by = "user_sample_id") %>% pull(Locality)

phylo_greenhouse@sam_data$Species %<>% str_replace("Spartina anglica", "S. anglicus")
phylo_greenhouse@sam_data$Species %<>% str_replace("Spartina maritima", "S. maritimus")
phylo_greenhouse@sam_data$Species %<>% str_replace("Spartina alternifllora", "S. alterniflorus")



phylo_greenhouse_prop <- transform_sample_counts(phylo_greenhouse, function(otu) otu/sum(otu))
ord.nmds.bray_elevation <- ordinate(phylo_greenhouse_prop, method="NMDS", distance="bray")


phylo_greenhouse_prop@sam_data$Species <- if_else(phylo_greenhouse_prop@sam_data$Species %in% c("S. alterniflorus", "S. maritimus", "S. anglicus"), glue("<i>{phylo_greenhouse_prop@sam_data$Species}</i>"), phylo_greenhouse_prop@sam_data$Species)

phylo_greenhouse_prop@sam_data$Species <- fct_relevel(phylo_greenhouse_prop@sam_data$Species, "Control", "<i>S. alterniflorus</i>", "<i>S. anglicus</i>", "<i>S. maritimus</i>")

  
  
greenhouse_ordination <- plot_ordination(phylo_greenhouse_prop, ord.nmds.bray_elevation, shape="compartment", color="Species", title="Bray NMDS") + 
  geom_point(size = 12) +
  theme(strip.text.x = element_text(size=30),
        axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=25),
        axis.title = element_text(size=30),
        legend.text = element_text(size=25),
        legend.title = element_blank(),
        plot.margin = margin(1,5,1,5, "cm"),
        legend.position = c(0.99, 0.2),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(colour = "black"),
        panel.background = element_blank()) +
  ggtitle("") +
  scale_colour_manual(values = c("grey60", "brown2", "palegreen3", "dodgerblue2")) +
  theme(legend.text = ggtext::element_markdown())


greenhouse_richness <- plot_richness(phylo_greenhouse_prop, x="compartment", measures=c("Shannon"), color="Species") +
  geom_point(size = 3.5) +
  theme(strip.text.x = element_text(size=25),
        axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size=25),
        legend.text = element_text(size=20),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none") +
  scale_colour_manual(values = c("grey71", "brown2", "palegreen3", "dodgerblue2"))



#################################
#        greenhouse PCA         #
#################################


png("FigureS4_greenhouse_PCA.png", width=1000, height=800)
greenhouse_ordination
dev.off()






