library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(patchwork)

read_delim("/Users/katieemelianova/Desktop/Spartina/Spartina2025/spartina_photosynq_results.txt") %>% colnames()

soil <- readxl::read_xlsx("/Users/katieemelianova/Desktop/Spartina/Spartina2025/spartina_soil_chemistry_results.xlsx")
soil <- read_delim("/Users/katieemelianova/Desktop/Spartina/Spartina2025/spartina_photosynq_results.txt") %>% 
  dplyr::select("Which sample number?", "Which species?", "Which locality?") %>%
  rename("sample" = `Which sample number?`, "species" = "Which species?", "locality" = "Which locality?") %>%
  mutate(species=case_when(species == "anglica" ~"S. anglicus",
                           species == "maritima" ~ "S. maritimus",
                           species == "alterniflora" ~ "S. alterniflorus")) %>%
  right_join(soil) %>%
  mutate(species=case_when(sample %in% c("C1", "C2", "C3") ~ "Control",
                           !(sample %in% c("C1", "C2", "C3")) ~ species)) %>% 
  #drop_na() %>%
  set_colnames(c("sample", "species", "locality", "labnum", "ph", "EC", "dryweight", "nitrogen", "carbon")) %>%
  filter(!(is.na(species)))

common_theme <- theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1),
                      legend.position="none",
                      axis.title.x = element_blank(),
                      axis.title = element_text(size=30),
                      axis.text = element_text(size=25))

nitrogen_plot <- soil %>% 
  ggplot(aes(x=species, y=nitrogen, fill=species)) + 
  geom_boxplot(width = 0.6) +
  common_theme +
  scale_fill_manual(values=c("grey89", "brown2", "palegreen3", "dodgerblue2")) +
  ylab("Total Nitrogen (g/kg-1)")
  
carbon_plot <- soil %>% 
  ggplot(aes(x=species, y=carbon, fill=species)) + 
  geom_boxplot(width = 0.6) +
  common_theme +
  scale_fill_manual(values=c("grey89", "brown2", "palegreen3", "dodgerblue2")) +
  ylab("Soil Organic Carbon (g/kg-1)")

cnratio_plot <- soil %>% 
  ggplot(aes(x=species, y=carbon/nitrogen, fill=species)) + 
  geom_boxplot(width = 0.6) +
  common_theme +
  scale_fill_manual(values=c("grey89", "brown2", "palegreen3", "dodgerblue2")) +
  ylab("Carbon/Nitrogen Ratio")

ph_plot <- soil %>% 
  ggplot(aes(x=species, y=ph, fill=species)) + 
  geom_boxplot(width = 0.6) +
  common_theme +
  scale_fill_manual(values=c("grey89", "brown2", "palegreen3", "dodgerblue2")) +
  ylab("pH")




#######################################
#              composite              #
#######################################


png("soil_chemistry_composite.png", width=1700, height=1000)
((nitrogen_plot | carbon_plot) / (cnratio_plot | ph_plot))  | plot_spacer() + (lef_boxplot + phi2_boxplot) +
  plot_layout(widths = c(1, 5))
dev.off()



