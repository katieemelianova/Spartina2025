library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyverse)
library(ggplot2)

soil <- readxl::read_xlsx("soil_analysis_ph_results.xlsx")
soil <- read_delim("/Users/katieemelianova/Desktop/Spartina/Spartina2025/spartina_photosynq_results.txt") %>% 
  dplyr::select("Which sample number?", "Which species?") %>%
  rename(`Which sample number?` = "sample", "Which species?" = "species") %>%
  mutate(species=case_when(species == "anglica" ~"Sporobolus anglicus",
                           species == "maritima" ~ "Sporobolus maritimus",
                           species == "alterniflora" ~ "Sporobolus alterniflorus")) %>%
  right_join(soil) %>%
  mutate(species=case_when(sample %in% c("C1", "C2", "C3") ~ "Control",
                           !(sample %in% c("C1", "C2", "C3")) ~ species)) %>% 
  drop_na()


soil %>% 
  ggplot(aes(x=species, y=`pH (1:10; w:v; 10 mM CaCl2)`, fill=species)) + 
  geom_boxplot(width = 0.6) +
  scale_fill_manual(values=c("grey89", "brown2", "palegreen3", "dodgerblue2"))
  


