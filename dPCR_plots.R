
library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(janitor)





samplespecies <- read_delim("~/Desktop/Spartina/Spartina2025/Rennes_Sampling.tsv") %>% 
  dplyr::select(`Sample number`, "Species") %>%
  set_colnames(c("User sample ID", "species"))


test <- read.csv("~/Desktop/burn_after_reading/DOME-SD-SPARTINA16S-RUN1_analysis_31_03_2026_08_33_15_UTC+02_00.csv", skip = 2)


test$Conc...cp.µL...undiluted.sample. %<>% as.numeric()


test2 <- readxl::read_xlsx("~/Desktop/burn_after_reading/KatiedPCR.xlsx")
jmf <- readxl::read_xlsx("~/Desktop/burn_after_reading/JMF-2508-06-1.xlsx")
test2$`Sample.NTC.Control` <- test2$`JMF sample ID` %>% str_sub(-2,-1) %>% as.integer() %>% as.character()

concentrations <- inner_join(test, test2, by="Sample.NTC.Control") %>% 
  dplyr::select("User sample ID", `Conc...cp.µL...undiluted.sample.`) %>% 
  mutate(`User sample ID` = str_remove(`User sample ID`, " r"))




results <- inner_join(concentrations, samplespecies, by="User sample ID")


png("16S_copies_barplot.png", height=500, width=500)
ggplot(results, aes(x=species, y=(`Conc...cp.µL...undiluted.sample.`), fill=species)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("brown2", "palegreen3", "dodgerblue2")) +
  theme(axis.title = element_text(size=20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=20),
        strip.background =element_rect(fill="skyblue2"),
        strip.text = element_text(20),
        legend.title = element_blank(),
        legend.text = element_text(size=15)) +
  ylab("16S Copies per Microlitre")
dev.off()


