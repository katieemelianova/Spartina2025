
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


results <- inner_join(concentrations, sample_mapping, by="User sample ID")


png("FigureS2_16S_copies_barplot.png", height=500, width=500)
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

png("FigureS2_16S_copies_barplot.png", height=500, width=500)
ggplot(results, aes(x=species, y=(`Conc...cp.µL...undiluted.sample.`), fill=species)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("brown2", "palegreen3", "dodgerblue2")) +
  theme(axis.title = element_text(size=20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=20),
        strip.background =element_rect(fill="skyblue2"),
        strip.text = element_text(20))
        #legend.title = element_blank(),
        #legend.text = element_text(size=15)) +
  #ylab("16S Copies per Microlitre")
dev.off()

