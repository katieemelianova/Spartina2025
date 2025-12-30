library(dada2)
library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyverse)
library(phyloseq)


# I took this from the JMF code, it just alters the filenames
convert_to_filtered_file_name <- function(x, reads_dir) {
  x %>%
    basename() %>%
    str_replace("\\.fastq.gz$", ".filtered.fastq.gz") %>% 
    file.path(reads_dir, .)
}

################################################################################################
######             this bit doesnt take too long and is fine to run locally             ########
##### I need the resultig output to combine with taxonomy info to make final phyloseq object.  #
################################################################################################

run_jr_analysis <- function(){
  # provide fastq path, read in fastqs and set filtered filenames
  reads_directory <- ("/Users/katieemelianova/Desktop/Spartina/spartina_amplicon_code/jr_amplicons")
  #reads_directory <- ("/lisc/scratch/spartina/JR_amplicons")
  read_1_fastq_files <- list.files(reads_directory, pattern = "_1\\.fastq.gz$", full.names = TRUE)
  read_2_fastq_files <- list.files(reads_directory, pattern = "_2\\.fastq.gz$", full.names = TRUE)
  read_1_filtered_fastq_files <- read_1_fastq_files %>% convert_to_filtered_file_name(reads_directory)
  read_2_filtered_fastq_files <- read_2_fastq_files %>% convert_to_filtered_file_name(reads_directory)
  
  out <- filterAndTrim(read_1_fastq_files, read_1_filtered_fastq_files, read_2_fastq_files, read_2_filtered_fastq_files, truncLen=c(220,150),
                       maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                       compress=TRUE, multithread=FALSE) 
  
  # read them in again to get existing files only, I guess some files do not pass filtering and dont get created
  read_1_filtered_fastq_files <- list.files(reads_directory, pattern = "_1\\.filtered.fastq.gz$", full.names = TRUE)
  read_2_filtered_fastq_files <- list.files(reads_directory, pattern = "_2\\.filtered.fastq.gz$", full.names = TRUE)
  
  # dereplicate
  derep1s <- derepFastq(read_1_filtered_fastq_files, verbose=TRUE)
  derep2s <- derepFastq(read_2_filtered_fastq_files, verbose=TRUE)
  
  # error estimation
  err1 <- learnErrors(derep1s, multithread=TRUE)
  err2 <- learnErrors(derep2s, multithread=TRUE)
  plotErrors(err1, nominalQ=TRUE)
  
  # run dada inference
  dada1 <- dada(derep1s, err=err1, multithread=TRUE)
  dada2 <- dada(derep2s, err=err2, multithread=TRUE)
  
  # merge paired reads
  merged_jr <- mergePairs(
    dadaF = dada1,
    dadaR = dada2,
    derepF = derep1s,
    derepR = derep2s,
    justConcatenate = TRUE)
  
  # make an ASV table
  seqtab_jr <- makeSequenceTable(merged_jr)
  
  # remove chimeras
  seqtab.nochim_jr <- removeBimeraDenovo(seqtab_jr, method="consensus", multithread=TRUE, verbose=TRUE)
  
  # way to subset the seqtable
  #getSequences(seqtab.nochim)[1:5]
  
  return(seqtab.nochim_jr)
}

########################################################################################################################################
###### this bit takes ages and I ran it on lisc and moved the resulting taxonomy_species_america to local once it had run              #
##### therefore do not run this function while locally, rather run it once on lisc, move across to local and then read in and proceed. #
########################################################################################################################################
assign_taxonomy_jr <- function(seqtab) {
  taxonomy_levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  taxa_jr <- assignTaxonomy(
    seqs = seqtab,
    refFasta = "/lisc/scratch/spartina/spartina_amplicon_code/silva_nr99_v138.2_toGenus_trainset.fa.gz",
    taxLevels = taxonomy_levels,
    multithread=12
  )
  taxa_jr <- addSpecies(taxa_jr, "/lisc/scratch/spartina/spartina_amplicon_code/silva_v138.2_assignSpecies.fa.gz")
  write.table(taxa_jr, "/Users/katieemelianova/Desktop/Spartina/spartina_amplicon_code/taxonomy_species_america", quote=FALSE)
}



sequence_table <- run_jr_analysis()
# not to be run locally
#assign_taxonomy_jr(sequence_table)



taxonomy_species_america <- read_delim("/Users/katieemelianova/Desktop/Spartina/spartina_amplicon_code/taxonomy_species_america") %>% set_colnames(c("sequence", "Domain", "Phylum", "Class", "Order", "Family")) %>% column_to_rownames("sequence") %>% separate(Family, sep=" ", c("Family", "Genus"), extra="merge")
samplesheet <- read_delim("/Users/katieemelianova/Desktop/Spartina/spartina_amplicon_code/JR_metadata.csv")

metadata_jr <- samplesheet %>% dplyr::select(Run, HOST, isolation_source) %>% 
  mutate(depth=case_when(grepl("5-10", isolation_source) ~ "5-10",
                         grepl("0-5", isolation_source) ~ "0-5")) %>%
  mutate(height=case_when(grepl("Tall", isolation_source) ~ "Tall Phenotype",
                          grepl("Short", isolation_source) ~ "Short Phenotype",
                          grepl("Medium", isolation_source) ~ "Medium")) %>%
  mutate(compartment=case_when(grepl("Endosphere", isolation_source) ~ "Endosphere",
                               grepl("Rhizosphere", isolation_source) ~ "Rhizosphere",
                               grepl("sediment", isolation_source) ~ "sediment")) %>%
  mutate(sample_id=Run) %>%
  dplyr::select(-isolation_source) %>%
  column_to_rownames("Run")




OTU_JR <- otu_table(sequence_table, taxa_are_rows=FALSE)
TAX_JR <- tax_table(taxonomy_species_america) %>% set_rownames(rownames(taxonomy_species_america)) %>% set_colnames(c("Domain", "Phylum", "Class", "Order", "Family", "Genus"))
sample_names(OTU_JR) <- str_remove(sample_names(OTU_JR), "_1.filtered.fastq.gz")

ps_jr <- phyloseq(OTU_JR, 
                  sample_data(metadata_jr), 
                  TAX_JR)


#######################################################################################################################
#       after running all above code, I save RDS and this can e reloaded and can start analysis from there            #
#######################################################################################################################

saveRDS(ps_jr, "/Users/katieemelianova/Desktop/Spartina/Spartina2025/JR_amplicons.rds")

##############################################################
##############################################################
#   Begin analysis here unless you want to rerun something   #
##############################################################
##############################################################


phylo_rolando <- readRDS("/Users/katieemelianova/Desktop/Spartina/Spartina2025/JR_amplicons.rds")
phylo_rolando <- tax_filter(
  phylo_rolando,
  min_prevalence = 10,
  prev_detection_threshold = 5,
  min_total_abundance = 5,
  min_sample_abundance = 5)


phylo_rolando@sam_data$study <- "usa"
phylo_rennes@sam_data$study <- "eur"

phylo_rennes_prop <- transform_sample_counts(phylo_rennes, function(otu) otu/sum(otu))
phylo_rolando_prop <- transform_sample_counts(phylo_rolando, function(otu) otu/sum(otu))
phylo_both <- merge_phyloseq(phylo_rolando_prop, phylo_rennes_prop)
colnames(phylo_both@sam_data) <- colnames(phylo_both@sam_data) %>% str_replace("HOST", "Species")

phylo_rennes_prop %>% psmelt() %>% filter(Genus == "Candidatus Thiodiazotropha") %>% nrow()
phylo_rolando_prop %>% psmelt() %>% filter(Genus == "Candidatus Thiodiazotropha") %>% nrow()

phylo_rennes_prop %>% psmelt() %>% filter(Genus == "Sedimenticola") %>% nrow()
phylo_rolando_prop %>% psmelt() %>% filter(Genus == "Sedimenticola") %>% nrow()





ord.nmds.phylo_both <- ordinate(phylo_both_prop, method="NMDS", distance="bray")
phylo_rolando_prop@sam_data$

subset_taxa(phylo_both, Family == "Desulfosarcinaceae") %>% 
  subset_samples(compartment %in% c("root", "Endosphere")) %>%
  subset_taxa(!(is.na(Family))) %>% 
  tax_glom("Family") %>%
  plot_bar(fill="Family") + facet_wrap(~sample_Species, scales="free_x", ncol=3)

phylo_both_melted <- phylo_both %>% psmelt()

phylo_both_melted %>% group_by(study, Genus) %>%
  summarise(nFeatures = n()) %>% 
  group_by(Genus) %>% 
  mutate(sumfeats = sum(nFeatures)) %>%
  mutate(percent_total = (nFeatures/sumfeats) * 100) %>%
  data.frame

phylo_both_melted %>% filter(Genus == "Candidatus Thiodiazotropha") %>% group_by(study) %>% summarise(nFeatures = n())


phylo_both_melted %>%#
  filter(Genus == "Candidatus Thiodiazotropha") %>%
  group_by(study) %>%
  summarise(nFeatures = n()) 
  
phylo_both_melted %>% filter(study == "eur" & Genus == "Candidatus Thiodiazotropha") %>% nrow()
phylo_both_melted %>% filter(study == "usa" & Genus == "Candidatus Thiodiazotropha") %>% nrow()

phylo_both_melted %>% filter(study == "eur" & Genus == "Sedimenticola") %>% nrow()
phylo_both_melted %>% filter(study == "usa" & Genus == "Sedimenticola") %>% nrow()



phylo_both_melted %>% filter(Genus == "Candidatus Thiodiazotropha" & study == "eur")

# plot ordination plot
plot_ordination(phylo_both_prop, ord.nmds.phylo_both, color="study", title="Bray NMDS", shape="compartment") + 
  geom_point(size = 7) +
  theme(strip.text.x = element_text(size=25),
        axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size=25),
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        panel.background = element_blank()) +
  ggtitle("")

