library(microViz)



phylo_rennes <- readRDS("/Users/katieemelianova/Desktop/Spartina/JMF_results/JMF-2508-06_16S_raw_phyloseq.rds")
sample_info <- read_tsv("/Users/katieemelianova/Desktop/Spartina/Spartina2025//Rennes_Sampling.tsv") %>% 
  dplyr::select(Locality, Species, `Sample number`)


phylo_rennes@sam_data$User_sample_ID_number <- phylo_rennes@sam_data %>% data.frame() %>% pull(User_sample_ID) %>% str_split_i(" ", 1)
phylo_rennes@sam_data$compartment <- phylo_rennes@sam_data %>% data.frame() %>% pull(User_sample_ID) %>% str_split_i(" ", 2) %>% 
  str_replace("root", "r") %>%
  str_replace("rhizome", "z") %>%
  str_replace("soil", "s") %>%
  str_replace("r", "root") %>%
  str_replace("z", "rhizome") %>%
  str_replace("s", "soil") %>% 
  replace_na('Unknown')

phylo_rennes@sam_data$Locality <- phylo_rennes@sam_data %>% data.frame() %>% left_join(sample_info, by = c("User_sample_ID_number" = "Sample number")) %>% pull(Locality) %>% replace_na('Unknown')
phylo_rennes@sam_data$Species <- phylo_rennes@sam_data %>% data.frame() %>% left_join(sample_info, by = c("User_sample_ID_number" = "Sample number")) %>% pull(Species) %>% replace_na('Unknown')

phylo_rennes = subset_samples(phylo_rennes, Species != "Unknown")



plot_richness(phylo_rennes, x="compartment", measures=c("Shannon", "Simpson"), color="Species") +
  geom_point(size = 3.5) +
  theme(strip.text.x = element_text(size=25),
        axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size=25),
        legend.text = element_text(size=20),
        axis.title.x = element_blank(),
        legend.title = element_blank()) 


phylo_rennes_prop <- transform_sample_counts(phylo_rennes, function(otu) otu/sum(otu))
ord.nmds.bray_jr <- ordinate(phylo_rennes_prop, method="NMDS", distance="bray")

plot_ordination(phylo_rennes_prop, ord.nmds.bray_jr, color="Species", title="Bray NMDS", shape="compartment") + 
  geom_point(size = 7) +
  theme(strip.text.x = element_text(size=25),
        axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size=25),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20)) +
  ggtitle("")












# geometric means function
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# convert to deseq object specifying the 
phylo_rennes_deseq <- phyloseq_to_deseq2(phylo_rennes, ~Species)

# row-wise geometric mean of ASV counts
phylo_rennes_deseq_geomeans <- apply(counts(phylo_rennes_deseq), 1, gm_mean)

# estimate size factors and run DESeq
phylo_rennes_deseq = estimateSizeFactors(phylo_rennes_deseq, geoMeans = phylo_rennes_deseq_geomeans)
phylo_rennes_deseq = DESeq(phylo_rennes_deseq, fitType="local")




res = res[order(res$padj, na.last=NA), ]
alpha = 0.001

# get significant results sorted by adjusted pvalue and
res_alterniflora_maritima = results(phylo_rennes_deseq, contrast=c("Species", "Spartina alternifllora","Spartina maritima"))
res_anglica_maritima = results(phylo_rennes_deseq, contrast=c("Species", "Spartina anglica","Spartina maritima"))
res_anglica_alterniflora = results(phylo_rennes_deseq, contrast=c("Species", "Spartina alternifllora","Spartina anglica"))

res_alterniflora_maritima = res_alterniflora_maritima[order(res_alterniflora_maritima$padj, na.last=NA), ]
res_anglica_maritima = res_anglica_maritima[order(res_anglica_maritima$padj, na.last=NA), ]
res_anglica_alterniflora = res_anglica_alterniflora[order(res_anglica_alterniflora$padj, na.last=NA), ]

res_alterniflora_maritima = res_alterniflora_maritima[(res_alterniflora_maritima$padj < alpha), ]
res_anglica_maritima = res_anglica_maritima[(res_anglica_maritima$padj < alpha), ]
res_anglica_alterniflora = res_anglica_alterniflora[(res_anglica_alterniflora$padj < alpha), ]

res_alterniflora_maritima = cbind(data.frame(res_alterniflora_maritima), as.matrix(tax_table(phylo_rennes)[rownames(res_alterniflora_maritima), ]))
res_anglica_maritima = cbind(data.frame(res_anglica_maritima), as.matrix(tax_table(phylo_rennes)[rownames(res_anglica_maritima), ]))
res_anglica_alterniflora = cbind(data.frame(res_anglica_alterniflora), as.matrix(tax_table(phylo_rennes)[rownames(res_anglica_alterniflora), ]))


res_alterniflora_maritima %>% filter(Family != "Mitochondria") %>% pull(Family) %>% table() %>% sort(decreasing = TRUE) %>% head(10) %>% data.frame()
res_anglica_maritima %>% filter(Family != "Mitochondria") %>% pull(Family) %>% table() %>% sort(decreasing = TRUE) %>% head(10) %>% data.frame()
res_anglica_alterniflora %>% filter(Family != "Mitochondria") %>% pull(Family) %>% table() %>% sort(decreasing = TRUE) %>% head(10) %>% data.frame()


phylo_rennes_prop %<>% subset_taxa(!(Family %in% c("Incertae Sedis", "Mitochondria")))

phylo_rennes_prop %>%
  tax_fix() %>%
  tax_transform("clr", rank = "Family") %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc() %>%
  ord_plot(color = "Species", shape = "compartment", size = 2) +
  scale_colour_brewer(palette = "Dark2")

phylo_rennes_prop %>%
  tax_fix() %>%
  #tax_transform("clr", rank = "Genus") %>%
  dist_calc(dist = "bray") %>%
  ord_calc(method="NMDS") %>%
  #ord_plot_iris(tax_level = "Genus", ord_plot = "above", anno_colour = "Species")
  ord_plot(color = "Species", shape = "compartment", size = 2)

