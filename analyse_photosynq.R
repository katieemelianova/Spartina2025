library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(car)
library(lme4)
library(lsmeans)

read_delim("/Users/katieemelianova/Desktop/Spartina/Spartina2025/spartina_photosynq_results.txt") %>% colnames()

psynq <- read_delim("/Users/katieemelianova/Desktop/Spartina/Spartina2025/spartina_photosynq_results.txt") %>%
  dplyr::select(LEF, NPQt, PhiNPQ, Phi2, "Light Intensity (PAR)", PhiNO, `Which sample number?`, `Which species?`, `Which locality?`, `which number measurement?`) %>%
  dplyr::rename(sample_id = `Which sample number?`,
         species = `Which species?`,
         locality = `Which locality?`,
         PAR = `Light Intensity (PAR)`,
         measurement=`which number measurement?`)

psynq$species %<>% 
  str_replace("anglica", "Sporobolus anglicus") %>%
  str_replace("alterniflora", "Sporobolus alterniflorus") %>%
  str_replace("maritima", "Sporobolus maritimus")
  
  

# take the average of repeat measurements to get a single value per plant per metric
mean_psynq <- psynq %>% 
  group_by(sample_id) %>%
  summarize(LEF=mean(LEF),
            NPQt=mean(NPQt),
            PhiNPQ=mean(PhiNPQ),
            Phi2=mean(Phi2),
            PAR=mean(PAR),
            PhiNO=mean(PhiNO)) %>%
  left_join((psynq %>% dplyr::select(sample_id, species, locality)), 
            by="sample_id") %>% 
  unique()



####################################
#       phinpq vs ph2.             #
####################################
# linear inverse realtionship
# test per species of the pearson correlation coefficient. it looks to be about the same for all three species
# once you have established a good linear relationship for each species, you can use phi2 to perform a t-test for each species pair
# need to first transform dta to meet normality assumptions

png("phi2_phinpq.png", height=500, width=500)
mean_psynq %>% 
  ggplot(aes(x=PhiNPQ, y=Phi2, colour=species)) +
  geom_point(size=4) +
  ylim(0.4, 0.9) +
  xlim(-0.1, 0.4) +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=17))
dev.off()

mean_psynq %>%
  drop_na() %>%
  summarise(correl = cor(Phi2, PhiNPQ))
# -0.93 correlation pearson for all species combined
# because there is a good inveerse relationship between Phi2 and PhiNPQ we can 
# use Phi2 to test for any significant differences between species



##################################################################################
#          test for significanct difference betwen Phi2 between species pairs    #
##################################################################################

# first use the shapiro test to see if Phi2 is significantly different from normal distribution, removing outliers
mean_psynq %>% drop_na() %>% filter(Phi2 > 0 & Phi2 < 1) %>% pull(Phi2) %>% shapiro.test
#mean_psynq %>% drop_na() %>% filter(Phi2 > 0 & Phi2 < 1) %>% pull(Phi2) %>% densityPlot()

# yes it is sig diff from normal distribution, so need to transform
# reflect and sqrt normalise phi2 to make it normal
mean_psynq %>% drop_na() %>% filter(Phi2 > 0 & Phi2 < 1) %>% mutate(transformed_phi2=sqrt(max(Phi2) - Phi2)) %>% pull(transformed_phi2) %>% shapiro.test 
#mean_psynq %>% drop_na() %>% filter(Phi2 > 0 & Phi2 < 1) %>% mutate(transformed_phi2=sqrt(max(Phi2) - Phi2)) %>% pull(transformed_phi2) %>% densityPlot()

# lock in the new transformed values
mean_psynq %<>% drop_na() %>% filter(Phi2 > 0 & Phi2 < 1) %>% mutate(transformed_phi2=sqrt(max(Phi2) - Phi2))

# do a t test for difference between species - pairwise

# maritima - alterniflora
t.test(x = mean_psynq %>% filter(species == "Sporobolus maritimus") %>% pull(transformed_phi2), 
       y = mean_psynq %>% filter(species == "Sporobolus alterniflorus") %>% pull(transformed_phi2), 
       alternative = c("two.sided"))

# maritima - anglica
t.test(x = mean_psynq %>% filter(species == "Sporobolus maritimus") %>% pull(transformed_phi2), 
       y = mean_psynq %>% filter(species == "Sporobolus anglicus") %>% pull(transformed_phi2), 
       alternative = c("two.sided"))


# alterniflora - anglica
t.test(x = mean_psynq %>% filter(species == "Sporobolus anglicus") %>% pull(transformed_phi2), 
       y = mean_psynq %>% filter(species == "Sporobolus alterniflorus") %>% pull(transformed_phi2), 
       alternative = c("two.sided"))


# conclusion: parents have significantly different Phi2 but their hybriod is not different to either one
# plot could be the scatter plot on left and thern boxplot of values per species with significance
# phi2 is the amount of light going towards photosynthesis, PhiNPQ is the amount going towards non photochemica quenching i.e. wasted in stress
# alterniflora has less going to photosynthesis and maritima has more going to photosynthesis
# a significant difference betrween parents, and the hybrid has an intermediate efficiency of potosynthesis

#png("phi2_per_species_boxplot.png", height=1500, width=1000)
phi2_boxplot <- mean_psynq %>%
  mutate(species = case_when(species == "Sporobolus alterniflorus" ~ "S. alterniflorus",
                             species == "Sporobolus anglicus" ~"S. anglicus",
                             species == "Sporobolus maritimus" ~ "S. maritimus")) %>%
  ggplot(aes(x=species, y=Phi2, fill=species)) + 
  geom_boxplot(lwd=0.8) +
  scale_fill_manual(values=c("brown2", "palegreen3", "dodgerblue2")) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1),
        axis.title = element_text(size=30),
        axis.text = element_text(size=25),
        legend.position="none",
        axis.title.x = element_blank())
#dev.off()

##########################
#       LEF              #
##########################

# I had thought that PAR was variable in this protocol but it turns out that its just the ambient light level
# thats why there is such low variation between the values
# i had analysed LEF vs PAR in the past but now that there doesnt seem to be much variation (measured or expected) in PAR its not worth us including it
# instead we can just repeat the analyses performed for the Phi2 - compare LEF between species


# first use the shapiro test to see if LEF is significantly different from normal distribution, removing outliers
mean_psynq %>% drop_na() %>% pull(LEF) %>% shapiro.test

# different so reflect transform, check again
mean_psynq %>% drop_na() %>% mutate(transformed_LEF=sqrt(max(LEF) - LEF)) %>% pull(transformed_LEF) %>% shapiro.test 

# now shapiro test satisfied, lock in values
mean_psynq %<>% drop_na() %>% mutate(transformed_LEF=sqrt(max(LEF) - LEF))

# maritima - alterniflora
t.test(x = mean_psynq %>% filter(species == "Sporobolus maritimus") %>% pull(transformed_LEF), 
       y = mean_psynq %>% filter(species == "Sporobolus alterniflorus") %>% pull(transformed_LEF), 
       alternative = c("two.sided"))

# maritima - anglica
t.test(x = mean_psynq %>% filter(species == "Sporobolus maritimus") %>% pull(transformed_LEF), 
       y = mean_psynq %>% filter(species == "Sporobolus anglicus") %>% pull(transformed_LEF), 
       alternative = c("two.sided"))


# alterniflora - anglica
t.test(x = mean_psynq %>% filter(species == "Sporobolus anglicus") %>% pull(transformed_LEF), 
       y = mean_psynq %>% filter(species == "Sporobolus alterniflorus") %>% pull(transformed_LEF), 
       alternative = c("two.sided"))


lef_boxplot <- mean_psynq %>%
  mutate(species = case_when(species == "Sporobolus alterniflorus" ~ "S. alterniflorus",
                             species == "Sporobolus anglicus" ~"S. anglicus",
                             species == "Sporobolus maritimus" ~ "S. maritimus")) %>%
  ggplot(aes(x=species, y=LEF, fill=species)) + 
  geom_boxplot(lwd=0.8) +
  scale_fill_manual(values=c("brown2", "palegreen3", "dodgerblue2")) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1),
        axis.title = element_text(size=30),
        axis.text = element_text(size=25),
        legend.position="none",
        axis.title.x = element_blank())


#png("lef_ph2_boxplot.png", height=1000, width=800)
(lef_boxplot + phi2_boxplot)
#dev.off()
