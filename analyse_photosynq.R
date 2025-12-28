library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyverse)
library(ggplot2)



psynq <- read_delim("/Users/katieemelianova/Desktop/Spartina/France2025/photosynq_results.txt") %>%
  dplyr::select(LEF, NPQt, PhiNPQ, Phi2, "Light Intensity (PAR)", PhiNO, `Which sample number?`, `Which species?`, `Which locality?`) %>%
  rename(sample_id = `Which sample number?`,
         species = `Which species?`,
         locality = `Which locality?`,
         PAR = `Light Intensity (PAR)`)

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

mean_psynq %>% 
  ggplot(aes(x=PhiNPQ, y=Phi2, colour=species)) +
  geom_point(size=4) +
  ylim(0.4, 0.9) +
  xlim(-0.1, 0.4)

mean_psynq %>%
  summarise(correl = cor(Phi2, PhiNPQ))
# -0.93 correlation pearson for all species combined
# because there is a good inveerse relationship between Phi2 and PhiNPQ we can 
# use Phi2 to test for any significant differences between species



##################################################################################
#          test for significanct difference betwen Phi2 between species pairs    #
##################################################################################

# first use the shapiro test to see if Phi2 is significantly different from normal distribution, removing outliers
mean_psynq %>% drop_na() %>% filter(Phi2 > 0 & Phi2 < 1) %>% pull(Phi2) %>% shapiro.test
mean_psynq %>% drop_na() %>% filter(Phi2 > 0 & Phi2 < 1) %>% pull(Phi2) %>% densityPlot()

# yes it is, so need to transform
# reflect and sqrt normalise phi2 to make it normal
mean_psynq %>% drop_na() %>% filter(Phi2 > 0 & Phi2 < 1) %>% mutate(transformed_phi2=sqrt(max(Phi2) - Phi2)) %>% pull(transformed_phi2) %>% shapiro.test 
mean_psynq %>% drop_na() %>% filter(Phi2 > 0 & Phi2 < 1) %>% mutate(transformed_phi2=sqrt(max(Phi2) - Phi2)) %>% pull(transformed_phi2) %>% densityPlot()

# lock in the new transformed values
mean_psynq %<>% drop_na() %>% filter(Phi2 > 0 & Phi2 < 1) %>% mutate(transformed_phi2=sqrt(max(Phi2) - Phi2))

# do a t test for difference between species - pairwise

# maritima - alterniflora
t.test(x = mean_psynq %>% filter(species == "maritima") %>% pull(transformed_phi2), 
       y = mean_psynq %>% filter(species == "alterniflora") %>% pull(transformed_phi2), 
       alternative = c("two.sided"))

# maritima - anglica
t.test(x = mean_psynq %>% filter(species == "maritima") %>% pull(transformed_phi2), 
       y = mean_psynq %>% filter(species == "anglica") %>% pull(transformed_phi2), 
       alternative = c("two.sided"))


# alterniflora - anglica
t.test(x = mean_psynq %>% filter(species == "anglica") %>% pull(transformed_phi2), 
       y = mean_psynq %>% filter(species == "alterniflora") %>% pull(transformed_phi2), 
       alternative = c("two.sided"))


####################################
#       LEF vs PAR             #
####################################


mean_psynq %>% 
  ggplot(aes(x=PAR, y=LEF, colour=species)) +
  geom_point(size=4) +
  ylim(30, 60) +
  xlim(110, 180)


mean_psynq %>%
  group_by(species) %>%
  summarise(correl = cor(LEF, PAR))







mean_psynq %>% filter(species == "maritima") %>% pull(Phi2) %>% shapiro.test

mean_psynq %>%
  ggplot(aes(x=species, y=PhiNPQ)) + 
  geom_boxplot()





