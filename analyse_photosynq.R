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
            PAR=mean(`Light Intensity (PAR)`),
            PhiNO=mean(PhiNO)) %>%
  left_join((psynq %>% dplyr::select(sample_id, species, locality)), 
            by="sample_id") %>% 
  unique()


plot(mean_psynq$LEF, mean_psynq$PAR)



####################################
#       phinpq vs ph2.             #
####################################
# linear inverse realtionship
# test per species of the pearson correlation coefficient. it looks to be about the same for all three species
# once you have established a good linear relationship for each species, you can use phi2 to perform an ANOVA
# to perform anova though you need to check for normality etc and assumptions of test are satisfied

mean_psynq %>% 
  ggplot(aes(x=PhiNPQ, y=Phi2, colour=species)) +
  geom_point(size=4) +
  ylim(0.4, 0.9) +
  xlim(-0.1, 0.4)


###############################################################################################
#          test for significanct difference betwen Phi2 between alterniflora and maritima    #
###############################################################################################



shapiro.test(mean_psynq$Phi2)
# reflect and sqrt normalise phi2 to make it normal
# I used the shapiro test to test for normality and did a density plot to look for normal distribution
mean_psynq %<>% drop_na() %>% filter(Phi2 > 0) %>% mutate(reflect=sqrt(max(Phi2) - Phi2))

mean_psynq %>% filter(species == "maritima") %>% pull(reflect) %>% length()
mean_psynq %>% filter(species == "alterniflora") %>% pull(reflect) %>% length()
mean_psynq %>% filter(species == "anglica") %>% pull(reflect) %>% length()

# do a t test for difference between species - pairwise
t.test(x = mean_psynq %>% filter(species == "maritima") %>% pull(reflect), 
       y = mean_psynq %>% filter(species == "alterniflora") %>% pull(reflect), 
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





