library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(car)
library(lme4)
library(lsmeans)


psynq <- read_delim("/Users/katieemelianova/Desktop/Spartina/France2025/photosynq_results.txt") %>%
  dplyr::select(LEF, NPQt, PhiNPQ, Phi2, "Light Intensity (PAR)", PhiNO, `Which sample number?`, `Which species?`, `Which locality?`, `which number measurement?`) %>%
  rename(sample_id = `Which sample number?`,
         species = `Which species?`,
         locality = `Which locality?`,
         PAR = `Light Intensity (PAR)`,
         measurement=`which number measurement?`)

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


# conclusion: parents have significantly different Phi2 but their hybriod is not different to either one
# plot could be the scatter plot on left and thern boxplot of values per species with significance
# phi2 is the amount of light going towards photosynthesis, PhiNPQ is the amount going towards non photochemica quenching i.e. wasted in stress
# alterniflora has less going to photosynthesis and maritima has more going to photosynthesis
# a significant difference betrween parents, and the hybrid has an intermediate efficiency of potosynthesis

mean_psynq %>%
  ggplot(aes(x=species, y=Phi2)) + 
  geom_boxplot()


####################################
#       LEF vs PAR             #
####################################


# PAR does not go above 140 or so for maritima but goes higher for the other two species
# this could point to the leaves being more crispy, it doesnt seem that this is a consistent problem but one which is specific to maritima
# light might have an issue getting through the leaves of this species
psynq %>% 
  filter(PAR > 100 & PAR < 180 & LEF > 30 & LEF < 180) %>%
  ggplot(aes(x=PAR, y=LEF, colour=species)) +
  geom_point(size=4) +
  geom_smooth(method="lm", aes(fill = species), alpha=0.2) + 
  theme(panel.background = element_blank(),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=12))
  



# I am using linear mixed effects model to account for the multiple measurementsd taken of the same plant
# I am removing outliers
# then I am using lstrends to get and compare the slopes of the dots of each species
# if I use random effects for measurement and it has a singular warning, can I then agree that there is no need to include it?
m <- lmer(LEF ~ PAR * species + (1|measurement), data = (psynq %>% filter(PAR > 100 & PAR < 180 & LEF > 30 & LEF < 180)))
mls <- lstrends(m, "species", var="PAR")
pairs(mls)

test

qqnorm(resid(m))
qqline(resid(m))

# no significant difference in slopes

plot(resid(m))
plot(m)







