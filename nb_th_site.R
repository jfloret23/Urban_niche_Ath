library(ggplot2)
library(readr)
library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(FSA)
library(multcompView)

tri.to.squ<-function(x)
{
  rn <- row.names(x)
  cn <- colnames(x)
  an <- unique(c(cn,rn))
  myval <-  x[!is.na(x)]
  mymat <-  matrix(1,nrow=length(an),ncol=length(an),dimnames=list(an,an))
  for(ext in 1:length(cn))
  {
    for(int in 1:length(rn))
    {
      if(is.na(x[row.names(x)==rn[int],colnames(x)==cn[ext]])) next
      mymat[row.names(mymat)==rn[int],colnames(mymat)==cn[ext]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
      mymat[row.names(mymat)==cn[ext],colnames(mymat)==rn[int]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
    }
    
  }
  return(mymat)
}

setwd("/home/justine/Documents/These/City_experiment")
getwd()

Nb_th_site <- read_csv("Nb_th_site.csv")

# Destroyed site
table(Nb_th_site$type) # tot number of site sowed
table(Nb_th_site$type, Nb_th_site$Destroyed) # site destroyed

# Site with thaliana : presence / absence 
site_with_th <- Nb_th_site %>% dplyr::filter(Presence == "Yes")
table(site_with_th$type)

# Removing the destroyed site
Nb_th_site <- Nb_th_site %>% dplyr::filter(Destroyed == "Exist")


# Preliminary graph 
#----------------------------
# Number of Thaliana per site with only Thaliana site
site_with_th <- Nb_th_site %>% dplyr::filter(Presence == "Yes")
ggplot(site_with_th, aes(x = type, y= Nb_Thaliana_germ, fill = type)) + geom_boxplot() + theme_minimal() + geom_jitter(shape=16, position=position_jitter(0.2), size = 2.5)

# Number of Thaliana with all site
# site without Thaliana, number is 0 
Nb_th_site2 <- Nb_th_site %>% replace(is.na(.), 0)
ggplot(Nb_th_site2, aes(x = type, y= Nb_Thaliana_germ, fill = type)) + geom_boxplot() + theme_minimal() + geom_jitter(shape=16, position=position_jitter(0.2), size = 2.5)



# Testing establishment per environment
# ==================================================

# Parametric test 
# ----------------------------
# Tukey & anova
establishment <- aov(lm(Nb_th_site2$Nb_Thaliana_germ ~ factor(Nb_th_site2$type)))
summary.lm(establishment)
estbl <- TukeyHSD(establishment)

# T test + no correction
pp <- pairwise.t.test(Nb_th_site2$Nb_Thaliana_germ, Nb_th_site2$type, p.adjust.method = "none")
mat <- tri.to.squ(pp$p.value)
myletters <- multcompLetters(mat,compare="<=",threshold=0.05,Letters=letters)
myletters_df <- data.frame(group=names(myletters$Letters),letter = myletters$Letters )

# T test + Bonferroni correction
pp <- pairwise.t.test(Nb_th_site2$Nb_Thaliana_germ, Nb_th_site2$type, p.adjust.method = "bonf")
mat <- tri.to.squ(pp$p.value)
myletters <- multcompLetters(mat,compare="<=",threshold=0.05,Letters=letters)
myletters_df <- data.frame(group=names(myletters$Letters),letter = myletters$Letters )

# T test + Holm correction
pp <- pairwise.t.test(Nb_th_site2$Nb_Thaliana_germ, Nb_th_site2$type, p.adjust.method = "holm")
mat <- tri.to.squ(pp$p.value)
myletters <- multcompLetters(mat,compare="<=",threshold=0.05,Letters=letters)
myletters_df <- data.frame(group=names(myletters$Letters),letter = myletters$Letters )


# Non parametric test
# ----------------------------
# Kruskal and Wallis test 
result <- kruskal.test(Nb_Thaliana_germ ~ factor(type), data = Nb_th_site2)
dunnTest(Nb_Thaliana_germ ~ factor(type), data = Nb_th_site2)

# Wilcoxon pairwise test + Bonferroni correction
pp <- pairwise.wilcox.test(Nb_th_site2$Nb_Thaliana_germ, Nb_th_site2$type, p.adjust.method = "bonf")
mat <- tri.to.squ(pp$p.value)
myletters <- multcompLetters(mat,compare="<=",threshold=0.01,Letters=letters)
myletters_df <- data.frame(group=names(myletters$Letters),letter = myletters$Letters )


# Presence/absence A. thaliana model 
# ==================================================

Nb_th_site2 <- Nb_th_site2 %>% mutate(Presence = replace(Presence, Presence == "Yes", 1))
Nb_th_site2 <- Nb_th_site2 %>% mutate(Presence = replace(Presence, Presence == "No", 0))
Nb_th_site2$Presence <- as.numeric(Nb_th_site2$Presence)

modnul <- glm(Presence ~ 1, data = Nb_th_site2 ,family = binomial)
mod1 <- glm(Presence ~ type, data = Nb_th_site2 ,family = binomial)
mod2 <- glm(Presence ~ type, data = Nb_th_site2 ,family = quasibinomial)
mod3 <- glm(Presence ~ type, data = Nb_th_site2 ,family = binomial(link = "logit"))


summary(mod1)
summary(mod2)
summary(mod3)
summary(modnul)

aov(mod1)
aov(mod2)
aov(mod3)
aov(modnul)

qqnorm(resid(mod1))
qqnorm(resid(mod2))
qqnorm(resid(mod3))
qqnorm(resid(modnul))


pp <- pairwise.wilcox.test(Nb_th_site2$Presence, Nb_th_site2$type, p.adjust.method = "bonf")
mat <- tri.to.squ(pp$p.value)
myletters <- multcompLetters(mat,compare="<=",threshold=0.05,Letters=letters)
myletters_df <- data.frame(group=names(myletters$Letters),letter = myletters$Letters )
