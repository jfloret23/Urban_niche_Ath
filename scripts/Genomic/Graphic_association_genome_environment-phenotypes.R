#' ---
#'title: "QTL regions"
#'author: "Justine Floret"
#'date: "2026-02-24"
#'output: 
#'  html_document: 
#'    toc: yes
#'    toc_float: yes
#'    theme: cerulean
#'    highlight: pygments
#'    df_print: paged
#' ---

#+ include=FALSE
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE
)

#+ message=FALSE, warning=FALSE
library(readr)
library(dplyr)
library(tidyr)
library(data.table)
library(gtable)
library(ggplot2)
library(ggrepel)
library(knitr)
library(stringr) #remove space in string
library(viridis)



#' # I. Preparing the data

#' QTL positions
#+ message=FALSE, warning=FALSE
dormancy <- read_csv("/home/justine/sciebo/Genomic/assays/QTL_position/Dormancy_position_quantile.csv")
dormancy <- dormancy %>% mutate(across(c(2:4), as.numeric)) 
dormancy <- dormancy %>% mutate(across(c(1,9), as.factor)) 
dormancy <- dormancy[,-c(5:8)]
dormancy$Factor <- "Dormancy"
dormancy <- dormancy %>% relocate(Factor, .after = Pos)

#+ message=FALSE, warning=FALSE
exsituFT <- read_csv("/home/justine/sciebo/Genomic/assays/QTL_position/meanExsituFT_quantile_associations.csv")
exsituFT <- exsituFT %>% mutate(across(c(2:4), as.numeric)) 
exsituFT <- exsituFT %>% mutate(across(c(1,9), as.factor)) 
exsituFT <- exsituFT[,-c(5:8)]
exsituFT$Factor <- "exsituFT"
exsituFT <- exsituFT %>% relocate(Factor, .after = Pos)

#+ message=FALSE, warning=FALSE
popsize <- read_csv("/home/justine/sciebo/Genomic/assays/QTL_position/Popsize23-22_position_quantile.csv")
popsize <- popsize %>% filter(Census_year == "23")
popsize <- popsize %>% mutate(across(c(2:5), as.numeric)) 
popsize <- popsize %>% mutate(across(c(1,6), as.factor)) 
popsize <- popsize[,-c(3,7)]
popsize$Factor <- "Popsize"
popsize <- popsize %>% relocate(Factor, .after = Pos)

#+ message=FALSE, warning=FALSE
envIndexes <- read_csv("/home/justine/sciebo/Genomic/assays/QTL_position/EnvIndexes_sum_effects_quantile.csv")
envIndexes <- envIndexes %>% mutate(across(c(5,6), as.numeric)) 
envIndexes <- envIndexes %>% mutate(across(c(3,4), as.factor)) 
envIndexes <- envIndexes %>% relocate(pvalue, .after = Factors)
envIndexes <- envIndexes %>% relocate(R, .after = pvalue)
colnames(envIndexes) <- c("Chrm", "Pos", "Factor", "pvalue", "R", "Effect" )

#+ message=FALSE, warning=FALSE
envFactors_subdata <- read_csv("/home/justine/sciebo/Genomic/assays/QTL_position/EnvFactors_effects_quantile.csv")
envFactors_subdata <- envFactors_subdata %>% mutate(across(c(5,6), as.numeric)) 
envFactors_subdata <- envFactors_subdata %>% mutate(across(c(3,4), as.factor)) 
envFactors_subdata <- envFactors_subdata %>% relocate(pvalue, .after = Factors)
envFactors_subdata <- envFactors_subdata %>% relocate(R, .after = pvalue)
colnames(envFactors_subdata) <- c("Chrm", "Pos", "Factor", "pvalue", "R", "Effect" )
envFactors_subdata %>% distinct(Factor)%>% print(n=25)

#' Factors to keep
#+ message=FALSE, warning=FALSE
envFactors_subdata <- envFactors_subdata %>% filter(Factor %in% c("C.N.ratio", "DD", "NO3", "PO4", "SO4", "Cl","F",  "Tmax_04", "Tmean_04", "Tmin_04"))


#' Merge the dataset QTL positions
#+ message=FALSE, warning=FALSE
QTL_position <- rbind(dormancy, exsituFT)
QTL_position <- rbind(QTL_position, popsize)
QTL_position <- rbind(QTL_position, envIndexes)
QTL_position <- rbind(QTL_position, envFactors_subdata)


#' p_value treshold 0.01
#+ message=FALSE, warning=FALSE
QTL_position$Effect <- as.factor(QTL_position$Effect)
QTL_position$pvalue <- as.numeric(QTL_position$pvalue)
QTL_position <- QTL_position %>% mutate(Effect = ifelse(pvalue >= 0.01, "Neutral", Effect))
QTL_position <- QTL_position %>% mutate(Effect = ifelse(Effect == 1, "Negatif", Effect))
QTL_position <- QTL_position %>% mutate(Effect = ifelse(Effect == 3, "Positif", Effect))
QTL_position$Effect <- as.factor(QTL_position$Effect)


#' Threshold with FDR
#' threshold -log10(p) > 2.5
#+ message=FALSE, warning=FALSE
QTL_position$Effect <- as.factor(QTL_position$Effect)
QTL_position$pvalue <- as.numeric(QTL_position$pvalue)
QTL_position$logpvalue <- -log10(QTL_position$pvalue)
QTL_position <- QTL_position %>% mutate(Effect = ifelse(logpvalue < 2.5, "Neutral", Effect))
QTL_position <- QTL_position %>% mutate(Effect = ifelse(Effect == 1, "Negatif", Effect))
QTL_position <- QTL_position %>% mutate(Effect = ifelse(Effect == 3, "Positif", Effect))
QTL_position$Effect <- as.factor(QTL_position$Effect)
QTL_position <- QTL_position %>% mutate(logpvalue = ifelse(R >0, logpvalue, logpvalue*-1))

#' threshold -log10(p) > 2.9
#+ message=FALSE, warning=FALSE
envFactors_subdata$Effect <- as.factor(envFactors_subdata$Effect)
envFactors_subdata$pvalue <- as.numeric(envFactors_subdata$pvalue)
envFactors_subdata$logpvalue <- -log10(envFactors_subdata$pvalue)
envFactors_subdata <- envFactors_subdata %>% mutate(Effect = ifelse(logpvalue < 2.9, "Neutral", Effect))
envFactors_subdata <- envFactors_subdata %>% mutate(Effect = ifelse(Effect == 1, "Negatif", Effect))
envFactors_subdata <- envFactors_subdata %>% mutate(Effect = ifelse(Effect == 3, "Positif", Effect))
envFactors_subdata$Effect <- as.factor(envFactors_subdata$Effect)
envFactors_subdata <- envFactors_subdata %>% mutate(logpvalue = ifelse(R >0, logpvalue, logpvalue*-1))

QTL_position <- rbind(QTL_position, envFactors_subdata)


#' Add color and effect columns for graphic
#+ message=FALSE, warning=FALSE
QTL_position <- QTL_position %>% mutate(Effect_color = paste(Factor, Effect, sep = "_"))
QTL_position <- QTL_position %>% mutate(Effect_color = if_else(Effect == "Neutral", "Neutral", Factor))
QTL_position$Effect_color <- as.factor(QTL_position$Effect_color)
summary(QTL_position)

phenotype <- c("Dormancy", "exsituFT", "Popsize" )
Indexes <- c("Moss", "Type", "Shannon", "Temperature", "Moisture", "Reaction", "Nutrients", "Salinity", "Disturbance", "Sealed", "Vege_cover","Light")
Factors <- c("C.N.ratio", "DD", "NO3", "PO4", "SO4", "Cl","F", "Tmax_04", "Tmean_04", "Tmin_04")
QTL_position <- QTL_position %>% mutate(Category = ifelse(Factor %in% phenotype, "Phenotype", "Indexes" ))
QTL_position <- QTL_position %>% mutate(Category = ifelse(Effect == "Neutral", "Neutral", Category ))
QTL_position <- QTL_position %>% mutate(Category = ifelse(Factor %in% Factors, "Factors", Category ))

QTL_position$logpvalue <- -log10(QTL_position$pvalue)
QTL_position <- QTL_position %>% mutate(logpvalue = ifelse(R >0, logpvalue, logpvalue*-1))

QTL_position %>% distinct(Factor) %>% print(n=25)
write.csv(QTL_position, file="QTL_position.csv")



#' # Graphic

#' ## Color palet
#+ message=FALSE, warning=FALSE
color_all <- c("Dormancy" ="#68228B", 
               "exsituFT" = "#EE3A8C", 
               "Popsize" ="#EEA9B8", 
               "Neutral"= "#D1D1D1",
               
               "Type" = "#EE9A49", 
               "Moss" ="#008B45",
               "Shannon" = "#6CA6CD",
               "Temperature" = "#0000FF", 
               "Reaction" = "#76EEC6", 
               "Nutrients" = "#AB82FF", 
               "Salinity" = "#8B864E", 
               "Disturbance" = "#EE4000", 
               "Sealed"= "#363636", 
               "Vege_cover"= "#66CD00",
               "Light" = "#EEC900",
               "Moisture" = "#B0E2FF",
               
               "C.N.ratio" = "#FF7F00",
               "DD" = "#00CD00", 
               "NO3" = "#00CED1",
               "PO4" ="#1C86EE",
               "SO4" = "#EEC900",
               "Cl" = "#AB82FF",
               "F" = "#4A708B",
               "Tmax_04" = "#8B008B",
               "Tmean_04" = "#FFAEB9")



#' # Graphic Phenotype
#=============================================================================================================

#+ message=FALSE, warning=FALSE
QTL_position  %>% 
  filter(Factor %in% phenotype ) %>%
  ggplot(aes(x=Pos, y=-log10(pvalue), color = Effect_color, alpha = Category, size = Category)) + 
  geom_point() +  
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = color_all)+
  scale_alpha_manual(values = c(0.2, 0.5, 0.5))+
  scale_size_manual(values = c(1,2,2))+
  ggtitle(paste("QTL positions of phenotype, median AF per 50kb windows (+10 SNP per window)"))


QTL_position  %>% 
  filter(Factor %in% phenotype ) %>%
  ggplot(aes(x=Pos, y=logpvalue, color = Effect_color, alpha = Effect, size = Effect)) + 
  geom_point() +  
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = color_all)+
  scale_alpha_manual(values = c(0.8, 0.1, 0.8))+
  scale_size_manual(values = c(2,1.5,2))+
  ggtitle(paste("QTL positions of phenotype, median AF per 50kb windows (+10 SNP per window)"))
#ggsave("QTL positions of phenotype - median AF (+10 SNP per 50kb windows)_effect.svg", plot =  last_plot(),dpi = 300, bg = "white", scale = 2 )





#' # Graphic Environmental Indexes
#=============================================================================================================

#+ message=FALSE, warning=FALSE
QTL_position  %>% 
  filter(Factor %in% Indexes ) %>%
  ggplot(aes(x=Pos, y=-log10(pvalue), color = Effect_color, alpha = Category, size = Category)) + 
  geom_point() +  
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = color_all)+
  scale_alpha_manual(values = c(0.5, 0.2, 0.5))+
  scale_size_manual(values = c(2,1,2))+
  ggtitle(paste("QTL positions of environmental indexes, median AF per 50kb windows (+10 SNP per window)"))

QTL_position  %>% 
  filter(Factor %in% Indexes ) %>%
  ggplot(aes(x=Pos, y=logpvalue, color = Effect_color,alpha = Category, size = Category)) + 
  geom_point(size = 2.5) +  
  geom_hline(yintercept = 0, colour = "darkgrey") +
  ylim(-4,4) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = color_all)+
  scale_alpha_manual(values = c(0.5, 0.2, 0.5))+
  scale_size_manual(values = c(2,1,2))+
  ggtitle(paste("QTL positions of environmental indexes, median AF per 50kb windows (+10 SNP per window)"))






#' ## Graphic Environmental Factors subdata
#=============================================================================================================

#+ message=FALSE, warning=FALSE
QTL_position  %>% 
  filter(Factor %in% Factors ) %>%
  ggplot(aes(x=Pos, y=logpvalue, color = Effect_color,alpha = Effect, size = Effect)) + 
  geom_point() +  
  geom_hline(yintercept = 0, colour = "darkgrey") +
  ylim(-4,4) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = color_all)+
  scale_alpha_manual(values = c(0.8, 0.1, 0.8))+
  scale_size_manual(values = c(2,1.5,2))+
  ggtitle(paste("QTL positions of environmental indexes, median AF per 50kb windows (+10 SNP per window)"))




#' # Graphic all QTL
#=============================================================================================================

#+ message=FALSE, warning=FALSE
QTL_position  %>% 
  ggplot(aes(x=Pos, y=-log10(pvalue), color = Effect_color, alpha = Category, size = Category)) + 
  geom_point() +  
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = color_all)+
  scale_alpha_manual(values = c(0.5, 0.5, 0.2, 0.5))+
  scale_size_manual(values = c(2, 2,1,2))+
  ggtitle(paste("QTL positions median AF per 50kb windows (+10 SNP per window)"))


QTL_position  %>% 
  ggplot(aes(x=Pos, y=logpvalue, color = Effect_color, alpha = Effect, size = Category, shape = Category)) + 
  geom_point() +  
  geom_hline(yintercept = 0, colour = "darkgrey") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = color_all)+
  scale_alpha_manual(values = c(0.8, 0.1, 0.8))+
  scale_size_manual(values = c(2,2,2,1.5))+
  scale_shape_manual(values = c(8,19,19,5))+
  ggtitle(paste("QTL positions median AF per 50kb windows (+10 SNP per window)"))
#ggsave("QTL positions - median AF (+10 SNP per 50kb windows)_effect.svg", plot =  last_plot(),dpi = 120, bg = "white", scale = 2 )




#' # Per chrm
#=============================================================================================================

#+ message=FALSE, warning=FALSE
QTL_position  %>% 
  filter(Chrm == "Chr1" ) %>%
  ggplot(aes(x=Pos, y=logpvalue, color = Effect_color, alpha = Category, shape = Category)) + 
  geom_point(size = 2) +  
  geom_jitter( height = 0.75) +
  geom_hline(yintercept = 0, colour = "darkgrey") +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = color_all)+
  scale_alpha_manual(values = c(0.7,0.7, 0.2, 0.7))+
  scale_shape_manual(values= c(19,19, 19 ,5))+
  ggtitle(paste("Chr1 positions of phenotype, median AF per 50kb windows (+10 SNP per window)"))
#ggsave("Chr1 QTL positions - median AF (+10 SNP per 50kb windows).svg", plot =  last_plot(),dpi = 300, bg = "white", scale = 2 )


QTL_position  %>% 
  filter(Chrm == "Chr2" ) %>%
  ggplot(aes(x=Pos, y=logpvalue, color = Effect_color, alpha = Category, shape = Category)) + 
  geom_point(size = 2) +  
  geom_jitter( height = 0.75) +
  geom_hline(yintercept = 0, colour = "darkgrey") +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = color_all)+
  scale_alpha_manual(values = c(0.7,0.7, 0.2, 0.7))+
  scale_shape_manual(values= c(19,19, 19 ,5))+
  ggtitle(paste("Chr2 positions of phenotype, median AF per 50kb windows (+10 SNP per window)"))
#ggsave("Chr2 QTL positions - median AF (+10 SNP per 50kb windows).svg", plot =  last_plot(),dpi = 300, bg = "white", scale = 2 )

QTL_position  %>% 
  filter(Chrm == "Chr2" ) %>%
  filter(Factor %in% Factors | Factor %in% phenotype ) %>%
  ggplot(aes(x=Pos, y=logpvalue, color = Effect_color, alpha = Category, shape = Category)) + 
  geom_point(size = 2) +  
  geom_jitter( height = 0.75) +
  geom_hline(yintercept = 0, colour = "darkgrey") +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = color_all)+
  scale_alpha_manual(values = c(0.7,0.7, 0.2, 0.7))+
  scale_shape_manual(values= c(19,19, 19 ,5))+
  ggtitle(paste("Chr2 positions of phenotype, median AF per 50kb windows (+10 SNP per window)"))
#ggsave("Chr4 QTL positions - median AF (+10 SNP per 50kb windows).svg", plot =  last_plot(),dpi = 300, bg = "white", scale = 2 )



QTL_position  %>% 
  filter(Chrm == "Chr3" ) %>%
  ggplot(aes(x=Pos, y=logpvalue, color = Effect_color, alpha = Category, shape = Category)) + 
  geom_point(size = 2) +  
  geom_jitter( height = 0.75) +
  geom_hline(yintercept = 0, colour = "darkgrey") +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = color_all)+
  scale_alpha_manual(values = c(0.7,0.7, 0.2, 0.7))+
  scale_shape_manual(values= c(19,19, 19 ,5))+
  ggtitle(paste("Chr3 positions of phenotype, median AF per 50kb windows (+10 SNP per window)"))
#ggsave("Chr3 QTL positions - median AF (+10 SNP per 50kb windows).svg", plot =  last_plot(),dpi = 300, bg = "white", scale = 2 )


QTL_position  %>% 
  filter(Chrm == "Chr4" ) %>%
  filter(logpvalue > -5) %>%
  ggplot(aes(x=Pos, y=logpvalue, color = Effect_color, alpha = Category, shape = Category)) + 
  geom_point(size = 2) +  
  geom_jitter( height = 0.75) +
  geom_hline(yintercept = 0, colour = "darkgrey") +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = color_all)+
  scale_alpha_manual(values = c(0.7,0.7, 0.2, 0.7))+
  scale_shape_manual(values= c(19,19, 19 ,5))+
  ggtitle(paste("Chr4 positions of phenotype, median AF per 50kb windows (+10 SNP per window)"))
#ggsave("Chr4 QTL positions - median AF (+10 SNP per 50kb windows).svg", plot =  last_plot(),dpi = 300, bg = "white", scale = 2 )

QTL_position  %>% 
  filter(Chrm == "Chr4" ) %>%
  ggplot(aes(x=Pos, y=logpvalue, color = Effect_color, alpha = Category, shape = Category)) + 
  geom_point(size = 2) +  
  geom_jitter( height = 0.75) +
  geom_hline(yintercept = 0, colour = "darkgrey") +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = color_all)+
  scale_alpha_manual(values = c(0.7,0.7, 0.2, 0.7))+
  scale_shape_manual(values= c(19,19, 19 ,5))+
  ggtitle(paste("Chr4 positions of phenotype, median AF per 50kb windows (+10 SNP per window)"))
#ggsave("Chr4 QTL positions - median AF (+10 SNP per 50kb windows).svg", plot =  last_plot(),dpi = 300, bg = "white", scale = 2 )


QTL_position  %>% 
  filter(Chrm == "Chr4" ) %>%
  filter(Factor %in% Factors | Factor %in% phenotype ) %>%
  ggplot(aes(x=Pos, y=logpvalue, color = Effect_color, alpha = Category, shape = Category)) + 
  geom_point(size = 2) +  
  geom_jitter( height = 0.75) +
  geom_hline(yintercept = 0, colour = "darkgrey") +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = color_all)+
  scale_alpha_manual(values = c(0.7,0.7, 0.2, 0.7))+
  scale_shape_manual(values= c(19,19, 19 ,5))+
  ggtitle(paste("Chr4 positions of phenotype, median AF per 50kb windows (+10 SNP per window)"))
#ggsave("Chr4 QTL positions - median AF (+10 SNP per 50kb windows).svg", plot =  last_plot(),dpi = 300, bg = "white", scale = 2 )



QTL_position  %>% 
  filter(Chrm == "Chr5" ) %>%
  ggplot(aes(x=Pos, y=logpvalue, color = Effect_color, alpha = Category, shape = Category)) + 
  geom_point(size = 2) +  
  geom_jitter( height = 0.75) +
  geom_hline(yintercept = 0, colour = "darkgrey") +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = color_all)+
  scale_alpha_manual(values = c(0.7,0.7, 0.2, 0.7))+
  scale_shape_manual(values= c(19,19, 19 ,5))+
  ggtitle(paste("Chr5 positions of phenotype, median AF per 50kb windows (+10 SNP per window)"))
#ggsave("Chr5 QTL positions - median AF (+10 SNP per 50kb windows).svg", plot =  last_plot(),dpi = 300, bg = "white", scale = 2 )


