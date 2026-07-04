#' ---
#'title: "Allelic frequency per environmental indexes sum"
#'author: "Justine Floret"
#'date: "2026-04-11"
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
library(forcats) # re order hte x axis
library(stringr) #remove space in string
# Stat. libray
library(lmerTest) #glht
library(lsmeans) #pairwise comparison
library(emmeans)
library(glmmTMB) # negative binomial model
library(vegan) #RDA
library(missMDA) #Missing data imputation, pca
library("Hmisc")  # correlation
library("corrplot")  # correlation
library("plot3D")



#' # I. Preparing the data

#' ## Allelic frequency data
#+ message=FALSE, warning=FALSE
Median_AF_50kb_noContamination <- read_csv("/home/justine/sciebo/Genomic/assays/F5_pools_medianAF_50kb_noContamination.csv")

Median_AF_50kb_noContamination$Chrm <- as.factor(Median_AF_50kb_noContamination$Chrm)
Median_AF_50kb_noContamination$Pos_mean <- as.numeric(Median_AF_50kb_noContamination$Pos_mean)
Median_AF_50kb_noContamination$Sample_ID <- as.factor(Median_AF_50kb_noContamination$Sample_ID)
Median_AF_50kb_noContamination$median_AF <- as.numeric(Median_AF_50kb_noContamination$median_AF)
Median_AF_50kb_noContamination$SNP_tot <- as.numeric(Median_AF_50kb_noContamination$SNP_tot)

summary(Median_AF_50kb_noContamination)


#' ## Environmental data
#+ message=FALSE, warning=FALSE
Env_Index <- read_csv("/home/justine/sciebo/2_Ecological_parameters/2_assays/Ecological_parameters.csv")

Env_Index$Site_ID <- as.factor(Env_Index$Site_ID)
Env_Index$Type <- as.factor(Env_Index$Type)

Env_Index <- Env_Index %>%
  mutate(Type = recode(Type, "grass" = "Meadow",
                       "pavement" = "Pavement",
                       "Tree bed" = "Tree_bed", 
                       "wall"="Wall"))

Env_Index$Light_intensity <- -1* Env_Index$Shading
summary(Env_Index)

#' Keeping only the ecological variables associating with the establishment. 
#+ message=FALSE, warning=FALSE
factor_estblishment <- c("Site_ID", "Type", "Moss_abundance", "Shannon_index", "Light_intensity", 
                 "sum_Disturbance.Severity", "sum_Moisture", "sum_Temperature", "sum_Reaction", "sum_Nutrients", "sum_Salinity",
                 "Vegetation_cover_05_22", "Sealed_Surface_05_22")

Env_Index <- Env_Index[,(names(Env_Index) %in% factor_estblishment)] 

Env_Index_standardized <- decostand(Env_Index[,c(3:13)], method = "standardize") # Scale and center variables
Env_Index_standardized <- cbind(Env_Index[,c(1,2)], Env_Index_standardized)

#' ## Merged SNP data and the ecological factors
#+ message=FALSE, warning=FALSE
median_AF_Env_Index_standardized <- inner_join(Median_AF_50kb_noContamination, Env_Index_standardized, join_by("Sample_ID"== "Site_ID"))
median_AF_Env_Index_standardized <- median_AF_Env_Index_standardized %>% filter(SNP_tot>10)
summary(median_AF_Env_Index_standardized)
median_AF_Env_Index_standardized %>% distinct(Sample_ID) %>% count()

median_AF_Env_Index <- inner_join(Median_AF_50kb_noContamination, Env_Index, join_by("Sample_ID"== "Site_ID"))
median_AF_Env_Index <- median_AF_Env_Index %>% filter(SNP_tot>10)




#' # II. Association AF - Environmental indexes

#+ message=FALSE, warning=FALSE
for (chr in chrm_list) {
  data_chr <- median_AF_Env_Index_standardized %>% filter(Chrm == chr)
  positions <- unique(data_chr$Pos_mean)
  for (pos in positions) {
    data_pos <- data_chr %>% filter(Pos_mean == pos)
    if (nrow(data_pos) == 0) next
    corresults_Type <- cor.test(data_pos$median_AF, data_pos$Type,  method = "spearman")
    corresults_Moss <- cor.test(data_pos$median_AF, data_pos$Moss_abundance,  method = "spearman")
    corresults_Shannon <- cor.test(data_pos$median_AF, data_pos$Shannon_index,  method = "spearman")
    corresults_Temperature <- cor.test(data_pos$median_AF, data_pos$sum_Temperature,  method = "spearman")
    corresults_Moisture <- cor.test(data_pos$median_AF, data_pos$sum_Moisture,  method = "spearman")
    corresults_Reaction <- cor.test(data_pos$median_AF, data_pos$sum_Reaction,  method = "spearman")
    corresults_Nutrients <- cor.test(data_pos$median_AF, data_pos$sum_Nutrients,  method = "spearman")
    corresults_Salinity <- cor.test(data_pos$median_AF, data_pos$sum_Salinity,  method = "spearman")
    corresults_Disturbance <- cor.test(data_pos$median_AF, data_pos$sum_Disturbance.Severity,  method = "spearman")
    corresults_Sealed <- cor.test(data_pos$median_AF, data_pos$Sealed_Surface_05_22,  method = "spearman")
    corresults_Vegetation <- cor.test(data_pos$median_AF, data_pos$Vegetation_cover_05_22,  method = "spearman")
    corresults_Light <- cor.test(data_pos$median_AF, data_pos$Light_intensity,  method = "spearman")
    
    tmp <- c(chr,pos, 
             corresults_Type[["p.value"]], corresults_Type[["estimate"]], 
             corresults_Moss[["p.value"]], corresults_Moss[["estimate"]], 
             corresults_Shannon[["p.value"]], corresults_Shannon[["estimate"]], 
             corresults_Temperature[["p.value"]], corresults_Temperature[["estimate"]], 
             corresults_Moisture[["p.value"]], corresults_Moisture[["estimate"]], 
             corresults_Reaction[["p.value"]], corresults_Reaction[["estimate"]], 
             corresults_Nutrients[["p.value"]], corresults_Nutrients[["estimate"]], 
             corresults_Salinity[["p.value"]], corresults_Salinity[["estimate"]], 
             corresults_Disturbance[["p.value"]], corresults_Disturbance[["estimate"]], 
             corresults_Sealed[["p.value"]], corresults_Sealed[["estimate"]], 
             corresults_Vegetation[["p.value"]], corresults_Vegetation[["estimate"]], 
             corresults_Light[["p.value"]], corresults_Light[["estimate"]])
    
    Env_model <- rbind(Env_model, tmp)
  }
}

Env_model <- Env_model[-1,] 
Env_model <- Env_model %>% mutate(across(c(2:26), as.numeric)) 




#' # III. Randomization x 1000 

#+ eval=F, echo=T
EnvIndex_rand_df <- vector("list", 1000)  

#+ eval=F, echo=T
for (i in 1:1000) {
  # 1) Randomize 
  pool_map <- median_AF_Env_Index_standardized %>%
    distinct(Sample_ID,.keep_all = TRUE) %>%
    mutate(Randomized_Type = sample(Type),
           Randomized_Moss = sample(Moss_abundance),
           Randomized_Shannon = sample(Shannon_index),
           Randomized_Temperature = sample(sum_Temperature),
           Randomized_Moisture = sample(sum_Moisture),
           Randomized_Reaction = sample(sum_Reaction),
           Randomized_Nutrients = sample(sum_Nutrients),
           Randomized_Salinity = sample(sum_Salinity),
           Randomized_Disturbance = sample(sum_Disturbance.Severity),
           Randomized_Sealed = sample(Sealed_Surface_05_22),
           Randomized_Vegetation = sample(Vegetation_cover_05_22),
           Randomized_Light = sample(Light_intensity))
  data_random <- left_join(median_AF_Env_Index_standardized[,1:5], pool_map[, c(3, 18:29)], by = "Sample_ID")
  
  # 2) models for all Chr × Pos
  chrm_list <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
  EnvIndex_rand_model <- tibble("Chrm" = "1", "Pos" = 1,
                                Type_pvalue=1, Type_R=1, 
                                Moss_pvalue = 1, Moss_R=1,
                                Shannon_pvalue = 1, Shannon_R=1,
                                Temperature_pvalue = 1, Temperature_R=1,
                                Moisture_pvalue = 1, Moisture_R=1,
                                Reaction_pvalue = 1, Reaction_R=1,
                                Nutrients_pvalue = 1, Nutrients_R=1,
                                Salinity_pvalue = 1, Salinity_R=1,
                                Disturbance_pvalue = 1, Disturbance_R=1,
                                Sealed_pvalue = 1, Sealed_R=1, 
                                Vege_cover_pvalue = 1,  Vege_cover_R=1,
                                Light_pvalue = 1,  Light_R=1)
  
  for (chr in chrm_list) {
    data_chr <- data_random %>% filter(Chrm == chr)
    positions <- unique(data_chr$Pos_mean)
    for (pos in positions) {
      data_pos <- data_chr %>% filter(Pos_mean == pos)
      if (nrow(data_pos) == 0) next
      corresults_Type <- cor.test(data_pos$median_AF, data_pos$Randomized_Type,  method = "spearman")
      corresults_Moss <- cor.test(data_pos$median_AF, data_pos$Randomized_Moss,  method = "spearman")
      corresults_Shannon <- cor.test(data_pos$median_AF, data_pos$Randomized_Shannon,  method = "spearman")
      corresults_Temperature <- cor.test(data_pos$median_AF, data_pos$Randomized_Temperature,  method = "spearman")
      corresults_Moisture <- cor.test(data_pos$median_AF, data_pos$Randomized_Moisture,  method = "spearman")
      corresults_Reaction <- cor.test(data_pos$median_AF, data_pos$Randomized_Reaction,  method = "spearman")
      corresults_Nutrients <- cor.test(data_pos$median_AF, data_pos$Randomized_Nutrients,  method = "spearman")
      corresults_Salinity <- cor.test(data_pos$median_AF, data_pos$Randomized_Salinity,  method = "spearman")
      corresults_Disturbance <- cor.test(data_pos$median_AF, data_pos$Randomized_Disturbance,  method = "spearman")
      corresults_Sealed <- cor.test(data_pos$median_AF, data_pos$Randomized_Sealed,  method = "spearman")
      corresults_Vegetation <- cor.test(data_pos$median_AF, data_pos$Randomized_Vegetation,  method = "spearman")
      corresults_Light <- cor.test(data_pos$median_AF, data_pos$Randomized_Light,  method = "spearman")
      
    
      tmp <- c(chr,pos, 
               corresults_Type[["p.value"]], corresults_Type[["estimate"]], 
               corresults_Moss[["p.value"]], corresults_Moss[["estimate"]], 
               corresults_Shannon[["p.value"]], corresults_Shannon[["estimate"]], 
               corresults_Temperature[["p.value"]], corresults_Temperature[["estimate"]], 
               corresults_Moisture[["p.value"]], corresults_Moisture[["estimate"]], 
               corresults_Reaction[["p.value"]], corresults_Reaction[["estimate"]], 
               corresults_Nutrients[["p.value"]], corresults_Nutrients[["estimate"]], 
               corresults_Salinity[["p.value"]], corresults_Salinity[["estimate"]], 
               corresults_Disturbance[["p.value"]], corresults_Disturbance[["estimate"]], 
               corresults_Sealed[["p.value"]], corresults_Sealed[["estimate"]], 
               corresults_Vegetation[["p.value"]], corresults_Vegetation[["estimate"]], 
               corresults_Light[["p.value"]], corresults_Light[["estimate"]])
      
      EnvIndex_rand_model <- rbind(EnvIndex_rand_model, tmp)
    }
  }
      
  # 3) Post-processing
  EnvIndex_rand_model <- EnvIndex_rand_model[-1, ]
  
  # save result of this i
  EnvIndex_rand_df[[i]] <- EnvIndex_rand_model
  print(paste("randomization", i, "done"))
}


#' Create one data frame
#+ eval=F, echo=T
#saveRDS(EnvIndexes_rand_models, "EnvIndexes_rand_models" )
EnvIndex_rand_models_df <- purrr::map2_df(
  EnvIndex_rand_df,
  1:1000,
  ~ dplyr::mutate(.x, iter = .y))


#' Merge the data set  
#+ eval=F, echo=T
Env_model$iter <- "Observed"
EnvIndex_rand_models_df <- rbind(EnvIndex_rand_models_df,Env_model)
EnvIndex_rand_models_df <- EnvIndex_rand_models_df %>% mutate(across(c(2:26), as.numeric)) 




#' # IV. Analysis observed vs randomized data

#+ message=FALSE, warning=FALSE
Env_model <- EnvIndex_rand_models_df %>% filter(iter == "Observed")



#' # Type

#' ## Estimates
#+ message=FALSE, warning=FALSE
EnvIndex_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Type_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Type model estimate mod(Alt_frq ~ Type) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- EnvIndex_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Type_R, probs = c(.95)),
                                    Q05 = quantile(Type_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Env_model[,c(1,2,3,4)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Type_R>Q95 & Type_pvalue <0.05 | Type_R<Q05 & Type_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Type = ifelse(Estimate_5 == "TRUE" & Type_R >0, "Positif",
                                                           ifelse(Estimate_5 == "TRUE" & Type_R <0, "Negatif", "Neutral")))
EnvFactors_effects <- Quantile_df[,c(1,2,5,6,8)]

#' ## Graphic p-value & estimate
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Type_pvalue), color = Effect_Type, alpha = Effect_Type)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c( "red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Type - median AF (+10 SNP per 50kb windows)"))


#' ## Correlation allelic frequency variable
Quantile_df %>% filter(Estimate_5 ==TRUE) %>% slice_min(Type_pvalue, n= 1)

type_df <- median_AF_Env_Index %>% filter(Chrm == "Chr3" & Pos_mean == "2325000")
type_df <- type_df %>% mutate(across(c(2:17), as.numeric)) 
type_df <- type_df %>%
  mutate(Type = recode(Type, "1" = "Meadow",
                       "2" = "Pavement",
                       "3" = "Tree_bed", 
                       "4"="Wall"))
type_df$Type <- as.factor(type_df$Type)
colors <- c("#6ece58", "#22a884", "#355f8d","#440154")


type_chr3 <- glmmTMB(median_AF ~ Type, 
                     data = type_df %>% filter(Chrm== "Chr3"),
                     family = gaussian)
summary(type_chr3)
car::Anova(type_chr3, type ="II")

type_emm <- emmeans::ref_grid(type_chr3)
emm <- emmeans::emmeans(type_emm, ~ Type)
pairs(emm )


type_df %>% filter(Chrm== "Chr3") %>%
  mutate(Type = fct_relevel(Type, "Meadow", "Tree_bed", "Pavement","Wall")) %>%
  ggplot(aes(x=Type, y=median_AF, fill = Type)) + 
  geom_boxplot() +
  scale_fill_manual(values = colors)+
  theme_minimal()



#' # Moss abundance

#' ## Estimates
#+ message=FALSE, warning=FALSE
EnvIndex_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Moss_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Moss abundace model estimate mod(Alt_frq ~ Moss) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- EnvIndex_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Moss_R, probs = c(.95)),
                                    Q05 = quantile(Moss_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Env_model[,c(1,2,5,6)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Moss_R>Q95 & Moss_pvalue <0.05 | Moss_R<Q05 & Moss_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Moss = ifelse(Estimate_5 == "TRUE" & Moss_R >0, "Positif",
                                                           ifelse(Estimate_5 == "TRUE" & Moss_R <0, "Negatif", "Neutral")))
EnvIndexes_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], , by = c("Chrm", "Pos"))

#' ## Graphic p-value & estimate
#+ message=FALSE, warning=FALSE
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Moss_pvalue), color = Effect_Moss, alpha = Effect_Moss)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c( "red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Moss abundance - median AF (+10 SNP per 50kb windows)"))


#' ## Correlation allelic frequency variable
Quantile_df %>% filter(Estimate_5 ==TRUE) %>% slice_min(Moss_pvalue, n= 1)
moss_df <- median_AF_Env_Index %>% filter(Chrm == "Chr1" & Pos_mean == "5275000")
moss_df <- moss_df %>% mutate(across(c(2:17), as.numeric)) 

allelic_frq <- as.numeric(unlist(moss_df[moss_df$Chrm == "Chr1","median_AF"]))
variable <- as.numeric(unlist(moss_df[moss_df$Chrm == "Chr1","Moss_abundance"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation
moss_df %>%
  ggplot(aes(x=median_AF, y=Moss_abundance)) +
  geom_point(size=4, color="#008B45") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.75, y=0.8, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("Moss abundance ~ allelic frequency Chr1"))





#' # Shannon Index

#' ## Estimates
#+ message=FALSE, warning=FALSE
EnvIndex_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Shannon_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Shannon Index model estimate mod(Alt_frq ~ Shannon) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- EnvIndex_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Shannon_R, probs = c(.95)),
                                    Q05 = quantile(Shannon_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Env_model[,c(1,2,7,8)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Shannon_R>Q95 & Shannon_pvalue <0.05 | Shannon_R<Q05 & Shannon_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Shannon = ifelse(Estimate_5 == "TRUE" & Shannon_R >0, "Positif",
                                                              ifelse(Estimate_5 == "TRUE" & Shannon_R <0, "Negatif", "Neutral")))

EnvIndexes_effects <- left_join(EnvIndexes_effects, Quantile_df[,c(1,2,5,6,8)], , by = c("Chrm", "Pos"))

#' ## Graphic p-value & estimate
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Shannon_pvalue), color = Effect_Shannon, alpha = Effect_Shannon)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c( "red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Shannon Index - median AF (+10 SNP per 50kb windows)"))


#' ## Correlation allelic frequency variable
Quantile_df %>% filter(Estimate_5 ==TRUE) %>% slice_min(Shannon_pvalue, n= 1)
shannon_df <- median_AF_Env_Index %>% filter(Chrm == "Chr2" & Pos_mean == "12625000")
shannon_df <- shannon_df %>% mutate(across(c(2:17), as.numeric)) 

allelic_frq <- as.numeric(unlist(shannon_df[shannon_df$Chrm == "Chr2","median_AF"]))
variable <- as.numeric(unlist(shannon_df[shannon_df$Chrm == "Chr2","Shannon_index"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation
shannon_df %>% 
  ggplot(aes(x=median_AF, y=Shannon_index)) +
  geom_point(size=4, color="#6CA6CD") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.75, y=2.5, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("Shannon_index ~ allelic frequency Chr2"))




#' # Temperature 

#' ## Estimates
#+ message=FALSE, warning=FALSE
EnvIndex_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Temperature_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Temperature model estimate mod(Alt_frq ~ Temperature) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- EnvIndex_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Temperature_R, probs = c(.95)),
                                    Q05 = quantile(Temperature_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Env_model[,c(1,2,9,10)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Temperature_R>Q95 & Temperature_pvalue <0.05 | Temperature_R<Q05 & Temperature_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Temperature = ifelse(Estimate_5 == "TRUE" & Temperature_R >0, "Positif",
                                                           ifelse(Estimate_5 == "TRUE" & Temperature_R <0, "Negatif", "Neutral")))
EnvIndexes_effects <- left_join(EnvIndexes_effects, Quantile_df[,c(1,2,5,6,8)], , by = c("Chrm", "Pos"))


#' ## Graphic p-value & estimate
#+ message=FALSE, warning=FALSE
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Temperature_pvalue), color = Effect_Temperature, alpha = Effect_Temperature)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c("red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Temperature - median AF (+10 SNP per 50kb windows)"))


#' ## Correlation allelic frequency variable
Quantile_df %>% filter(Estimate_5 ==TRUE) %>% slice_min(Temperature_pvalue, n= 1)
temperature_df <- median_AF_Env_Index %>% filter( Chrm == "Chr1" & Pos_mean == "28525000" | 
                                               Chrm == "Chr2" & Pos_mean == "12175000" |
                                               Chrm == "Chr5" & Pos_mean == "4675000")
temperature_df <- temperature_df %>% mutate(across(c(2:17), as.numeric)) 


#' Chr1
allelic_frq <- as.numeric(unlist(temperature_df[temperature_df$Chrm == "Chr1","median_AF"]))
variable <- as.numeric(unlist(temperature_df[temperature_df$Chrm == "Chr1","sum_Temperature"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation
temperature_df %>% filter(Chrm== "Chr1") %>%
  ggplot(aes(x=median_AF, y=sum_Temperature)) +
  geom_point(size=4, color="#0000FF") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.75, y=7, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("Temperature index ~ allelic frequency Chr1"))


#' Chr2
allelic_frq <- as.numeric(unlist(temperature_df[temperature_df$Chrm == "Chr2","median_AF"]))
variable <- as.numeric(unlist(temperature_df[temperature_df$Chrm == "Chr2","sum_Temperature"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation
temperature_df %>% filter(Chrm== "Chr2") %>%
  ggplot(aes(x=median_AF, y=sum_Temperature)) +
  geom_point(size=4, color="#0000FF") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.25, y=7, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("Temperature index ~ allelic frequency Chr2"))


#' Chr5
allelic_frq <- as.numeric(unlist(temperature_df[temperature_df$Chrm == "Chr5","median_AF"]))
variable <- as.numeric(unlist(temperature_df[temperature_df$Chrm == "Chr5","sum_Temperature"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation
temperature_df %>% filter(Chrm== "Chr5") %>%
  ggplot(aes(x=median_AF, y=sum_Temperature)) +
  geom_point(size=4, color="#0000FF") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.65, y=7, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("Temperature index ~ allelic frequency Chr5"))




#' # Moisture

#' ## Estimates
#+ message=FALSE, warning=FALSE
EnvIndex_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Moisture_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Moisture model estimate mod(Alt_frq ~ Moisture) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- EnvIndex_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Moisture_R, probs = c(.95)),
                                    Q05 = quantile(Moisture_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Env_model[,c(1,2,11,12)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Moisture_R>Q95 & Moisture_pvalue <0.05 | Moisture_R<Q05 & Moisture_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Moisture = ifelse(Estimate_5 == "TRUE" & Moisture_R >0, "Positif",
                                                              ifelse(Estimate_5 == "TRUE" & Moisture_R <0, "Negatif", "Neutral")))

EnvIndexes_effects <- left_join(EnvIndexes_effects, Quantile_df[,c(1,2,5,6,8)], , by = c("Chrm", "Pos"))

#' ## Graphic p-value & estimate
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Moisture_pvalue), color = Effect_Moisture, alpha = Effect_Moisture)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c( "red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Moisture - median AF (+10 SNP per 50kb windows)"))

#' ## Correlation allelic frequency variable
Quantile_df %>% filter(Estimate_5 ==TRUE) %>% slice_min(Moisture_pvalue, n= 1)
moisture_df <- median_AF_Env_Index %>% filter(Chrm == "Chr2" & Pos_mean == "18925000")
moisture_df <- moisture_df %>% mutate(across(c(2:17), as.numeric)) 

allelic_frq <- as.numeric(unlist(moisture_df[moisture_df$Chrm == "Chr2","median_AF"]))
variable <- as.numeric(unlist(moisture_df[moisture_df$Chrm == "Chr2","sum_Moisture"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation
moisture_df %>% 
  ggplot(aes(x=median_AF, y=sum_Moisture)) +
  geom_point(size=4, color="#B0E2FF") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.75, y=4.5, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("Moisture index ~ allelic frequency Chr2"))




#' # Reaction

#' ## Estimates
#+ message=FALSE, warning=FALSE
EnvIndex_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Reaction_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Reaction model estimate mod(Alt_frq ~ Reaction) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- EnvIndex_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Reaction_R, probs = c(.95)),
                                    Q05 = quantile(Reaction_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Env_model[,c(1,2,13,14)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Reaction_R>Q95 & Reaction_pvalue <0.05 | Reaction_R<Q05 & Reaction_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Reaction = ifelse(Estimate_5 == "TRUE" & Reaction_R >0, "Positif",
                                                               ifelse(Estimate_5 == "TRUE" & Reaction_R <0, "Negatif", "Neutral")))

EnvIndexes_effects <- left_join(EnvIndexes_effects, Quantile_df[,c(1,2,5,6,8)], , by = c("Chrm", "Pos"))

#' ## Graphic p-value & estimate
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Reaction_pvalue), color = Effect_Reaction, alpha = Effect_Reaction)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c( "red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Reaction - median AF (+10 SNP per 50kb windows)"))


#' ## Correlation allelic frequency variable
Quantile_df %>% filter(Estimate_5 ==TRUE) %>% slice_min(Reaction_pvalue, n= 1)
reaction_df <- median_AF_Env_Index %>% filter( Chrm == "Chr4" & Pos_mean == "1525000" | 
                                                    Chrm == "Chr5" & Pos_mean == "20875000")
reaction_df <- reaction_df %>% mutate(across(c(2:17), as.numeric)) 


#' Chr4
allelic_frq <- as.numeric(unlist(reaction_df[reaction_df$Chrm == "Chr4","median_AF"]))
variable <- as.numeric(unlist(reaction_df[reaction_df$Chrm == "Chr4","sum_Reaction"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation
reaction_df %>% filter(Chrm== "Chr4") %>%
  ggplot(aes(x=median_AF, y=sum_Reaction)) +
  geom_point(size=4, color="#76EEC6") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.75, y=7, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("Reaction index ~ allelic frequency Chr4"))
#ggsave("sum Reaction ~ allelic frequency Chr4.svg", plot =  last_plot(),dpi = 120, bg = "white", scale = 1)


#' Chr5
allelic_frq <- as.numeric(unlist(reaction_df[reaction_df$Chrm == "Chr5","median_AF"]))
variable <- as.numeric(unlist(reaction_df[reaction_df$Chrm == "Chr5","sum_Reaction"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation
reaction_df %>% filter(Chrm== "Chr5") %>%
  ggplot(aes(x=median_AF, y=sum_Reaction)) +
  geom_point(size=4, color="#76EEC6") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.25, y=5, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("Reaction index ~ allelic frequency Chr5"))




#' # Nutrient

#' ## Estimates
#+ message=FALSE, warning=FALSE
EnvIndex_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Nutrients_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Nutrient model estimate mod(Alt_frq ~ Nutrient) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- EnvIndex_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Nutrients_R, probs = c(.95)),
                                    Q05 = quantile(Nutrients_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Env_model[,c(1,2,15,16)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Nutrients_R>Q95 & Nutrients_pvalue <0.05 | Nutrients_R<Q05 & Nutrients_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Nutrients = ifelse(Estimate_5 == "TRUE" & Nutrients_R >0, "Positif",
                                                               ifelse(Estimate_5 == "TRUE" & Nutrients_R <0, "Negatif", "Neutral")))

EnvIndexes_effects <- left_join(EnvIndexes_effects, Quantile_df[,c(1,2,5,6,8)], , by = c("Chrm", "Pos"))

#' ## Graphic p-value & estimate
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Nutrients_pvalue), color = Effect_Nutrients, alpha = Effect_Nutrients)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c("red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Nutrient - median AF (+10 SNP per 50kb windows)"))





#' # Salinity

#' ## Estimates
#+ message=FALSE, warning=FALSE
EnvIndex_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Salinity_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Salinity model estimate mod(Alt_frq ~ Salinity) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- EnvIndex_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Salinity_R, probs = c(.95)),
                                    Q05 = quantile(Salinity_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Env_model[,c(1,2,17,18)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Salinity_R>Q95 & Salinity_pvalue <0.05 | Salinity_R<Q05 & Salinity_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Salinity = ifelse(Estimate_5 == "TRUE" & Salinity_R >0, "Positif",
                                                               ifelse(Estimate_5 == "TRUE" & Salinity_R <0, "Negatif", "Neutral")))

EnvIndexes_effects <- left_join(EnvIndexes_effects, Quantile_df[,c(1,2,5,6,8)], , by = c("Chrm", "Pos"))

#' ## Graphic p-value & estimate
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Salinity_pvalue), color = Effect_Salinity, alpha = Effect_Salinity)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c("red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Salinity - median AF (+10 SNP per 50kb windows)"))



#' ## Correlation allelic frequency variable
Quantile_df %>% filter(Estimate_5 ==TRUE) %>% slice_min(Salinity_pvalue, n= 1)
salinity_df <- median_AF_Env_Index %>% filter(Chrm == "Chr2" & Pos_mean == "19275000")
salinity_df <- salinity_df %>% mutate(across(c(2:17), as.numeric)) 

allelic_frq <- as.numeric(unlist(salinity_df[salinity_df$Chrm == "Chr2","median_AF"]))
variable <- as.numeric(unlist(salinity_df[salinity_df$Chrm == "Chr2","sum_Salinity"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation
salinity_df %>% 
  ggplot(aes(x=median_AF, y=sum_Salinity)) +
  geom_point(size=4, color="#8B864E") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.25, y=0.85, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("Salinity index ~ allelic frequency Chr2"))




#' # Disturbance

#' ## Estimates
#+ message=FALSE, warning=FALSE
EnvIndex_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Disturbance_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Disturbance model estimate mod(Alt_frq ~ Disturbance) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#+ message=FALSE, warning=FALSE
Quantile_df <- EnvIndex_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Disturbance_R, probs = c(.95)),
                                    Q05 = quantile(Disturbance_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Env_model[,c(1,2,19,20)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Disturbance_R>Q95 & Disturbance_pvalue <0.01 | Disturbance_R<Q05 & Disturbance_pvalue <0.01 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Disturbance = ifelse(Estimate_5 == "TRUE" & Disturbance_R >0, "Positif",
                                                                  ifelse(Estimate_5 == "TRUE" & Disturbance_R <0, "Negatif", "Neutral")))
EnvIndexes_effects <- left_join(EnvIndexes_effects, Quantile_df[,c(1,2,5,6,8)], , by = c("Chrm", "Pos"))

#' ## Graphic p-value & estimate
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Disturbance_pvalue), color = Effect_Disturbance, alpha = Effect_Disturbance)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c("lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Disturbance - median AF (+10 SNP per 50kb windows)"))



#' ## Correlation allelic frequency variable
Quantile_df %>% filter(Estimate_5 ==TRUE) %>% slice_min(Disturbance_pvalue, n= 1)
disturbance_df <- median_AF_Env_Index %>% filter( Chrm == "Chr1" & Pos_mean == "29575000" | 
                                                Chrm == "Chr2" & Pos_mean == "8875000")
disturbance_df <- disturbance_df %>% mutate(across(c(2:17), as.numeric)) 


#' Chr1
allelic_frq <- as.numeric(unlist(disturbance_df[disturbance_df$Chrm == "Chr1","median_AF"]))
variable <- as.numeric(unlist(disturbance_df[disturbance_df$Chrm == "Chr1","sum_Disturbance.Severity"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation
disturbance_df %>% filter(Chrm== "Chr1") %>%
  ggplot(aes(x=median_AF, y=sum_Disturbance.Severity)) +
  geom_point(size=4, color="#EE4000") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.25, y=0.75, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("Disturbance severity index ~ allelic frequency Chr1"))


#' Chr2
allelic_frq <- as.numeric(unlist(disturbance_df[disturbance_df$Chrm == "Chr2","median_AF"]))
variable <- as.numeric(unlist(disturbance_df[disturbance_df$Chrm == "Chr2","sum_Disturbance.Severity"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation
disturbance_df %>% filter(Chrm== "Chr2") %>%
  ggplot(aes(x=median_AF, y=sum_Disturbance.Severity)) +
  geom_point(size=4, color="#EE4000") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.25, y=0.75, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("Disturbance severity index ~ allelic frequency Chr2"))




#' # Sealed surface  

#' ## Estimates
#+ message=FALSE, warning=FALSE
EnvIndex_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Sealed_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Sealed surface model estimate mod(Alt_frq ~ Sealed) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#+ message=FALSE, warning=FALSE
Quantile_df <- EnvIndex_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Sealed_R, probs = c(.95)),
                                    Q05 = quantile(Sealed_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Env_model[,c(1,2,21,22)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Sealed_R>Q95 & Sealed_pvalue <0.05 | Sealed_R<Q05 & Sealed_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Sealed = ifelse(Estimate_5 == "TRUE" & Sealed_R >0, "Positif",
                                                             ifelse(Estimate_5 == "TRUE" & Sealed_R <0, "Negatif", "Neutral")))
EnvIndexes_effects <- left_join(EnvIndexes_effects, Quantile_df[,c(1,2,5,6,8)], , by = c("Chrm", "Pos"))

#' ## Graphic p-value & estimate
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Sealed_pvalue), color = Effect_Sealed, alpha = Effect_Sealed)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c( "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Sealed surface - median AF (+10 SNP per 50kb windows)"))


#' ## Correlation allelic frequency variable
Quantile_df %>% filter(Estimate_5 ==TRUE) %>% slice_min(Sealed_pvalue, n= 1)
sealed_df <- median_AF_Env_Index %>% filter(Chrm == "Chr2" & Pos_mean == "675000")
sealed_df <- sealed_df %>% mutate(across(c(2:17), as.numeric)) 

allelic_frq <- as.numeric(unlist(sealed_df[sealed_df$Chrm == "Chr2","median_AF"]))
variable <- as.numeric(unlist(sealed_df[sealed_df$Chrm == "Chr2","Sealed_Surface_05_22"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation
sealed_df %>% 
  ggplot(aes(x=median_AF, y=Sealed_Surface_05_22)) +
  geom_point(size=4, color="#363636") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.25, y=0.95, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("Sealed surface ~ allelic frequency Chr2"))




#' # Vegetation cover

#' ## Estimates
#+ message=FALSE, warning=FALSE
EnvIndex_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Vege_cover_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Vegetation cover model estimate mod(Alt_frq ~ Vegetation cover) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- EnvIndex_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Vege_cover_R, probs = c(.95)),
                                    Q05 = quantile(Vege_cover_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Env_model[,c(1,2,23,24)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Vege_cover_R>Q95 & Vege_cover_pvalue <0.05 | Vege_cover_R<Q05 & Vege_cover_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Vege_cover = ifelse(Estimate_5 == "TRUE" & Vege_cover_R >0, "Positif",
                                                                 ifelse(Estimate_5 == "TRUE" & Vege_cover_R <0, "Negatif", "Neutral")))
EnvIndexes_effects <- left_join(EnvIndexes_effects, Quantile_df[,c(1,2,5,6,8)], , by = c("Chrm", "Pos"))

#' ## Graphic p-value & estimate
#+ message=FALSE, warning=FALSE
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Vege_cover_pvalue), color = Effect_Vege_cover, alpha = Effect_Vege_cover)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c( "red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Vegetation cover - median AF (+10 SNP per 50kb windows)"))



#' ## Correlation allelic frequency variable
Quantile_df %>% filter(Estimate_5 ==TRUE) %>% slice_min(Vege_cover_pvalue, n= 1)
vege_cover_df <- median_AF_Env_Index %>% filter( Chrm == "Chr1" & Pos_mean == "29325000" | 
                                                    Chrm == "Chr3" & Pos_mean == "175000" |
                                                    Chrm == "Chr5" & Pos_mean == "25000" )
vege_cover_df <- vege_cover_df %>% mutate(across(c(2:17), as.numeric)) 


#' Chr1
allelic_frq <- as.numeric(unlist(vege_cover_df[vege_cover_df$Chrm == "Chr1","median_AF"]))
variable <- as.numeric(unlist(vege_cover_df[vege_cover_df$Chrm == "Chr1","Vegetation_cover_05_22"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation
vege_cover_df %>% filter(Chrm== "Chr1") %>%
  ggplot(aes(x=median_AF, y=Vegetation_cover_05_22)) +
  geom_point(size=4, color="#66CD00") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.25, y=0.35, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("Vegetation cover index ~ allelic frequency Chr1"))


#' Chr3
allelic_frq <- as.numeric(unlist(vege_cover_df[vege_cover_df$Chrm == "Chr3","median_AF"]))
variable <- as.numeric(unlist(vege_cover_df[vege_cover_df$Chrm == "Chr3","Vegetation_cover_05_22"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation
vege_cover_df %>% filter(Chrm== "Chr3") %>%
  ggplot(aes(x=median_AF, y=Vegetation_cover_05_22)) +
  geom_point(size=4, color="#66CD00") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.25, y=0.35, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("Vegetation cover index ~ allelic frequency Chr3"))
ggsave("Vegetation cover ~ allelic frequency Chr3.svg", plot =  last_plot(),dpi = 120, bg = "white", scale = 1)


#' Chr5
allelic_frq <- as.numeric(unlist(vege_cover_df[vege_cover_df$Chrm == "Chr5","median_AF"]))
variable <- as.numeric(unlist(vege_cover_df[vege_cover_df$Chrm == "Chr5","Vegetation_cover_05_22"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation
vege_cover_df %>% filter(Chrm== "Chr5") %>%
  ggplot(aes(x=median_AF, y=Vegetation_cover_05_22)) +
  geom_point(size=4, color="#66CD00") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.25, y=0.35, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("Vegetation cover index ~ allelic frequency Chr5"))




#' # Light Intensity

#' ## Estimates
#+ message=FALSE, warning=FALSE
EnvIndex_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Light_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Light intensity model estimate mod(Alt_frq ~ Light) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- EnvIndex_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Light_R, probs = c(.95)),
                                    Q05 = quantile(Light_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Env_model[,c(1,2,25,26)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Light_R>Q95 & Light_pvalue <0.05 | Light_R<Q05 & Light_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Light = ifelse(Estimate_5 == "TRUE" & Light_R >0, "Positif",
                                                            ifelse(Estimate_5 == "TRUE" & Light_R <0, "Negatif", "Neutral")))
EnvIndexes_effects <- left_join(EnvIndexes_effects, Quantile_df[,c(1,2,5,6,8)], , by = c("Chrm", "Pos"))

#' ## Graphic p-value & estimate
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Light_pvalue), color = Effect_Light, alpha = Effect_Light)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values =  c( "red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Light intensity - median AF (+10 SNP per 50kb windows)"))




#' # Combined effect of the factor

#' Pivot and merged the files
#+ message=FALSE, warning=FALSE
effect <- EnvIndexes_effects %>% select(Chrm, Pos, contains("Effect")) %>%
  pivot_longer(
    cols = -c(Chrm, Pos),
    names_to = "Factors",
    names_prefix = "Effect_",
    values_to = "Effect")

pvalue <- EnvIndexes_effects %>% select(Chrm, Pos, contains("pvalue")) %>%
  pivot_longer(
    cols = -c(Chrm, Pos),
    names_to = "Factors",
    names_prefix = "_pvalue",
    values_to = "pvalue")
pvalue$Factors <- gsub("_pvalue", "", pvalue$Factors)

estimate <- EnvIndexes_effects %>% select(Chrm, Pos, contains("_R")) %>%
  select(!"Effect_Reaction") %>%
  pivot_longer(
    cols = -c(Chrm, Pos),
    names_to = "Factors",
    names_prefix = "_R",
    values_to = "R")
estimate$Factors <- gsub("_R", "", estimate$Factors)

EnvIndexes_effects2 <- left_join(effect, estimate , by = c("Chrm", "Pos", "Factors"))
EnvIndexes_effects2 <- left_join(EnvIndexes_effects2, pvalue , by = c("Chrm", "Pos", "Factors"))





