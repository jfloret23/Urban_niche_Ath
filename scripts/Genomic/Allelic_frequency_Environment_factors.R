#' ---
#'title: "Allelic frequency per factors limiting the environment"
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
Env_factors <- read_csv("/home/justine/sciebo/2_Ecological_parameters/2_assays/Ecological_parameters.csv")

Env_factors$Site_ID <- as.factor(Env_factors$Site_ID)
Env_factors$Type <- as.factor(Env_factors$Type)

Env_factors <- Env_factors %>%
  mutate(Type = recode(Type, "grass" = "Meadow",
                       "pavement" = "Pavement",
                       "Tree bed" = "Tree_bed", 
                       "wall"="Wall"))
summary(Env_factors)
colnames(Env_factors)

#' Keeping only the ecological variables associating with the establishment. 
#+ message=FALSE, warning=FALSE
factor_estblishment <- c("Site_ID", "Type", "Moss_abundance", "Shannon_index", "Sealed_Surface_09_21", "Vegetation_cover_09_21", 
                         "Light_intensity", "DD_04.22","Tmin_03.22", "res_Tmin", "Tmax_03.22", "res_Tmax", "Tmean_03.22", "res_Tmean",
                         "Hmean_03.22", "res_Hmean", "Root_Volume_estimate", "Sk_Vol", "pH", "Bulk_density", "PO4_estimate", 
                         "NO3_estimate", "SO4_estimate", "Cl_estimate", "F_estimate", "C.N.ratio")

Env_factors <- Env_factors[,(names(Env_factors) %in% factor_estblishment)] 
Env_factors <- Env_factors %>% drop_na()

Env_factors_standardized <- decostand(Env_factors[,c(3:26)], method = "standardize") # Scale and center variables
Env_factors_standardized <- cbind(Env_factors[,1:2], Env_factors_standardized)

#' ## Merged SNP data and the ecological factors
#+ message=FALSE, warning=FALSE
median_AF_Env_factors_standardized <- inner_join(Median_AF_50kb_noContamination, Env_factors_standardized, join_by("Sample_ID"== "Site_ID"))
median_AF_Env_factors_standardized <- median_AF_Env_factors_standardized %>% filter(SNP_tot>10)

summary(median_AF_Env_factors_standardized)
median_AF_Env_factors_standardized %>% distinct(Sample_ID) %>% count()
colnames(median_AF_Env_factors_standardized)

median_AF_Env_factors <- inner_join(Median_AF_50kb_noContamination, Env_factors, join_by("Sample_ID"== "Site_ID"))
median_AF_Env_factors <- median_AF_Env_factors %>% filter(SNP_tot>10)


#' # II. Association AF - Environmental indexes


#+ message=FALSE, warning=FALSE
for (chr in chrm_list) {
  data_chr <- median_AF_Env_factors_standardized %>% filter(Chrm == chr)
  positions <- unique(data_chr$Pos_mean)
  for (pos in positions) {
    data_pos <- data_chr %>% filter(Pos_mean == pos)
    if (nrow(data_pos) == 0) next
    corresults_Type <- cor.test(data_pos$median_AF, data_pos$Type,  method = "spearman")
    corresults_Moss <- cor.test(data_pos$median_AF, data_pos$Moss_abundance,  method = "spearman")
    corresults_Shannon <- cor.test(data_pos$median_AF, data_pos$Shannon_index,  method = "spearman")
    corresults_Light <- cor.test(data_pos$median_AF, data_pos$Light_intensity,  method = "spearman")
    corresults_DD <- cor.test(data_pos$median_AF, data_pos$DD_04.22,  method = "spearman")
    corresults_Tmin_03 <- cor.test(data_pos$median_AF, data_pos$Tmin_03.22,  method = "spearman")
    corresults_Tmax_03 <- cor.test(data_pos$median_AF, data_pos$Tmax_03.22,  method = "spearman")
    corresults_Tmean_03 <- cor.test(data_pos$median_AF, data_pos$Tmean_03.22,  method = "spearman")
    corresults_Hmean_03 <- cor.test(data_pos$median_AF, data_pos$Hmean_03.22,  method = "spearman")
    corresults_Root <- cor.test(data_pos$median_AF, data_pos$Root_Volume_estimate,  method = "spearman")
    corresults_Sk <- cor.test(data_pos$median_AF, data_pos$Sk_Vol,  method = "spearman")
    corresults_pH <- cor.test(data_pos$median_AF, data_pos$pH,  method = "spearman")
    corresults_Bulk <- cor.test(data_pos$median_AF, data_pos$Bulk_density,  method = "spearman")
    corresults_PO4 <- cor.test(data_pos$median_AF, data_pos$PO4_estimate,  method = "spearman")
    corresults_NO3 <- cor.test(data_pos$median_AF, data_pos$NO3_estimate,  method = "spearman")
    corresults_SO4 <- cor.test(data_pos$median_AF, data_pos$SO4_estimate,  method = "spearman")
    corresults_Cl <- cor.test(data_pos$median_AF, data_pos$Cl_estimate,  method = "spearman")
    corresults_F <- cor.test(data_pos$median_AF, data_pos$F_estimate,  method = "spearman")
    corresults_Sealed <- cor.test(data_pos$median_AF, data_pos$Sealed_Surface_09_21,  method = "spearman")
    corresults_Vegetation <- cor.test(data_pos$median_AF, data_pos$Vegetation_cover_09_21,  method = "spearman")
    corresults_Tmin_04 <- cor.test(data_pos$median_AF, data_pos$res_Tmin,  method = "spearman")
    corresults_Tmax_04 <- cor.test(data_pos$median_AF, data_pos$res_Tmax,  method = "spearman")
    corresults_Tmean_04 <- cor.test(data_pos$median_AF, data_pos$res_Tmean,  method = "spearman")
    corresults_Hmean_04 <- cor.test(data_pos$median_AF, data_pos$res_Hmean,  method = "spearman")
    corresults_C.N.ratio <- cor.test(data_pos$median_AF, data_pos$C.N.ratio,  method = "spearman")

    
    tmp <- c(chr,pos, 
             corresults_Type[["p.value"]], corresults_Type[["estimate"]], 
             corresults_Moss[["p.value"]],  corresults_Moss[["estimate"]],
             corresults_Shannon[["p.value"]],  corresults_Shannon[["estimate"]],
             corresults_Light[["p.value"]], corresults_Light[["estimate"]],
             corresults_DD[["p.value"]], corresults_DD[["estimate"]], 
             corresults_Tmin_03[["p.value"]], corresults_Tmin_03[["estimate"]], 
             corresults_Tmax_03[["p.value"]], corresults_Tmax_03[["estimate"]], 
             corresults_Tmean_03[["p.value"]], corresults_Tmean_03[["estimate"]], 
             corresults_Hmean_03[["p.value"]], corresults_Hmean_03[["estimate"]], 
             corresults_Root[["p.value"]], corresults_Root[["estimate"]], 
             corresults_Sk[["p.value"]], corresults_Sk[["estimate"]], 
             corresults_pH[["p.value"]], corresults_pH[["estimate"]], 
             corresults_Bulk[["p.value"]], corresults_Bulk[["estimate"]], 
             corresults_PO4[["p.value"]], corresults_PO4[["estimate"]], 
             corresults_NO3[["p.value"]], corresults_NO3[["estimate"]], 
             corresults_SO4[["p.value"]], corresults_SO4[["estimate"]], 
             corresults_Cl[["p.value"]], corresults_Cl[["estimate"]], 
             corresults_F[["p.value"]], corresults_F[["estimate"]], 
             corresults_Sealed[["p.value"]], corresults_Sealed[["estimate"]], 
             corresults_Vegetation[["p.value"]],  corresults_Vegetation[["estimate"]],
             corresults_Tmin_04[["p.value"]], corresults_Tmin_04[["estimate"]], 
             corresults_Tmax_04[["p.value"]], corresults_Tmax_04[["estimate"]], 
             corresults_Tmean_04[["p.value"]], corresults_Tmean_04[["estimate"]], 
             corresults_Hmean_04[["p.value"]], corresults_Hmean_04[["estimate"]], 
             corresults_C.N.ratio[["p.value"]], corresults_C.N.ratio[["estimate"]])
    
    Env_model <- rbind(Env_model, tmp)
  }
}

Env_model <- Env_model[-1,] 
Env_model <- Env_model %>% mutate(across(c(2:52), as.numeric)) 




#' # III. Randomization x 1000 

#' ## Model with randomization: 10 times
#+ eval=F, echo=T
Envfactors_rand_models <- vector("list", 1000)  

#+ eval=F, echo=T
for (i in 1:1000) {
  # 1) Randomize 
  pool_map <- median_AF_Env_factors_standardized %>%
    distinct(Sample_ID,.keep_all = TRUE) %>%
    mutate(Randomized_Type = sample(Type),
           Randomized_Moss = sample(Moss_abundance),
           Randomized_Shannon = sample(Shannon_index),
           Randomized_Light = sample(Light_intensity),
           Randomized_DD = sample(DD_04.22),
           Randomized_Tmin_03 = sample(Tmin_03.22),
           Randomized_Tmax_03 = sample(Tmax_03.22),
           Randomized_Tmean_03 = sample(Tmean_03.22),
           Randomized_Hmean_03 = sample(Hmean_03.22),
           Randomized_Root = sample(Root_Volume_estimate),
           Randomized_Sk = sample(Sk_Vol),
           Randomized_pH = sample(pH),
           Randomized_Bulk = sample(Bulk_density),
           Randomized_PO4 = sample(PO4_estimate),
           Randomized_NO3 = sample(NO3_estimate),
           Randomized_SO4 = sample(SO4_estimate),
           Randomized_Cl = sample(Cl_estimate),
           Randomized_F = sample(F_estimate),
           Randomized_Sealed = sample(Sealed_Surface_09_21),
           Randomized_Vegetation = sample(Vegetation_cover_09_21),
           Randomized_Tmin_04 = sample(res_Tmin),
           Randomized_Tmax_04 = sample(res_Tmax),
           Randomized_Tmean_04 = sample(res_Tmean),
           Randomized_Hmean_04 = sample(res_Hmean),
           Randomized_C.N = sample(C.N.ratio))
  data_random <- left_join(median_AF_Env_factors_standardized[,1:5], pool_map[, c(3, 31:55)], by = "Sample_ID")
  
  # 2) models for all Chr × Pos
  chrm_list <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
  Envfactors_rand_dfmodel <- tibble("Chrm" = "1", "Pos" = 1,
                                    Type_pvalue=1, Type_R=1, 
                                    Moss_pvalue = 1, Moss_R=1,
                                    Shannon_pvalue = 1, Shannon_R=1,
                                    Light_pvalue = 1,  Light_R=1,
                                    DD_pvalue = 1, DD_R=1,
                                    Tmin_03_pvalue = 1, Tmin_03_R=1,
                                    Tmax_03_pvalue = 1, Tmax_03_R=1,
                                    Tmean_03_pvalue = 1, Tmean_03_R=1,
                                    Hmean_03_pvalue = 1, Hmean_03_R=1,
                                    Root_pvalue = 1, Root_R=1,
                                    Sk_pvalue = 1, Sk_R=1,
                                    pH_pvalue = 1, pH_R=1,
                                    Bulk_pvalue = 1, Bulk_R=1,
                                    PO4_pvalue = 1, PO4_R=1,
                                    NO3_pvalue = 1, NO3_R=1,
                                    SO4_pvalue = 1, SO4_R=1,
                                    Cl_pvalue = 1, Cl_R=1,
                                    F_pvalue = 1, F_R=1,
                                    Sealed_pvalue = 1, Sealed_R=1, 
                                    Vege_cover_pvalue = 1,  Vege_cover_R=1,
                                    Tmin_04_pvalue = 1, Tmin_04_R=1,
                                    Tmax_04_pvalue = 1, Tmax_04_R=1,
                                    Tmean_04_pvalue = 1, Tmean_04_R=1,
                                    Hmean_04_pvalue = 1, Hmean_04_R=1,
                                    C.N.ratio_pvalue = 1, C.N.ratio_R=1)
  
  for (chr in chrm_list) {
    data_chr <- data_random %>% filter(Chrm == chr)
    positions <- unique(data_chr$Pos_mean)
    for (pos in positions) {
      data_pos <- data_chr %>% filter(Pos_mean == pos)
      if (nrow(data_pos) == 0) next
      corresults_Type <- cor.test(data_pos$median_AF, data_pos$Randomized_Type,  method = "spearman")
      corresults_Moss <- cor.test(data_pos$median_AF, data_pos$Randomized_Moss,  method = "spearman")
      corresults_Shannon <- cor.test(data_pos$median_AF, data_pos$Randomized_Shannon,  method = "spearman")
      corresults_Light <- cor.test(data_pos$median_AF, data_pos$Randomized_Light,  method = "spearman")
      corresults_DD <- cor.test(data_pos$median_AF, data_pos$Randomized_DD,  method = "spearman")
      corresults_Tmin_03 <- cor.test(data_pos$median_AF, data_pos$Randomized_Tmin_03,  method = "spearman")
      corresults_Tmax_03 <- cor.test(data_pos$median_AF, data_pos$Randomized_Tmax_03,  method = "spearman")
      corresults_Tmean_03 <- cor.test(data_pos$median_AF, data_pos$Randomized_Tmean_03,  method = "spearman")
      corresults_Hmean_03 <- cor.test(data_pos$median_AF, data_pos$Randomized_Hmean_03,  method = "spearman")
      corresults_Root <- cor.test(data_pos$median_AF, data_pos$Randomized_Root,  method = "spearman")
      corresults_Sk <- cor.test(data_pos$median_AF, data_pos$Randomized_Sk,  method = "spearman")
      corresults_pH <- cor.test(data_pos$median_AF, data_pos$Randomized_pH,  method = "spearman")
      corresults_Bulk <- cor.test(data_pos$median_AF, data_pos$Randomized_Bulk,  method = "spearman")
      corresults_PO4 <- cor.test(data_pos$median_AF, data_pos$Randomized_PO4,  method = "spearman")
      corresults_NO3 <- cor.test(data_pos$median_AF, data_pos$Randomized_NO3,  method = "spearman")
      corresults_SO4 <- cor.test(data_pos$median_AF, data_pos$Randomized_SO4,  method = "spearman")
      corresults_Cl <- cor.test(data_pos$median_AF, data_pos$Randomized_Cl,  method = "spearman")
      corresults_F <- cor.test(data_pos$median_AF, data_pos$Randomized_F,  method = "spearman")
      corresults_Sealed <- cor.test(data_pos$median_AF, data_pos$Randomized_Sealed,  method = "spearman")
      corresults_Vegetation <- cor.test(data_pos$median_AF, data_pos$Randomized_Vegetation,  method = "spearman")
      corresults_Tmin_04 <- cor.test(data_pos$median_AF, data_pos$Randomized_Tmin_04,  method = "spearman")
      corresults_Tmax_04 <- cor.test(data_pos$median_AF, data_pos$Randomized_Tmax_04,  method = "spearman")
      corresults_Tmean_04 <- cor.test(data_pos$median_AF, data_pos$Randomized_Tmean_04,  method = "spearman")
      corresults_Hmean_04 <- cor.test(data_pos$median_AF, data_pos$Randomized_Hmean_04,  method = "spearman")
      corresults_C.N.ratio <- cor.test(data_pos$median_AF, data_pos$Randomized_C.N, method = "spearman")
      
      
      
      tmp <- c(chr,pos, 
               corresults_Type[["p.value"]], corresults_Type[["estimate"]], 
               corresults_Moss[["p.value"]],  corresults_Moss[["estimate"]],
               corresults_Shannon[["p.value"]],  corresults_Shannon[["estimate"]],
               corresults_Light[["p.value"]], corresults_Light[["estimate"]],
               corresults_DD[["p.value"]], corresults_DD[["estimate"]], 
               corresults_Tmin_03[["p.value"]], corresults_Tmin_03[["estimate"]], 
               corresults_Tmax_03[["p.value"]], corresults_Tmax_03[["estimate"]], 
               corresults_Tmean_03[["p.value"]], corresults_Tmean_03[["estimate"]], 
               corresults_Hmean_03[["p.value"]], corresults_Hmean_03[["estimate"]], 
               corresults_Root[["p.value"]], corresults_Root[["estimate"]], 
               corresults_Sk[["p.value"]], corresults_Sk[["estimate"]], 
               corresults_pH[["p.value"]], corresults_pH[["estimate"]], 
               corresults_Bulk[["p.value"]], corresults_Bulk[["estimate"]], 
               corresults_PO4[["p.value"]], corresults_PO4[["estimate"]], 
               corresults_NO3[["p.value"]], corresults_NO3[["estimate"]], 
               corresults_SO4[["p.value"]], corresults_SO4[["estimate"]], 
               corresults_Cl[["p.value"]], corresults_Cl[["estimate"]], 
               corresults_F[["p.value"]], corresults_F[["estimate"]], 
               corresults_Sealed[["p.value"]], corresults_Sealed[["estimate"]], 
               corresults_Vegetation[["p.value"]],  corresults_Vegetation[["estimate"]],
               corresults_Tmin_04[["p.value"]], corresults_Tmin_04[["estimate"]], 
               corresults_Tmax_04[["p.value"]], corresults_Tmax_04[["estimate"]], 
               corresults_Tmean_04[["p.value"]], corresults_Tmean_04[["estimate"]], 
               corresults_Hmean_04[["p.value"]], corresults_Hmean_04[["estimate"]], 
               corresults_C.N.ratio[["p.value"]], corresults_C.N.ratio[["estimate"]])
      
      Envfactors_rand_dfmodel <- rbind(Envfactors_rand_dfmodel, tmp)
    }
  }
      
  # 3) Post-processing
  Envfactors_rand_dfmodel <- Envfactors_rand_dfmodel[-1, ]
  
  # save result of this i
  Envfactors_rand_models[[i]] <- Envfactors_rand_dfmodel
  print(paste("randomization", i, "done"))
}


#' Create one data frame
#+ eval=F, echo=T
#saveRDS(EnvIndexes_rand_models, "EnvIndexes_rand_models" )
Envfactors_rand_models_df <- purrr::map2_df(
  Envfactors_rand_models,
  1:1000,
  ~ dplyr::mutate(.x, iter = .y))


#' Merge the data set  
#+ eval=F, echo=T
Env_model$iter <- "Observed"
Envfactors_rand_models_df <- rbind(Envfactors_rand_models_df,Env_model)
Envfactors_rand_models_df <- Envfactors_rand_models_df %>% mutate(across(c(2:34), as.numeric)) 






#' # IV. Analysis observed vs randomized data

#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df <- Envfactors_rand_models_df %>% mutate("Data" = if_else(iter != "Observed", "Randomized", "Observed"))
Envfactors_model <- Envfactors_rand_models_df %>% filter(iter == "Observed")

colnames(Envfactors_rand_models_df)




#' # Type
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Type_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Type model estimate mod(Alt_frq ~ Type) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Type_R, probs = c(.95)),
                                    Q05 = quantile(Type_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,3,4)], by = c("Chrm", "Pos"))
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




#' # Moss abundance
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Moss_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Moss abundace model estimate mod(Alt_frq ~ Moss) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Moss_R, probs = c(.95)),
                                    Q05 = quantile(Moss_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,5,6)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Moss_R>Q95 & Moss_pvalue <0.05 | Moss_R<Q05 & Moss_pvalue <0.05, "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Moss = ifelse(Estimate_5 == "TRUE" & Moss_R >0, "Positif",
                                                           ifelse(Estimate_5 == "TRUE" & Moss_R <0, "Negatif", "Neutral")))
EnvFactors_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], by = c("Chrm", "Pos"))

#' ## Graphic p-value & estimate
#+ message=FALSE, warning=FALSE
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Moss_pvalue), color = Effect_Moss, alpha = Effect_Moss)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c( "red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Moss abundance - median AF (+10 SNP per 50kb windows)"))




#' # Shannon Index
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Shannon_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Shannon Index model estimate mod(Alt_frq ~ Shannon) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Shannon_R, probs = c(.95)),
                                    Q05 = quantile(Shannon_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,7,8)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Shannon_R>Q95 & Shannon_pvalue <0.05 | Shannon_R<Q05 & Shannon_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Shannon = ifelse(Estimate_5 == "TRUE" & Shannon_R >0, "Positif",
                                                              ifelse(Estimate_5 == "TRUE" & Shannon_R <0, "Negatif", "Neutral")))

EnvFactors_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], , by = c("Chrm", "Pos"))

#' ## Graphic p-value & estimate
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Shannon_pvalue), color = Effect_Shannon, alpha = Effect_Shannon)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c( "red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Shannon Index - median AF (+10 SNP per 50kb windows)"))





#' # Light Intensity
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Light_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Light intensity model estimate mod(Alt_frq ~ Light) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Light_R, probs = c(.95)),
                                    Q05 = quantile(Light_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,9,10)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Light_R>Q95 & Light_pvalue <0.05 | Light_R<Q05 & Light_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Light = ifelse(Estimate_5 == "TRUE" & Light_R >0, "Positif",
                                                            ifelse(Estimate_5 == "TRUE" & Light_R <0, "Negatif", "Neutral")))
EnvFactors_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], , by = c("Chrm", "Pos"))

#' ## Graphic p-value & estimate
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Light_pvalue), color = Effect_Light, alpha = Effect_Light)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c("red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Light intensity - median AF (+10 SNP per 50kb windows)"))





#' # Degree Day 
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(DD_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of DD model estimate mod(Alt_frq ~ DD) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(DD_R, probs = c(.95)),
                                    Q05 = quantile(DD_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,11,12)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(DD_R>Q95 & DD_pvalue <0.05 | DD_R<Q05 & DD_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_DD = ifelse(Estimate_5 == "TRUE" & DD_R >0, "Positif",
                                                           ifelse(Estimate_5 == "TRUE" & DD_R <0, "Negatif", "Neutral")))

EnvFactors_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], by = c("Chrm", "Pos"))


#' ## Graphic p-value & estimate
#+ message=FALSE, warning=FALSE
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(DD_pvalue), color = Effect_DD, alpha = Effect_DD)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c("lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association DD - median AF (+10 SNP per 50kb windows)"))


#' ## Correlation allelic frequency variable
Quantile_df %>% filter(Estimate_5 ==TRUE) %>% slice_min(DD_pvalue, n= 1)

DD_df <- median_AF_Env_factors %>% filter(Chrm == "Chr4" & Pos_mean == "13225000")
DD_df <- DD_df %>% mutate(across(c(2:30), as.numeric)) 

allelic_frq <- as.numeric(unlist(DD_df[DD_df$Chrm == "Chr4","median_AF"]))
variable <- as.numeric(unlist(DD_df[DD_df$Chrm == "Chr4","DD_04.22"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation

DD_df %>% 
  ggplot(aes(x=median_AF, y=DD_04.22)) +
  geom_point(size=4, color="#00CD00") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.35, y=39380, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("DD index ~ allelic frequency Chr4"))
#ggsave("DD index ~ allelic frequency Chr4.svg", plot =  last_plot(),dpi = 120, bg = "white", scale = 1)




#' # March min. temperature 
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Tmin_03_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Tmin_03 model estimate mod(Alt_frq ~ Tmin_03) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Tmin_03_R, probs = c(.95)),
                                    Q05 = quantile(Tmin_03_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,13,14)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Tmin_03_R>Q95 & Tmin_03_pvalue <0.05 | Tmin_03_R<Q05 & Tmin_03_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Tmin_03 = ifelse(Estimate_5 == "TRUE" & Tmin_03_R >0, "Positif",
                                                         ifelse(Estimate_5 == "TRUE" & Tmin_03_R <0, "Negatif", "Neutral")))

#' ## Graphic p-value & estimate
#+ message=FALSE, warning=FALSE
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Tmin_03_pvalue), color = Effect_Tmin_03, alpha = Effect_Tmin_03)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c("lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Tmin_03 - median AF (+10 SNP per 50kb windows)"))

EnvFactors_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], by = c("Chrm", "Pos"))






#' # March max. temperature 
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Tmax_03_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Tmax_03_R model estimate mod(Alt_frq ~ Tmax_03) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Tmax_03_R, probs = c(.95)),
                                    Q05 = quantile(Tmax_03_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,15,16)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Tmax_03_R>Q95 & Tmax_03_pvalue <0.05 | Tmax_03_R<Q05 & Tmax_03_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Tmax_03 = ifelse(Estimate_5 == "TRUE" & Tmax_03_R >0, "Positif",
                                                              ifelse(Estimate_5 == "TRUE" & Tmax_03_R <0, "Negatif", "Neutral")))

#' ## Graphic p-value & estimate
#+ message=FALSE, warning=FALSE
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Tmax_03_pvalue), color = Effect_Tmax_03, alpha = Effect_Tmax_03)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c("red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Tmax_03 - median AF (+10 SNP per 50kb windows)"))

EnvFactors_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], by = c("Chrm", "Pos"))





#' # March mean temperature 
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Tmean_03_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Tmean_03_R model estimate mod(Alt_frq ~ Tmean_03) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Tmean_03_R, probs = c(.95)),
                                    Q05 = quantile(Tmean_03_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,17,18)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Tmean_03_R>Q95 & Tmean_03_pvalue <0.05 | Tmean_03_R<Q05 & Tmean_03_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Tmean_03 = ifelse(Estimate_5 == "TRUE" & Tmean_03_R >0, "Positif",
                                                              ifelse(Estimate_5 == "TRUE" & Tmean_03_R <0, "Negatif", "Neutral")))

#' ## Graphic p-value & estimate
#+ message=FALSE, warning=FALSE
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Tmean_03_pvalue), color = Effect_Tmean_03, alpha = Effect_Tmean_03)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c("red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Tmean_03 - median AF (+10 SNP per 50kb windows)"))

EnvFactors_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], by = c("Chrm", "Pos"))






#' # March mean humidity 
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Hmean_03_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Hmean_03_R model estimate mod(Alt_frq ~ Hmean_03) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Hmean_03_R, probs = c(.95)),
                                    Q05 = quantile(Hmean_03_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,19,20)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Hmean_03_R>Q95 & Hmean_03_pvalue <0.05 | Hmean_03_R<Q05 & Hmean_03_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Hmean_03 = ifelse(Estimate_5 == "TRUE" & Hmean_03_R >0, "Positif",
                                                               ifelse(Estimate_5 == "TRUE" & Hmean_03_R <0, "Negatif", "Neutral")))

#' ## Graphic p-value & estimate
#+ message=FALSE, warning=FALSE
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Hmean_03_pvalue), color = Effect_Hmean_03, alpha = Effect_Hmean_03)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c("red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Hmean_03 - median AF (+10 SNP per 50kb windows)"))

EnvFactors_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], by = c("Chrm", "Pos"))





#' # Root volume
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Root_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Root volume model estimate mod(Alt_frq ~ Root) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Root_R, probs = c(.95)),
                                    Q05 = quantile(Root_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,21,22)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Root_R>Q95 & Root_pvalue <0.05 | Root_R<Q05 & Root_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Root = ifelse(Estimate_5 == "TRUE" & Root_R >0, "Positif",
                                                               ifelse(Estimate_5 == "TRUE" & Root_R <0, "Negatif", "Neutral")))

#' ## Graphic p-value & estimate
#+ message=FALSE, warning=FALSE
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Root_pvalue), color = Effect_Root, alpha = Effect_Root)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c("red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Root Volume - median AF (+10 SNP per 50kb windows)"))

EnvFactors_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], by = c("Chrm", "Pos"))





#' # Skeleton volume
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Sk_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Skeleton volume model estimate mod(Alt_frq ~ Skeleton) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Sk_R, probs = c(.95)),
                                    Q05 = quantile(Sk_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,23,24)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Sk_R>Q95 & Sk_pvalue <0.05 | Sk_R<Q05 & Sk_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Sk = ifelse(Estimate_5 == "TRUE" & Sk_R >0, "Positif",
                                                           ifelse(Estimate_5 == "TRUE" & Sk_R <0, "Negatif", "Neutral")))

#' ## Graphic p-value & estimate
#+ message=FALSE, warning=FALSE
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Sk_pvalue), color = Effect_Sk, alpha = Effect_Sk)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c("red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Sk Volume - median AF (+10 SNP per 50kb windows)"))

EnvFactors_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], by = c("Chrm", "Pos"))







#' # pH
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(pH_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of pH model estimate mod(Alt_frq ~ pH) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(pH_R, probs = c(.95)),
                                    Q05 = quantile(pH_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,25,26)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(pH_R>Q95 & pH_pvalue <0.05 | pH_R<Q05 & pH_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_pH = ifelse(Estimate_5 == "TRUE" & pH_R >0, "Positif",
                                                         ifelse(Estimate_5 == "TRUE" & pH_R <0, "Negatif", "Neutral")))

#' ## Graphic p-value & estimate
#+ message=FALSE, warning=FALSE
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(pH_pvalue), color = Effect_pH, alpha = Effect_pH)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c("lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association pH - median AF (+10 SNP per 50kb windows)"))

EnvFactors_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], by = c("Chrm", "Pos"))






#' # Bulk density
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Bulk_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Bulk density model estimate mod(Alt_frq ~ Bulk density) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Bulk_R, probs = c(.95)),
                                    Q05 = quantile(Bulk_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,27,28)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Bulk_R>Q95 & Bulk_pvalue <0.05 | Bulk_R<Q05 & Bulk_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Bulk = ifelse(Estimate_5 == "TRUE" & Bulk_R >0, "Positif",
                                                         ifelse(Estimate_5 == "TRUE" & Bulk_R <0, "Negatif", "Neutral")))

#' ## Graphic p-value & estimate
#+ message=FALSE, warning=FALSE
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Bulk_pvalue), color = Effect_Bulk, alpha = Effect_Bulk)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c("red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Bulk density - median AF (+10 SNP per 50kb windows)"))

EnvFactors_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], by = c("Chrm", "Pos"))






#' # PO4
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(PO4_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of PO4 model estimate mod(Alt_frq ~ PO4) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(PO4_R, probs = c(.95)),
                                    Q05 = quantile(PO4_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,29,30)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(PO4_R>Q95 & PO4_pvalue <0.05 | PO4_R<Q05 & PO4_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_PO4 = ifelse(Estimate_5 == "TRUE" & PO4_R >0, "Positif",
                                                           ifelse(Estimate_5 == "TRUE" & PO4_R <0, "Negatif", "Neutral")))

#' ## Graphic p-value & estimate
#+ message=FALSE, warning=FALSE
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(PO4_pvalue), color = Effect_PO4, alpha = Effect_PO4)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c("lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association PO4 - median AF (+10 SNP per 50kb windows)"))

EnvFactors_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], by = c("Chrm", "Pos"))


#' ## Correlation allelic frequency variable
Quantile_df %>% filter(Estimate_5 ==TRUE) %>% slice_min(PO4_pvalue, n= 1)

PO4_df <- median_AF_Env_factors %>% filter(Chrm == "Chr3" & Pos_mean == "8125000"|
                                            Chrm == "Chr4" & Pos_mean == "13325000")
PO4_df <- PO4_df %>% mutate(across(c(2:30), as.numeric)) 

#' Chr3
allelic_frq <- as.numeric(unlist(PO4_df[PO4_df$Chrm == "Chr3","median_AF"]))
variable <- as.numeric(unlist(PO4_df[PO4_df$Chrm == "Chr3","PO4_estimate"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation

PO4_df %>% filter(Chrm == "Chr3") %>%
  ggplot(aes(x=median_AF, y=PO4_estimate)) +
  geom_point(size=4, color="#1C86EE") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.35, y=1300, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("PO4 index ~ allelic frequency Chr3"))
#ggsave("PO4 index ~ allelic frequency Chr3.svg", plot =  last_plot(),dpi = 120, bg = "white", scale = 1)


#' Chr4
allelic_frq <- as.numeric(unlist(PO4_df[PO4_df$Chrm == "Chr4","median_AF"]))
variable <- as.numeric(unlist(PO4_df[PO4_df$Chrm == "Chr4","PO4_estimate"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation

PO4_df %>% filter(Chrm == "Chr4") %>%
  ggplot(aes(x=median_AF, y=PO4_estimate)) +
  geom_point(size=4, color="#1C86EE") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.35, y=1300, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("PO4 index ~ allelic frequency Chr4"))
#ggsave("PO4 index ~ allelic frequency Chr4.svg", plot =  last_plot(),dpi = 120, bg = "white", scale = 1)




#' # NO3
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(NO3_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of NO3 model estimate mod(Alt_frq ~ NO3) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(NO3_R, probs = c(.95)),
                                    Q05 = quantile(NO3_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,31,32)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(NO3_R>Q95 & NO3_pvalue <0.05 | NO3_R<Q05 & NO3_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_NO3 = ifelse(Estimate_5 == "TRUE" & NO3_R >0, "Positif",
                                                          ifelse(Estimate_5 == "TRUE" & NO3_R <0, "Negatif", "Neutral")))

#' ## Graphic p-value & estimate
#+ message=FALSE, warning=FALSE
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(NO3_pvalue), color = Effect_NO3, alpha = Effect_NO3)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c("red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association NO3 - median AF (+10 SNP per 50kb windows)"))

EnvFactors_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], by = c("Chrm", "Pos"))


#' ## Correlation allelic frequency variable
Quantile_df %>% filter(Estimate_5 ==TRUE) %>% slice_min(NO3_pvalue, n= 1)

NO3_df <- median_AF_Env_factors %>% filter(Chrm == "Chr1" & Pos_mean == "7725000"|
                                             Chrm == "Chr3" & Pos_mean == "525000"|
                                             Chrm == "Chr4" & Pos_mean == "18075000")
NO3_df <- NO3_df %>% mutate(across(c(2:30), as.numeric)) 

#' Chr1
allelic_frq <- as.numeric(unlist(NO3_df[NO3_df$Chrm == "Chr1","median_AF"]))
variable <- as.numeric(unlist(NO3_df[NO3_df$Chrm == "Chr1","NO3_estimate"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation

NO3_df %>% filter(Chrm == "Chr1") %>%
  ggplot(aes(x=median_AF, y=NO3_estimate)) +
  geom_point(size=4, color="#00CED1") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.65, y=5000, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("NO3 index ~ allelic frequency Chr1"))
#ggsave("NO3 index ~ allelic frequency Chr1.svg", plot =  last_plot(),dpi = 120, bg = "white", scale = 1)

#' Chr3
allelic_frq <- as.numeric(unlist(NO3_df[NO3_df$Chrm == "Chr3","median_AF"]))
variable <- as.numeric(unlist(NO3_df[NO3_df$Chrm == "Chr3","NO3_estimate"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation

NO3_df %>% filter(Chrm == "Chr3") %>%
  ggplot(aes(x=median_AF, y=NO3_estimate)) +
  geom_point(size=4, color="#00CED1") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.65, y=5000, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("NO3 index ~ allelic frequency Chr3"))
#gsave("NO3 index ~ allelic frequency Chr3.svg", plot =  last_plot(),dpi = 120, bg = "white", scale = 1)



#' Chr4
allelic_frq <- as.numeric(unlist(NO3_df[NO3_df$Chrm == "Chr4","median_AF"]))
variable <- as.numeric(unlist(NO3_df[NO3_df$Chrm == "Chr4","NO3_estimate"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation

NO3_df %>% filter(Chrm == "Chr4") %>%
  ggplot(aes(x=median_AF, y=NO3_estimate)) +
  geom_point(size=4, color="#00CED1") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.35, y=3000, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("NO3 index ~ allelic frequency Chr4"))
#ggsave("NO3 index ~ allelic frequency Chr4.svg", plot =  last_plot(),dpi = 120, bg = "white", scale = 1)






#' # SO4
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(SO4_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of SO4 model estimate mod(Alt_frq ~ SO4) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(SO4_R, probs = c(.95)),
                                    Q05 = quantile(SO4_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,33,34)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(SO4_R>Q95 & SO4_pvalue <0.05 | SO4_R<Q05 & SO4_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_SO4 = ifelse(Estimate_5 == "TRUE" & SO4_R >0, "Positif",
                                                          ifelse(Estimate_5 == "TRUE" & SO4_R <0, "Negatif", "Neutral")))

#' ## Graphic p-value & estimate
#+ message=FALSE, warning=FALSE
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(SO4_pvalue), color = Effect_SO4, alpha = Effect_SO4)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c("red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association SO4 - median AF (+10 SNP per 50kb windows)"))

EnvFactors_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], by = c("Chrm", "Pos"))


#' ## Correlation allelic frequency variable
Quantile_df %>% filter(Estimate_5 ==TRUE) %>% slice_min(SO4_pvalue, n= 1)

SO4_df <- median_AF_Env_factors %>% filter(Chrm == "Chr2" & Pos_mean == "14875000"|
                                             Chrm == "Chr4" & Pos_mean == "13225000")
SO4_df <- SO4_df %>% mutate(across(c(2:30), as.numeric)) 

#' Chr2
allelic_frq <- as.numeric(unlist(SO4_df[SO4_df$Chrm == "Chr2","median_AF"]))
variable <- as.numeric(unlist(SO4_df[SO4_df$Chrm == "Chr2","SO4_estimate"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation

SO4_df %>% filter(Chrm == "Chr2") %>%
  ggplot(aes(x=median_AF, y=SO4_estimate)) +
  geom_point(size=4, color="#EEC900") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.25, y=700, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("SO4 index ~ allelic frequency Chr2"))
#ggsave("SO4 index ~ allelic frequency Chr2.svg", plot =  last_plot(),dpi = 120, bg = "white", scale = 1)


#' Chr4
allelic_frq <- as.numeric(unlist(SO4_df[SO4_df$Chrm == "Chr4","median_AF"]))
variable <- as.numeric(unlist(SO4_df[SO4_df$Chrm == "Chr4","SO4_estimate"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation

SO4_df %>% filter(Chrm == "Chr4") %>%
  ggplot(aes(x=median_AF, y=SO4_estimate)) +
  geom_point(size=4, color="#EEC900") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.45, y=700, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("SO4 index ~ allelic frequency Chr4"))
#ggsave("SO4 index ~ allelic frequency Chr4.svg", plot =  last_plot(),dpi = 120, bg = "white", scale = 1)





#' # Cl
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Cl_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Cl model estimate mod(Alt_frq ~ Cl) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Cl_R, probs = c(.95)),
                                    Q05 = quantile(Cl_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,35,36)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Cl_R>Q95 & Cl_pvalue <0.05 | Cl_R<Q05 & Cl_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Cl = ifelse(Estimate_5 == "TRUE" & Cl_R >0, "Positif",
                                                          ifelse(Estimate_5 == "TRUE" & Cl_R <0, "Negatif", "Neutral")))

#' ## Graphic p-value & estimate
#+ message=FALSE, warning=FALSE
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Cl_pvalue), color = Effect_Cl, alpha = Effect_Cl)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c("red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Cl - median AF (+10 SNP per 50kb windows)"))

EnvFactors_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], by = c("Chrm", "Pos"))


#' ## Correlation allelic frequency variable
Quantile_df %>% filter(Estimate_5 ==TRUE) %>% slice_min(Cl_pvalue, n= 1)

Cl_df <- median_AF_Env_factors %>% filter(Chrm == "Chr2" & Pos_mean == "14875000")
Cl_df <- Cl_df %>% mutate(across(c(2:30), as.numeric)) 

#' Chr2
allelic_frq <- as.numeric(unlist(Cl_df[Cl_df$Chrm == "Chr2","median_AF"]))
variable <- as.numeric(unlist(Cl_df[Cl_df$Chrm == "Chr2","Cl_estimate"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation

Cl_df %>% filter(Chrm == "Chr2") %>%
  ggplot(aes(x=median_AF, y=Cl_estimate)) +
  geom_point(size=4, color="#AB82FF") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.25, y=3500, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("Cl index ~ allelic frequency Chr2"))
#ggsave("Cl index ~ allelic frequency Chr2.svg", plot =  last_plot(),dpi = 120, bg = "white", scale = 1)







#' # F
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(F_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of F model estimate mod(Alt_frq ~ F) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(F_R, probs = c(.95)),
                                    Q05 = quantile(F_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,37,38)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(F_R>Q95 & F_pvalue <0.05 | F_R<Q05 & F_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_F = ifelse(Estimate_5 == "TRUE" & F_R >0, "Positif",
                                                         ifelse(Estimate_5 == "TRUE" & F_R <0, "Negatif", "Neutral")))

#' ## Graphic p-value & estimate
#+ message=FALSE, warning=FALSE
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(F_pvalue), color = Effect_F, alpha = Effect_F)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c("red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association F - median AF (+10 SNP per 50kb windows)"))

EnvFactors_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], by = c("Chrm", "Pos"))



#' ## Correlation allelic frequency variable
Quantile_df %>% filter(Estimate_5 ==TRUE) %>% slice_min(F_pvalue, n= 1)

F_df <- median_AF_Env_factors %>% filter(Chrm == "Chr2" & Pos_mean == "19625000")
F_df <- F_df %>% mutate(across(c(2:30), as.numeric)) 

#' Chr2
allelic_frq <- as.numeric(unlist(F_df[F_df$Chrm == "Chr2","median_AF"]))
variable <- as.numeric(unlist(F_df[F_df$Chrm == "Chr2","F_estimate"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation

F_df %>% filter(Chrm == "Chr2") %>%
  ggplot(aes(x=median_AF, y=F_estimate)) +
  geom_point(size=4, color="#4A708B") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.25, y=3500, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("F index ~ allelic frequency Chr2"))
#ggsave("F index ~ allelic frequency Chr2.svg", plot =  last_plot(),dpi = 120, bg = "white", scale = 1)






#' # Sealed surface  
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Sealed_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Sealed surface model estimate mod(Alt_frq ~ Sealed) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Sealed_R, probs = c(.95)),
                                    Q05 = quantile(Sealed_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,39,40)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Sealed_R>Q95 & Sealed_pvalue <0.05 | Sealed_R<Q05 & Sealed_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Sealed = ifelse(Estimate_5 == "TRUE" & Sealed_R >0, "Positif",
                                                             ifelse(Estimate_5 == "TRUE" & Sealed_R <0, "Negatif", "Neutral")))

#' ## Graphic p-value & estimate
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Sealed_pvalue), color = Effect_Sealed, alpha = Effect_Sealed)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c( "red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Sealed surface - median AF (+10 SNP per 50kb windows)"))

EnvFactors_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], , by = c("Chrm", "Pos"))





#' # Vegetation cover
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Vege_cover_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Vegetation cover model estimate mod(Alt_frq ~ Vegetation) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Vege_cover_R, probs = c(.95)),
                                    Q05 = quantile(Vege_cover_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,41,42)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Vege_cover_R>Q95 & Vege_cover_pvalue <0.05 | Vege_cover_R<Q05 & Vege_cover_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Vege_cover = ifelse(Estimate_5 == "TRUE" & Vege_cover_R >0, "Positif",
                                                                 ifelse(Estimate_5 == "TRUE" & Vege_cover_R <0, "Negatif", "Neutral")))

#' ## Graphic p-value & estimate
#+ message=FALSE, warning=FALSE
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Vege_cover_pvalue), color = Effect_Vege_cover, alpha = Effect_Vege_cover)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c( "red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Vegetation cover - median AF (+10 SNP per 50kb windows)"))

EnvFactors_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], , by = c("Chrm", "Pos"))






#' # April min. temperature 
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Tmin_04_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Tmin_04 model estimate mod(Alt_frq ~ Tmin_04) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Tmin_04_R, probs = c(.95)),
                                    Q05 = quantile(Tmin_04_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,43,44)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Tmin_04_R>Q95 & Tmin_04_pvalue <0.05 | Tmin_04_R<Q05 & Tmin_04_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Tmin_04 = ifelse(Estimate_5 == "TRUE" & Tmin_04_R >0, "Positif",
                                                              ifelse(Estimate_5 == "TRUE" & Tmin_04_R <0, "Negatif", "Neutral")))

#' ## Graphic p-value & estimate
#+ message=FALSE, warning=FALSE
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Tmin_04_pvalue), color = Effect_Tmin_04, alpha = Effect_Tmin_04)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c( "red" , "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Tmin_04 - median AF (+10 SNP per 50kb windows)"))

EnvFactors_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], by = c("Chrm", "Pos"))






#' # April max. temperature 
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Tmax_04_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Tmax_04_R model estimate mod(Alt_frq ~ Tmax_04) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Tmax_04_R, probs = c(.95)),
                                    Q05 = quantile(Tmax_04_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,45,46)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Tmax_04_R>Q95 & Tmax_04_pvalue <0.05 | Tmax_04_R<Q05 & Tmax_04_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Tmax_04 = ifelse(Estimate_5 == "TRUE" & Tmax_04_R >0, "Positif",
                                                              ifelse(Estimate_5 == "TRUE" & Tmax_04_R <0, "Negatif", "Neutral")))

#' ## Graphic p-value & estimate
#+ message=FALSE, warning=FALSE
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Tmax_04_pvalue), color = Effect_Tmax_04, alpha = Effect_Tmax_04)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c("red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Tmax_04 - median AF (+10 SNP per 50kb windows)"))

EnvFactors_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], by = c("Chrm", "Pos"))



#' ## Correlation allelic frequency variable
Quantile_df %>% filter(Estimate_5 ==TRUE) %>% slice_min(Tmax_04_pvalue, n= 1)

Tmax_04 <- median_AF_Env_factors %>% filter(Chrm == "Chr3" & Pos_mean == "3625000"|
                                              Chrm == "Chr4" & Pos_mean == "8725000")
Tmax_04 <- Tmax_04 %>% mutate(across(c(2:30), as.numeric)) 

#' Chr3
allelic_frq <- as.numeric(unlist(Tmax_04[Tmax_04$Chrm == "Chr3","median_AF"]))
variable <- as.numeric(unlist(Tmax_04[Tmax_04$Chrm == "Chr3","res_Tmax"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation

Tmax_04 %>% filter(Chrm == "Chr3") %>%
  ggplot(aes(x=median_AF, y=res_Tmax)) +
  geom_point(size=4, color="#EE6363") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.25, y=-1.5, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("Tmax_04 index ~ allelic frequency Chr3"))
#ggsave("Tmax_04 index ~ allelic frequency Chr3.svg", plot =  last_plot(),dpi = 120, bg = "white", scale = 1)


#' Chr4
allelic_frq <- as.numeric(unlist(Tmax_04[Tmax_04$Chrm == "Chr4","median_AF"]))
variable <- as.numeric(unlist(Tmax_04[Tmax_04$Chrm == "Chr4","res_Tmax"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation

Tmax_04 %>% filter(Chrm == "Chr4") %>%
  ggplot(aes(x=median_AF, y=res_Tmax)) +
  geom_point(size=4, color="#EE6363") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.55, y=-1.5, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("Tmax_04 index ~ allelic frequency Chr4"))
#ggsave("Tmax_04 index ~ allelic frequency Chr4.svg", plot =  last_plot(),dpi = 120, bg = "white", scale = 1)



#' ## Correlation allelic frequency variable
Quantile_df %>% filter(Estimate_5 ==TRUE) %>% slice_min(Tmax_04_pvalue, n= 1)

Tmax_04 <- median_AF_Env_factors %>% filter(Chrm == "Chr3" & Pos_mean == "3625000"|
                                              Chrm == "Chr4" & Pos_mean == "8725000")
Tmax_04 <- Tmax_04 %>% mutate(across(c(2:30), as.numeric)) 

#' Chr3
allelic_frq <- as.numeric(unlist(Tmax_04[Tmax_04$Chrm == "Chr3","median_AF"]))
variable <- as.numeric(unlist(Tmax_04[Tmax_04$Chrm == "Chr3","res_Tmax"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation

Tmax_04 %>% filter(Chrm == "Chr3") %>%
  ggplot(aes(x=median_AF, y=res_Tmax)) +
  geom_point(size=4, color="#EE6363") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.25, y=-1.5, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("Tmax_04 index ~ allelic frequency Chr3"))
#ggsave("Tmax_04 index ~ allelic frequency Chr3.svg", plot =  last_plot(),dpi = 120, bg = "white", scale = 1)


#' Chr4
allelic_frq <- as.numeric(unlist(Tmax_04[Tmax_04$Chrm == "Chr4","median_AF"]))
variable <- as.numeric(unlist(Tmax_04[Tmax_04$Chrm == "Chr4","res_Tmax"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation

Tmax_04 %>% filter(Chrm == "Chr4") %>%
  ggplot(aes(x=median_AF, y=res_Tmax)) +
  geom_point(size=4, color="#EE6363") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.55, y=-1.5, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("Tmax_04 index ~ allelic frequency Chr4"))
#ggsave("Tmax_04 index ~ allelic frequency Chr4.svg", plot =  last_plot(),dpi = 120, bg = "white", scale = 1)




#' # April mean temperature 
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Tmean_04_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Tmean_04_R model estimate mod(Alt_frq ~ Tmean_04) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Tmean_04_R, probs = c(.95)),
                                    Q05 = quantile(Tmean_04_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,47,48)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Tmean_04_R>Q95 & Tmean_04_pvalue <0.05 | Tmean_04_R<Q05 & Tmean_04_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Tmean_04 = ifelse(Estimate_5 == "TRUE" & Tmean_04_R >0, "Positif",
                                                               ifelse(Estimate_5 == "TRUE" & Tmean_04_R <0, "Negatif", "Neutral")))

#' ## Graphic p-value & estimate
#+ message=FALSE, warning=FALSE
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Tmean_04_pvalue), color = Effect_Tmean_04, alpha = Effect_Tmean_04)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c("red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Tmean_04 - median AF (+10 SNP per 50kb windows)"))

EnvFactors_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], by = c("Chrm", "Pos"))


#' ## Correlation allelic frequency variable
Quantile_df %>% filter(Estimate_5 ==TRUE) %>% slice_min(Tmean_04_pvalue, n= 1)

Tmean_04 <- median_AF_Env_factors %>% filter(Chrm == "Chr4" & Pos_mean == "18475000")
Tmean_04 <- Tmean_04 %>% mutate(across(c(2:30), as.numeric)) 

#' Chr4
allelic_frq <- as.numeric(unlist(Tmean_04[Tmean_04$Chrm == "Chr4","median_AF"]))
variable <- as.numeric(unlist(Tmean_04[Tmean_04$Chrm == "Chr4","res_Tmean"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation

Tmean_04 %>% filter(Chrm == "Chr4") %>%
  ggplot(aes(x=median_AF, y=res_Tmean)) +
  geom_point(size=4, color="#FFAEB9") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.25, y=-1.5, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("Tmean_04 index ~ allelic frequency Chr4"))
#ggsave("Tmean_04 index ~ allelic frequency Chr4.svg", plot =  last_plot(),dpi = 120, bg = "white", scale = 1)




#' # April mean humidity 
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(Hmean_04_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of Hmean_04_R model estimate mod(Alt_frq ~ Hmean_04) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(Hmean_04_R, probs = c(.95)),
                                    Q05 = quantile(Hmean_04_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,49,50)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(Hmean_04_R>Q95 & Hmean_04_pvalue <0.05 | Hmean_04_R<Q05 & Hmean_04_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_Hmean_04 = ifelse(Estimate_5 == "TRUE" & Hmean_04_R >0, "Positif",
                                                               ifelse(Estimate_5 == "TRUE" & Hmean_04_R <0, "Negatif", "Neutral")))

#' ## Graphic p-value & estimate
#+ message=FALSE, warning=FALSE
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(Hmean_04_pvalue), color = Effect_Hmean_04, alpha = Effect_Hmean_04)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c("red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association Hmean_04 - median AF (+10 SNP per 50kb windows)"))

EnvFactors_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], by = c("Chrm", "Pos"))


#' ## Correlation allelic frequency variable
Quantile_df %>% filter(Estimate_5 ==TRUE) %>% slice_min(Hmean_04_pvalue, n= 1)







#' # C/N ratio
#--------------------------------------------------------------------------------------

#' ## Estimates
#+ message=FALSE, warning=FALSE
Envfactors_rand_models_df %>%
  ggplot(aes(x=Pos, y=abs(C.N.ratio_R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of C.N.ratio_R model estimate mod(Alt_frq ~ C.N.ratio) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))

#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Envfactors_rand_models_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(C.N.ratio_R, probs = c(.95)),
                                    Q05 = quantile(C.N.ratio_R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join( Quantile_df, Envfactors_model[,c(1,2,51,52)], by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(C.N.ratio_R>Q95 & C.N.ratio_pvalue <0.05 | C.N.ratio_R<Q05 & C.N.ratio_pvalue <0.05 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect_C.N.ratio = ifelse(Estimate_5 == "TRUE" & C.N.ratio_R >0, "Positif",
                                                               ifelse(Estimate_5 == "TRUE" & C.N.ratio_R <0, "Negatif", "Neutral")))

#' ## Graphic p-value & estimate
#+ message=FALSE, warning=FALSE
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(C.N.ratio_pvalue), color = Effect_C.N.ratio, alpha = Effect_C.N.ratio)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c("red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Association C.N.ratio - median AF (+10 SNP per 50kb windows)"))

EnvFactors_effects <- left_join(EnvFactors_effects, Quantile_df[,c(1,2,5,6,8)], by = c("Chrm", "Pos"))



#' ## Correlation allelic frequency variable
Quantile_df %>% filter(Estimate_5 ==TRUE) %>% slice_min(C.N.ratio_pvalue, n= 1)

C.N_04 <- median_AF_Env_factors %>% filter(Chrm == "Chr1" & Pos_mean == "825000"|
                                               Chrm == "Chr3" & Pos_mean == "75000"|
                                             Chrm == "Chr4" & Pos_mean == "10225000"|
                                             Chrm == "Chr5" & Pos_mean == "25000")
C.N_04 <- C.N_04 %>% mutate(across(c(2:30), as.numeric)) 

#' Chr1
allelic_frq <- as.numeric(unlist(C.N_04[C.N_04$Chrm == "Chr1","median_AF"]))
variable <- as.numeric(unlist(C.N_04[C.N_04$Chrm == "Chr1","C.N.ratio"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation

C.N_04 %>% filter(Chrm == "Chr1") %>%
  ggplot(aes(x=median_AF, y=C.N.ratio)) +
  geom_point(size=4, color="#FF7F00") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.25, y=7, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("C.N_04 index ~ allelic frequency Chr1"))
#ggsave("C.N_04 index ~ allelic frequency Chr1.svg", plot =  last_plot(),dpi = 120, bg = "white", scale = 1)


#' Chr3
allelic_frq <- as.numeric(unlist(C.N_04[C.N_04$Chrm == "Chr3","median_AF"]))
variable <- as.numeric(unlist(C.N_04[C.N_04$Chrm == "Chr3","C.N.ratio"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation

C.N_04 %>% filter(Chrm == "Chr3") %>%
  ggplot(aes(x=median_AF, y=C.N.ratio)) +
  geom_point(size=4, color="#FF7F00") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.25, y=7, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("C.N_04 index ~ allelic frequency Chr3"))
#ggsave("C.N_04 index ~ allelic frequency Chr3.svg", plot =  last_plot(),dpi = 120, bg = "white", scale = 1)


#' Chr4
allelic_frq <- as.numeric(unlist(C.N_04[C.N_04$Chrm == "Chr4","median_AF"]))
variable <- as.numeric(unlist(C.N_04[C.N_04$Chrm == "Chr4","C.N.ratio"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation

C.N_04 %>% filter(Chrm == "Chr4") %>%
  ggplot(aes(x=median_AF, y=C.N.ratio)) +
  geom_point(size=4, color="#FF7F00") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.35, y=20, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("C.N_04 index ~ allelic frequency Chr4"))
#ggsave("C.N_04 index ~ allelic frequency Chr4.svg", plot =  last_plot(),dpi = 120, bg = "white", scale = 1)


#' Chr5
allelic_frq <- as.numeric(unlist(C.N_04[C.N_04$Chrm == "Chr5","median_AF"]))
variable <- as.numeric(unlist(C.N_04[C.N_04$Chrm == "Chr5","C.N.ratio"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
correlation

C.N_04 %>% filter(Chrm == "Chr5") %>%
  ggplot(aes(x=median_AF, y=C.N.ratio)) +
  geom_point(size=4, color="#FF7F00") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.35, y=20, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("C.N_04 index ~ allelic frequency Chr5"))
#ggsave("C.N_04 index ~ allelic frequency Chr5.svg", plot =  last_plot(),dpi = 120, bg = "white", scale = 1)




#' # Combined effect of the factor

#' Pivot and merged the files
#+ message=FALSE, warning=FALSE
effect <- EnvFactors_effects %>% select(Chrm, Pos, contains("Effect")) %>%
  pivot_longer(
    cols = -c(Chrm, Pos),
    names_to = "Factors",
    names_prefix = "Effect_",
    values_to = "Effect")

pvalue <- EnvFactors_effects %>% select(Chrm, Pos, contains("pvalue")) %>%
  pivot_longer(
    cols = -c(Chrm, Pos),
    names_to = "Factors",
    names_prefix = "_pvalue",
    values_to = "pvalue")
pvalue$Factors <- gsub("_pvalue", "", pvalue$Factors)

estimate <- EnvFactors_effects %>% select(Chrm, Pos, contains("_R")) %>%
  select(!"Effect_Root") %>%
  pivot_longer(
    cols = -c(Chrm, Pos),
    names_to = "Factors",
    names_prefix = "_R",
    values_to = "R")
estimate$Factors <- gsub("_R", "", estimate$Factors)

EnvFactors_effects2 <- left_join(effect, estimate , by = c("Chrm", "Pos", "Factors"))
EnvFactors_effects2 <- left_join(EnvFactors_effects2, pvalue , by = c("Chrm", "Pos", "Factors"))
#write_csv(EnvFactors_effects2, "EnvFactors_effects_quantile.csv")




