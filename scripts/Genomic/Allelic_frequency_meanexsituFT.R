#' ---
#'title: "Allelic frequency & mean ex situ FT"
#'author: "Justine Floret"
#'date: "2026-01-07"
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
# Stat. libray
library(lmerTest) #glht
library(lsmeans) #pairwise comparison
library(emmeans)
library(glmmTMB) # negative binomial model
library(vegan) #RDA
library(missMDA) #Missing data imputation, pca
library(ggpubr) #correlation line in plot



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


#' ## Ex situ flowering time
#+ message=FALSE, warning=FALSE
data_FT_F6 <- read_csv("/home/justine/sciebo/3_Exsitu_phenotype/Flowering_time/assay/F6_tplant_FT_2025.csv")
data_FT_F6 <- data_FT_F6 %>% select(ID_F5_parent, Pop, Individual, Tray, Column, Row, FT, Comments, brought_to_greenhouse )
data_FT_F6$ID_F5_parent <- as.factor(data_FT_F6$ID_F5_parent)
data_FT_F6 <- data_FT_F6 %>% dplyr::mutate(Comments= replace_na(Comments, "Alive"))
data_FT_F6$Comments <- as.factor(data_FT_F6$Comments)
data_FT_F6$Tray <- as.factor(data_FT_F6$Tray)
data_FT_F6 <- data_FT_F6 %>% filter(Comments == "Alive")

#' Ensure date columns are Date type
data_FT_F6$brought_to_greenhouse <- as.Date(data_FT_F6$brought_to_greenhouse, format = "%Y/%m/%d")
data_FT_F6$FT <- as.Date(data_FT_F6$FT , format = "%m/%d/%y")

#' Create FT column: number of days between the sowing and the flowering day, if the plant did not flower during the experiment a late flowering time is put (120 days)
data_FT_F6$FT_day <- as.integer(data_FT_F6$FT - data_FT_F6$brought_to_greenhouse)
data_FT_F6$FT_day <- as.numeric(data_FT_F6$FT_day)


#' Keep the the mean FT per population
#+ message=FALSE, warning=FALSE
data_FT_F6 <- data_FT_F6 %>% group_by(Pop) %>% mutate("mean_FT" = mean(FT_day)) %>% distinct(Pop, .keep_all=TRUE)
data_FT_F6$Pop <- as.factor(data_FT_F6$Pop)

#' ## Merged SNP data, FT ex situ
#+ message=FALSE, warning=FALSE
median_AF_existuFT <- left_join(Median_AF_50kb_noContamination, data_FT_F6[,c(2,11)], join_by("Sample_ID"== "Pop"))
median_AF_existuFT <- median_AF_existuFT %>% filter(SNP_tot>10)
summary(median_AF_existuFT)

#' ## Number of populations 
median_AF_existuFT %>% filter(Chrm == "Chr1" & Pos_mean == 75000)  %>% count()




#' # II. Association AF - exsitu FT 


#' ## Mean FT per population association
#+ message=FALSE, warning=FALSE
chrm_list <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
ExsituFT_meanpop_model <- tibble("Chrm" = "1", "Pos" = 1, pvalue = 1, R = 1)

#' Run the model
#+ message=FALSE, warning=FALSE
for (chr in chrm_list) {
  data_chr <- median_AF_existuFT %>% filter(Chrm == chr)
  positions <- unique(data_chr$Pos_mean)
  for (pos in positions) {
    data <- data_chr %>% filter(Pos_mean == pos)
    if (nrow(data) == 0) next
    corresults <- cor.test(data$median_AF, data$mean_FT,  method = "spearman")
    tmp <- c(chr,pos, corresults[["p.value"]], corresults[["estimate"]])
    ExsituFT_meanpop_model <- rbind(ExsituFT_meanpop_model, tmp)}}

#+ message=FALSE, warning=FALSE
ExsituFT_meanpop_model <- ExsituFT_meanpop_model %>% mutate(across(c(2:4), as.numeric)) 
ExsituFT_meanpop_model <- ExsituFT_meanpop_model[-1,]



#' # III. Randomization x 1000 


#' ## Mean FT per population association
#+ eval=F, echo=T
ExsituFT_meanpop_rand_models <- vector("list", 1000)  
for (i in 1:1000) {
  # 1) Randomize type
  pool_map <- median_AF_existuFT %>%
    distinct(Sample_ID,mean_FT) %>%
    mutate(Randomized_exsitu = sample(mean_FT))
  data_random <- median_AF_existuFT[,c(1:5)]
  data_random <- left_join(median_AF_existuFT, pool_map[, c(1, 3)], by = "Sample_ID")
  
  # 2) models for all Chr × Pos
  chrm_list <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
  ExsituFT_meanpop_rand_model <- tibble(Chrm = "1", Pos = 1, pvalue = 1, R = 1)
  
  for (chr in chrm_list) {
    data_chr <- data_random %>% filter(Chrm == chr)
    positions <- unique(data_chr$Pos_mean)
    for (pos in positions) {
      data <- data_chr %>% filter(Pos_mean == pos)
      if (nrow(data) == 0) next
      corresults <- cor.test(data$median_AF, data$Randomized_exsitu,  method = "spearman")
      tmp <- c(chr,pos, corresults[["p.value"]], corresults[["estimate"]])
      ExsituFT_meanpop_rand_model <- rbind(ExsituFT_meanpop_rand_model, tmp)
    }
  }
  
  # 3) Post-processing
  ExsituFT_meanpop_rand_model <- ExsituFT_meanpop_rand_model[-1, ]
  
  # save result of this i
  ExsituFT_meanpop_rand_models[[i]] <- ExsituFT_meanpop_rand_model
  print(paste("randomization", i, "done"))
}


#' Create one data frame
#+ eval=F, echo=T
ExsituFT_meanpop_rand_model_df <- purrr::map2_df(
  ExsituFT_meanpop_rand_models,
  1:1000,
  ~ dplyr::mutate(.x, iter = .y))


#' Merge the data set  
#+ eval=F, echo=T
ExsituFT_meanpop_model$iter <- "Observed"
ExsituFT_meanpop_rand_model_df <- rbind(ExsituFT_meanpop_rand_model_df,ExsituFT_meanpop_model)
ExsituFT_meanpop_rand_model_df <- ExsituFT_meanpop_rand_model_df %>% mutate(across(c(2:4), as.numeric)) 



#' # IV.  Analysis observed vs randomized data -  Correlation with mean FT per populations

#+ message=FALSE, warning=FALSE
ExsituFT_meanpop_rand_model_df <- ExsituFT_meanpop_rand_model_df %>% mutate("Data" = if_else(iter != "Observed", "Randomized", "Observed"))


#' ##  R estimates  
#+ message=FALSE, warning=FALSE
ExsituFT_meanpop_rand_model_df %>%
  ggplot(aes(x=Pos, y=abs(R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of R spearman correlation(Alt_frq ~ exsitu_FT) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))


#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- ExsituFT_meanpop_rand_model_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(R, probs = c(.95)),
                                    Q05 = quantile(R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
Quantile_df <-left_join(ExsituFT_meanpop_model, Quantile_df, by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(R>Q95 & pvalue <0.01 | R<Q05 & pvalue <0.01 , "TRUE","FALSE"))


Quantile_df <- Quantile_df %>% mutate(across(c(2:6), as.numeric))
Quantile_df <- Quantile_df %>% mutate(Effect = ifelse(Estimate_5 == "TRUE" & R >0, "Positif",
                                                      ifelse(Estimate_5 == "TRUE" & R <0, "Negatif", "Neutral")))


#' ## Graphic of the p-value in model and in the 5% distribution of the randomized estimate

Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(pvalue), color = Effect, alpha = Effect)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c( "red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Spearman correlation mean exsitu FT - median AF (+10 SNP per 50kb windows) with 1000 permutations"))



#' ## Correlation at QTL pic
Quantile_df %>% filter(Estimate_5 ==TRUE) %>% slice_min(pvalue, n=10)

exsituFT_df <- median_AF_existuFT %>% filter(Chrm == "Chr4" & Pos_mean == "275000" |
                                                  Chrm == "Chr2" & Pos_mean == "14725000")
exsituFT_df <- exsituFT_df %>% mutate(across(c(2, 4:6), as.numeric)) 

allelic_frq <- as.numeric(unlist(exsituFT_df[exsituFT_df$Chrm == "Chr2","median_AF"]))
variable <- as.numeric(unlist(exsituFT_df[exsituFT_df$Chrm == "Chr2","mean_FT"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
exsituFT_df %>% filter(Chrm== "Chr2") %>%
  ggplot(aes(x=median_AF, y=mean_FT)) +
  geom_point(size=4, color="#EE3A8C") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.75, y=1.69, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("Exsitu FT ~ allelic frequency Chr2"))


allelic_frq <- as.numeric(unlist(exsituFT_df[exsituFT_df$Chrm == "Chr4","median_AF"]))
variable <- as.numeric(unlist(exsituFT_df[exsituFT_df$Chrm == "Chr4","mean_FT"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
exsituFT_df %>% filter(Chrm== "Chr4") %>%
  ggplot(aes(x=median_AF, y=mean_FT)) +
  geom_point(size=4, color="#EE3A8C") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.75, y=1.69, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("Exsitu FT ~ allelic frequency Chr4"))


