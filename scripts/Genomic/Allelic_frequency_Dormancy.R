#' ---
#'title: "Allelic frequency and Dormancy"
#'author: "Justine Floret"
#'date: "2026-01-09"
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
library(qqman) # Manhattan plot
library(knitr)
library(stringr) #remove space in string
# Stat. libray
library(lmerTest) #glht
library(lsmeans) #pairwise comparison
library(emmeans)
library(glmmTMB) # negative binomial model
library(vegan) #RDA
library(missMDA) #Missing data imputation, pca


#' # I. Preparing the data

#' ## 1. Allelic frequency data
#+ message=FALSE, warning=FALSE
Median_AF_50kb_noContamination <- read_csv("/home/justine/sciebo/Genomic/assays/F5_pools_medianAF_50kb_noContamination.csv")

Median_AF_50kb_noContamination$Chrm <- as.factor(Median_AF_50kb_noContamination$Chrm)
Median_AF_50kb_noContamination$Pos_mean <- as.numeric(Median_AF_50kb_noContamination$Pos_mean)
Median_AF_50kb_noContamination$Sample_ID <- as.factor(Median_AF_50kb_noContamination$Sample_ID)
Median_AF_50kb_noContamination$median_AF <- as.numeric(Median_AF_50kb_noContamination$median_AF)
Median_AF_50kb_noContamination$SNP_tot <- as.numeric(Median_AF_50kb_noContamination$SNP_tot)

summary(Median_AF_50kb_noContamination)


#' ## 2. Ex situ Dormancy
#+ message=FALSE, warning=FALSE
dormancy_data <- read_delim("/home/justine/sciebo/3_Exsitu_phenotype/Dormancy/2_assays/dormancy_transplanted_Indv_2025.csv")
dormancy_data$ID_F5_parent <- as.factor(dormancy_data$ID_F5_parent)
dormancy_data$Counting_Rep <- as.factor(dormancy_data$Counting_Rep)
dormancy_data$run <- as.factor(dormancy_data$run)
dormancy_data$seed_count_germinated <- as.numeric(dormancy_data$seed_count_germinated)
dormancy_data$seed_count_total <- as.numeric(as.character(dormancy_data$seed_count_total))

#' Transforming the data set for analysis
dormancy_data <- dormancy_data %>% mutate(seed_count_total = replace_na(seed_count_total,0))
dormancy_data <- dormancy_data %>% mutate(
  ID_F5_parent = case_when(
    ID_F5_parent == "165-4" ~ "156-4",
    ID_F5_parent == "165-8" ~ "156-8",
    ID_F5_parent == "165-17" ~ "156-17",
    TRUE ~ ID_F5_parent))

dormancy_data$ID_parent <- dormancy_data$ID_F5_parent
dormancy_data <-dormancy_data %>% separate_wider_delim(ID_F5_parent, "-", names = c("ID_pop", "ID_pop_indv"))
dormancy_data$ID_pop <- as.factor(dormancy_data$ID_pop)
dormancy_data$ID_parent <- as.factor(dormancy_data$ID_parent)
dormancy_data <- dormancy_data %>% filter(Counting_Rep ==1)

#' Replacing the run by the number of week since bagging the first plant
dormancy_data <- dormancy_data %>%
  mutate(Week_afterbagging = case_when(run == 1 ~ 5, run == 2 ~ 8, run == 3 ~ 11, run == 4 ~ 13,) |> as.factor(),
         run = as.factor(run),
         Block = as.factor(Block))

#' When total seed count = 0 or NA, replace by 30
dormancy_data <- dormancy_data %>%
  mutate(seed_count_total = replace(seed_count_total, seed_count_total ==0, 30))
dormancy_data <- dormancy_data %>% dplyr::mutate(seed_count_total = replace_na(seed_count_total, 30))
summary(dormancy_data)

#' Filter out individual with NA data 
subset(dormancy_data,
       seed_count_germinated < 0 |
         seed_count_total < 0 |
         seed_count_germinated > seed_count_total |
         is.na(seed_count_germinated) |
         is.na(seed_count_total))
dormancy_data <- dormancy_data %>% filter(seed_count_germinated <= seed_count_total)


#' ### Extract the germination slope  

#' I need at least three observations (ideally 4+) per individual, i.e. more than 3 germination run.  
#' Some individual did no have enough seeds for three run of germination.  
#+ message=FALSE, warning=FALSE
dormancy_data_filtered <- dormancy_data %>% filter(Week_afterbagging != 11 | Block !=1)
dormancy_data_filtered$Week_afterbagging <- as.numeric(dormancy_data_filtered$Week_afterbagging)
dormancy_data_filtered <- dormancy_data_filtered %>%
  group_by(ID_parent) %>%
  filter(n_distinct(Week_afterbagging) >= 3) %>%
  ungroup()

model_slope_nested <- glmer(
  cbind(seed_count_germinated, seed_count_total - seed_count_germinated) ~ Week_afterbagging +
    (Week_afterbagging | ID_pop) +
    (Week_afterbagging | ID_pop:ID_parent),
  data = dormancy_data_filtered,
  family = binomial
)
summary(model_slope_nested)


#' Extract just the slope for Week_afterbagging per individual
#+ message=FALSE, warning=FALSE
population_slopes <- coef(model_slope_nested)$ID_pop
population_slopes <- as.data.frame(population_slopes)
population_slopes$ID_pop <- rownames(population_slopes)
colnames(population_slopes)[2] <- "Germ_slope"

ggplot(population_slopes, aes(x=Germ_slope)) + geom_histogram()

#' ## Merged SNP data, Dormancy
#+ message=FALSE, warning=FALSE
median_AF_Dormancy <- left_join(Median_AF_50kb_noContamination, population_slopes[,c(2,3)], join_by("Sample_ID"== "ID_pop"))
median_AF_Dormancy <- median_AF_Dormancy %>% filter(SNP_tot> 10)
summary(median_AF_Dormancy)




#' # II Association AF - Germination slope 

#' Run the model
#+ message=FALSE, warning=FALSE
chrm_list <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
germination_slope_model <- tibble("Chrm" = "1", "Pos" = 1, pvalue = 1, R = 1)

for (chr in chrm_list) {
  data_chr <- median_AF_Dormancy %>% filter(Chrm == chr)
  positions <- unique(data_chr$Pos_mean)
  for (pos in positions) {
    data <- data_chr %>% filter(Pos_mean == pos)
    if (nrow(data) == 0) next
    corresults <- cor.test(data$median_AF, data$Germ_slope,  method = "spearman")
    tmp <- c(chr,pos, corresults[["p.value"]], corresults[["estimate"]])
    germination_slope_model <- rbind(germination_slope_model, tmp)
  }
}

#' Process the data
#+ message=FALSE, warning=FALSE
germination_slope_model <- germination_slope_model[-1,] 
germination_slope_model <- germination_slope_model %>% mutate(across(c(2:4), as.numeric)) 



#' # III. Randomization x 1000 

#' Randomization: 1000 times
#+ eval=F, echo=T
Dormancy_rand_models <- vector("list", 1000)  

#+ eval=F, echo=T
for (i in 1:1000) {
  # 1) Randomize type
  pool_map <- median_AF_Dormancy %>%
    distinct(Sample_ID,Germ_slope) %>%
    mutate(Randomized_germslope = sample(Germ_slope))
  data_random <- left_join(median_AF_Dormancy, pool_map[, c(1, 3)], by = "Sample_ID")
  
  # 2) models for all Chr × Pos
  chrm_list <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
  germination_slope_model <- tibble("Chrm" = "1", "Pos" = 1, pvalue = 1, R = 1)
  
  for (chr in chrm_list) {
    data_chr <- data_random %>% filter(Chrm == chr)
    positions <- unique(data_chr$Pos_mean)
    for (pos in positions) {
      data <- data_chr %>% filter(Pos_mean == pos)
      if (nrow(data) == 0) next
      corresults <- cor.test(data$median_AF, data$Randomized_germslope,  method = "spearman")
      tmp <- c(chr,pos, corresults[["p.value"]], corresults[["estimate"]])
      germination_slope_model <- rbind(germination_slope_model, tmp)
    }
  }
  
  # 3) Post-processing
  germination_slope_model <- germination_slope_model[-1, ]
  
  # save result of this i
  Dormancy_rand_models[[i]] <- germination_slope_model
  print(paste("randomization", i, "done"))
}


#' Create one data frame
#+ eval=F, echo=T
saveRDS(Dormancy_rand_models, "Dormancy_rand_models" )
Dormancy_rand_model_df <- purrr::map2_df(
  Dormancy_rand_models,
  1:1000,
  ~ dplyr::mutate(.x, iter = .y))


#' Merge the data set  
#+ eval=F, echo=T
germination_slope_model$iter <- "Observed"
Dormancy_rand_model_df <- rbind(Dormancy_rand_model_df,germination_slope_model)
Dormancy_rand_model_df <- Dormancy_rand_model_df %>% mutate(across(c(2:4), as.numeric)) 



#' # IV. Analysis observed vs randomized data

#+ message=FALSE, warning=FALSE
Dormancy_rand_model_df <- Dormancy_rand_model_df %>% mutate("Data" = if_else(iter != "Observed", "Randomized", "Observed"))


#' ## Quantile distribution - extraction
#+ message=FALSE, warning=FALSE
Quantile_df <- Dormancy_rand_model_df %>% filter(iter != "Observed") %>%
  group_by(Chrm, Pos) %>% summarise(Q95 = quantile(R, probs = c(.95)),
                                    Q05 = quantile(R, probs = c(.05))) 
Quantile_df$Pos <- as.numeric(Quantile_df$Pos)
germination_slope_model$Pos <- as.numeric(germination_slope_model$Pos)
Quantile_df <-left_join(germination_slope_model, Quantile_df, by = c("Chrm", "Pos"))
Quantile_df <- Quantile_df %>% mutate(across(c(2:6), as.numeric))

Quantile_df <- Quantile_df %>% mutate(Estimate_5 = ifelse(R>Q95 & pvalue <0.01 | R<Q05 & pvalue <0.01 , "TRUE","FALSE"))
Quantile_df <- Quantile_df %>% mutate(Effect = ifelse(Estimate_5 == "TRUE" & R >0, "Positif",
                                                      ifelse(Estimate_5 == "TRUE" & R <0, "Negatif", "Neutral")))



#' ## All estimate   
#+ message=FALSE, warning=FALSE
Dormancy_rand_model_df %>%
  ggplot(aes(x=Pos, y=abs(R), color = iter, alpha = Data)) + 
  geom_line()  +
  theme_minimal() +
  labs(title = "Absolute value of R spearman correlation(Alt_frq ~ germination slope) per position with 1000 permutations") +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  scale_color_manual(values = c("Observed" = "red"))+ 
  scale_alpha_manual(values = c(1, 0.05))


#' ## Graphic of the p-value in model and in the 5% distribution of the randomized estimates + Effect
Quantile_df  %>% 
  ggplot(aes(x=Pos, y=-log10(pvalue), color = Effect, alpha = Effect)) + 
  geom_point(alpha=0.5)  +
  scale_color_manual(values = c("red", "lightgrey", "#00BFFF")) +
  theme_minimal() +
  facet_wrap(~Chrm , strip.position="bottom", scales = "free_x") +
  ggtitle(paste("Spearman correlation(Alt_frq ~ germination slope) Germination slope - median AF (+10 SNP per 50kb windows) with 1000 permutations"))


#' ## Correlation at QTL pic
Quantile_df %>% filter(Estimate_5 ==TRUE) %>% slice_min(pvalue)

germination_df <- median_AF_Dormancy %>% filter(Chrm == "Chr2" & Pos_mean == "14425000" |
                                              Chrm == "Chr5" & Pos_mean == "10025000")
germination_df <- germination_df %>% mutate(across(c(2, 4:6), as.numeric)) 

allelic_frq <- as.numeric(unlist(germination_df[germination_df$Chrm == "Chr2","median_AF"]))
variable <- as.numeric(unlist(germination_df[germination_df$Chrm == "Chr2","Germ_slope"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
germination_df %>% filter(Chrm== "Chr2") %>%
  ggplot(aes(x=median_AF, y=Germ_slope)) +
  geom_point(size=4, color ="#68228B") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  theme_minimal() +
  ggtitle(paste("Germ_slope ~ allelic frequency Chr2"))


allelic_frq <- as.numeric(unlist(germination_df[germination_df$Chrm == "Chr5","median_AF"]))
variable <- as.numeric(unlist(germination_df[germination_df$Chrm == "Chr5","Germ_slope"]))
correlation <- cor.test(allelic_frq,variable ,method = "spearman")
germination_df %>% filter(Chrm== "Chr5") %>%
  ggplot(aes(x=median_AF, y=Germ_slope)) +
  geom_point(size=4, color ="#68228B") +
  geom_smooth(method=lm , color="darkgrey", se=TRUE) +
  geom_text(x=0.75, y=1.69, 
            label=paste("pvalue =",round(correlation[["p.value"]], 3),
                        "and R =",round(correlation[["estimate"]],3))) +
  theme_minimal() +
  ggtitle(paste("Germ_slope ~ allelic frequency Chr5"))
