# Swan Portalier 07/2026
# Perform Spearman's rank correlation test


args <- commandArgs(trailingOnly = TRUE)
### Command-line arguments ###

#input paths
#path of the file with list of population with non 0 population sizes (site_list_non0size.txt)
path_nonempty_sitelist <- args[1]
#path of file with populations sizes
path_popsizes <- args[2]
#path of file with environmental parameter
path_env_param <- args[3]
#path of the file with frequency outputed by simulation (windowed_replicate_{replicate}.txt)
path_simu_freqs <- args[4]

#output path
#path of the file to save the Spearman correlation Rho statistic
path_estimate <- args[5]
#path of the file to save the Spearman correlation p_value
path_p_value <- args[6]

#replicate index
replicate <- args[7]
###############################################


site_list = read.csv(glue::glue("{path_nonempty_sitelist}/site_list_non0size.txt"), header = FALSE)
site_list <- site_list$V1

#population sizes (frequency matrix populations are in this order)
pop_sizes = read.csv(glue::glue("{path_popsizes}/Population_size_JF.csv"))
pop_names_freq <- pop_sizes$Site_ID
pop_names_freq <- pop_names_freq[pop_names_freq %in% site_list]

#dataframe with environmental values
env = read.csv(glue::glue("{path_env_param}/Biodiversity_index_JF.csv")
rownames(env) <- env$Site_ID
env <- env["Shannon_index"]
env <- env[as.integer(rownames(env)) %in% site_list, , drop = FALSE]



#frequencies (rows are in the order of pop_names_freq)
#the populations in frequencies are ALL in site_list
x <- read.csv(glue::glue("{path_simu_freqs}/windowed_replicate_{replicate}.txt"), sep = " ", header = FALSE)
rownames(x) <- pop_names_freq

#keeping only populations in common between frequencies and environmental variable
comm <- merge (env, x, by = 'row.names', all.x = FALSE, all.y = FALSE)
comm$Row.names <- as.integer (comm$Row.names)
comm <- comm[order(comm$Row.names), ]
rownames(comm) <- comm$Row.names
env = comm["Shannon_index"]

x = as.matrix(comm[, c(-1, -2)])

p_vals_SNPs = c()
estimates = c()

# Performing the test
for (snp in 1:ncol(x)){
    
    freq_SNP = x[, snp]

    p_vals_SNPs <- c(p_vals_SNPs, cor.test(as.vector(env)[[1]], freq_SNP,  method = "spearman")$p.value)
    estimates <- c(estimates, cor.test(as.vector(env)[[1]], freq_SNP,  method = "spearman")$estimate)
    if (snp == 1){
        print(cor.test(as.vector(env)[[1]], freq_SNP,  method = "spearman")$p.value)
    }
    if (snp %% 1000 == 0){
        print("Number of SNP to analyse:")
        print(ncol(x))
        print("current SNP:")
        print(snp)
    }
}


# Saving the p values and estimate from Spearman test
write.table(p_vals_SNPs, file = glue::glue("{path_p_value}/pvals_{replicate}.csv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(estimates, file = glue::glue("{path_estimate}/estimates_{replicate}.csv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
