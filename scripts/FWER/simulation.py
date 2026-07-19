# Swan Portalier 07/2026
# Running the simulation of the experimental setup, assuming no selection
import configparser
import numpy as np
import pandas as pd
import sys
import importlib
import functions as fun
importlib.reload(fun)
import math

### Command-line arguments ###

replicate = sys.argv[1]
population_number = sys.argv[2]
#input paths
# path of site list "site_list_{population_number}.txt"
pop_list_path = sys.argv[3]

# path of populations sizes "Population_size_JF.csv"
pop_sizes_path = sys.argv[4]

#path to simulation parameters file "parameters.ini"
param_path = sys.argv[5]

#path for SNP recombination maps "rec_rates_{chromosome}.npy"
path_maps = sys.argv[6]


# output paths 
# path of list of sites with only those in which there are more than 0 indiviuals (average over April sampling over 3 years) (site_list_non0size.txt)
site_list_non0size = sys.argv[7]

# path of frequencies obtained from simulated F5 "replicate_{replicate}.txt"
path_freqs_out = sys.argv[8]
###################################################









#list of sites
pop_list = pd.read_csv(f"{pop_list_path}/site_list.txt", header = None).iloc[:, 0]
pop_list = np.array(pop_list.tolist())

#getting population sizes 
pop_sizes = pd.read_csv(f"{pop_sizes_path}/Population_size_JF.csv", dtype={6: pd.Int64Dtype()})
pop_sizes.index = pop_sizes["Site_ID"]
pop_sizes = pop_sizes.loc[np.isin(pop_sizes["Site_ID"], pop_list), :]
aprils = pop_sizes[["April_22", "April_23", "April_24"]]
pop_sizes["April_mean"] = aprils.mean(axis = 1)
#One has a pop size of 0, at all April time points
pop_sizes = pop_sizes.loc[pop_sizes["April_mean"]!=0, :]
pop_sizes["Site_ID"].to_csv(f"{site_list_non0size}/site_list_non0size.txt", index = False, header=False)


#loading parameters
config = configparser.ConfigParser()
config.read(f"{param_path}/parameters.ini")
N_ini = int(config["general"]["N_ini"])
mu = float(config["general"]["mu"])
output_path = config["general"]["output_path"]
genome = config["general"]["genome"]
genome = np.array(genome.split(","), dtype = int).tolist()

#seeds
seed = np.random.SeedSequence().entropy
rng = np.random.default_rng(seed)



#total number of nucleotides
L = np.sum(genome)


#CROSSING OVER rates between locus i and locus i+1 (indices are the same as in all_effects)
#it is a CROSSING OVER , not recombination rate. There can be 2n CROSSING OVER between 2 loci , and then no recombination there
#length has to be chromosome size-1 for each rec_rates_{chromosome}.npy
#organized as a dictionary {chromosome : list of crossing over rates between each locus}
interSNP_genetic_dist = {} 
#in the same loop, dictionary with boundaries of chromosome is created (boundaries in index of genotypes matrix), 
#and also dict of loci of each chromosome (in index of genotypes matrix)
bounds_chr = {} #index of first and last nucleotide of each chromosome (python index starting at 0)
loci_chr = {}
genome_0 = np.cumsum([0] + genome)
for chromosome in range (len(genome)):
    interSNP_genetic_dist_chr = np.load(f"{path_maps}/rec_rates_{chromosome}.npy")
    interSNP_genetic_dist_chr[np.isnan(interSNP_genetic_dist_chr)] = 0
    interSNP_genetic_dist_chr[interSNP_genetic_dist_chr < 0] = 0
    interSNP_genetic_dist_chr = interSNP_genetic_dist_chr
    interSNP_genetic_dist[chromosome] = interSNP_genetic_dist_chr
    bounds_chr[chromosome] = [int(genome_0 [chromosome]), int(genome_0 [chromosome + 1]) - 1]
    loci_chr[chromosome] = list(range(bounds_chr[chromosome][0], bounds_chr[chromosome][1] + 1))


for i in range(len(genome)):
    if len(interSNP_genetic_dist[i]) != genome[i] - 1 :
        print("WARNING, len of co rates objects must be length of chromosomes - 1", flush = True)

#to save frequencies in F5
final_freq_pops = []


for site_index, site_popsize in enumerate(pop_sizes["April_mean"]) :
    N_post = math.ceil(site_popsize)
    print(f"{site_index}th site, with popsize {N_post}", flush = True)
    print(f"({len(pop_sizes)} populations in total)", flush = True)
    print("test1", flush = True)
    #building genotypes of an F1
    # (genotype of individual i is chromosome at row i and chromosome at row i + 1 (all in python index))
    base_individual = np.concat([np.repeat(0, L).reshape(1, L), np.repeat(1, L).reshape(1, L)], axis = 0)

    #genotyeps of all F1s
    genotypes = np.tile(base_individual,(N_ini,1))


    cc0 = []
    c = 0
    ################################################################################
    # generating F2: (no Bottleneck)
    ################################################################################

    #fitness computation (here neutral, every plant has the same fitness)
    fitnesses = np.repeat(1/(genotypes.shape[0]//2), (genotypes.shape[0]//2))

    #reproduction 
    genotypes, cc0_, c_ = fun.reproduction (genotypes, fitnesses, N_ini, interSNP_genetic_dist, genome, bounds_chr, rng)

    ################################################################################
    # 3 generations 
    ################################################################################
    for i in range(3):
        #bottleneck (with only N_post plants that grow)
        if i == 0:
            #individuals that survive the bottleneck
            inc_bott = rng.choice(np.arange(0, N_ini), N_post, replace=False)
            #genotypes index that survive the bottleneck
            g_bott = []
            for k in inc_bott:
                g_bott.extend ([int(k)*2, int(k*2 + 1)])
            genotypes = genotypes[g_bott, :]
        print(f"F{i + 3} in progress", flush = True)
        fitnesses = np.repeat(1/(genotypes.shape[0]//2), (genotypes.shape[0]//2))
        genotypes, cc0_, c_ = fun.reproduction (genotypes, fitnesses, N_post, interSNP_genetic_dist, genome, bounds_chr, rng)
        cc0.extend(cc0_)
        c = c + c_
        print(f"F{i + 3} done", flush = True)

    #current population frequency
    curr_freq = np.sum(genotypes, axis=0)/2/N_post
    final_freq_pops.append(curr_freq)


output = np.array(final_freq_pops)
np.save(f"{output_path}/replicate_{replicate}.txt", output)
