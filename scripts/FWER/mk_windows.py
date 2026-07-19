# Swan Portalier 07/2026
# This script takes the output of simulation in input and group SNPs into windows of 50 kb, with frequency computed as the median frequency of SNPs in the window




import numpy as np
import pandas as pd
import sys


### Command-line arguments ###

#input paths
#SNPs positions (F5_vcf_SNP_pos_{population_number}.npy)
snp_pos_path = sys.argv[1]
#path to frequencies from simulation output (replicate_{replicate}.txt.npy) (shape = population X windows)
freqs_path = sys.argv[2]




#output paths
#paths of the file to save boundaries of windows (windows_boundaries_chr{chr}.npy)
path_windows_bounds = sys.argv[3]
#path for file to save windows which actually SNPs (not every window of 50 kb do):
path_withSNP_windows = sys.argv[4]
#path for the file to save the median frequency of each window, in each population:
path_windowed_frequency = sys.argv[5]

population_number = sys.argv[6]
#############################################""




#path = os.path.dirname(os.path.abspath(__file__))
#arc_path = "/projects/ag-demeaux/Swan/ARCs_Justine_simulations/ARC_non_vectorized"
replicates = np.arange(1, 500)







#loading SNP positions
snp_pos = np.load(f"F5_vcf_SNP_pos_{population_number}.npy")
chromosomes = np.unique (snp_pos[0, :])


#list of simulations which did not finished due to time limit
unfinished = []
c = 0
for replicate in replicates:
    all_chr_medfreqs = np.array([]).reshape(68, 0)
    #loading frequencies
    for chr in chromosomes:

        #loading SNP positions
        snp_pos_chr = snp_pos

        
        mask_current_chromosome = snp_pos_chr[0, :] == chr
        snp_pos_chr = snp_pos_chr[:, mask_current_chromosome]
        snp_pos_chr = np.astype(snp_pos_chr[1, :], int)


        #position of boundaries between windows
        windows_bound = np.arange(0, snp_pos_chr[-1], 50000)[1:]




        #finding window each SNP belongs to
        w_of_eachSNPs = np.searchsorted(windows_bound, snp_pos_chr)
        #saving boundaries of windows 
        np.save(f"{path_windows_bounds}/windows_boundaries_chr{chr}.npy", windows_bound)



        




        frequencies_replicate = np.load(f"{freqs_path}/replicate_{replicate}.txt.npy")

        #keeping only the focal chromosome
        #keeping only those included in the set of SNP on which the test is performed in real data
        frequencies_replicate = frequencies_replicate[:, mask_current_chromosome]

        #compute median inside each window
        df_frequencies = pd.DataFrame(frequencies_replicate.T)
        df_windows = pd.DataFrame({"window" : w_of_eachSNPs})
        df_windows = pd.concat([df_windows, df_frequencies], axis = 1)
        #saving lists with windows that contain SNPs
        np.save(f"{path_withSNP_windows}/windows_with_SNPs{chr}.npy", w_of_eachSNPs)
        df_windows = df_windows.groupby(["window"]).median()
        

        df_windows = np.array(df_windows).T

        if df_windows.shape[0] != 68:
            break 

        all_chr_medfreqs = np.concat([all_chr_medfreqs, df_windows], axis = 1)
    c +=1
    print("done_simulations:", c, flush=True)

    if df_windows.shape[0] != 68:
        unfinished.append (replicate)
        continue 
    

    #saving median frequency in each window for the replicate
    np.savetxt(f"{path_windowed_frequency}/windowed_replicate_{replicate}.txt", all_chr_medfreqs)
    print(replicate, flush=True)



print("unfinished: \n", unfinished, flush=True)

