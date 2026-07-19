# Swan Portalier 07/2026
# Data parsing and input data preparation


import numpy as np
import pandas as pd
import sys
import gzip


### Command-line arguments ###

########################################
# Paths
########################################
#set of analysed populations (66 or 32 populations)
population_number = sys.argv[1] #must be equal to 32 or 66

#Path to vcf
vcf_F5 = sys.argv[2]

#This script uses thaliana genetic map from Loehmuller lab to find genetic distances between SNPs from VCF
#the map of each chromosome must be named arab_chr{chromosome}_map_loess.txt
#path for maps from Loehmuller lab
maps_path = sys.argv[3]

#path for maps of SNPs (output)
output_SNPmap = sys.argv[4]

#path to output the SNP positions
snps_positions_output = sys.argv[5]

#path to output list of populations
path_site_list = sys.argv[6]

########################################


if population_number == 66:
  ######## This chunk outputs a txt with only accessions for which we have the genotype. Here for the set of 66 populations. 
  #for now there are 69 populations in the output list, but only 66 will be used because we don't have Shannon index for the others
  vcf = "F5_ParentalFilter_GeneRegions_DPclean_segregatingSNP_polarizedRef_noContamination.vcf.gz"
  #output file
  site_list = f"{path_site_list}/site_list.txt"
  with gzip.open(vcf, "rt") as input, open(site_list, "w") as output:
      for line in input:
          if line.startswith("#CHROM"):
              print(line)
              sites = line.strip()
              sites = sites.split("\t")[9:]
              print(sites)
              break
      for site in sites[:-1]:
          output.write(site + "\n")
      output.write(sites[-1])


else:
  ######## Here is the same, for the 32 populations set of genotyped populations (this set is in subdata_32pop_JF)
  site_list = "site_list.txt"
  pops_32 = pd.read_csv("subdata_32pop_JF.csv")
  pops_32 = pops_32["Sample_ID"].tolist()
  with open(site_list, "w") as output:
      for i in pops_32[:-1]:
          output.write(f'{i}' + "\n")
      output.write(f"{pops_32[-1]}")




###############################################################
# Getting chromosome and positions of each SNP in F5
###############################################################
chrom_pos = [[], []]
with gzip.open(vcf_F5, "rt") as vcf:
    c = 0
    for line in vcf:
        if line.startswith("#CHROM"):
            col_names = line.strip().split("\t")
            continue
        if line.startswith("#"):
            continue
        line = line.strip().split("\t")
        chrom_pos[0].append(line[0])
        chrom_pos[1].append(line[1])
        c += 1
        if c%1000 == 0:
            print(c)
chrom_pos = np.array(chrom_pos)

#saving in numpy format
np.save(f"{snps_positions_output}/F5_vcf_SNP_pos_{population_number}.npy", chrom_pos)


###############################################################
# Building genetic map of SNPs
###############################################################
chrom_pos = np.load (snps_positions_output)

for chr in [1, 2, 3, 4, 5]:
    map = pd.read_csv(f"{maps_path}/arab_chr{chr}_map_loess.txt", sep = " ")
    map.rename(columns = {'Rate(cM/Mb)' : 'Rate(M/b)'}, inplace=True)
    map["Rate(M/b)"] = map["Rate(M/b)"] / 100 / 1e6


    #distances between each interval
    map_only = (np.array(map["Rate(M/b)"])[:-1] + np.array(map["Rate(M/b)"])[1:]) / 2
    map_only = map_only.tolist()
    map_only = [float(np.array(map["Rate(M/b)"])[0])] + map_only + [float(np.array(map["Rate(M/b)"])[-1])]
    map_only = np.array(map_only)


    chrom_pos_curr = chrom_pos
    chrom_pos_curr = chrom_pos_curr[1, chrom_pos_curr[0, :] == f"Chr{chr}"]
    chrom_pos_curr = np.array([int(i) for i in chrom_pos_curr])

    distances = []

    #looping over positions of SNPs (from vcf)
    for idx, snp in enumerate(chrom_pos_curr[:-1]):
        #physical distance between current SNP and the next one
        current_physical_dist = chrom_pos_curr[idx + 1] - snp
        print(current_physical_dist)
        print(idx, snp, len(chrom_pos_curr))
        snp_bins = np.searchsorted ( map["Position(bp)"], [snp, chrom_pos_curr[idx + 1]])
        #in case one boundary of the interval is one of the limits in the downloaded map
        if snp in np.array(map["Position(bp)"]) :
            distances.append(((map.loc[map["Position(bp)"] == snp, "Rate(M/b)"]).tolist()[0]) * current_physical_dist)
            continue
        if chrom_pos_curr[idx + 1] in np.array(map["Position(bp)"]):
            distances.append(((map.loc[map["Position(bp)"] == chrom_pos_curr[idx + 1], "Rate(M/b)"]).tolist()[0]) * current_physical_dist)
            continue

        if snp_bins[0] == snp_bins[1] :
            distances.append(map_only[snp_bins[0]] * current_physical_dist)
        else:
            distances.append(map.loc[snp_bins[0], "Rate(M/b)"] * current_physical_dist)
    distances = np.array(distances)
    np.save(f"{output_SNPmap}/SNP_recombination_map_chr{chr}.npy", distances)