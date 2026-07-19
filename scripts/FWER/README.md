# Deternining the p-value threshold to control the familiy-wise error rate 
The following pipeline has been run for a set of N=32 and N=66 sampling sites (or population hereafter) separately, for a target FWER of 0.01, 0.05 and 0.25.

## Software requirements
The following are required:
- Python 3.12.3
- numpy 2.4.2
- pandas 3.0.1
- R 4.4.2

## Step 1: Data parsing and input data preparation: Getting list of analysed populations Calculate the genetic distances between SNPs Getting list of analysed SNPs
Input files:
- VCF file
F5_ParentalFilter_GeneRegions_DPclean_segregatingSNP_polarizedRef_noContamination.vcf.gz
- Recombination maps from LohmuellerLab

<h5><strong><code>input_preparation.py</code></strong></h5>


```bash
module load python

#number of populations (must be equal to 66 or 32)
population_number=66

#path to the folder with recombination maps from LohmuellerLab. Those are the files arab_chr{chromosome_index}_map_loess.txt from the LohmuellerLab's github (https://github.com/LohmuellerLab/arabidopsis_recomb_maps)
maps_path=path/folder_maps

#vcf of the F5 population
vcf_F5=path/F5_ParentalFilter_GeneRegions_DPclean_segregatingSNP_polarizedRef_noContamination.vcf

#path for maps of SNPs (output) (SNP_recombination_map_chr{chromosomes}.npy)
output_SNPmap=path/folder_output_SNPmap

#path to output the SNP positions (output)
snps_positions_output=path/folder_positions_output

#path to output list of populations
path_site_list=path/folder_populations_output

python input_preparation.py $population_number $vcf_F5 $maps_path $vcf_F5 $output_SNPmap $snps_positions_output
```

Output files:
- List of sampling sites included: `site_list_${population_number}.txt`
- List of SNPs in the input to be simulated: `F5_vcf_SNP_pos_${population_number}.npy`
- Genetic distances between the SNPs to be simulated: `SNP_recombination_map_chr${chromosomes}.npy`

</br>

## Step 2: Perform the simulations
Input files:
- File with population sizes: `population_size_JF.cov`
- List of sampling sites included: `site_list_${population_number}.txt`
- Genetic distances between the SNPs to be simulated: `rec_rates_${chromosomes}.npy`
- Simulation parameters: `parameters.ini`

<h5><strong><code>simulation.py</code></strong></h5>

```bash
module load python

#simulation replicate
replicate=1

#size of the set of populations of the analysis (must be 66 or 32)
population_number=66

#input paths
# path of site list "site_list_{population_number}.txt"
pop_list_path=path/folder_site_list

# path of populations sizes "Population_size_JF.csv"
pop_sizes_path=path/folder_Population_size_JF

#path to simulation parameters file "parameters.ini"
param_path=path/folder_parameters

#path for SNP recombination maps "rec_rates_{chromosome}.npy"
path_maps=path/folder_rec_rates


# output paths 
# path of list of sites with only those in which there are more than 0 indiviuals (average over April sampling over 3 years) 
site_list_non0size=path/folder_site_list

# path of frequencies obtained from simulated F5 "replicate_{replicate}.txt.npy"
path_freqs_out=path/folder_freqs_out

python simulation.py $replicate $population_number $pop_list_path $pop_sizes_path $param_path $path_maps $site_list_non0size $path_freqs_out
```


We recommend this script to be run in parallel with a job array.


Output files: 
- Array with the frequencies from the F5 simulated population (population*SNPs): `replicate_${replicate}.txt.npy`
- Path of list of sites with only those in which there are more than 0 indiviuals `site_list_non0size.txt`
</br>

## Step 3: Reconstruct the SNP windows as used in `/Genomics`

Input files:
- Simulation outputs: `windowed_${replicate}.txt.npy`
- List of SNPs in the input to be simulated: `F5_vcf_SNP_pos_${population_number}.npy`

<h5 a><strong><code>mk_windows.py</code></strong></h5>

```bash
module load python

#input paths
#path to SNPs positions (F5_vcf_SNP_pos_{population_number}.npy)
snp_pos_path=path/folder_positions_output

#path to frequencies from simulation output (replicate_{replicate}.txt.npy) (shape = population X windows)
freqs_path=path/folder_freqs_out



#output paths
#paths of the file to save boundaries of windows (windows_boundaries_chr{chr}.npy)
path_windows_bounds=path/folder_windowsBoundaries

#path for file to save windows which actually SNPs (not every window of 50 kb do):
path_withSNP_windows=path/folder_windows_WithSNPs

#path for the file to save the median frequency of each window, in each population (windowed_replicate_{replicate}.txt):
path_windowed_frequency=path/folder_windowed_freqs

#number of populations (66 or 32)
population_number=66

python mk_windows.py $snp_pos_path $freqs_path $path_windows_bounds $path_withSNP_windows $path_windowed_frequency $population_number
```

Output files: 
- Median frequency of simulation output for each window: `windowed_replicate_${replicate}.txt`
- Boundaries of windows: `windows_boundaries_chr${chr}.npy`
- List of windows that actually contain SNPs: `windows_with_SNPs${chr}.npy`

</br>

## Step 4: Perform Spearman's rank correlation test
Input files: 
- Simulated SNP windows: `windowed_replicate_${replicate}.txt`
- File with biodiverity index measure for each population: `Biodiversity_index_J.csv`
- File with population sizes: `Population_size.py`
- List of sampling sites included: `pop_list_${population_number}.txt`

<h5 a><strong><code>read_corr_envVSfreq.R</code></strong></h5>

```bash
module load R

#input paths
#path of the file with list of population with non 0 population sizes (site_list_non0size.txt)
path_nonempty_sitelist=path/folder_site_list

#path of file with populations sizes (Population_size_JF.csv)
path_popsizes=path/folder_Population_size_JF

#path of file with environmental parameter (Biodiversity_index_JF.csv)
path_env_param=path/folder_Biodiversity_index

#path of the file with frequency outputed by simulation (windowed_replicate_{replicate}.txt)
path_simu_freqs=path/folder_windowed_freqs

#output path
#path of the file to save the Spearman correlation Rho statistic (estimates_{replicate}.csv)
path_estimate=path/folder_estimates

#path of the file to save the Spearman correlation p_value (pvals_{replicate}.csv)
path_p_value=path/folder_pvalues

#replicate index
replicate=1

Rscript read_corr_envVSfreq.R $path_nonempty_sitelist $path_popsizes $path_env_param $path_simu_freqs $path_estimate $path_p_value $replicate


```

Output files: 
- Spearman correlation $\rho$ statistics: `estimates_${replicate}.csv`
- Spearman correlation p_values: `pvals_${replicate}.csv`


</br>


## Step 5: Determine the p-value threshold for the target FWER
Input files: 
- Spearman correlation $\rho$ statistics: `estimates_{replicate}.csv`
- Spearman correlation p_values: `pvals_${replicate}.csv`


<h5 a><strong><code>thresh_finder.py</code></strong></h5>

```bash
module load python

#path of the files with pvalues calculated from simulation output (pvals_{replicate}.csv)
p_vals_path=path/folder_pvalues

# FWER for which we want to find the corresponding p value threshold
target_fwer=0.05


python thresh_finder.py $p_vals_path $target_fwer > results_${population_number}_FWER${target_fwer}.txt

```
We recommend this script to be run in parallel with a job array.


