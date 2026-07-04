# Dozens of genetic variants sustain adaptation to urban spatial heterogeneity in Arabidopsis thaliana

This repository contains the data and analysis scripts associated with the pape:  

doi: [https://doi.org/10.64898/2026.06.27.734931](https://doi.org/10.64898/2026.06.27.734931)

## Abstract

Urbanization creates mosaics of microhabitats that sustain diverse plant communities and offer untapped opportunities
to understand contemporary ecological and evolutionary processes. Using a city-wide colonization experiment in Cologne
(Germany), we identified the environmental factors limiting establishment and persistence of the ruderal species
Arabidopsis thaliana and the genetic variants enabling adaptation. We show that standing genetic variation is essential
for realizing the species’ urban niche, enabling rapid adaptation to fine-scale gradients in disturbance, vegetation,
and soil conditions. Despite originating from only two parental genotypes, populations evolved substantial adaptive
differentiation within three generations, revealing a highly polygenic basis of local adaptation. More broadly, our
approach establishes cities as powerful open-air laboratories for uncovering and fostering ecological and evolutionary
processes that generate and sustain biodiversity in human-dominated landscapes.

## Overview

This repository provides the data, code, and documentation needed to reproduce the analyses presented in the paper. It is intended to support transparency, reuse, and reproducibility.

## Repository contents

- `data/` — Input data used in the analyses.
- `scripts/` — R scripts or other code used for data processing and statistical analyses.
- `README.md` — This file.

## Data description

The repository includes the following data files:

**Dormancy_exsitu.ods** 
Release of primary dormancy across time in transplanted A. thaliana populations. Germination rate (number of germinated
seeds relative to total seeds) was measured on the progeny of 359 transplanted individuals from 87 populations at
different time points after fruit ripening (number of weeks).

**Ecological_parameters.ods**
Environmental factors characterizing the urban habitats. The ecological factors were measured at different key
phenological stages of A. thaliana between September 2021 and May 2022 across 393 sites in Cologne. 

**F5_pools_medianAF_50kb_noContamination.csv**	
Median allele frequencies of F5 generation, calculated in 50 kb windows across the genome. Sequencing reads were mapped
to the A. thaliana Col0 reference genome, and alleles were polarized so that the MIL2 genotype was considered as the
reference. 

**FT_exsitu.ods**
Ex situ flowering time data from a total of 958 progeny from 138 transplanted populations (17 meadows, 38 pavements, 23 tree beds and 34 walls) were sown under longday conditions.

**Population_biology_insitu.ods**
Date of the habitat site monitoring. Each year, sites were monitored at several critical stages of the A. thaliana life
cycle: germination, flowering, and fruiting. Sites where no population was established were further monitored at least
twice per year (October and April) to ensure we did not miss a population that would have germinated late. These
observations enabled quantification of establishment and persistence rates, population size, and onset of flowering:

## Scripts



## Reproducing the analyses

To reproduce the results in the paper:

1. Clone this repository.
2. Install the required software and packages.
3. Run the scripts in the order indicated in the `scripts/` folder.
4. Generated outputs will be saved in the `output/` folder.


## Citation

If you use this repository, please cite the associated paper:

Dozens of genetic variants sustain adaptation to urban spatial heterogeneity in Arabidopsis thaliana
Justine Floret, Anja Linstädter, Huiyao Zhang, Vera Hesen, Swan Portalier, Lee Weinand, Gaelle Bustarret, Kirsten Bell,
Fabrice Roux, Tahir Ali, Margarita Takou, Gregor Schmitz, Stanislav Kopriva, Juliette de Meaux
bioRxiv 2026.06.27.734931; doi: https://doi.org/10.64898/2026.06.27.734931

## Funder Information Declared

Deutsche Forschungsgemeinschaft, 
https://ror.org/018mejw64, TRR341, EXC 2048/1-390686111, CRC1644

## Contact

For questions, please contact:  
Justine Floret: justine.marie.floret@gmail.com
Juliette de Meaux: jdemeaux@uni-koeln.de
