
import numpy as np



def reproduction (genotypes, fitnesses, new_N, interSNP_genetic_dist, genome, bounds_chr, rng):
    cc0_ = []
    c = 0
    fitnesses = fitnesses / np.sum(fitnesses)
    previous_N = len(fitnesses)
    #choice of parents
    parents = rng.choice (list(range(previous_N)), new_N, p = fitnesses)
    new_genotypes = np.empty(shape = (2 * new_N, np.sum(genome)), dtype = int)
    for index_parent, parent in enumerate(parents):
        prev_genotype = genotypes[[2* parent, 2* parent + 1], :]
        new_genotype = np.empty(shape = (2, np.sum(genome)), dtype = int)
        #two gametes
        for g in [0, 1] :
            gamete = np.repeat(None, np.sum(genome))
            for chr in range(len(genome)):
                c +=1
                bounds_chr_current = bounds_chr[chr]
                prev_chr = prev_genotype [:, bounds_chr_current[0]:bounds_chr_current[1] + 1]
                #genetic map 
                gen_map = np.cumsum(interSNP_genetic_dist[chr])
                
                #total genetic distance of the current chromosome
                tot_dist = gen_map[-1]
                
                #choosing crossing over location
                co_number = rng.poisson(tot_dist)
                if co_number == 0:
                    gamete[bounds_chr_current[0]: bounds_chr_current[1] + 1] = prev_chr[rng.binomial(1, 1/2, 1)[0], :]
                    cc0_.append(0)
                    continue



                #locations of crossing over
                co_loc = rng.uniform (0, tot_dist, co_number)
                co_loc = np.searchsorted(gen_map, co_loc) # if i in co_loc then co between i and i + 1 (python index)

                #to check the expected number of crossing over is correct
                if chr == 0:
                    cc0_.append(np.sum(co_loc == 100000))
                
                #locations of recombinations 
                rec_loc = np.unique(co_loc, return_counts=True)


                rec_loc = rec_loc[0][ rec_loc[1]%2 == 1]
                if len(rec_loc) == 0 :
                    gamete[bounds_chr_current[0]: bounds_chr_current[1] + 1] = prev_chr[rng.binomial(1, 1/2, 1)[0], :]
                    continue
                rec_loc = np.sort(rec_loc)
                #number of SNPs between successive recombinations
                rec_loc_dist = np.concatenate([[rec_loc[0] + 1], rec_loc[1:] - rec_loc[:-1]])

                #adding number of sites after last recombination
                rec_loc_dist = np.concatenate([rec_loc_dist, np.array([genome[chr]- sum (rec_loc_dist)])])
                n_rec = len(rec_loc_dist)
                if rng.uniform (0, 1) < 1/2:
                    gamete_chr_index = np.repeat(([0, 1]*(n_rec//2 + 1))[:n_rec], rec_loc_dist)
                else :
                    gamete_chr_index = np.repeat(([1, 0]*(n_rec//2 + 1))[:n_rec], rec_loc_dist)
                gamete_new_chromosome = prev_chr[gamete_chr_index, list(range(genome[chr]))]

                gamete[bounds_chr_current[0]: bounds_chr_current[1] + 1] = gamete_new_chromosome

            new_genotype[g, :] = gamete
        new_genotypes[[2*index_parent, index_parent * 2 + 1], :] = new_genotype
        
    return new_genotypes, cc0_, c

