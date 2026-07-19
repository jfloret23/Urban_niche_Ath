# Swan Portalier 07/2026
# Determine the p-value threshold for the target FWER

import numpy as np
import pandas as pd
import sys

### Command-line arguments ###
#paths
#input paths
#path of the files with pvalues calculated from simulation output (pvals_{replicate}.csv)
p_vals_path = sys.argv[1]

target_fwer = sys.argv[2]
####################################

replicates = np.arange(1, 500)

finished_replicates = []
p_val_array = []
for replicate in replicates: 
    p_val_replicate = pd.read_csv(f"{p_vals_path}/pvals_{replicate}.csv", header = None)[0].tolist()
    p_val_array.append(p_val_replicate)
    finished_replicates.append (replicate)


p_val_array = np.array(p_val_array)



#getting distribution of mins p-values
mins_array = np.min(p_val_array, axis = 1)

#threshold that yields the target FWER
thr = np.quantile(a = mins_array, q = target_fwer)
print(f"The threshold that yields a FWER of {target_fwer} is {thr}")

#-log10(found threshold)
thrlog = -np.log10(np.quantile(a = mins_array, q = target_fwer))
print(f"-log10(threshold that yields a FWER of {target_fwer}) is {thrlog}")

