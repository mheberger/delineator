"""
Used this quick script to get a stratified sample of the gages in a region.
I wanted to make sure to get samples across a range of watershed areas, 
so I split the population of gages into 10 groups or strata and drew an 
equal number of samples from each. 
"""

import numpy as np


def random_stratified_sample_indices(n, q, samples_per_stratum):
    strata_size = n // q
    strata_indices = np.arange(0, n, strata_size)

    indices = []
    for i in range(q):
        stratum_start = strata_indices[i]
        stratum_end = strata_indices[i + 1] if i < q - 1 else n

        stratum_indices = np.arange(stratum_start, stratum_end)
        stratum_sample_indices = np.random.choice(stratum_indices, samples_per_stratum, replace=False)
        indices.extend(stratum_sample_indices)

    return indices


# Usage: set the size of the population, number of strata, and run
population_size = 214  # Size of the population
num_strata = 10  # Number of strata
samples_per_stratum = 2  # Number of samples to draw from each stratum

# Get the indices to use in lookups -- this is your stratified random sample
sample_indices = random_stratified_sample_indices(population_size, num_strata, samples_per_stratum)
print("Sample Indices:", sample_indices)
