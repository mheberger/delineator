"""
Draw random samples based on quintiles.

Assumes the population is sorted. 
Divide the population into 5 strata (quintiles)
and draw an equal number of samples from each. 
"""


import random
import math

# Parameters
population_size = 882
n_quintiles = 5
samples_per_quintile = 4

# Calculate quintile sizes
base_size = population_size // n_quintiles
remainder = population_size % n_quintiles

# Create ranges for each quintile
ranges = []
start = 1
for i in range(n_quintiles):
    # Add an extra item to early quintiles if we have remainder
    quintile_size = base_size + (1 if i < remainder else 0)
    end = start + quintile_size - 1
    ranges.append((start, end))
    start = end + 1

# Sample from each range
for i, (start, end) in enumerate(ranges, 1):
    samples = random.sample(range(start, end + 1), samples_per_quintile)
    print(f"Quintile {i}: {sorted(samples)}")