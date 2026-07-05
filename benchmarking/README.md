# Validating and Benchmarking the `delineator` Python package

Matthew Heberger, matt@mghydro.com
July 2026

These scripts and files are for testing the speed and accuracy 
of the Python package `delineator` and comparison to other
watershed delineation methods. This is described in my 2026 manuscript
submitted to the Journal *Computers & Geosciences*. 

My goal is to provide all the information here to repeat these experiments and to
verify the results my validation and benchmarking.

These files are organized into four main experiments. Each folder contains
a README with more information, and the scripts and data needed to repeat the 
experiments.	 

## Data Preparation TauDEM

To use use TauDEM, the MERIT-Hydro data required some preparation. See the folder 
`data preparation for Taudem` for the scripts used for this. 


## Computer specs

I ran all of the benchmarking experiments on my laptop:

```
Dell Latitude 5510 laptop
Intel(R) Core(TM) i7-10610U CPU @ 1.80GHz
4 Core(s), 8 Logical Processor(s)
32 GB RAM
Windows 64-bit operating system, x64-based processor
```

## Overview of the Experiments

### Experiment 1: Validating the correctness of the output

Watershed delienation with a flow-direction raster is deterministic, so 
any software package should give the same result when it is using the 
same input data, or flow-direction raster. 

Experiment 1 tests whether the watersheds created by `delineator` are the 
same as those created with `TauDEM` and `pysheds`. 

The MERIT-Hydro raster data required some processing for use with TauDEM. 
See the folder `data preparation for Taudem`.


### Experiment 2: Checking the accuracy of delineated watersheds

This experiment was to determine whether the watersheds created with 
`delineator` are reasonable and accurate, i.e. whether the boundaries are true, 
by comparing them to US Geological Survey watershed boundaries for river gaging
stations in the continental United States. 

In one sense, this experiment is a test of the MERIT-Hydro dataset, which has a
relatively coarse resolution (~90 m) compared to publicly-available US datasets
like the 10 m resolution National Elevation Dataset.


### Experiment 3: Benchmarking delineation time

For this experiment, I timed the watershed delineation process
using 3 software packages:

 - TauDEM v5, https://hydrology.usu.edu/taudem/taudem5/
 - pysheds v0.4, https://github.com/pysheds/pysheds 
 - delineator v2.2.0 

The results for experiments 3 and 4 are summarized in the Excel workbook 
`Delineator Benchmarking Results.xlsx`. 


### Experiment 4:  Benchmarking delineation for varying raster input size 

We saw in Experiment 3 that the delineation time for `pysheds` and `TauDEM` 
varies depending on the size of the flow direction raster. Here, I tested
the effect of input file size more methodically. 


Note: For the benchmarking experiments, I had to decide what to wrap in the timer. 
I am only including the "marginal" processing time for delineating a watershed, 
after the preliminaries are taken care of. For example, with TauDEM, the first step
of snapping pour points is done in batch mode for many outlets at once, and outputs
a shapefile which is the input for the next step. I have not included this in the 
benchmarking. 
