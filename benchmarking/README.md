# Validating and Benchmarking the `delineator` Python package

Matthew Heberger, matt@mghydro.com

July 2026

These scripts and files are for testing the speed and accuracy 
of the Python package `delineator` and comparison to other
watershed delineation methods. This is described in my 2026 manuscript
submitted to the Journal *Computers & Geosciences*. This document describes
the files in the `delineator` GitHub repository: https://github.com/mheberger/delineator.

My goal is to provide all the information here to repeat these experiments and to
verify the results of my validation and benchmarking.

These files are organized into four main experiments. Each folder contains
a README with more information, and the scripts and data needed to repeat the 
experiments. The folder structure is as follows:

```
📁 delineator\
└─ 📁 benchmarking\
	├─📁 data preparation for Taudem\
	├─📁 Experiment 1 - validate watershed boundaries\
	├─📁 Experiment 2 - check accuracy\
	├─📁 Experiment 3 - benchmark speed, vary basin area\
	├─📁 Experiment 4 - benchmark speed, vary input size\
	├─Benchmarking Results.xlsx
	└─README.md
```

The purpose of **Experiment 1** was to validate the correctness of the watershed
boundaries output by `delineator`. Watershed delienation with a flow-direction
raster is deterministic, so any software package should give the same result when 
it is using the same input data or flow-direction raster. Experiment 1 tests whether the watersheds created by `delineator` are the 
same as those created with two industry-standard software packages, `TauDEM` and `pysheds`. 
The MERIT-Hydro raster data required some processing before it could be used with TauDEM. 
Data processing scripts are in folder `data preparation for Taudem`.


In **Experiment 2**, the purpose was to determine whether the watersheds created with 
`delineator` are reasonable and accurate, i.e., whether the boundaries are true, 
by comparing them to US Geological Survey watershed boundaries for river gaging
stations in the continental United States. 
In one sense, this experiment is a test of the MERIT-Hydro dataset, which has a
relatively coarse resolution (~90 m) compared to publicly available US datasets
like the 10 m resolution National Elevation Dataset.


**Experiment 3** tests the claim that `delineator` is faster than other software 
packages. Here, I used Python scripts to benchmark delineation time for
three software packages:

 - TauDEM v5, https://hydrology.usu.edu/taudem/taudem5/
 - pysheds v0.4, https://github.com/pysheds/pysheds 
 - delineator v2.2.0 


We saw in Experiment 3 that the delineation time for `pysheds` and `TauDEM` 
varies depending on the size of the flow direction raster. In **Experiment 4**,
 I tested the effect of input file size more methodically. 

The results for experiments 3 and 4 are summarized in the Excel workbook 
`Delineator Benchmarking Results.xlsx`. 

Note: For the benchmarking experiments, I had to decide what to wrap in the timer. 
I am only including the "marginal" processing time for delineating a watershed, 
after the preliminaries are taken care of. For example, with TauDEM, the first step
of snapping pour points is done in batch mode for many outlets at once and outputs
a shapefile which is the input for the next step. I have not included this in the 
benchmarking. 

## Computer specs

I ran all the benchmarking experiments on my circa 2020 laptop:

```
Dell Latitude 5510 laptop
Intel(R) Core(TM) i7-10610U CPU @ 1.80GHz
4 Core(s), 8 Logical Processor(s)
32 GB RAM
Windows 64-bit operating system, x64-based processor
```
