# Benchmarking Experiment 3: Comparing delineation time to other software packages

For this experiment, I delineated a set of watersheds of varying sizes. For the
locations, I used the catalog of the GRDC, the Global Runoff Data Center. 

## Sampling

In working with various watershed delineation software, I noticed that it 
generally takes longer when you are dealing with large files. I wanted to 
control for this somewhat, so I decided to draw samples from 3 of the major
basins, which I call megabasins. I set out to draw 30 samples from each, but
the Iceland basin (#27), does not have that many gages, so I ended up with 22 
for that basin and a total sample size of 82 watershed outlets.

I also wished to create watersheds of various sizes. So I stratified the 
gages into 10 groups by area, with the script `sampler.py`. 

I also wished to test the performance on very large watersheds. Based on my
experience running the Global Watersheds web app, these are very popular, 
with many requests for the Amazon, Nile, and Mississippi every day. So I 
manually selected the 3 gages in the Amazon basin (megabasin 62) for inclusion. 
I also added one ungaged point, the Amazon River at its mouth, to include the
world's largest watershed.


## Watershed Delineation Scripts

The delineation process is controlled with these scripts:

- `delineator_timer.py`
- `pysheds_timer.py`
- `taudem_timer.py`


TauDEM is a Windows executable. The script `taudem_timer.py` is a test harness
that invokes TauDEM using the `subprocess` module. The delineation time that
I recorded and report here does not include the time to read input geodata, 
snap pour points, or to write output data. 

One of the features of TauDEM is the ability to run it in parallel. To make 
the benchmarking comparison fair, I tried to run TauDEM as fast as possible
on my laptop. My understanding is that, when running TauDEM with `mpiexec`, 
you should set `-n` to the number of physical CPU cores on your machine 
(not logical/hyperthreaded cores). My mid-range Dell laptop shows (from running `msinfo32`):

Processor: Intel(R) Core(TM) i7-1061OU CPU @ 1.80 GHz, 2304 Mhz, 4 Core(s), 8 Logical Processor(s)

Similarly, the script `pysheds_timer.py` benchmarks the time required
to delineate individual watersheds using `pysheds`. Here, the timer
does not include the time to read the flow direction and flow 
accumulation raster data files, and reports the time required to
create the watershed polygon, but not saving it to disk. 

## Results

I compiled the results of this experiment and made plots for the manuscript
in `Benchmarking Results.xlsx`.

