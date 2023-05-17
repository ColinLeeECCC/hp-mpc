# Massively Parallel Hierarchical Agglomerative Clustering

This program performs the hierarchical agglomerative clustering algorithm on large datasets by spreading it over many nodes, allowing access to much larger memory. The code can also store periodic checkpoints, to allow for the program to be stopped and restarted in case the problem size is too big to be completed in a single run with the maximum allowed wallclock time. The algorithm is divided into two parts - first, the dissimilarity matrix is computed (using either the $R^2$, Euclidean distance (EuD), or the product of the two EuD$\times R^2$ as the metric) and this is stored to the filesystem. The matrix is then read in and the PQueues (which allow for quickly finding the lowest value in the matrix) are populated with this data and then the clustering proceeds stepwise, reducing the dissimilarity matrix using the chosen linkage function.

All parameters are specified in the `arglist.lst` file (by default; the name of the `.lst` file can be specified on the command line with the `f` option. The other option is `a` which specifies the data is in the aircraft format as opposed to the gridded data format. Both formats (currently) require NetCDF input files. The command line options do not have dashes because it breaks ord_soumet :upside_down_face:

To run the algorithm, make sure your NetCDF data is in the correct format and your `arglist.lst` file is in order.

## Arglist parameters file

### Lst file lines for gridded model data
1. input directory
2. forecast hour
3. start timestep
4. start date `YYYY MM DD`
5. end date `YYYY MM DD`
6. metric (0: $R^2$, 1: EuD, 2: EuD$\times R^2$
7. linkage (0-7)
8. grid selection
9. field
10. human-readable description
11. output directory
12. stop after calcluating dissimilarity matrix ([F]/T)
13. checkpoint frequency
14. stop after checkpointing

### Lst file lines for observational data
1. input netcdf file
6. metric (0: $R^2$, 1: EuD, 2: EuD$\times R^2$
7. linkage (0-7)
8. filter selection
9. field
11. output directory
12. stop after calcluating dissimilarity matrix ([F]/T)
13. checkpoint frequency
14. stop after checkpointing

## Input Data

All input data must be in NetCDF format.

### For gridded model output data

Use the fst2nc utility to convert your model output files directly to NetCDF (to save space, only include the fields you want to cluster), using the same filename as GEM-MACH output the fst files but adding the extension `.netcdf4.compressed`. At present, the program is searching for the files as `YYYYMMDDFF_SSSSSSp` and skiping by 60 SSSSSS each file, as this is how the 2.5km GEM-MACH outputs data.

### For observations

Use the included utitly csv_to_netcdf. First put your input observation data into a CSV file where the first row contains only the name of the field you want in your NetCDF file. Subsequent rows represent the datapoints you want to cluster and the observations are the columns along those rows. All metadata should be left out, do not include any kind of row or column headers or labels, just the measurements.

For example, if you want to cluster a monitoring network by NO$_2$ timeseries, each row will (from 2 on) will be the NO2 concentration timeseries for a single station.

If you want to do factor analysis instead, each row will represent a single time-point and the columns will be the concentrations of different species or the mass-spec signals.

## Compiling

On ppp5 or 6, simply run `./compile_csv2netcdf.sh` and `./compile.sh` which will generate the binaries you need.

## Running

You can simply call `./run_mpi.sh` which will use ord_soumet to submit a job to the compute cluster with reasonable defaults. You will probably want to tinker with the wallclock time, and the MPI / OpenMP configuration for your individual problem size. Command line parameters to the program can be specified after the `--`

For example, to get 4.5 hours of wallclock time on 8 nodes with 10 MPI tasks and 8 OpenMP threads per node:
```./run_mpi.sh -wtime 270 -nodesx 8 -nodesy 10 -ompnum 8 -- my-arglist.lst```

For observation files, you will also need to specify the `a` flag and an arglist name.
```./run_mpi.sh -wtime 270 -nodesx 200 -nodesy 80 -ompnum 1 -- a arglist_obs.lst```

# For more help

Contact Colin Lee <colin.lee@ec.gc.ca>