# RTE Stability and Noise-Assessment Workflow

This repository provides a reproducible Python workflow for evaluating the numerical stability and high-frequency noise behaviour of the Reduction-to-Equator (RTE) transformation applied to aeromagnetic grids.

The workflow is intended for low-latitude magnetic data, where Fourier-domain directional transformations can be sensitive to magnetic inclination, declination, grid sampling, and short-wavelength noise. It is written as a general guide that can be adapted to any aeromagnetic dataset supplied as regular XYZ grids.

## Purpose

The repository provides scripts to:

1. Regrid total magnetic intensity (TMI) point data to a regular 50 m grid using a continuous-curvature/minimum-curvature-type workflow.
2. Compute Fourier-domain RTE transfer-function gain statistics.
3. Estimate strict and robust RTE condition numbers.
4. Compare high-frequency variability before and after RTE.
5. Produce CSV summary tables for reproducible reporting.

## Repository structure

```text
benue-rte-stability-assessment/
|
├── README.md
├── LICENSE
├── requirements.txt
├── environment.yml
├── .gitignore
|
├── scripts/
│   ├── 01_grid_TMI_50m_minimum_curvature.py
│   ├── 02_rte_condition_number.py
│   └── 03_rte_noise_table.py
|
└── outputs/
    └── RTE_stability_summary_values.csv
```

## Input data format

The scripts expect XYZ CSV files with the following column headings:

```text
X,Y,Z
```

where:

- `X` = projected easting coordinate, preferably in metres
- `Y` = projected northing coordinate, preferably in metres
- `Z` = magnetic anomaly value in nT

For consistency, input grids should be in a projected coordinate reference system such as UTM. The scripts convert non-numeric values in `Z` to `NoData`, which is useful when exported grids contain dummy values.

## Data availability note

The raw aeromagnetic grids are not included in this repository because geophysical survey data may be subject to institutional, governmental, or commercial restrictions. Users should supply their own TMI and RTE grids in XYZ CSV format.

This repository includes only:

- processing scripts
- environment files
- reproducibility instructions
- derived summary values

## Software environment

Create the Conda environment using:

```bash
conda env create -f environment.yml
conda activate maggrid
```

Alternatively, install the required packages manually:

```bash
conda create -n maggrid python=3.11 -y
conda activate maggrid
conda install -c conda-forge pandas numpy scipy pygmt netcdf4 h5netcdf -y
```

For the condition-number and noise-assessment scripts only, the minimum required packages are:

```bash
conda install -c conda-forge pandas numpy scipy -y
```

`pyGMT` is required only for the optional regridding step.

## Workflow

### 1. Optional regridding of TMI data

Run:

```bash
python scripts/01_grid_TMI_50m_minimum_curvature.py --input C:/data/TMI.csv --out-csv C:/data/TMI_50m_minimum_curvature.csv --out-nc C:/data/TMI_50m_minimum_curvature.nc --cell-size 50
```

This script applies `pyGMT.blockmean` followed by `pyGMT.surface` to produce a regular grid. The procedure provides a continuous-curvature/minimum-curvature-type grid from irregular or denser point data.

### 2. RTE transfer-function condition number

Run:

```bash
python scripts/02_rte_condition_number.py --input C:/data/TMI_50m_minimum_curvature.csv --output C:/data/RTE_condition_number_report.csv --cell-size 50 --inclination -8.94 --declination 0.21
```

This script computes the RTE transfer-function gain in the Fourier domain over the non-zero wavenumber domain.

The robust condition number is computed as:

```text
kappa_robust = max(G) / P1(G)
```

where:

- `G` is the RTE transfer-function gain
- `P1(G)` is the first-percentile gain value

The robust condition number is less sensitive to isolated near-zero directional gains than the strict maximum/minimum condition number.

### 3. High-frequency noise comparison

Run:

```bash
python scripts/03_rte_noise_table.py --tmi C:/data/TMI_50m_minimum_curvature.csv --rte C:/data/RTE_50m.csv --output C:/data/RTE_noise_assessment_table.csv --window-size 15
```

High-frequency variability is estimated as the standard deviation of the residual obtained after subtracting a moving-window local mean.

For a 50 m grid, the default 15-cell window corresponds to:

```text
15 × 50 m = 750 m
```

The RTE/TMI high-frequency noise ratio is computed as:

```text
noise_ratio = HFSD_RTE / HFSD_TMI
```

where:

- `HFSD_TMI` is the high-frequency standard deviation before RTE
- `HFSD_RTE` is the high-frequency standard deviation after RTE

A ratio greater than 1 indicates relative amplification of high-frequency variability after RTE. A ratio less than 1 indicates no high-frequency noise amplification under this diagnostic.

## Example summary values

The example output included in this repository reports the following values for one aeromagnetic grid from the Benue Trough, Nigeria:

| Parameter | Value |
|---|---:|
| Magnetic inclination | -8.94° |
| Magnetic declination | 0.21° |
| Maximum RTE gain | 1.00 |
| Median RTE gain | 0.988 |
| Robust condition number, max/P1 | 7.84 |
| Valid overlapping cells | 17,777,198 |
| Valid overlap | 97.03% |
| TMI high-frequency standard deviation | 15.69 nT |
| RTE high-frequency standard deviation | 14.95 nT |
| RTE/TMI high-frequency noise ratio | 0.953 |

These values are provided only as an example of how the workflow output can be summarized. Users should recompute the metrics using their own magnetic grids, RTE parameters, and grid spacing.

## Suggested reporting items

When reporting results from this workflow, include:

- magnetic inclination and declination
- target inclination and declination used for RTE
- grid cell size
- number and percentage of valid overlapping cells
- RTE gain minimum, percentiles, median, and maximum
- strict and robust condition numbers
- TMI and RTE high-frequency standard deviations
- RTE/TMI high-frequency noise ratio

## License

This repository is released under the license specified in the `LICENSE` file.
