# Surface Atoms Analysis

This *python* script analyzes which atoms of a molecule (in xyz file, for instance) are exposed to the vacuum and is exposed area. 

![example](https://raw.githubusercontent.com/johnatanmucelini/surface_analysis_standalone/master/figure.png?token=AF3ZDJ3EMVHYIY3HY4G5773AQCTUS)

The atoms exposed to the vacuum could be considered *surface atoms* and often play important roles in the physical-chemical properties of small clusters/nanoclusters/nanostructures.
I strongly recomend this analysis for the investigations of structures with more than 20 atoms. 

For a detailed description of the methodology, I suggest the Section V of [this document](https://pubs.acs.org/doi/suppl/10.1021/acs.jpcc.9b09561/suppl_file/jp9b09561_si_001.pdf).

It also calculates the effective coordination number: **J. Appl. Phys. 2011, 109, 023502**.

## For newbies:

To run this script you need *python* and some *python packages* that are detailed below. 

If you don't have python, I recommend [this tutorial](https://varhowto.com/install-miniconda-ubuntu-20-04/).


## Prerequisites:

The following python packages are prerequisites:
- **numpy**
- **scipy**
- **atomic simulation environment**
- **scikit-learn**

The **pandas** package is necessary to  save the data of the atomic structures analyzed, which I strongly recommend.

If you employ Anaconda package management, you can install the packages with the following commands:
```bash 
conda install numpy scipy pandas scikit-learn 
conda install -c conda-forge ase
```


## Usage

Script help messange:

```bash
usage: surface_analysis_std.py [-h] --mol MOL [--r_adatom R_ADATOM] [--r_atoms R_ATOMS] [--ssamples SSAMPLES] [--sp_file SP_FILE] [--save_json SAVE_JSON]

This script calculates the atoms' exposition to the vacuum.

optional arguments:
  -h, --help            show this help message and exit
  --mol MOL             the molecules (xyz, geometry.in, etc) to analyze.
  --r_adatom R_ADATOM   the radius for the adatom radius. (Default=1.1)
  --r_atoms R_ATOMS     the radius for the mol atoms. (Default=dav/2)
  --ssamples SSAMPLES   the number of points to distribute over each atom. (Default=1000)
  --sp_file SP_FILE     if defined (ex: surf.xyz), the surface points found will be printed in this file.
  --json_file SAVE_JSON if defined (ex: dados.json), all the data will be saved in a json file.
```


## Examples (and tests)

To analyze a single molecule:

```
$ python surface_analysis.py --mol test_case/Ag_30_g_283_.xyz > output.txt
```

To analyse all molecules and save the data in a json file:

```bash
$ python surface_analysis.py --mol "$(ls test_case)" --json_file full_analysis.json > output2.txt
```

## Cite us, please! :)

If you employed this methodology, please, cite **J. Phys. Chem. C 2020, 124, 1, 1158–1164**.
