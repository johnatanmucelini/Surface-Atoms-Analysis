# Surface Atoms Analysis

This *python* script analyzes which atoms of a molecule (in xyz file, for instance) are exposed to the vacuum and is exposed area. 

![example](https://raw.githubusercontent.com/johnatanmucelini/surface_analysis_standalone/master/figure.png?token=AF3ZDJ3EMVHYIY3HY4G5773AQCTUS)

The atoms exposed to the vacuum could be considered *surface atoms* and often play important roles in the physical-chemical properties of small clusters/nanoclusters/nanostructures.
I strongly recomend this analysis for the investigations of structures with more than 20 atoms. 

For a detailed description of the methodology, I suggest the Section V of [**this document**](https://pubs.acs.org/doi/suppl/10.1021/acs.jpcc.9b09561/suppl_file/jp9b09561_si_001.pdf).

This scripy also calculates the effective coordination number: **J. Appl. Phys. 2011, 109, 023502**.

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

```
usage: surface_analysis.py [-h] --mol mol1 [mol2 ...] [--r_adatom val] [--ssamples val] [--r_atoms [val1 [val2 ...]]] [--save_surf file.xyz] [--save_json file.json]

This script calculate the atoms exposition to the vacuum.

required arguments:
  --mol mol1 [mol2 ...]  One or more molecular files (xyz, geometry.in, etc) to analyze.

optional arguments:
  --r_adatom val         The radius of the adatom. (Default=1.1)
  --ssamples val         The (approximately) number of points to distribute over each atom. (Default=1000)
  --r_atoms [val1 [val2 ...]]
                         This flag controls the radii of the atoms in molecular files: If not defined, the atomic radii will defined as half of the average bond distance. If a single float value was
                         priveded, it will be the radius for every atom on the molecular files. If N float values were provided, were N is the number of atoms in the molecular files, they will be the
                         radius for each atom following the sequence of atoms in the molecular file. (Default=dav/2)
  --save_surf file.xyz   If defined, the position of the surface points found are writen in this xyz file as H atoms.
  --save_json file.json  If defined, all the collected data are writen in this json file.
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
