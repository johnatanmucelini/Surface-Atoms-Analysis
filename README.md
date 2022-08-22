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


## The code

#### Required packages

- numpy
- scipy
- atomic simulation environment
- scikit-learn

If you employ Anaconda package management, you can install the packages with the following commands:
```bash 
conda install numpy scipy pandas scikit-learn 
conda install -c conda-forge ase
```

#### Running the code

To analyze a single molecule and save the script output:

```
$ python surface_analysis.py --mol test_case/Ag_30_g_283_.xyz > output.txt
```

To analyse all molecules and save the data in a json file:

```bash
$ python surface_analysis.py --mol ./test_case/* --save_json full_analysis.json
```

Please, check the code help messange to see all parameters and their descriptions.

## Cite me, please! :)

If you employed this methodology, please, cite **J. Phys. Chem. C 2020, 124, 1, 1158–1164**.
