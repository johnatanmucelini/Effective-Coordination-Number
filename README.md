# Effective Coordination Number

This *python* script calculates the effective coordination number and average bond distance from almost any atomic structure file type. It work with any atomic structuresfile readble by ase.io.read(), including structures with periodic boundary conditions.

See this reference: **J. Appl. Phys. 2011, 109, 023502**.

## For newbies:

To run this script you need *python* and some *python packages* that are detailed below. 

If you don't have python, I recommend [this tutorial](https://varhowto.com/install-miniconda-ubuntu-20-04/).

## Prerequisites:

The following python packages are prerequisites:
- **numpy**
- **scipy**
- **atomic simulation environment**

The **pandas** package is necessary to save the data of the atomic structures analyzed, which I strongly recommend.

If you employ Anaconda package management, you can install the packages with the following commands:

```bash
conda install numpy scipy pandas
conda install -c conda-forge ase
```

## Usage

Script help messange:

```
usage: ecn.py [-h] --files file1 [file2 ...] [--save_json file.json]

The script mensure the Effective Coordination Number (ECN) and the average bound distance (dav), 
for one or more atomic structures. It work with any atomic structuresfile readble by the 
ase.io.read function, including structures with periodic boundary conditions.

required arguments:
  --files file1 [file2 ...]
                        One ore more structure files (xyz, geometry.in, POSCAR, etc) to analyze.

optional arguments:
  --save_json file.json
                        If defined, the data collecte will be writen in this json file.
```


## Examples (and tests)

To analyze all the example structures, including molecules (xyz files) and a crystal (vasp-POSCAR file), and save the output:

```
$ python ecn.py --files test_case/* > output2.txt
```

To analyse all the example structures and save data found in a json file:

```bash
$ python ecn.py --files test_case/* --json_file full_analysis.json 
```

## Cite:

If you employed this methodology, please, cite **J. Appl. Phys. 2011, 109, 023502**.
