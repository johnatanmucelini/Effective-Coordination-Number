# Effective-Coordination-Number

This *python* script calculates the effective coordination number and average bond distance from almost any atomic structure file type. It work with any atomic structuresfile readble by ase.io.read(), including structures with periodic boundary conditions.

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
conda install -c conda-forge ase
```

## Usage

Script help messange:

```bash
usage: ecn.py [-h] [--files FILES [FILES ...]] [--json_file file.json]

The script mensure the Effective Coordination Number (ECN) and the average bound distance (dav), for one or more atomic structures. It work with any atomic structuresfile readble by ase.io.read(), including structures with periodic boundary conditions.

optional arguments:
  -h, --help                 show this help message and exit
  --files FILES [FILES ...]  the molecules (xyz, geometry.in, etc) to analyze.
  --json_file file.json      the name of a json file to save the data.
```


## Examples (and tests)

To analyze all the example structures, including molecules (xyz files) and a crystal (vasp-POSCAR file), and save the output:

```
$ python ecn.py --files test_case/* > output2.txt
``

To analyse all the example structures and save data found in a json file:

```bash
$ python ecn.py --files test_case/* --json_file full_analysis.json 
```

## Cite:

If you employed this methodology, please, cite **J. Appl. Phys. 2011, 109, 023502**.
