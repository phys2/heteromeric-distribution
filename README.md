# Homo- and heteromeric distribution of TRPC channels

Some proteins like TRPC1, TRPC4 and TRPC5 exist in various homomeric and heteromeric 
configurations. 

This script calculates the frequency of each configuration by building and 
solving systems of linear equations from protein abundances taken from an 
AP-MS dataset (where APs targeting one of C1/C4/C5 were chained in multiple 
combinations).

## Sample input data

There is a sample abundance data file in the `sample-data` directory. 
The file structure is explained in the file itself.

## How to run

Simply pass the abundance file:

    heteromers.py sample-data/TRPC-aundances.tsv

You can save both the results and the linear equations to a file:

    heteromers.py --out results.tsv --out-systems systems.tsv sample-data/TRPC-aundances.tsv

More options are available, see help with:

    heteromers.py --help
