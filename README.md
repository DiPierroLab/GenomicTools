# GenomicTools

This is a set of tools for comparative genomics, written in python. The workflow involves running OrthoFinder on a set of genome annotations, and generating homology matrices/dot plots representing pairwise comparisons of genomes. These homology matrices are stored as files with rows representing pairs of genes, one from each of the two genomes compared, that are in the same orthogroup (i.e. are homologous). We have developed a simple algorithm, APES, for identifying microsynteny using perfectly conserved gene order (which we call nanosynteny) as "anchors".

See our paper (https://doi.org/10.1101/2025.04.03.647042) for details on how we've used these tools.

## Installation:

Just clone this repository into your Python site packages directory. For me (using MacOS), this is ~/.local/lib/python3.9/site-packages

### Dependencies:

This package uses Python 3.9, and we used Anaconda 3 to manage package dependencies. These dependencies are:
* NumPy 1.26.4
* SciPy 1.11.4
* Matplotlib 3.5.1
* NetworkX 2.7.1
* Biopython 1.79

## Examples:

See GenomicTools/examples/ for an example of running APES to infer microsynteny between the human and mouse genomes.

## Data formats:

Our workflow depends on a few standardized data structures, which we describe here.

### Species data:

Each genome annotation is summarized by an N X 12 "species data" array, which is saved in .csv (comma separated) files as follows:

Species data (comma separated):
* Row 1: # Genus_species
* Row 2: # column labels
	- Column 1: chromosome name
	- Column 2: chromosome number
	- Column 3: relative index (the gene index, starting with 1 for each chromosome)
	- Column 4: absolute index (the gene index, unique to each gene in the genome)
	- Column 5: orthogroup (an integer)
	- Column 6: gene name
	- Column 7: MAKER gene label
	- Column 8: empty
	- Column 9: empty
	- Column 10: empty
	- Column 11: empty
	- Column 12: empty
* Rows n (n >= 3): Data 

### Chromosome information:

This is a dictionary for each species, generated from the species data summarizing the size of each chromosome.

Chromosome information (dictionary):
- 'Chromosome name': {'number': X, 'size': Y}

To explain, the number of keys in "first layer" of the dictionary is equal to the number of chromosomes in the genome. The value associated to each chromosome name is another dictionary, where 'number' is the chromosome index, and 'size' is the number of genes within the chromosome.

### Dot plot / homology matrix:

Each comparison of a pair of genomes is summarized by an N X 9 "dot plot" array. We use "homology matrix" synonymously with "dot plot" - "dot plot" seems to imply an actual plot rather than a data structure, but we use it anyway in this context. These arrays are saved as follows in .csv (comma separated) files:
	
Dot plot files (comma separated):
* Row 1: # Genus_species A, Genus_species B
* Row 2: # column labels
	- Column 1: chromosome name A
	- Column 2: relative index A
	- Column 3: chromosome name B
	- Column 4: relative index B
	- Column 5: empty
	- Column 6: empty
	- Column 7: empty
	- Column 8: empty
	- Column 9: empty
* Rows n (n >= 3): Data

### Synteny blocks:

From a dot plot, we use APES to estimate a set of "nanosynteny" (perfectly conserved gene order) or microsynteny blocks. We represent synteny blocks as a list of N X 9 arrays, where N varies from block to block. These arrays are saved as follows in .csv files:
	
Dot plot files (comma separated):
* Row 1: # Genus_species A, Genus_species B
* Row 2: # column labels
	- Column 1: chromosome name A
	- Column 2: relative index A
	- Column 3: chromosome name B
	- Column 4: relative index B
	- Column 5: empty
	- Column 6: empty
	- Column 7: empty
	- Column 8: empty
	- Column 9: empty
* Row 2n+1 (n >= 1): # (indicates a new synteny block)
* Row 2n+2 (n >= 1): synteny block data
