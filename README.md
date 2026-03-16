# GenomicTools

This is a set of tools for comparative genomics, written in python. The workflow involves running OrthoFinder on a set of genome annotations, and generating homology matrices/dot plots representing pairwise comparisons of genomes. These homology matrices are stored as files with rows representing pairs of genes, one from each of the two genomes compared, that are in the same orthogroup (i.e. are homologous). We have developed a simple algorithm, APES, for identifying microsynteny using perfectly conserved gene order (which we call nanosynteny) as "anchors".

See our paper (https://doi.org/10.1101/2025.04.03.647042) for details on how we've used these tools.

## Examples:

See GenomicTools/examples/ for an example of running APES to infer microsynteny between the human and mouse genomes.

## Data formats:

Our workflow depends on a few standardized data structures, which we describe here.

### Species data:

Each genome annotation is summarized by an N X 12 "species data" array, which is arranged as follows:

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
* Rows 3-: Data 

### Chromosome information:

This is a dictionary for each species, generated from the species data summarizing the size of each chromosome.

Chromosome information (dictionary):
- 'Chromosome name': {'number': X, 'size': Y}

To explain, the number of keys in "first layer" of the dictionary is equal to the number of chromosomes in the genome. The value associated to each chromosome name is another dictionary, where 'number' is the chromosome index, and 'size' is the number of genes within the chromosome.

### Dot plot / homology matrix:

Each comparison of a pair of genomes is summarized by an N X 9 "dot plot" array. We use "homology matrix" synonymously with "dot plot" - "dot plot" seems to imply an actual plot rather than a data structure, but we use it anyway in this context. These arrays are arranged as follows:
	
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
* Rows 3-: Data
