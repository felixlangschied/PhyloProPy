# PhyloProPy

`PhyloProPy` is a Python package designed to parse and manipulate phylogenetic profiles generated with [PhyloProfile](https://github.com/BIONF/PhyloProfile). It allows for the loading, processing, and visualization of phylogenetic profiles as 2D-matrices.

## Features

- [**PhyloProfile Python Class**](#phyloprofile-class):
  - **Filtering**: Filter the PhyloProfile stored in the object based on a list of genes or taxonomic IDs.
  - **Slicing**: Return specific portions of the PhyloProfile as a Pandas DataFrame.
  - **Lineage slices**: Extract a slice of the PhyloProfile containing members of a lineage based on information in the NCBI Taxonomy.
- [**Phylo-tSNE**](#phylo-tsne):
  - **Visualization**: Project (large) phylogenetic profiles into 2D space.
  - **Color lineages**: Label datapoints given a taxonomic level

## PhyloProfile Class

### Initialization

Create a `PhyloProfile` object by providing a path to your phyloprofile file along with other optional parameters.

```
# stores FAS foreward scores per default
pp = PhyloProfile(path='/path/to/phyloprofile')

# binary profile
pp = PhyloProfile(path='/path/to/phyloprofile', style='binary')

# order profile according to taxonomic distance to a reference species (left to right)
pp = PhyloProfile(path='/path/to/phyloprofile', reference=9606)
pp = PhyloProfile(path='/path/to/phyloprofile', reference='Mus musculus')

```

### Filtering and Slicing

Filter or slice the phyloprofile based on genes or taxa.

```
# Filtering 
pp.filter_pp(genes=['gene1', 'gene2'], taxa=['9606', '10090'])

# Slicing
slice_dataframe = pp.slice(genes=['gene1', 'gene2'])

# Use slicing to generate a dataframe copy of the profile stored in the PhyloProfile class
slice_dataframe = pp.slice()
```

### Lineage Analysis

Extract a slice of the PhyloProfile based on a specific lineage

```
lineage_slice = pp.lineage_slice('Metazoa')
```

### Binary Transformation

Convert the FAS scores in the phyloprofile matrix to binary values.

```
pp.fas_to_binary()
```

### Writing Output

Write the processed phyloprofile to a file.

```
pp.write(path='./output.phyloprofile')
```

## Phylo-tSNE








