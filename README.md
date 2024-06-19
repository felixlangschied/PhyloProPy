# PhyloProPy

`PhyloProPy` is a Python package designed to parse and manipulate phylogenetic profiles generated with [PhyloProfile](https://github.com/BIONF/PhyloProfile). It allows for the loading, processing, and visualization of phylogenetic profiles as 2D-matrices.

## Features

- [**PhyloProfile Python Class**](#phyloprofile-class):
  - **I/O**: Read and write PhyloProfile files
  - **Taxonomic order:** Order the PhyloProfile matrix according to the taxonomic distance to a reference species
  - **Filtering**: Filter the PhyloProfile stored in the object based on a list of genes or taxonomic IDs
  - **Slicing**: Return specific portions of the PhyloProfile as a Pandas DataFrame
  - **Lineage slices**: Extract a slice of the PhyloProfile containing members of a lineage based on information in the NCBI Taxonomy
- [**Phylo-tSNE**](#standalone-phylo-tsne):
  - **Visualization**: Project (large) phylogenetic profiles into 2D space
  - **Color lineages**: Label datapoints according to a taxonomic level

## PhyloProfile Class

### Initialization

Create a `PhyloProfile` object by providing a path to your phyloprofile file along with other optional parameters.
```
# stores FAS foreward scores per default
pp = PhyloProfile(path='/path/to/PhyloProPy/data/medium.phyloprofile')

# binary profile
pp = PhyloProfile(path='/path/to/PhyloProPy/data/medium.phyloprofile', style='binary')

# store entries from the OrthoID column
pp = PhyloProfile(path='/path/to/PhyloProPy/data/medium.phyloprofile', style='orthoid')

# order profile according to taxonomic distance to a reference species (left to right)
pp = PhyloProfile(path='/path/to/PhyloProPy/data/medium.phyloprofile', reference=9606)
pp = PhyloProfile(path='/path/to/PhyloProPy/data/medium.phyloprofile', reference='Mus musculus')
pp.set_reference('Homo_sapiens')

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

### Visualization

Project large Phyloprofiles to 2D using UMAP 

```
# return plotly figure
fig = pp.two_d_plot(orient='species', taxlevel='kingdom')
fig.show()

# return dataframe
umap_df = pp.two_d_plot(orient='genes', return_as='dataframe')
```

### Binary Transformation

Convert the FAS scores in the phyloprofile matrix to binary values.
```
pp.fas_to_binary()
```

### Working with the NCBI Taxonomy 

PhyloProPy uses the [NCBI Taxonomy functionality of the ETE3 toolkit](http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html) under the hood. Use it for even more control over your PhyloProfile object.

```
ncbi = pp.ncbi
name2taxid = ncbi.get_name_translator(['Homo sapiens', 'primates'])
```

### Writing Output

Write the processed phyloprofile to a file.
```
pp.write(path='./output.phyloprofile')
```









