# PhyloProPy

`PhyloProPy` is a Python class designed to parse and manipulate phylogenetic profile data generated with [https://github.com/BIONF/PhyloProfile](https://github.com/BIONF/PhyloProfile). It allows for the loading, processing, and visualization of phyloprofile matrices.

## PhyloProfile Class

### Features

*   **Initialization**: Parse a PhyloProfile file and store it as a Pandas DataFrame in the PhyloProfile.
*   **Filtering**: Filter the PhyloProfile stored in the object based on a list of genes or taxonomic IDs.
*   **Slicing**: Return specific portions of the PhyloProfile as a Pandas DataFrame.
*   **Lineage slices**: Extract a slice of the PhyloProfile containing members of a lineage based on information in the NCBI Taxonomy.
*   **Binary Transformation**: Convert FAS scores to binary values for simplified analysis.
*   **Output**: Write the processed PhyloProfile to a file.

### Initialization

Create a `PhyloProfile` object by providing a path to your phyloprofile file along with other optional parameters.

```
# stores FAS foreward scores per default
pp = PhyloProfile(path='/path/to/phyloprofile')
# binary profile
pp = PhyloProfile(path='/path/to/phyloprofile', style='binary')
# order profile according to taxonomic distance to a reference species (left to right)
pp = PhyloProfile(path='/path/to/phyloprofile', reference=9606)
pp = PhyloProfile(path='/path/to/phyloprofile', reference='Mus)
```

### Filtering and Slicing

Filter or slice the phyloprofile based on genes or taxa.

```
# Filtering 
pp.filter_pp(genes=['gene1', 'gene2'], taxa=['9606', '10090'])
# Slicing
slice = pp.slice(genes=['gene1', 'gene2'])
```

### Binary Transformation

Convert the FAS scores in the phyloprofile matrix to binary values.

pythonCopy code

`pp.fas_to_binary()`

### Writing Output

Write the processed phyloprofile to a file.

pythonCopy code

`pp.write(path='./output.phyloprofile')`

### Visualization

Plot the phyloprofile as a heatmap.

pythonCopy code

`pp.plot()`

### Lineage Analysis

Extract a slice of the PhyloProfile based on a specific lineage.

pythonCopy code

`lineage_slice = pp.lineage_slice('Mammalia')`

Installation
------------

To use `PhyloProfile`, ensure you have the required dependencies installed:

bashCopy code

`pip install pandas seaborn ete3`

Clone or download this repository to your local machine, and you're ready to incorporate `PhyloProfile` into your project.

Contribution
------------

Contributions to improve `PhyloProfile` are welcome. Please fork the repository and submit a pull request with your enhancements.

License
-------

This project is licensed under the MIT License - see the LICENSE file for details.
