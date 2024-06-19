import pandas as pd
import os
from ete3 import NCBITaxa
from PhyloProPy.load_phyloprofile import phyloprofile2matrix, sort_phyloprofile
from PhyloProPy.plotting_tools import plot_tsne, phylo_heatmap, dimension_reduced_phyloprofile
from PhyloProPy.logger import phyloprofile_logger
from PhyloProPy.mapping import check_taxonomy_input
import logging


###################################################################################################
# ToDo: Use lazy loading for some packages, eg. plotly, seaborn, matplotlib
###################################################################################################

class PhyloProfile():
    """
    Parse a PhyloProfile file and store it as a Pandas DataFrame.
    """
    def __init__(
        self, path='', style='fasf', from_custom=False, fasF_filter=0.0, fasB_filter=0.0, fillna=0, resolve_coorthologs=True, reference='', debug=False, silent=False, 
    ):
        """
        style: ['fasf', 'fasb', 'binary', 'orthoid', 'ncRNA'] -> How to fill cells of phyloprofile matrix, (In case of co-orthologs: Maxmimum FAS-score, List of orthoIDs)
        from_custom: bool -> Switch to True if your phyloprofile file was exported from PhyloProfile and contains a "%Spec" column
        fasF_filter: float -> Orthologs with a lower FAS-Foreward score will not be loaded into the matrix 
        fasB_filter: float -> Orthologs with a lower FAS-Backward score will not be loaded into the matrix
        fillna: char -> Fill cells without orthologs with fillna
        resolve_coorthologs: bool -> If True maintain only maximum score, if False fill cells with list of scores
        reference: int/str -> NCBI Taxonomy ID or Species name of the Seed species (of the fDOG analysis). Re-orders the columns of the matrix so that the seed species is left and the most distantly related target species is right.
        debug: bool -> More verbose
        silent: bool -> Less verbose
        """
        logger = phyloprofile_logger(debug=debug, silent=silent)
        # taxonomy
        logger.info('Reading NCBI Taxonomy')
        self.ncbi = NCBITaxa()
        #data
        if not path:  # load example data
            logger.info('No path specified. Loading example phyloprofile')
            path  = os.path.dirname(__file__) + '/data/medium.phyloprofile'
        self.matrix, self.outmatrix = phyloprofile2matrix(path, self.ncbi, style, from_custom, fasF_filter, fasB_filter, fillna, resolve_coorthologs, reference)
        self.style = style

    def to_binary(self):
        if self.style == 'fasf' or self.style == 'fasb':
            self.matrix = self.matrix.applymap(lambda x: 1 if x > 0 else 0)
        elif self.style == 'orthoid':
            self.matrix = self.matrix.astype(str).applymap(lambda x: 0 if x == '0' else 1)
        else:
            self.matrix = self.matrix
        self.style == 'binary'

    def write_csv(self, path='./output.phyloprofile'):
        with open(path, 'w') as of:
            of.write('geneID\tncbiID\torthoID\tFAS_F\tFAS_B\n')
            for geneid, row in self.outmatrix.iterrows():
                for taxid, entrylist in row.items():
                    if entrylist == 0:
                        continue
                    for entry in entrylist:
                        orthoid, fasf, fasb = entry
                        of.write(f'{geneid}\t{taxid}\t{orthoid}\t{fasf}\t{fasb}\n')

    def filter_profile(self, genes=None, taxa=None):
        """Filter the the PhyloProfile based on a list of genes or taxids. Irreversible but can be used for writing."""
        if genes:
            self.matrix = self.matrix.filter(genes, axis='index')
            self.outmatrix = self.outmatrix.filter(genes, axis='index')
        if taxa:
            taxa = [taxon if str(taxon).startswith('ncbi') else f'ncbi{taxon}' for taxon in taxa]
            self.matrix = self.matrix.filter(taxa, axis='columns')
            self.outmatrix = self.outmatrix.filter(taxa, axis='columns')

    def slice(self, genes=None, taxa=None):
        """Return a DataFrame slice of a PhyloProfile."""
        if genes:
            return self.matrix.filter(genes, axis='index')
        if taxa:
            taxa = [taxon if str(taxon).startswith('ncbi') else f'ncbi{taxon}' for taxon in taxa]
            return self.matrix.filter(taxa, axis='columns')
        return self.matrix

    def set_reference(self, reference):
        _, order = sort_phyloprofile(self.matrix, self.ncbi, reference)
        self.matrix = self.matrix[order]
        self.outmatrix = self.outmatrix[order]
         
    def print(self):
        """Print the phyloenetic profile dataframe"""
        print(self.matrix)

    def display(self):
        """Display the phyloenetic profile dataframe"""
        display(self.matrix)

    def plot(self):
        """Plot the phylogenetic profile as a simple heatmap"""
        sns.heatmap(self.matrix)

    def lineage_slice(self, lineage):
        """Return a DataFrame containing only taxa in lineage"""
        input = check_taxonomy_input(lineage, self.ncbi)
        if not input:
            logger = logging.getLogger('phyloprofile')
            logger.error(f'Could not find "{lineage}" in the NCBI Taxonomy')
            return None
     
        descendants = self.ncbi.get_descendant_taxa(input)
        return self.matrix.filter([f'ncbi{taxid}' for taxid in descendants])

    def two_d_plot(self, orient='species', taxlevel='species', update_taxonomy=False, seed=42, jitter=0.0, method='umap', scaler='None', return_as='figure', **kwargs):
        """
        method: ['umap', 'PCA', 'tSNE', 'MDS']
        
        Project phylogenetic profile into 2D space and scatterplot.
        Accepts **kwargs of plotly.express.scatter
        Returns as a plotly express "figure" or as the dimension-reduced "dataframe"
        """
        logger = logging.getLogger('phyloprofile')
        if orient == 'species':
            transpose = True
        elif orient == 'genes':
            transpose = False
        else:
            raise ValueError(f'Unknown orientation "{orient}". Choose "species" or "genes".')
        
        # reduce dimension
        logger.info(f'Reducing dimensions')
        red_df = dimension_reduced_phyloprofile(
            self.matrix, taxlevel, self.ncbi, 
            update_taxonomy=update_taxonomy, method=method, jitter=jitter, scaler=scaler, transpose=transpose, seed=seed, 
            **kwargs
        )
        if return_as == 'dataframe':
            return red_df
        elif return_as == 'figure':
            logger.info(f'Generating plot')
            fig = plot_tsne(red_df, **kwargs)
            logger.info(f'Done')
            return fig
        else:
            raise ValueError(f'Cannot return result as "{return_as}". Choose "figure" or "dataframe"')

    def plot(self, clustermethod='average', names=True, **kwargs):
        """Plot phylogenetic profile as simple heatmap."""
        if names:
            taxids = [int(taxid.replace('ncbi', '')) for taxid in self.matrix.columns]
            taxid2name = self.ncbi.get_taxid_translator(taxids)
            taxid2name = {f'ncbi{taxid}': name for taxid, name in taxid2name.items()}
            return phylo_heatmap(self.matrix.rename(columns=taxid2name), clustermethod, **kwargs)
        else:
            return phylo_heatmap(self.matrix, clustermethod, **kwargs)

    def genes(self):
        return self.matrix.index

    def taxa(self, return_as='int'):
        if return_as == 'int':
            return [int(taxon.replace('ncbi', '')) for taxon in self.matrix.columns]
        # elif return_as == 'names':
            
        # return self.matrix.columns


    def write_orthoxml():
        pass
        


