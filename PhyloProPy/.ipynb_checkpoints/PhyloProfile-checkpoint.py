import pandas as pd
import seaborn as sns
from ete3 import NCBITaxa
from PhyloProPy.load_phyloprofile import phyloprofile2matrix
from PhyloProPy.logger import phyloprofile_logger
from PhyloProPy.mapping import check_taxonomy_input
import logging


class PhyloProfile():
    """
    Parse a PhyloProfile file and store it as a Pandas DataFrame.
    """
    def __init__(
        self, path, style='fasf', from_custom=False, fasF_filter=0.0, fasB_filter=0.0, reference='', debug=False, silent=False
    ):
        """
        style: ['fasf', 'fasb', 'binary', 'orthoid'] -> How to fill cells of phyloprofile matrix
        from_custom: bool -> Switch to True if your phyloprofile file was exported from PhyloProfile and contains a "%Spec" column
        fasF_filter: float -> Orthologs with a lower FAS-Foreward score will not be loaded into the matrix
        fasB_filter: float -> Orthologs with a lower FAS-Backward score will not be loaded into the matrix
        reference: int/str -> NCBI Taxonomy ID or Species name of the Seed species (of the fDOG analysis). Re-orders the columns of the matrix so that the seed species is left and the most distantly related target species is right.
        debug: bool -> More verbose
        silent: bool -> Less verbose
        """
        logger = phyloprofile_logger(debug=debug, silent=silent)
        logger.info('Reading NCBI Taxonomy')
        self.ncbi = NCBITaxa()
        self.matrix, self.outmatrix = phyloprofile2matrix(path, self.ncbi, style, from_custom, fasF_filter, fasB_filter, reference)

    def fas_to_binary(self):
        # ToDO: Check if matrix numeric
        self.matrix = self.matrix.applymap(lambda x: 1 if x > 0 else 0)
        # for orthoid matrix
        # self.matrix = self.matrix.groupdf.applymap(lambda x: False if x == '0' else True)

    def write(self, path='./output.phyloprofile'):
        with open(path, 'w') as of:
            of.write('geneID\tncbiID\torthoID\tFAS_F\tFAS_B\n')
            for geneid, row in self.outmatrix.iterrows():
                for taxid, entry in row.items():
                    if entry == 0:
                        continue
                    orthoid, fasf, fasb = entry
                    of.write(f'{geneid}\t{taxid}\t{orthoid}\t{fasf}\t{fasb}\n')

    def filter_pp(self, genes=None, taxa=None):
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
        


def main():
    testpath = '/home/felixl/PycharmProjects/cellulases/data/metazoa_cellulase.phyloprofile'
    pp = PhyloProfile(testpath, from_custom=True, style='binary')
    metazoa = pp.lineage_slice('Horst')
    print(metazoa)
    #pp.print()

if __name__ == "__main__":
    main()
