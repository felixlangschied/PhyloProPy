import pandas as pd
import logging
from PhyloProPy.mapping import check_taxonomy_input


def order_taxa(tree, reference):
    """Order taxa according to a tree object from the ete3 package."""
    def initial_list(t, reference):
        copylist = []
        node = t.search_nodes(name=reference)[0]
        while node:
            for subnode in node:
                copylist.append(subnode.name)
            node = node.up
        return copylist


    def order_species(ilist):
        specorder = []
        for name in ilist:
            if name not in specorder:
                specorder.append(name)
        return specorder
    
    
    def parse_mappingfile(path):
        n2t = {}
        with open(path) as sh:
            for line in sh:
                taxid, name = line.strip().split()
                n2t[name] = f'ncbi{taxid}'
                if len(name.split('_')) > 2:
                    shortname = '_'.join(name.split('_')[:2])
                    n2t[shortname] = f'ncbi{taxid}'
        return n2t

    ##################################################################################################
    initial = initial_list(tree, reference)
    specorder = order_species(initial)
    ordernames = [f'ncbi{taxid}' for taxid in specorder]

    return ordernames


def sort_phyloprofile(df, ncbi, reference):
    # logging
    logger = logging.getLogger('phyloprofile')
    logger.info(f'Reordering matrix according to "{reference}"')

    # find all taxids
    taxids = [int(taxid.replace('ncbi', '')) for taxid in df.columns]
    tree = ncbi.get_topology(taxids)
    
    # check format of reference
    reference = check_taxonomy_input(reference, ncbi)
    if not reference:
        logger.warning(f'Could not map {reference} to exactly one NCBI taxonomy ID. Skipping Ordering')
        logger.warning(f'Received: {name2taxid}')

    # check that reference is valid
    if not reference in taxids:
        logger.warning(f'Could not find {reference} in the taxonomy IDs of your PhyloProfile file. Skipping ordering..')
        return df
    
    # retrieve order
    order = order_taxa(tree, str(reference))

    return df[order], order


def phyloprofile2matrix(path, ncbi, style, from_custom, fasF_filter, fasB_filter, reference):
    """
    Convert a phyloprofile file into a 2D matrix.
    Creates a copy of matrix containing the forward and backward FAS scores for writing phyloprofile output files.
    """


    def initialize_phyloprofile_df(path, gene_idx, taxa_idx):
        taxa = set()
        genes = set()
        with open(path) as fh:
            header = next(fh)
            for line in fh:
                dl = line.strip().split('\t')[:5]
                taxa.add(dl[taxa_idx])
                genes.add(dl[gene_idx])
        if not all(s.startswith('ncbi') for s in taxa):
            raise ValueError(f'Taxids in PhyloProfile file do not start with "ncbi". Alternatively, you might need to set "from_custom" to True.')
        return pd.DataFrame(index=list(genes), columns=list(taxa))

    def fill_phyloprofile_dataframe(df, path, style, gene_idx, taxa_idx, ortho_idx, fasf_idx, fasb_idx):
        outdf = df.copy()
        with open(path) as fh:
            header = next(fh)
            for line in fh:
                dl = line.strip().split('\t')
                gene = dl[gene_idx]
                orthoid = dl[ortho_idx]
                taxid = dl[taxa_idx]
                fasf = dl[fasf_idx]
                fasb = dl[fasb_idx]

                if float(fasf) <= fasF_filter or float(fasb) <= fasB_filter:
                    continue

                if style == 'orthoid':
                    df.loc[gene, taxid] = orthoid
                elif style == 'fasf':
                    df.loc[gene, taxid] = float(fasf)
                elif style == 'fasb':
                    df.loc[gene, taxid] = float(fasb)
                elif style == 'binary':
                    df.loc[gene, taxid] = 1
                else:
                    raise ValueErro(f'Cannot fill matrix in style "{style}". Choose "orthoid", "fasf", "fasb" or "binary"')
                outdf.loc[gene, taxid] = (orthoid, fasf, fasb)

        return df.fillna(0), outdf.fillna(0)

    ##################################################################################################
    logger = logging.getLogger('phyloprofile')
    
    if from_custom:
        gene_idx, taxa_idx, ortho_idx, fasf_idx, fasb_idx = 0, 3, 1, 5, 6
    else:
        gene_idx, taxa_idx, ortho_idx, fasf_idx, fasb_idx = 0, 1, 2, 3, 4

            
    logger.info(f'Initializing PhyloProfile matrix')
    df = initialize_phyloprofile_df(path, gene_idx, taxa_idx)
    if reference: 
        df, _ = sort_phyloprofile(df, ncbi, reference)
    logger.info(f'Loading PhyloProfile matrix')
    df, outdf = fill_phyloprofile_dataframe(df, path, style, gene_idx, taxa_idx, ortho_idx, fasf_idx, fasb_idx)
    return df, outdf
