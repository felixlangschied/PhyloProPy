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
    order = [taxid for taxid in order if taxid in df.columns]

    return df[order], order


def phyloprofile2matrix(path, ncbi, style, from_custom, fasF_filter, fasB_filter, fillna, resolve_coorthologs, reference):
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
        """Fill the empty dataframe with a list of values to accomodate co-orthologs. Then resolve the lists."""
        def update_cell(df, gene, taxid, value):
            cell_value = df.at[gene, taxid]
            if isinstance(cell_value, float):
                df.loc[gene, taxid] = [value]
            elif isinstance(cell_value, list):
                df.loc[gene, taxid].append(value)
            else:
                raise ValueError('Wrong initialzation of DataFrame for loading PhyloProfile.')

        ##################################################################
        outdf = df.copy()
        with open(path) as fh:
            header = next(fh)
            for line in fh:
                dl = line.strip().split('\t')
                if len(dl) <= fasf_idx:
                    fasf, fasb = 1, 1
                    gene, orthoid, taxid = dl[gene_idx], dl[ortho_idx], dl[taxa_idx]
                elif len(dl) <= fasb_idx:
                    fasb = 1
                    gene, orthoid, taxid, fasf = dl[gene_idx], dl[ortho_idx], dl[taxa_idx], float(dl[fasf_idx])
                else:
                    fasf, fasb = dl[fasf_idx], dl[fasb_idx]
                    if fasf == 'NA':
                        continue
                    gene, orthoid, taxid, fasf, fasb = dl[gene_idx], dl[ortho_idx], dl[taxa_idx], float(dl[fasf_idx]), float(dl[fasb_idx])
                
                # apply filter
                if fasf < fasF_filter or fasb < fasB_filter:
                    continue

                # Fill with list
                if style == 'orthoid':
                    update_cell(df, gene, taxid, orthoid)
                elif style == 'ncRNA':
                    seedscore = fasf
                    if seedscore == 0.0:
                        seedscore = 0.5
                    update_cell(df, gene, taxid, seedscore)
                elif style in ['fasf', 'fasb']:
                    value = fasf if style == 'fasf' else fasb
                    update_cell(df, gene, taxid, value)
                elif style == 'binary':
                    df.loc[gene, taxid] = 1
                else:
                    raise ValueError(f'Cannot fill matrix in style "{style}". Choose "orthoid", "fasf", "fasb" or "binary"')
                
                # Update outdf 
                update_cell(outdf, gene, taxid, (orthoid, fasf, fasb))

        # resolve lists for scores
        df, outdf = df.fillna(fillna), outdf.fillna(fillna)
        if resolve_coorthologs and style in ['fasf', 'fasb', 'ncRNA']:
            df = df.applymap(lambda x: max(x) if isinstance(x, list) else x)
        
        return df, outdf

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
