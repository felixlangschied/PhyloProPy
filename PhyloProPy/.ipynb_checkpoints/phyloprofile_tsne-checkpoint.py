import argparse
import pandas as pd
from PhyloProPy.PhyloProfile import PhyloProfile
from ete3 import NCBITaxa
from sklearn.decomposition import PCA
import numpy as np
from sklearn.preprocessing import StandardScaler, RobustScaler, QuantileTransformer
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
from sklearn.manifold import TSNE
from collections import Counter
pd.set_option('mode.chained_assignment', None)
import logging


def taxlevel_from_lineage(ncbi, lineage, taxlevel):
    node_id_dict = ncbi.get_name_translator(lineage)
    for node_name, node_id in node_id_dict.items():
        node_rank = ncbi.get_rank(node_id)
        if list(node_rank.values())[0] == taxlevel:
            return node_name


def retrieve_taxa_mapping(taxids4download, taxlevel, update_taxonomy):
    ncbi = NCBITaxa()
    if update_taxonomy:
        ncbi.update_taxonomy_database()

    # initialize output
    taxid2name = ncbi.get_taxid_translator(taxids4download)
    taxid2lineage = {}
    taxid2levelname = {}

    # get lineage and levelname
    for taxid in taxids4download:
        nodelineage = ncbi.get_lineage(taxid)
        nodid2lineage = ncbi.get_taxid_translator(nodelineage)
        lineage = list(nodid2lineage.values())
        levelname = taxlevel_from_lineage(ncbi, lineage, taxlevel)

        taxid2lineage[taxid] = lineage
        taxid2levelname[taxid] = levelname

    return taxid2name, taxid2lineage, taxid2levelname


def dimension_reduced_phyloprofile(
    df, taxid2name, taxid2levelname,
    method='PCA', jitter='',
    n_components=2, scaler=StandardScaler(), transpose=True, seed=42
):
    """
    Take a 2D representation of a phylogenetic profile and apply dimensionality reduction
    """
    # Standardize the features
    if transpose:
        df = df.transpose()
    if scaler:
        scaled_data = scaler.fit_transform(df)
    else:
        scaled_data = df

    # reduce dimensions
    if method == 'PCA':
        pca = PCA(n_components=n_components)
        result = pca.fit_transform(scaled_data)
    elif method == 'tSNE':
        tsne = TSNE(n_components=2, random_state=seed)
        result = tsne.fit_transform(scaled_data)
    else:
        raise ValueError(f'Unknown method "{method}". Choose "PCA" or "tSNE"')
        
    # store result in dataframe
    red_df = pd.DataFrame(data=result, columns=[f'PC{i}' for i in range(1, n_components+1)])
    if all(s.startswith('ncbi') for s in df.index):
        red_df['taxid'] = [taxid.replace('ncbi', '') for taxid in df.index]
        red_df['species'] = red_df.taxid.apply(lambda x: taxid2name[int(x)])
        red_df['taxlevel'] = red_df.taxid.apply(lambda x: taxid2levelname[x])
        red_df['sum'] = df.sum(axis=1).values
    elif all(s.startswith('ncbi') for s in df.columns):
        red_df['gene'] = df.index
        red_df['sum'] = df.sum(axis=1).values

    # add jitter
    if jitter:
        x_jitter = np.random.normal(loc=0, scale=jitter, size=red_df['PC1'].size)
        y_jitter = np.random.normal(loc=0, scale=jitter, size=red_df['PC2'].size)
        red_df['PC1'] = red_df['PC1'] + x_jitter
        red_df['PC2'] = red_df['PC2'] + y_jitter

    return red_df


def parse_arguments():
    parser = argparse.ArgumentParser(description='Reduce dimensions of a phylogenetic profile and plot in 2D space')
    parser.add_argument('--path', type=str, required=True, help='Path to the phyloprofile file')
    parser.add_argument('--outpath', type=str, default='./phyloprofile_tsne.html', help='Path to the output HTML file')
    parser.add_argument('--taxlevel', type=str, default='species', choices=['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'], help='Taxonomic level for analysis')
    parser.add_argument('--orient', type=str, default='species', choices=['species', 'genes'], help='Orientation for analysis')
    parser.add_argument('--style', type=str, default='binary', choices=['binary', 'fasf', 'fasb'], help='Style of analysis')
    parser.add_argument('--scaler', type=str, default='StandardScaler', choices=['StandardScaler', 'RobustScaler', 'QuantileTransformer', 'None'], help='Scaler function to use')
    parser.add_argument('--method', type=str, default='tSNE', choices=['tSNE', 'PCA'], help='Dimensionality reduction method to use')
    parser.add_argument('--from_custom', action='store_true', help='Set if input file was exported from a higher-level rank PhyloProfile')
    parser.add_argument('--update_taxonomy', action='store_true', help='Set to update the NCBI taxonomy database')
    parser.add_argument('--seed', type=int, default=42, help='Seed for randomness')
    parser.add_argument('--jitter', type=float, default=0.3, help='Jitter for plot scatter points')
    parser.add_argument('--width', type=int, default=1000, help='Width of the scatter plot')
    parser.add_argument('--height', type=int, default=1000, help='Heigth of the scatter plot')
    return parser.parse_args()


def perform_tsne(
    df, taxid2name, taxid2levelname, method, jitter, n_components, scaler, transpose, seed, width, height
    
):
    red_df = dimension_reduced_phyloprofile(
        df, taxid2name, taxid2levelname, 
        method=method, jitter=jitter,
        n_components=2, scaler=scaler, transpose=transpose, seed=seed
    )
    # plot
    if 'taxid' in red_df.columns:
        fig = px.scatter(
            red_df, x='PC1', y='PC2', #title=f'{method} plot', 
            labels={'PC1': f'{method} 1', 'PC2': f'{method} 2'}, 
            hover_data={'species':True, 'taxlevel': True, 'PC1': False, 'PC2': False}, width=width, height=height,
            size='sum',
            color='taxlevel',
            #color_discrete_sequence=px.colors.qualitative.Vivid
        )
    else:
        fig = px.scatter(
            red_df, x='PC1', y='PC2', #title=f'{method} plot', 
            labels={'PC1': f'{method} 1', 'PC2': f'{method} 2'}, 
            hover_data={'gene' :True, 'PC1': False, 'PC2': False}, width=width, height=height,
            size='sum'
            #color_discrete_sequence=px.colors.qualitative.Vivid
        )
    return fig


def main():
    # arguments
    args = parse_arguments()
    logger = logging.getLogger('phyloprofile')
    
    # Convert string scaler to actual scaler object
    scaler_mapping = {
        'StandardScaler': StandardScaler(),
        'RobustScaler': RobustScaler(),
        'QuantileTransformer': QuantileTransformer(),
        'None': None
    }
    scaler = scaler_mapping[args.scaler]

    # load data    
    pp = PhyloProfile(args.path, style=args.style, from_custom=args.from_custom)

    # orientation
    if args.orient == 'species':
        transpose = True
        logger.info(f'Generating labels on "{args.taxlevel}" level')
        taxids4download = [taxid.replace('ncbi', '') for taxid in pp.matrix.columns]
        taxid2name, taxid2lineage, taxid2levelname = retrieve_taxa_mapping(taxids4download, args.taxlevel, args.update_taxonomy)
    elif args.orient == 'genes':
        transpose = False
        taxid2name, taxid2lineage, taxid2levelname = {}, {}, {}
    else:
        raise ValueError(f'Unknown orientation "{orient}". Choose "species" or "genes".')



    # tsne
    logger.info(f'Generating plot')
    fig = perform_tsne(
        pp.matrix, taxid2name, taxid2levelname, 
        method=args.method, jitter=args.jitter,
        n_components=2, scaler=scaler, transpose=transpose, seed=args.seed,
        width=args.width, height=args.height
    )
    
    # save
    logger.info(f'Writing plot to {args.outpath}')
    fig.write_html(args.outpath)
    logger.info(f'Done')

if __name__ == "__main__":
    main()