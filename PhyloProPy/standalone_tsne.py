import argparse
import pandas as pd
from PhyloProPy.PhyloProfile import PhyloProfile
from PhyloProPy.plotting_tools import retrieve_taxa_mapping, perform_tsne
from ete3 import NCBITaxa
from sklearn.decomposition import PCA
import numpy as np
from sklearn.preprocessing import StandardScaler, RobustScaler, QuantileTransformer
from sklearn.manifold import TSNE
from collections import Counter
pd.set_option('mode.chained_assignment', None)
import logging


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


def main():
    # arguments
    args = parse_arguments()
    logger = logging.getLogger('phyloprofile')
    


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
        n_components=2, scaler=args.scaler, transpose=transpose, seed=args.seed,
        width=args.width, height=args.height
    )
    
    # save
    logger.info(f'Writing plot to {args.outpath}')
    fig.write_html(args.outpath)
    logger.info(f'Done')

if __name__ == "__main__":
    main()