import argparse
import pandas as pd
import numpy as np
from collections import Counter
import logging


def phylo_heatmap(df, clustermethod, **kwargs):
    import plotly.express as px
    from plotly.colors import label_rgb
    from scipy.cluster.hierarchy import linkage, leaves_list
    
    color_map= [
        [0.0, 'white'],
        [0.00000001, label_rgb((242, 139, 13))],  # Blue
        [1.0, label_rgb((88, 117, 164))]  # orange
    ]

    # clustering
    if clustermethod:
        row_linkage = linkage(df, method=clustermethod)
        row_order = leaves_list(row_linkage)
        ordered_df = df.iloc[row_order, :]
    
    fig = px.imshow(df, color_continuous_scale=color_map, aspect="auto", **kwargs)

    return fig

def taxlevel_from_lineage(ncbi, lineage, taxlevel):
    node_id_dict = ncbi.get_name_translator(lineage)
    for node_name, node_id in node_id_dict.items():
        node_rank = ncbi.get_rank(node_id)
        if list(node_rank.values())[0] == taxlevel:
            return node_name


def retrieve_taxa_mapping(taxids4download, taxlevel, ncbi, update_taxonomy):
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
    df, taxlevel, ncbi,
    update_taxonomy, method, jitter, scaler, transpose, seed, n_components=2,
    **kwargs
):
    """
    Take a 2D representation of a phylogenetic profile and apply dimensionality reduction
    """
    
    
    from sklearn.preprocessing import StandardScaler, RobustScaler, QuantileTransformer
    
    logger = logging.getLogger('phyloprofile')
    # Convert string scaler to actual scaler object
    scaler_mapping = {
        'StandardScaler': StandardScaler(),
        'RobustScaler': RobustScaler(),
        'QuantileTransformer': QuantileTransformer(),
        'None': None
    }
    scaler = scaler_mapping[scaler]
    
    # Standardize the features
    if transpose:
        df = df.transpose()
    if scaler:
        scaled_data = scaler.fit_transform(df)
    else:
        scaled_data = df

    # reduce dimensions
    if method == 'PCA':
        from sklearn.decomposition import PCA
        pca = PCA(n_components=n_components)
        result = pca.fit_transform(scaled_data)
    elif method == 'tSNE':
        from sklearn.manifold import TSNE
        tsne = TSNE(n_components=2, random_state=seed)
        result = tsne.fit_transform(scaled_data)
    elif method == 'MDS':
        from sklearn.manifold import MDS
        mds = MDS(n_components=2, random_state=seed)
        result = mds.fit_transform(scaled_data)
    elif method == 'umap':
        import umap
        reducer = umap.UMAP(random_state=seed)
        result = reducer.fit_transform(scaled_data)
    else:
        raise ValueError(f'Unknown method "{method}". Choose "PCA" or "tSNE"')
        
    # store result in dataframe
    red_df = pd.DataFrame(data=result, columns=[f'PC{i}' for i in range(1, n_components+1)])
    if all(s.startswith('ncbi') for s in df.index):
        logger.info(f'Generating labels on "{taxlevel}" level')
        taxids4download = [taxid.replace('ncbi', '') for taxid in df.index]
        taxid2name, taxid2lineage, taxid2levelname = retrieve_taxa_mapping(taxids4download, taxlevel, ncbi, update_taxonomy)
        # assign labels
        red_df['taxid'] = [taxid.replace('ncbi', '') for taxid in df.index]
        red_df['species'] = red_df.taxid.apply(lambda x: taxid2name[int(x)])
        red_df['sum'] = df.sum(axis=1).values
        red_df['clade'] = red_df.taxid.apply(lambda x: taxid2levelname[x])
        red_df['clade'] = red_df.clade.apply(lambda x: 'NA' if x == None else x)
    elif all(s.startswith('ncbi') for s in df.columns):
        red_df['gene'] = df.index
        red_df['sum'] = df.sum(axis=1).values

    # add jitter
    if jitter:#
        np.random.seed(seed=seed)
        x_jitter = np.random.normal(loc=0, scale=jitter, size=red_df['PC1'].size)
        y_jitter = np.random.normal(loc=0, scale=jitter, size=red_df['PC2'].size)
        red_df['PC1'] = red_df['PC1'] + x_jitter
        red_df['PC2'] = red_df['PC2'] + y_jitter

    return red_df


def plot_tsne(
    red_df, method='tSNE', width=1000, height=1000, **kwargs
    
):
    import plotly.express as px
    # plot
    if 'taxid' in red_df.columns:
        fig = px.scatter(
            red_df, x='PC1', y='PC2', #title=f'{method} plot', 
            labels={'PC1': f'{method} 1', 'PC2': f'{method} 2'}, 
            hover_data={'species':True, 'clade': True, 'PC1': False, 'PC2': False}, width=width, height=height,
            size='sum',
            color='clade',
            **kwargs
            #color_discrete_sequence=px.colors.qualitative.Vivid
        )
    else:
        fig = px.scatter(
            red_df, x='PC1', y='PC2', #title=f'{method} plot', 
            labels={'PC1': f'{method} 1', 'PC2': f'{method} 2'}, 
            hover_data={'gene' :True, 'PC1': False, 'PC2': False}, width=width, height=height,
            size='sum',
            
            #color_discrete_sequence=px.colors.qualitative.Vivid
        )
    return fig

if __name__ == "__main__":
    main()