# Time-stamp: <2020-01-03 15:29:10 taoliu>

"""Description: SCKMER cluster cmd

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------

import os
import sys
import logging

from anndata import read_h5ad
import scanpy as sc
from sklearn.cluster import SpectralClustering
import hdbscan

# ------------------------------------
# own python modules
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def run( args ):
    """Cluster
    
    """
    # Parse options...
    options = args
    # end of parsing commandline options
    info = options.info
    warn = options.warn
    debug = options.debug
    error = options.error

    # take options
    h5ad_fname = options.afile
    out_fname = options.ofile

    method = options.method
    data_source = options.data_source
    n_neighbors = options.n_neighbors

    if data_source != "X":
        data_source = "X_" + data_source

    random_state = options.random_state
    dbscan_min_samples = options.min_samples
    dbscan_min_cluster_size = options.min_cluster_size
        
    # read h5ad
    adata = read_h5ad( h5ad_fname )

    # clustering
    if method == "louvain":
        adata = louvain_clustering ( adata, data_source, n_neighbors = n_neighbors )
    elif method == "leiden":
        adata = leiden_clustering ( adata, data_source, n_neighbors = n_neighbors )
    elif method == "dbscan":
        adata = dbscan_clustering ( adata, data_source, dbscan_min_samples = dbscan_min_samples, dbscan_min_cluster_size = dbscan_min_cluster_size )
    elif method == "spectral":
        adata = spectral_clustering ( adata, data_source, n_neighbors = n_neighbors )

    # write h5a
    adata.write_h5ad( filename = out_fname )
    return

def louvain_clustering ( adata, source, n_neighbors = 15, random_state = 40 ):
    adata = add_graph ( adata, source, n_neighbors = n_neighbors )
    sc.tl.louvain(adata, random_state=random_state, use_weights=True)
    return adata

def leiden_clustering ( adata, source, n_neighbors = 15, random_state = 40 ):
    adata = add_graph ( adata, source, n_neighbors = n_neighbors )
    sc.tl.leiden(adata, random_state=random_state, use_weights=True)
    return adata

def dbscan_clustering ( adata, source, dbscan_min_samples = 5, dbscan_min_cluster_size = 200 ):
    if source == "X":
        S = adata.X
    else:
        S = adata.obsm[ source ]
    clusterer = hdbscan.HDBSCAN( min_samples=dbscan_min_samples, min_cluster_size=dbscan_min_cluster_size ).fit( S )
    adata.obs['dbscan'] = clusterer.labels_
    adata.uns['dbscan'] = { 'params' : {'min_samples':dbscan_min_samples,
                                        'min_cluster_size':dbscan_min_cluster_size} }
    return adata

def spectral_clustering ( adata, source, spectral_n_neighbors = 15, random_state = 40 ):
    if source == "X":
        S = adata.X
    else:
        S = adata.obsm[ source ]
    clusterer = SpectralClustering( affinity="nearest_neighbors", n_neighbors = spectral_n_neighbors, random_state = random_state ).fit( S )
    adata.obs['spectral'] = clusterer.labels_
    adata.uns['spectral'] = { 'params' : {'n_neighbors': spectral_n_neighbors,
                                          'affinity' : "nearest_neighbors",
                                          'random_state': random_state} }
    return adata

def add_graph ( adata, source, n_neighbors ):
    sc.pp.neighbors( adata, use_rep = source, n_neighbors = n_neighbors, metric = 'euclidean' )
    return adata

#sc.tl.paga(adata)
#sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
#sc.tl.umap(adata, init_pos='paga')
#sc.pl.umap(adata, color="leiden")
#m = umap.UMAP(metric="euclidean", random_state=40, n_neighbors=50, min_dist=0.0, n_components=2).fit(adata.obsm['X_'+method])
#adata.obsm['X_umap'] = m.embedding_
