# Time-stamp: <2020-01-03 15:58:48 taoliu>

"""Description: SCKMER umap cmd

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
import umap

# ------------------------------------
# own python modules
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def run( args ):
    """UMAP
    
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

    data_source = options.data_source
    n_neighbors = options.n_neighbors
    min_dist = options.min_dist

    if data_source != "X":
        if data_source == "lsi":
            data_source = "lsa"
        data_source = "X_" + data_source

    random_state = options.random_state

    # read h5ad
    adata = read_h5ad( h5ad_fname )
    
    # UMAP
    m = umap.UMAP(metric="euclidean", init="spectral", random_state=random_state, n_neighbors=n_neighbors, min_dist=min_dist, n_components=2).fit(adata.obsm[ data_source ])
    adata.obsm['X_umap'] = m.embedding_
    adata.uns['umap'] = { 'params':{'metric':'euclidean',
                                    'init':'spectral',
                                    'random_state':random_state,
                                    'n_neighbors':n_neighbors,
                                    'min_dist':min_dist } }
    
    # write h5a
    adata.write_h5ad( filename = out_fname )
    return
