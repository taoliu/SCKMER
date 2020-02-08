# Time-stamp: <2020-02-06 11:28:09 taoliu>

"""Description: SCKMER view cmd

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

import matplotlib; matplotlib.use('Agg')

from anndata import read_h5ad
import scanpy as sc

# ------------------------------------
# own python modules
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def run( args ):
    """View
    
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
    out_dname = options.odir

    extra_keys = options.keys
    
    # read h5ad
    adata = read_h5ad( h5ad_fname )
    if extra_keys:
        sc.pl.umap(adata, color=extra_keys, save=out_dname, show=False)
    else:
        sc.pl.umap(adata, save=out_dname, show=False)

    return
