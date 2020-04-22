# Time-stamp: <2020-01-02 15:18:16 taoliu>

"""Description: SCKMER anndata cmd

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
import anndata as ad
import numpy as np

# ------------------------------------
# own python modules
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def run( args ):
    """Compile an AnnData object from a loom file
    
    """
    # Parse options...
    options = args
    # end of parsing commandline options
    info = options.info
    warn = options.warn
    debug = options.debug
    error = options.error

    # take options
    h5_fname = options.ifile
    out_fname = options.ofile
    
    # read
    adata = ad.read_loom( h5_fname, sparse=True )

    # add n_counts and n_kmers to adata.obs
    adata.obs['n_counts'] = np.sum(adata.X, axis=1).A1    
    adata.obs['n_kmers'] = np.sum(adata.X > 0, axis=1).A1

    # add n_cells to adata.var
    adata.var['n_cells'] = np.sum(adata.X > 0, axis=0).A1
    
    # save adata to h5ad
    adata.write_h5ad( filename = out_fname )
    return
    
