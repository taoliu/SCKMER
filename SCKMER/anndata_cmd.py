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

from scipy import sparse as sp
from anndata import AnnData
from pandas import read_table as p_read_table
import numpy as np

# ------------------------------------
# own python modules
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def run( args ):
    """Compile an AnnData object from a sparseMatrix (npz) file for count table, the barcode and the k-mer names files.
    
    """
    # Parse options...
    options = args
    # end of parsing commandline options
    info = options.info
    warn = options.warn
    debug = options.debug
    error = options.error

    # take options
    ct_fname = options.cfile
    barcode_fname = options.bfile
    kmer_fname = options.kfile
    out_fname = options.ofile
    
    # read the sparse matrix
    t = sp.load_npz( ct_fname )

    # read .var from kmer file
    var = load_var( kmer_fname )

    # read .obs from barcode_fname
    obs = load_obs( barcode_fname )

    # create anndata
    adata = AnnData( X = t, obs = obs, var = var )

    # add n_counts and n_kmers to adata.obs
    adata.obs['n_counts'] = np.sum(adata.X, axis=1).A1    
    adata.obs['n_kmers'] = np.sum(adata.X > 0, axis=1).A1

    # add n_cells to adata.var
    adata.var['n_cells'] = np.sum(adata.X > 0, axis=0).A1
    
    # save adata to h5ad
    adata.write_h5ad( filename = out_fname )

    return
    
def load_obs ( fname ):
    obs = p_read_table( fname )
    return obs


def load_var ( fname ):
    var = p_read_table( fname )
    return var
