# Time-stamp: <2020-01-02 15:50:16 taoliu>

"""Description: SCKMER tfidf cmd

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
from sklearn.feature_selection import SelectKBest

# ------------------------------------
# own python modules
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def run( args ):
    """Filter cells and kmers
    
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

    max_kmers = options.max_kmers
    max_counts = options.max_counts
    max_cells = options.max_cells
    
    min_kmers = options.min_kmers
    min_counts = options.min_counts
    min_cells = options.min_cells

    assert min_kmers>=0, "min_kmers must be >= 0"
    assert min_counts>=0, "min_counts must be >= 0"
    assert min_cells>=0, "min_cells must be >= 0"
    if max_kmers:
        assert max_kmers>=min_kmers, "If set, max_kmers must be >= min_kmers"
    if max_counts:
        assert max_counts>=min_counts, "If set, max_counts must be >= min_counts"
    if max_cells:
        assert max_cells>=min_cells, "If set, max_cells must be >= min_cells"

    # read h5ad
    adata = read_h5ad( h5ad_fname )

    # filtering upper bound
    if max_kmers:
        adata = adata[ adata.obs.n_kmers <= max_kmers, : ]
    if max_counts:
        adata = adata[ adata.obs.n_counts <= max_counts, : ]
    if max_cells:
        adata = adata[ :, adata.var.n_cells <= max_cells ]        

    # filtering lower bound
    if min_kmers != 0:
        adata = adata[ adata.obs.n_kmers >= min_kmers, : ]
    if min_counts != 0:
        adata = adata[ adata.obs.n_counts >= min_counts, : ]
    if min_cells != 0:
        adata = adata[ :, adata.var.n_cells >= min_cells ]

    # top top_variable_kmers kmers
    #S = adata.X
    #S = SelectKBest(chi2, k=top_variable_kmers).fit_transform(S)
    
    # write h5a
    adata.write_h5ad( filename = out_fname )

    return
