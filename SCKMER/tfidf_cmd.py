# Time-stamp: <2020-01-02 12:41:39 taoliu>

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

from scipy import sparse as sp
from anndata import read_h5ad
import sklearn
from sklearn.feature_extraction.text import TfidfTransformer

# ------------------------------------
# own python modules
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def run( args ):
    """Make TFIDF
    
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
    
    norm = options.norm
    smooth_idf = True #options.smooth_idf
    use_idf = True #options.use_idf
    sublinear_tf = options.sublinear_tf
    
    # read h5ad
    adata = read_h5ad( h5ad_fname )

    # TF/IDF re-weighting
    S = adata.X
    T = TfidfTransformer(norm = norm, smooth_idf = smooth_idf,
                         use_idf = use_idf, sublinear_tf = sublinear_tf).fit_transform(S)

    # replace
    adata.X = S
    
    # write h5a
    adata.write_h5ad( filename = out_fname )

    return
