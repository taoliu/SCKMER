# Time-stamp: <2020-01-02 17:35:32 taoliu>

"""Description: SCKMER decomp cmd

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
from time import time
from anndata import read_h5ad
import sklearn
from sklearn.decomposition import NMF, LatentDirichletAllocation, TruncatedSVD
from sklearn.preprocessing import Normalizer
from sklearn.pipeline import make_pipeline
from sklearn.feature_selection import VarianceThreshold

# ------------------------------------
# own python modules
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def run( args ):
    """Decompose
    
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
    if method == "lsi":
        method = "lsa"
    variance_threshold = options.variance_threshold
    n_components = options.n_components
    
    # read h5ad
    adata = read_h5ad( h5ad_fname )

    # decomp
    S = adata.X
    S = VarianceThreshold(variance_threshold).fit_transform(S)

    if method == "lsa":
        (m, T) = lsa_decomp ( S, n_components )
    elif method == "lda":
        (m, T) = lda_decomp ( S, n_components )        
    elif method == "nmf":
        (m, T) = nmf_decomp ( S, n_components )

    # replace
    adata.varm[method+"_components"] = m.components_.T
    adata.obsm["X_"+method] = T
    adata.uns[method] = {'params' : {'n_components' : n_components,
                                     'variance_threshold': variance_threshold}}
    #sc.pp.neighbors(adata, use_rep="X_"+method)
    #sc.tl.leiden(adata)
    #sc.tl.paga(adata)
    #sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
    #sc.tl.umap(adata, init_pos='paga')
    #sc.pl.umap(adata, color="leiden")
    #m = umap.UMAP(metric="euclidean", random_state=40, n_neighbors=50, min_dist=0.0, n_components=2).fit(adata.obsm['X_'+method])
    #adata.obsm['X_umap'] = m.embedding_
    # write h5a
    adata.write_h5ad( filename = out_fname )    


def nmf_decomp( t, n_components ):
    #print(f"Fit NMF with {n_components} components")
    nmf = NMF(n_components=n_components, init="nndsvd",random_state=1, alpha=.1, l1_ratio=.5, max_iter=100).fit( t )
    #print(f"Transform TD/IDF matrix with {n_components} components NMF")
    t_nmf = nmf.transform(t)
    #print("Reconstruction error (lower the better): %.3f" % nmf.reconstruction_err_ )
    #print("done in %0.3fs." % (time() - t0))
    return ( nmf, t_nmf )

def lda_decomp( t, n_components, learning_method="online", learning_offset=10.0, max_iter=20, random_state=1 ):
    #t0=time()    
    #print(f"Fit LDA with {n_components} components")
    lda = LatentDirichletAllocation(n_components=n_components, max_iter=max_iter, learning_method=learning_method, learning_offset=learning_offset, random_state=random_state).fit(t)
    #print(f"Transform TD/IDF matrix with {n_components} components LDA")    
    t_lda = lda.transform(t)
    score = lda.score(t)
    perplexity = lda.perplexity(t)
    #print("Approximate log likelihood score (higher the better): %.3f" % score)
    #print("Approximate perplexity (lower the better): %.3f" % perplexity)    
    #print("done in %0.3fs." % (time() - t0))
    return ( lda, t_lda )

def lsa_decomp( t, n_components, algorithm="arpack", random_state=1 ):
    #t0=time()
    #print(f"Fit LSA with {n_components} components")
    svd = TruncatedSVD(n_components, algorithm=algorithm, random_state=random_state)
    normalizer = Normalizer(copy=False)
    lsa = make_pipeline(svd, normalizer).fit(t)
    t_lsa = lsa.transform(t)
    explained_variance = svd.explained_variance_ratio_.sum()
    #print("Explained variance of the SVD step (higher the better): {}%".format(int(explained_variance * 100)))
    #print("done in %0.3fs." % (time() - t0))
    return ( lsa['truncatedsvd'], t_lsa )
