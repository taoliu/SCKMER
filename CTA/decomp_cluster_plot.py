#!/usr/bin/env python3

import scipy
import sys
import os
import umap
import numpy as np
import matplotlib
matplotlib.use('agg')
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import sklearn
from sklearn.decomposition import NMF
from sklearn.decomposition import LatentDirichletAllocation
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import Normalizer
from sklearn.pipeline import make_pipeline
import hdbscan
import seaborn as sns
from time import time
import pickle
import warnings
warnings.filterwarnings('ignore')

def main ():
    if len( sys.argv ) < 5:
        sys.stderr.write( "Data Decomp + UMAP on a given sparse matrix (.npz)\n" )
        sys.stderr.write( "need 4 parameter: <method> <n_components> <.npz> <known_labels_file>\n" )
        sys.stderr.write( "   method: Decomposition method can only be: NMF, LDA, or LSA\n" )
        sys.stderr.write( "   n_components: Number of components for decomposition method\n" )
        sys.stderr.write( "   .npz: Sparse Matrix (uniq_id x features) of counts saved by Scipy in NPZ format.\n" )
        sys.stderr.write( "   known_labels_file: tab-delimited file with 1st column as uniq_id and 2nd column as known label\n\n" )
        sys.stderr.write( "Note: if some backup files exist, decomposition or umap can be skipped:\n" )
        sys.stderr.write( "   transformed_{decomp_method}_{n_components}.sav: Backup file of decomposition transformed data\n" )
        sys.stderr.write( "   umap_embedding_{decomp_method}_{n_components}_cosine/euclidean_2.sav: Backup file of UMAP embedding data\n" )                
        sys.stderr.write( "\n" )
        exit(1)
    decomp_method = sys.argv[1]
    if not decomp_method in [ "NMF", "LDA", "LSA", "LSI" ]:
        sys.stderr.write(f"ERROR: Please choose method from: NMF, LDA, LSA, or LSI. Your input: {decomp_method}\n")
        exit(1)
    n_components = int( sys.argv[2] )
    if n_components < 1:
        sys.stderr.write(f"ERROR: n_components must be larger than 0. Your input: {n_components}\n")
        exit(1)        
    tfidf_filename = sys.argv[3]
    if not os.path.isfile( tfidf_filename ):
        sys.stderr.write(f"ERROR: .npz file {tfidf_filename} can't be found!\n")
        exit(1)        
    known_label_filename = sys.argv[4]
    if not os.path.isfile( known_label_filename ):
        sys.stderr.write(f"ERROR: known_label_filename file {known_label_filename} can't be found!\n")
        exit(1)        

    # load known labels dict
    uniq_ids, labels_dict = build_known_labels( known_label_filename )
    
    # decomposition
    t = decomp( tfidf_filename, decomp_method, n_components )

    # umap
    e_cosine = run_umap( t, decomp_method, n_components, "cosine", 2 )
    e_euclidean = run_umap( t, decomp_method, n_components, "euclidean", 2 )    

    # plot
    png_prefix = f"umap_fig_{decomp_method}_{n_components}_cosine"
    plot_umap_2d ( e_cosine, uniq_ids, labels_dict, png_prefix )

    png_prefix = f"umap_fig_{decomp_method}_{n_components}_euclidean"
    plot_umap_2d ( e_euclidean, uniq_ids, labels_dict, png_prefix )    

def build_known_labels ( label_file ):
    uniq_id_list = []
    labels_dict = {}
    with open( label_file,"r" ) as label_fhd:
        label_fhd.readline()      # skip first line
        for l in label_fhd:
            (uniq_id, label) = l.rstrip().split("\t")
            uniq_id_list.append( uniq_id )
            labels_dict[ uniq_id ] = label
    return ( uniq_id_list, labels_dict )
    
def decomp( tfidf_filename, decomp_method, n_components ):
    # check existing bk file first:
    tfidf_transformed_fn = f"transformed_{decomp_method}_{n_components}.sav"
    tfidf_model_fn = f"model_{decomp_method}_{n_components}.sav"    
    if os.path.isfile( tfidf_transformed_fn ):
        # load
        print( f"Load TF/IDF {decomp_method} transformed data" )
        t_trans = load_pickle_file( tfidf_transformed_fn )
        return t_trans
    else:
        # load TF/IDF reweighted count table
        # TF/IDF reweight
        t = scipy.sparse.load_npz( tfidf_filename )        
        # decompose
        if decomp_method == "NMF":
            ( model, t_trans ) = nmf_decomp( t, n_components )
        elif decomp_method == "LDA":
            ( model, t_trans ) = lda_decomp( t, n_components )
        elif decomp_method == "LSI" or decomp_method == "LSA":
            ( model, t_trans ) = lsa_decomp( t, n_components )
        else:
            exit(1)
        save_fit_trans_file( model, t_trans, tfidfmodel_fn, tfidf_transformed_fn )
        return t_trans

def load_pickle_file( fn ):
    with open( fn, "rb" ) as fhd:
        d = pickle.load( fhd )
    return d

def save_fit_trans_file( model, t, model_fn, trans_fn ):
    with open( model_fn, "wb" ) as fhd:
        pickle.dump( model, fhd )
    with open( trans_fn, "wb" ) as fhd:
        pickle.dump( t, fhd )        
    return

def nmf_decomp( t, n_components ):
    t0=time()
    print(f"Fit NMF with {n_components} components")
    nmf = NMF(n_components=n_components, init="nndsvd",random_state=1, alpha=.1, l1_ratio=.5, max_iter=100).fit( t )
    print(f"Transform TD/IDF matrix with {n_components} components NMF")
    t_nmf = nmf.transform(t)
    #print("Reconstruction error: %.3f" % nmf.reconstruction_err_ )
    print("done in %0.3fs." % (time() - t0))
    return ( nmf, t_nmf )

def lda_decomp( t, n_components ):
    t0=time()    
    print(f"Fit LDA with {n_components} components")
    lda = LatentDirichletAllocation(n_components=n_components, max_iter=20, learning_method='online', learning_offset=10., random_state=1).fit(t)
    print(f"Transform TD/IDF matrix with {n_components} components LDA")    
    t_lda = lda.transform(t)
    #score = lda.score(t)
    #print("Approximate likelihood score: %.3f" % score)
    print("done in %0.3fs." % (time() - t0))
    return ( lda, t_lda )

def lsa_decomp( t, n_components ):
    t0=time()
    print(f"Fit LSA with {n_components} components")
    svd = TruncatedSVD(n_components, algorithm="ARPACK", random_state=1)
    normalizer = Normalizer(copy=False)
    lsa = make_pipeline(svd, normalizer).fit(t)
    t_lsa = lsa.transform(t)
    explained_variance = svd.explained_variance_ratio_.sum()
    print("Explained variance of the SVD step: {}%".format(int(explained_variance * 100)))
    print("done in %0.3fs." % (time() - t0))
    return ( lsa, t_lsa )

def run_umap ( data, decomp_method, n_components, dist_name, n_dimensions ):
    # check existing bk file first:
    umap_model_fn = f"umap_model_{decomp_method}_{n_components}_{dist_name}_{n_dimensions}.sav"
    umap_embedding_fn = f"umap_embedding_{decomp_method}_{n_components}_{dist_name}_{n_dimensions}.sav"
    if os.path.isfile( umap_embedding_fn ):
        # load
        print( f"Load UMAP {dist_name} embedding" )
        e = load_pickle_file( umap_embedding_fn )
        return e
    else:
        print( f"UMAP ({dist_name})" )
        m = umap.UMAP(metric=dist_name, random_state=40, n_neighbors=50, min_dist=0.0, n_components=n_dimensions).fit(data)
        e = m.transform(data)
        save_fit_trans_file( m, e, umap_model_fn, umap_embedding_fn )
        return e

def plot_umap_2d ( e, uniq_ids, labels_dict, png_prefix ):
    # plot
    plt.rcParams.update({'font.size': 22})

    print("UMAP w/o labels")
    # umap
    plt.figure(figsize=(16, 16))
    plt.style.use('fast')
    plt.scatter(e[:, 0], e[:, 1], s=80, alpha=0.4)
    plt.savefig(f"{png_prefix}_wo_labels.png", format="png")

    # umap + known labels
    print("UMAP w/ known labels")
    classes = sorted(list(set(labels_dict.values())))
    map_class = { x:i for ( i, x ) in enumerate( classes )}
    known_labels = [ map_class[labels_dict[x]] for x in uniq_ids ]
    cm = ListedColormap(sns.color_palette("Paired",24))
    fig, ax = plt.subplots(1, figsize=(20, 16))
    sc = ax.scatter(e[:, 0], e[:, 1], s=80, c=known_labels, cmap=cm, alpha=0.4)
    cbar = fig.colorbar(sc, boundaries=np.arange( len(classes) + 1 )-0.5)
    cbar.set_ticks(np.arange( len(classes) ))
    cbar.set_ticklabels(classes)
    plt.title('UMAP with labels');
    plt.savefig(f"{png_prefix}_w_known_labels.png", format="png")

    # HDBSCAN
    print("UMAP w/ HDBSCAN clustering")
    hdbscan_clusterer = hdbscan.HDBSCAN( min_samples=50, min_cluster_size=100 ).fit(e)
    #cm = ListedColormap(sns.color_palette("Paired",12))    
    color_palette = sns.color_palette('Paired', 12)
    cluster_colors = [color_palette[x] if x >= 0 else (0.5, 0.5, 0.5) for x in hdbscan_clusterer.labels_]
    plt.figure(figsize=(16, 16))
    plt.style.use('fast')
    plt.scatter(e[:, 0], e[:, 1], s=80, c=cluster_colors, alpha=0.4)
    #plt.colorbar()
    plt.savefig(f"{png_prefix}_w_hdbscan_labels.png", format="png")    
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupted me! ;-) Bye!\n")
        sys.exit(0)
