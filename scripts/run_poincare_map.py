import sys
sys.path.append('/Users/cchu/Desktop/phd_work/hyperChromatin/src/PoincareMaps')

import scanpy as sc
import pandas as pd
import numpy as np
from refs import celltype_colors
import matplotlib.pyplot as plt
from IPython.display import Image
import multiprocessing as mp
import argparse
import tqdm
from main import *
from poincare_maps import *


rna_pca_name = '../results/02/rna_pca.celltype_label'
rna_simba_name = '../results/02/rna_simba.celltype_label'
atac_pca_name = '../results/02/atac_pca.celltype_label'
atac_std_pca_name = '../results/02/atac_std_pca.celltype_label'
atac_simba_name = '../results/02/atac_simba.celltype_label'
simba_multi_name = '../results/02/multi_simba.celltype_label'

features_names = {
    "RNA PCA": rna_pca_name, 
    "RNA Simba": rna_simba_name, 
    "ATAC PCA": atac_pca_name, 
    "ATAC PCA (standard)": atac_std_pca_name, 
    "ATAC Simba": atac_simba_name, 
    "Simba Multi": simba_multi_name
}

def run_poincare_map(file_name, output_name, epochs=1000, lr=0.1):
    print(f'Running {file_name}...')
    features, labels = prepare_data(file_name, 
                                with_labels=True, 
                                normalize=True, 
                                n_pca=0)
    
    poincare_coord, _ = compute_poincare_maps(features, labels,
                        output_name,
                        mode='features', k_neighbours=15, 
                        distlocal='minkowski', sigma=1.0, gamma=2.0,
                        color_dict=celltype_colors, epochs=epochs,
                        batchsize=32, lr=lr, earlystop=0.0001, cuda=0)
    
    poincare_coord_df = pd.DataFrame(np.concatenate([poincare_coord, labels.reshape(-1, 1)], axis=1))
    poincare_coord_df.columns = ['x', 'y', 'labels']
    poincare_coord_fn = f'{output_name}.labeled.csv'
    poincare_coord_df.to_csv(poincare_coord_fn, index=False, sep=',')
    
    model = PoincareMaps(poincare_coord)
    model.plot('ori', labels=labels, file_name=output_name, 
           title_name='Poincaré map',
           coldict=celltype_colors, 
           labels_order=None, 
           zoom=4, bbox=(1.1, 0.8), leg=False, ft='png')
    
    return poincare_coord_fn
    

if __name__ == '__main__':
    # Run from scripts folder

    parser = argparse.ArgumentParser()
    parser.add_argument('--epochs', type=int, default=5)
    parser.add_argument('--lr', type=float, default=0.1)
    args = parser.parse_args()

    file_names = [
        # rna_pca_name, 
        rna_simba_name, 
        # atac_std_pca_name,atac_simba_name, 
        # simba_multi_name
    ]
    epochs = args.epochs
    lr = args.lr
    output_names = [f'../results/03/{file_name.split("/")[-1]}.ep{epochs}_lr{lr}.poincare_coord' for file_name in file_names]

    poincare_coord_fns = {}

    for i, (file_name, output_name) in tqdm.tqdm(enumerate(zip(file_names, output_names)), total=len(file_names)):
        result = run_poincare_map(file_name, output_name, epochs, lr)
        poincare_coord_fns[file_names[i]] = result


    pd.DataFrame(poincare_coord_fns, index=[0]).to_csv('../results/03/poincare_coord_fns.csv', index=False)
