import os
import numpy as np
import pytorch_lightning as pl
import torch
from torch.utils.data import DataLoader
import anndata as ad
from typing import Optional, Dict, Any
from tqdm import tqdm

from nicheformer.models import Nicheformer
from nicheformer.data import NicheformerDataset

print(torch.__version__)
print(torch.version.cuda)
print(torch.cuda.is_available())
print(torch.cuda.device_count())
if torch.cuda.is_available():
    print(torch.cuda.get_device_name(0))
    

config = {
    'data_path': '/oak/stanford/groups/ccurtis2/users/bgsimon/SpatialWorkshop_Integration/objects/PATH-IC_CosMx_Analysis_coregistered/PATH-IC_CosMx_Analysis_Madrid1_with_HECoreID_checkedtxt_SingleR_niche_with_fibannot.h5ad', #'path/to/your/data.h5ad',  # Path to your AnnData file
    'technology_mean_path': '/oak/stanford/groups/ccurtis2/users/bgsimon/MultimodalSpatialIntegration/Nicheformer/nicheformer/data/model_means/cosmx_6k_panTMA_mean.npy', #'path/to/technology_mean.npy',  # Path to technology mean file
    'checkpoint_path': '/oak/stanford/groups/ccurtis2/users/bgsimon/MultimodalSpatialIntegration/Nicheformer/nicheformer.ckpt',  # Path to model checkpoint
    'output_path': '/oak/stanford/groups/ccurtis2/users/bgsimon/MultimodalSpatialIntegration/Nicheformer/nicheformer/my_outputs/5.12.26_3/MATMA1_with_NFembeddings_panMA.h5ad',  # Where to save the result, it is a new h5ad
    'output_dir': '/oak/stanford/groups/ccurtis2/users/bgsimon/MultimodalSpatialIntegration/Nicheformer/nicheformer/my_outputs/5.12.26_3',  # Directory for any intermediate outputs
    'batch_size': 32,
    'max_seq_len': 1500, 
    'aux_tokens': 30, 
    'chunk_size': 1000, # to prevent OOM
    'num_workers': 4,
    'precision': 32,
    'embedding_layer': -1,  # Which layer to extract embeddings from (-1 for last layer)
    'embedding_name': 'NFembeddings'  # Name suffix for the embedding key in adata.obsm
}
model = ad.read_h5ad('/oak/stanford/groups/ccurtis2/users/bgsimon/MultimodalSpatialIntegration/Nicheformer/nicheformer/data/model_means/model.h5ad')

import pandas as pd
import scipy.sparse as sp
import mygene
from pathlib import Path

# Set random seed
pl.seed_everything(42)

# Load data
adata_raw = ad.read_h5ad(config["data_path"])
technology_mean = np.load(config["technology_mean_path"])

print("Raw CosMx:", adata_raw.shape)
print("Model:", model.shape)
print("Technology mean:", technology_mean.shape)

assert technology_mean.shape[0] == model.n_vars

symbols = adata_raw.var_names.astype(str).tolist()

mg = mygene.MyGeneInfo()
query = mg.querymany(
    symbols,
    scopes="symbol",
    fields="ensembl.gene,symbol",
    species="human",
    as_dataframe=True,
    df_index=True,
    verbose=False,
)

symbol_to_ensembl = {}

for symbol in symbols:
    if symbol not in query.index:
        continue

    row = query.loc[symbol]

    if isinstance(row, pd.DataFrame):
        row = row.iloc[0]

    ens = row.get("ensembl.gene", None)

    if isinstance(ens, list):
        ens = ens[0]

    if pd.notna(ens):
        symbol_to_ensembl[symbol] = str(ens).split(".")[0]
        
model_gene_to_idx = {str(g): i for i, g in enumerate(model.var_names)}

common = []

for adata_idx, symbol in enumerate(symbols):
    ens = symbol_to_ensembl.get(symbol)

    if ens in model_gene_to_idx:
        common.append((adata_idx, symbol, ens, model_gene_to_idx[ens]))

print("Mapped genes:", len(common))
print("First mapped:", common[:10])

X_full = sp.lil_matrix((adata_raw.n_obs, model.n_vars), dtype=np.float32)

for adata_idx, symbol, ens, model_idx in tqdm(common, desc="Projecting CosMx to Nicheformer space"):
    X_full[:, model_idx] = adata_raw.X[:, adata_idx]

adata = ad.AnnData(
    X=X_full.tocsr(),
    obs=adata_raw.obs.copy(),
    var=model.var.copy()
)

adata.var_names = model.var_names.copy()

print(adata)

adata.obs["nicheformer_split"] = "train"
adata.obs["split"] = "train"

adata.obs["modality"] = 4  # spatial
adata.obs["specie"] = 5    # human
adata.obs["assay"] = 8     # CosMx

print(adata.obs["nicheformer_split"].value_counts())

assert adata.n_vars == model.n_vars
assert adata.n_vars == technology_mean.shape[0]
assert list(adata.var_names) == list(model.var_names)

print("adata shape:", adata.shape)
print("technology_mean shape:", technology_mean.shape)
print("mean non-NaN:", np.sum(~np.isnan(technology_mean)))

dataset = NicheformerDataset(
    adata=adata,
    technology_mean=technology_mean,
    split="train",
    max_seq_len=config["max_seq_len"],
    aux_tokens=config["aux_tokens"],
    chunk_size=config["chunk_size"],
    metadata_fields={"obs": ["modality", "specie", "assay"]}
)

dataloader = DataLoader(
    dataset,
    batch_size=config["batch_size"],
    shuffle=False,
    num_workers=config["num_workers"],
    pin_memory=True
)

# Load pre-trained model
model = Nicheformer.load_from_checkpoint(checkpoint_path=config['checkpoint_path'], strict=False)
model.eval()  # Set to evaluation mode

# Configure trainer
trainer = pl.Trainer(
    accelerator='gpu' if torch.cuda.is_available() else 'cpu',
    devices=1,
    default_root_dir=config['output_dir'],
    precision=config.get('precision', 32),
)

print("Extracting embeddings...")
embeddings = []
device = model.embeddings.weight.device

with torch.no_grad():
    for batch in tqdm(dataloader):
        # Move batch to device
        batch = {k: v.to(device) if isinstance(v, torch.Tensor) else v
                for k, v in batch.items()}

        # Get embeddings from the model
        emb = model.get_embeddings(
            batch=batch,
            layer=config.get('embedding_layer', -1)  # Default to last layer
        )
        embeddings.append(emb.cpu().numpy())


# Concatenate all embeddings
embeddings = np.concatenate(embeddings, axis=0)

# Store embeddings in AnnData object
embedding_key = f"X_niche_{config.get('embedding_name', 'embeddings')}"
adata.obsm[embedding_key] = embeddings

# Save updated AnnData
adata.write_h5ad(config['output_path'])

print(f"Embeddings saved to {config['output_path']} in obsm['{embedding_key}']")