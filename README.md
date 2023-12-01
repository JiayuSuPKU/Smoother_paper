# Scripts for "Smoother: a unified and modular framework for incorporating structural dependency in spatial omics data"
[![DOI](https://zenodo.org/badge/725407828.svg)](https://zenodo.org/doi/10.5281/zenodo.10242792)

## Description

This repository contains the scripts to reproduce results and figures for the Smoother paper:

Su, Jiayu, et al. "Smoother: A Unified and Modular Framework for Incorporating Structural Dependency in Spatial Omics Data." bioRxiv (2022): 2022-10.
https://www.biorxiv.org/content/10.1101/2022.10.25.513785v2.full

## Table of Contents

- [Dependencies](#dependencies)
- [Instructions](#instructions)
- [License](#license)

## Dependencies
To run the scripts in this repository, first make sure that all dependencies listed in the [environment.yml](environment.yml) file have been installed. 
We recommend creating a conda environment named `smoother_smoother` with:
```zsh
conda env create --file environment.yml
```
Now activate the environment and install the [Smoother package](https://github.com/JiayuSuPKU/Smoother), which is managed as a separate repository.
```zsh
conda activate smoother_smoother
pip install git+https://github.com/JiayuSuPKU/Smoother.git#egg=smoother
```

To run the full benchmark against existing methods, you may need to install additional packages, including:
1. (For cell-type deconvolution): [CARD (v1.1)](https://github.com/YMa-lab/CARD)
2. (For dimensionality reduction): [STAGATE (stagate-pyg v1.0.0)](https://github.com/QIFEIDKN/STAGATE_pyG), [SpaceFlow (v1.0.3)](https://github.com/hongleir/SpaceFlow). 

Please refer to the original repositories for installation instructions. Note `CARD` is an R package and its results were computed separately from the Python scripts/notebooks in this repository. 

### Troubleshooting
- `Smoother`'s dimensionality reduction module is built upon [scvi-tools](https://docs.scvi-tools.org/en/stable/index.html), which doesn't officially support Apple chips yet. To use `SCVI` and the corresponding `SpatialVAE` on Macs with Apple silicon, both [PyTorch](https://pytorch.org/) and [PyTorch Geometric](https://pytorch-geometric.readthedocs.io/en/latest/) (PyG) must be compiled with compatible wheel files. See [tips for installing PyG on M1 chips](https://github.com/rusty1s/pytorch_scatter/issues/241#issuecomment-1086887332). 
- `STAGATE` and `SpaceFlow` also require `PyG`. 
- Bug fix for `SpaceFlow`:
```zsh
git diff --unified=0 SpaceFlow/SpaceFlow.py

diff --git a/SpaceFlow/SpaceFlow.py b/SpaceFlow/SpaceFlow.py
index 170e4ff..de31b46 100644
--- a/SpaceFlow/SpaceFlow.py
+++ b/SpaceFlow/SpaceFlow.py
@@ -50 +50 @@ class SpaceFlow(object):
-            if gene_names:
+            if gene_names is not None:
@@ -52 +52 @@ class SpaceFlow(object):
-            if sample_names:
+            if sample_names is not None:
@@ -174 +174 @@ class SpaceFlow(object):
-        return nx.to_scipy_sparse_matrix(extended_graph, format='csr')
+        return nx.to_scipy_sparse_array(extended_graph, format='csr')
```

## Instructions

### Data and intermediate results
The `data/` folder contains the raw data used in this study. See [data/README.md](data/README.md) for details on the source of each dataset and links to download the data. 
The `results/` folder contains the results generated by the scripts. See [results/README.md](results/README.md) for details on the structure of the `results/` folder and where to download all intermediate results for reproducing the figures in the paper. Note that some intermediate results and figures in the notebooks may not be exactly the same as those in the paper due to the randomness in some analysis steps, but should be very close.

### Scripts
The `scripts/` folder contains all scripts used to generate the results and figures in the paper.
Here is a breakdown of all scripts in this repository by figure and analysis:

### Main figures
1. Fig 1: Overview of the Smoother framework. Cartoon illustration generated using [BioRender](https://app.biorender.com/user/signin).
2. Fig 2: Evaluation of spatial regularization effects on deconvolution accuracy using simulated data. See [scripts/synthetic_deconv/](scripts/synthetic_deconv/).
3. Fig 3: Smoother enhances cell-type deconvolution performance in various spatial omics data.
    1. Fig 3a-c: Breast cancer T cell infiltration (Visium). See [scripts/breast_cancer_infiltration/](scripts/breast_cancer_infiltration/).
    2. Fig 3d-g: Mouse brain neural subtypes (Visium). See [scripts/mouse_brain_visium/](scripts/mouse_brain_visium/).
    3. Fig 3h-k: Cross-modality deconvolution (RNA vs spatial CUT&Tag) of mouse embryonic data. See [scripts/embryo_cutntag/](scripts/embryo_cutntag/).
4. Fig 4: Smoother detects tumor-specific plasma cell subtypes in colorectal adenocarcinoma Stereo-seq slide. See [scripts/crc_stereo/](scripts/crc_stereo/).
5. Fig 5: Smoother enables spatially aware joint embeddings of single-cell and Slide-seqV2 data of human prostate and improves reference mapping accuracy. See [scripts/prostate_ref_mapping/](scripts/prostate_ref_mapping/).

### Supplementary figures
1. Fig S1: Spatial similarity patterns in different datasets. See [scripts/loss_design/figs1_real_data_similarity.ipynb](scripts/loss_design/figs1_real_data_similarity.ipynb).
2. Fig S2: Encoding boundary informations into the spatial prior. See [scripts/loss_design/figs2_swm_scaling.ipynb](scripts/loss_design/figs2_swm_scaling.ipynb).
3. Fig S3: Selecting covariance structures. See [scripts/loss_design/figs3_sp_loss_hyperparams.ipynb](scripts/loss_design/figs3_sp_loss_hyperparams.ipynb).
4. Fig S4-6: Data imputation and resolution enhancing on the DLPFC (Visium) and MBM (Slide-seqV2) datasets. See [scripts/imputation/](scripts/imputation/).
5. Fig S7-17: Benchmark deconvolution on synthetic data. See [scripts/synthetic_deconv/](scripts/synthetic_deconv/).
6. Fig S18: Breast cancer T cell infiltration (Visium). See [scripts/breast_cancer_infiltration/](scripts/breast_cancer_infiltration/).
7. Fig S19-22: Cross-modality deconvolution (RNA vs spatial CUT&Tag) of mouse embryonic data. See [scripts/embryo_cutntag/](scripts/embryo_cutntag/).
8. Fig S23-25: Smoother detects tumor-specific plasma cell subtypes in colorectal adenocarcinoma Stereo-seq slide. See [scripts/crc_stereo/](scripts/crc_stereo/).
9. Fig S26: Benchmark dimensionality reduction on the DLPFC (Visium) dataset. See [scripts/dlpfc_dr_benchmark/](scripts/dlpfc_dr_benchmark/).
10. Fig S27: SpatialVAE with single-cell and Slide-seqV2 data of human prostate. See [scripts/prostate_ref_mapping](scripts/prostate_ref_mapping).
11. Table S1: Computational cost analysis. See [scripts/runtime/tables1_check_runtime_memory_smoother.ipynb](scripts/runtime/tables1_check_runtime_memory_smoother.ipynb).


## License

This project is licensed under the [MIT License](LICENSE).
