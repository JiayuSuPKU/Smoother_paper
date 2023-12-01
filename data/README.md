# Data Directory Overview

## Folder structure
This folder contains the raw and processed data used in this study and is organized as follows:
```
data/
├── README.md
├── breast_cancer_visium
├── crc_stereo
├── cutntag_processed
├── dlpfc_sodb
├── ductal_carcinoma
├── mbm_slideseqv2
├── mouse_brain_visium
├── prostate_ref_mapping
├── sodb_samples
├── synthetic_deconv
```

All data files can be downloaded from [10.5281/zenodo.10223862 (DOI)](https://zenodo.org/records/10223862).

## Data sources
Descriptions and links to all public datasets analyzed in this study can be found in the corresponding Methods sections of the [Smoother paper](https://www.biorxiv.org/content/10.1101/2022.10.25.513785v2.full). These include:
1. The single-nucleus RNA-seq data and the 10x Visium data of mouse brain: `data/mouse_brain_visium/`
    - Downloaded from ArrayExpress (https://www.ebi.ac.uk/arrayexpress/experiments/) under accession IDs E-MTAB-11114 and E-MTAB-11115.
    - Kleshchevnikov V, Shmatko A, Dann E, Aivazidis A, King HW, Li T, et al. Cell2location maps fine-grained cell types in spatial transcriptomics. Nat Biotechnol. 2022;40(5):661-+.
2. The sci-Space mouse embryonic data: `data/synthethic_deconv/sci_space/`
    - Downloaded from GEO at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166692
    - Srivatsan SR, Regier MC, Barkan E, Franks JM, Packer JS, Grosjean P, et al. Embryo-scale, single-cell spatial transcriptomics. Science. 2021;373(6550):111-+.
3. The 10x Visium invasive ductal carcinoma (IDC) datasets:
    1. Visium data with DAPI and anti-CD3 staining in Fig 3a-c: `data/ductal_carcinoma/10x_visium/`
        - Downloaded from 10x Genomics at https://support.10xgenomics.com/spatial-gene-expression/datasets
        - Zhao E, Stone MR, Ren X, Guenthoer J, Smythe KS, Pulliam T, et al. Spatial transcriptomics at subspot resolution with BayesSpace. Nat Biotechnol. 2021;39(11):1375-+.
    2. Single-cell reference data for deconvolution in Fig 3a-c: `data/ductal_carcinoma/scref/`
        - See https://www.nature.com/articles/s41588-021-00911-1#data-availability.
        - Wu SZ, Al-Eryani G, Roden DL, Junankar S, Harvey K, Andersson A, et al. A single-cell and spatially resolved atlas of human breast cancers. Nat Genet. 2021;53(9):1334-47. 
    3. Visium data for histology analysis and boundary visualization in Fig S2: `data/breast_cancer_visium/`
        - Downloaded from https://doi.org/10.5281/zenodo.4739739
        - Wu SZ, Al-Eryani G, Roden DL, Junankar S, Harvey K, Andersson A, et al. A single-cell and spatially resolved atlas of human breast cancers. Nat Genet. 2021;53(9):1334-47. https://www.nature.com/articles/s41588-021-00911-1#data-availability.

4. The spatial-CUT&Tag mouse embryonic data: `data/cutntag_processed/`
    - Raw data downloaded from GEO at https://www.ncbi.xyz/geo/query/acc.cgi?acc=GSE165217 with preprocessing scripts from Zenodo at https://zenodo.org/record/5797109#.Y9wezOyZP0p. 
    - For space limit, here we only provide the CUT&Tag data with gene activity score after preprocessing.
    - Deng YX, Bartosovic M, Kukanja P, Zhang D, Liu Y, Su G, et al. Spatial-CUT&Tag: Spatially resolved chromatin modification profiling at the cellular level. Science. 2022;375(6581):681-+.
5. The CRC Stereo-seq and paired scRNA-seq data: `data/crc_stereo/`
    - Downloaded from CNGB Nucleotide Sequence Archive under accession ID CNP0002432 at https://db.cngb.org/search/project/CNP0002432/. The human prostate Slide-seqV2 data (48) were obtained from https://github.com/shenglinmei/ProstateCancerAnalysis.
    - Zhang R, Feng Y, Ma W, Guo Y, Luo M, Li Y, et al. Spatial transcriptome unveils a discontinuous inflammatory pattern in proficient mismatch repair colorectal adenocarcinoma. Fundamental Research. 2022.
6. The human prostate Slide-seqV2 data: `data/prostate_ref_mapping/`
    - ST data downloaded from https://github.com/shenglinmei/ProstateCancerAnalysis.
    - Pretrained SCVI model and Tabula Sapiens data downloaded from https://huggingface.co/scvi-tools/tabula-sapiens-prostate-scvi.
    - Hirz T, Mei S, Sarkar H, Kfoury Y, Wu S, Verhoeven BM, et al. Dissecting the immune suppressive human prostate tumor microenvironment via integrated single-cell and spatial transcriptomic analyses. Nat Commun. 2023;14(1):663.
7. The DLPFC/SpatialLIBD SODB data: `data/dlpfc_sodb/`
    - Preprocessed ST data downloaded from https://gene.ai.tencent.com/SpatialOmics/
    - Maynard KR, Collado-Torres L, Weber LM, Uytingco C, Barry BK, Williams SR, et al. Transcriptome-scale spatial gene expression in the human dorsolateral prefrontal cortex. Nat Neurosci. 2021;24(3):425-36.
8. Other SODB datasets inspected in Fig S1: `data/sodb_samples/`
    - Preprocessed ST data downloaded from https://gene.ai.tencent.com/SpatialOmics/
    - See [scripts/loss_design/figs1_real_data_similarity.ipynb](../scripts/loss_design/figs1_real_data_similarity.ipynb) for details on each inspected data.
9. The Melanoma Brain Metastasis (MBM) Slide-seqV2 data: `data/mbm_slideseqv2/`
    - Obtained from in-house collaboration with the Izar lab.
    - Biermann J, Melms JC, Amin AD, Wang YP, Caprio LA, Karz A, et al. Dissecting the treatment-naive ecosystem of human melanoma brain metastasis. Cell. 2022;185(14):2591-+.

