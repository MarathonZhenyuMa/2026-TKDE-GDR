# BGDR & CGDR for Multi-View Clustering

This repository contains the MATLAB implementation of the accepted paper:

Zhenyu Ma, Shengzhao Guo, Jingyu Wang, Feiping Nie and Xuelong Li, *Scalable Graph Discrete Reconstruction for Efficient Multi-View Clustering*, IEEE Transactions on Knowledge and Data Engineering (TKDE).

The code implements two efficient multi-view clustering methods:

- **BGDR**: Bipartite Graph Discrete Reconstruction
- **CGDR**: Compact Graph Discrete Reconstruction

Both methods directly learn discrete cluster indicators from multi-view data and use anchor-based graph construction for scalable clustering.

## Folder Structure

```text
.
+-- demo.m                  # Demo script for reproducing the paper settings
+-- MvC_BGDR.m              # Main function of BGDR
+-- MvC_CGDR.m              # Main function of CGDR
+-- Data/                   # Multi-view data sets in .mat format
+-- subfunc_BGDR/           # Shared utilities and BGDR subfunctions
`-- subfunc_CGDR/           # CGDR subfunctions
```

Each data file is expected to contain:

- `X`: a `1 x V` cell array, where `X{v}` is the feature matrix of the `v`-th view
- `label`: the ground-truth label vector, used only for evaluation


## Quick Start

Open MATLAB, enter this folder, and run:

```matlab
demo
```

The script runs BGDR and CGDR on the included data sets using the parameter settings reported in the paper, and prints clustering metrics in the command window:

```text
Dataset      Method  ACC     NMI    Purity   Precision   Recall  F-score  ARI   Time(s)
-----------------------------------------------------
MSRC_v1      BGDR    ...
MSRC_v1      CGDR    ...
```

## Usage

BGDR:

```matlab
[result,F,G,alpha,t,Obj,converge] = MvC_BGDR(X,label,k,h,initLabel,isNormal,maxIter);
```

CGDR:

```matlab
[result,F,G,alpha,t,Obj,converge] = MvC_CGDR(X,label,k,h,initLabel,isNormal,maxIter);
```

Arguments:

- `X`: multi-view data cell array
- `label`: ground-truth labels for evaluation
- `k`: number of bipartite graph neighbors, denoted as `r` in the paper
- `h`: hierarchy depth for anchor generation, with anchor number `m = 2^h`
- `initLabel`: initialization method, default is `'N2HI'`
- `isNormal`: whether to perform row-wise normalization, default is `1`
- `maxIter`: maximum number of iterations, default is `30`

Returned `result` contains seven clustering metrics:

```text
[ACC, NMI, Purity, Precision, Recall, F-score, ARI]
```

## Paper Parameter Settings

The demo script uses the following settings.

| Data set | Data file | h | m | BGDR k/r | CGDR k/r |
|---|---|---:|---:|---:|---:|
| MSRC_v1 | `MSRC_v1_data.mat` | 7 | 128 | 36 | 28 |
| Dermatology | `Dermatology_data.mat` | 8 | 256 | 38 | 32 |
| 100leaves | `100leaves_data.mat` | 10 | 1024 | 30 | 16 |
| mnist4 | `mnist4_data.mat` | 10 | 1024 | 14 | 10 |
| Digit4k | `Digit4k_data.mat` | 10 | 1024 | 40 | 40 |
| Hdigit | `Hdigit_data.mat` | 12 | 4096 | 40 | 32 |
| ALOI | `ALOI_data.mat` | 13 | 8192 | 8 | 8 |
| MNIST | `MNIST6w_data.mat` | 8 | 256 | 32 | 38 |

## Citation

If you find this code useful, please cite:

```bibtex
@ARTICLE{Ma2026GDR,
  author={Ma, Zhenyu and Guo, Shengzhao and Wang, Jingyu and Nie, Feiping and Li, Xuelong},
  journal={IEEE Transactions on Knowledge and Data Engineering}, 
  title={Scalable Graph Discrete Reconstruction for Efficient Multi-View Clustering}, 
  year={2026},
  volume={38},
  number={7},
  pages={4641-4657},
  doi={10.1109/TKDE.2026.3682510}
  }
```

## Contact

For questions about the code or experiments, please contact <u>zhenyu.ma@mail.nwpu.edu.cn</u>
