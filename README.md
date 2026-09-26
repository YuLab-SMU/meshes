# meshes: MeSH Enrichment and Semantic analyses

[![](https://img.shields.io/badge/release%20version-1.38.0-green.svg)](https://www.bioconductor.org/packages/meshes)
[![](https://img.shields.io/badge/devel%20version-1.39.2-green.svg)](https://github.com/YuLab-SMU/meshes)
[![Bioc](http://www.bioconductor.org/shields/years-in-bioc/meshes.svg)](https://www.bioconductor.org/packages/devel/bioc/html/meshes.html#since)
[![platform](http://www.bioconductor.org/shields/availability/devel/meshes.svg)](https://www.bioconductor.org/packages/devel/bioc/html/meshes.html#archives)
[![Build
Status](http://www.bioconductor.org/shields/build/devel/bioc/meshes.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/meshes/)

Provides tools for MeSH (Medical Subject Headings) enrichment and
semantic similarity analysis. The package uses MeSH term annotations
linked to Entrez Gene identifiers through gendoo, gene2pubmed and RBBH,
enabling over-representation analysis, gene set enrichment analysis and
semantic comparison of gene lists or ranked expression profiles. It
implements Resnik, Schlicker, Jiang, Lin and Wang similarity methods for
comparing MeSH terms, genes and gene groups across more than 70 species.

For details, please visit
<https://yulab-smu.top/biomedical-knowledge-mining-book/>.

## :arrow_double_down: Installation

Get the released version from Bioconductor:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("meshes")
```

Or install the development version from GitHub:

``` r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("YuLab-SMU/meshes")
```

Citation: ***G Yu***. Using meshes for MeSH term enrichment and semantic
analyses. ***Bioinformatics*** 2018, 34(21):3766-3767.
<https://dx.doi.org/10.1093/bioinformatics/bty410>.

## :writing_hand: Authors

Guangchuang YU <https://yulab-smu.top>

School of Basic Medical Sciences, Southern Medical University

## :sparkling_heart: Contributing

We welcome any contributions! By participating in this project you agree
to abide by the terms outlined in the [Contributor Code of
Conduct](CONDUCT.md).
