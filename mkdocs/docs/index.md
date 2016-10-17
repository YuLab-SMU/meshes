<!-- AddToAny BEGIN -->
<div class="a2a_kit a2a_kit_size_32 a2a_default_style">
<a class="a2a_dd" href="//www.addtoany.com/share"></a>
<a class="a2a_button_facebook"></a>
<a class="a2a_button_twitter"></a>
<a class="a2a_button_google_plus"></a>
<a class="a2a_button_pinterest"></a>
<a class="a2a_button_reddit"></a>
<a class="a2a_button_sina_weibo"></a>
<a class="a2a_button_wechat"></a>
<a class="a2a_button_douban"></a>
</div>
<script async src="//static.addtoany.com/menu/page.js"></script>
<!-- AddToAny END -->

<link rel="stylesheet" href="https://guangchuangyu.github.io/css/font-awesome.min.css">

[![releaseVersion](https://img.shields.io/badge/release%20version-1.4.20-blue.svg?style=flat)](https://bioconductor.org/packages/ggtree)
[![develVersion](https://img.shields.io/badge/devel%20version-0.99.7-blue.svg?style=flat)](https://github.com/GuangchuangYu/ggtree)
[![total](https://img.shields.io/badge/downloads-15283/total-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/ggtree)
[![month](https://img.shields.io/badge/downloads-1678/month-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/ggtree)


MeSH (Medical Subject Headings) is the NLM controlled vocabulary used to manually index articles for MEDLINE/PubMed.
MeSH terms were associated by Entrez Gene ID by three methods, `gendoo`, `gene2pubmed` and `RBBH`.
This association is fundamental for enrichment and semantic analyses.
`meshes` supports enrichment analysis (over-representation and gene set enrichment analysis) of gene list or whole expression profile.
The semantic comparisons of MeSH terms provide quantitative ways to
compute similarities between genes and gene groups. meshes implemented five
methods proposed by `Resnik`, `Schlicker`, `Jiang`, `Lin` and `Wang` respectively and
supports more than 70 species.

`meshes` is released within the [Bioconductor](https://bioconductor.org/packages/meshes/) project and the source code is hosted on <a href="https://github.com/GuangchuangYu/meshes"><i class="fa fa-github fa-lg"></i> GitHub</a>.


## <i class="fa fa-user"></i> Author

Guangchuang Yu, School of Public Health, The University of Hong Kong.


## <i class="fa fa-download"></i> Installation

Install `meshes` is easy, follow the guide on the [Bioconductor page](https://bioconductor.org/packages/meshes/):

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
## biocLite("BiocUpgrade") ## you may need this
biocLite("meshes")
```


## <i class="fa fa-cogs"></i> Overview

#### <i class="fa fa-angle-double-right"></i> Enrichment Analysis

+ Over-representation analysis
+ Gene set enrichment analysis

#### <i class="fa fa-angle-double-right"></i> Semantic Analysis

+ MeSH term semantic similarity
+ Gene semantic similarity

#### <i class="fa fa-angle-double-right"></i> Visualization

+ barplot
+ cnetplot
+ dotplot
+ enrichMap
+ gseaplot
+ upsetplot
+ simplot


<i class="fa fa-hand-o-right"></i> Find out details and examples on <i class="fa fa-book"></i> [Documentation](https://guangchuangyu.github.io/meshes/documentation/).

## <i class="fa fa-wrench"></i> Related Tools

<ul class="fa-ul">
	<li><i class="fa-li fa fa-angle-double-right"></i><a href="https://guangchuangyu.github.io/clusterProfiler">clusterProfiler</a> for Ontologies/pathways enrichment analyses</li>
	<li><i class="fa-li fa fa-angle-double-right"></i><a href="https://guangchuangyu.github.io/DOSE">DOSE</a> for Disease Ontology Semantic and Enrichment analyses</li>
	<li><i class="fa-li fa fa-angle-double-right"></i><a href="https://guangchuangyu.github.io/GOSemSim">GOSemSim</a> for GO semantic similarity measurement</li>
	<li><i class="fa-li fa fa-angle-double-right"></i><a href="https://guangchuangyu.github.io/ReactomePA">ReactomePA</a> for Reactome pathway analysis</li>
</ul>

<i class="fa fa-hand-o-right"></i> Find out more on [projects](https://guangchuangyu.github.io/#projects).


## <i class="fa fa-comment"></i> Feedback
<ul class="fa-ul">
	<li><i class="fa-li fa fa-hand-o-right"></i> Please make sure you [follow the guide](https://guangchuangyu.github.io/2016/07/how-to-bug-author/) before posting any issue/question</li>
	<li><i class="fa-li fa fa-bug"></i> For bugs or feature requests, please post to <i class="fa fa-github-alt"></i> [github issue](https://github.com/GuangchuangYu/meshes/issues)</li>
	<li><i class="fa-li fa fa-question"></i>  For user questions, please post to [Bioconductor support site](https://support.bioconductor.org/) and [Biostars](https://www.biostars.org/). We are following every post tagged with **meshes**</li>
	<li><i class="fa-li fa fa-commenting"></i> Join the group chat on <a href="https://twitter.com/hashtag/meshes"><i class="fa fa-twitter fa-lg"></i></a> and <a href="http://huati.weibo.com/k/meshes"><i class="fa fa-weibo fa-lg"></i></a></li>
</ul>

