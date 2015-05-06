Plant species traits and growth rates meta-analysis
--------

This repository contains all the code used in the manuscript:

* Title: "Plant species traits and growth rates: meta-analysis shows correlations change with plant size as predicted"
* Authors: Anais Gibert, Emma F. Gray, Mark Westoby,  Ian J. Wright, and Daniel S. Falster,
* Year of publication:
* doi:

Synopsis of the study
--------
The literature is inconsistent about empirical correlations between functional traits and plant growth rate (GR), casting doubt on the capacity of traits to predict growth.
Traits should influence growth in a size-dependent manner. We outline mechanisms and hypotheses based on new theory, and test these predictions for five traits using a meta-analysis of 108 studies (> 500 correlations).
Results were consistent with predictions. Specific leaf area was correlated with GR in small but not large plants. Correlations of GR with wood density and assimilation rate were not affected by size. Maximum height and seed mass were correlated with GR only in one plant size category.
We show that correlations between traits and GR change in a predictable way as a function of plant size. Our understanding of plant strategies should shift away from attributing slow vs fast growth to species throughout life, in favour of attributing growth trajectories.

Running the code
--------

Here we present the data and the code to perform the meta-analyses and all the figures from the paper. Once you have the required packages installed, run this command:

```
source("analysis.R")
```

Figures will be output to a directory `output`. Also saved in this file is a homogenized dataset, containing the data used in the analyses and included as supplementary material with the paper.

List of files available and explanation
--------

- `data/CompileData.csv`: raw data, needed to run the analyses
- `data/CompileData_meta.csv`:definition of columns in `data/ComplieData.csv`
- `R` directory containing functions used in analysis
- `analysis.R`: main script to run the analyses and generate all the figures and tables.
- `MS.tex`: manuscript in LaTex
- `references/complete.bib`: bibtex file with all references used in the meta-analyses and in the manuscript
- `references/meta-analyses.bib`: bibtex file with all references used in the meta-analyses
- `references/read.bib`: bibtex file with all references red to do the meta-analyses (all the studies used + studies red but discard from our meta-analyses)
- `ecol_let.bst`: latex style file used for formatting paper

Contributors
------------------------
Daniel Falster
Anais Gibert
