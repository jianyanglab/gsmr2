# GSMR2
 
GSMR2 (Generalised Summary-data-based Mendelian Randomisation v2) is an improved version of [GSMR](https://github.com/jianyanglab/gsmr), which uses GWAS summary statistics to test for a putative causal association between two phenotypes (e.g., a modifiable risk factor and a disease) based on a multi-SNP model

## Installation

### Packages dependency

Please install the following package first. 

```{r}
# survey
install.packages('survey')
```

### Installation of gsmr2

Please run the following to install the `gsmr2` package:
```
devtools::install_github("jianyanglab/gsmr2")
```

or download the package [here](https://yanglab.westlake.edu.cn/software/gsmr/static/gsmr_1.1.1.tar.gz) and install it using the following command
```
install.packages("~/Downloads/gsmr_1.1.1.tar.gz", type = "source", repo = NULL)
```

## Example analysis pipeline

Here is an example of the `gsmr2` workflow to get started:

```
Under construction
```

# GCTA version
We also implement the C++ version of GSMR/GSMR2 within the [GCTA](https://github.com/jianyangqt/gcta). 

# Citation

**Angli Xue**, Zhihong Zhu, Huanwei Wang, Longda Jiang, Peter M. Visscher, Jian Zeng, Jian Yang. Unravelling the complex causal effects of substance use behaviours on common diseases. *Under 2nd review*. [Preprint coming soon]

For questions, please email us at Jian Yang (jian.yang@westlake.edu.cn)

For bugs, please raise an issue in this repository and assign it to [@anglixue](https://github.com/anglixue).
