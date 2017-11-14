
# Variance Component Test for Longitudinal Microbiome Study

**VCmicrobiome.jl** is a Julia package for performing exact variance component tests in longitudinal microbiome study. It utilizes three types of exact tests from Tao-Hu's Julia package [VarianceComponentTest.jl](https://github.com/Tao-Hu/VarianceComponentTest.jl/blob/master/README.md)

* exact likelihood ratio test (eLRT)
* exact restricted likelihood ratio test (eRLRT)
* exact score test (eScore)

The input files for _VCmicrobiome.jl_ are microbiome `kernel distance matrix (.csv file)`, `covariates file (.csv file)` and `phenotype file (.csv file)`. You should have these input files prepared before running our program. The `output file (.out file)` is a simple comma-delimited file containing the _p_-values for each phenotype and each group of OTUs under a certain testing scheme (eLRT, eRLRT or eScore).

To use **VCmicrobiome.jl** , you need to call the `MicrobiomeVCTest()` function.

Please see the detailed documents in [Read the Docs](http://vcmicrobiomejl.readthedocs.io/en/latest/)

# Installation

To install _VCMicrobiome.jl_, open up Julia and then type

```julia
julia> Pkg.clone("https://github.com/JingZhai63/VCmicrobiome.jl.git")
```
