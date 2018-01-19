# Examples

This example illustrates the usages of different options by analyzing a simulated data set obtained from Pulmonary Microbiome Study. This demo data set contains 2964 OTUs and 50 individuals. The simulated microbiome data measured at three time points. The required/optional input files are

* Kernel matrix file: kernel.csv, kernel_baseline.csv
* Covariates file: covariates.csv, covariates_baseline.csv (3 covariates)
* Phenotye file: y.csv, y_baseline.csv

These data files come with our package, and they are available at PATH/TO/PACKAGE/DOCS/EXAMPLES.

# Basic usage

Three types of exact tests can be performed. Open up a Julia session and type

* exact likelihood ratio test (eLRT)
```julia
julia> using VCmicrobiome
julia> microvctest(kernelFile = "kernel.csv", covFile = "covariates.csv", responseFile = "y.csv", test = "eLRT", out = true)
```

* exact restricted likelihood ratio test (eRLRT)
```julia
julia> using VCmicrobiome
julia> microvctest(kernelFile = "kernel.csv", covFile = "covariates.csv", responseFile = "y.csv", test = "eRLRT", yIdx = 3)
```
* exact score test (eScore)
```julia
julia> using VCmicrobiome
julia> microvctest(kernelFile = "kernel.csv", covFile = "covariates.csv", responseFile = "y.csv", test = "eScore", yIdx = 3)
```
Note: if the input files are not at the current directory, you should specify the paths correctly.

You can also call **microvctest** from command line. For example, to perform eRLRT

```command
$ julia -E 'using VCmicrobiome; microvctest(kernelFile = "kernel.csv", covFile = "covariates.csv", responseFile = "y.csv", test = "eScore", yIdx = 3)
```

# Option `longitudinal`
If the study is not longitudinal designed or you want to analyze data for one time point, then
```julia
julia> using VCmicrobiome
julia> microvctest(kernelFile = "kernel_baseline.csv", covFile = "covariates_baseline.csv", responseFile = "y_baseline.csv", longtitudinal = false, test = "eRLRT", yIdx = 3)
```

# Option `fine`
This option is specified for localizing fine microbiome cluster effects. If you want to adjust for effect contributed by related cluster, then

```julia
julia> using VCmicrobiome
julia> microvctest(kernelFile = "kernel.csv", kadjFile = "kernel_adj.csv", fine = true ,covFile = "covariates.csv", responseFile = "y.csv", test = "eScore", yIdx = 3)
```

If the rank of V1 has high or full rank, the package will evoke the low rank approximation with default `lowRank = 0.4`. For example, if `fine = true`, the _microvctest_ will perform low rank approximation since microbiome kernel matrix usually has full rank. 

# Option `pvalueComputing`
Chi squared approximation is recommended (though you don't have to write it out specifically)
```julia
julia> using VCmicrobiome
julia> microvctest(kernelFile = "kernel_baseline.csv", covFile = "covariates_baseline.csv", responseFile = "y_baseline.csv", longitudinal = false, test = "eRLRT", pvalueComputing = "chi2")
```

Note: Option pvalueComputing is only for eLRT and eRLRT. For eScore, it employs the method of inverting characteristics function to compute the _p_-value.

# Option `out` and `outFile`

Users can specify if the test result will be written to ".out" file and the name of the outfile. The default of `out` is false, which the _microvctest_ return the value without writing the out file. The output will be written to "y-eLRT.out" at the working directory by default.
