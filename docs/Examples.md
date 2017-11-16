#Examples

This example illustrates the usages of different options by analyzing a simulated data set obtained from Pulmonary Microbiome Study. This demo data set contains 2964 OTUs and 50 individuals. The simulated microbiome data measured at three time points. The required/optional input files are

* Microbiome tree file: OTUs_tree.txt
* Microbiome count table: count.csv
* Kernel matrix file: kernel.csv, kernel_baseline.csv, MultiKerName.csv
* Covariates file: covariates.csv, covariates_baseline.csv which contains 3 covariates
* Phenotye file: y.csv, y_baseline.csv

These data files come with our package, and they are available at here.

# Basic usage

Three types of exact tests can be performed. Open up a Julia session and type

* exact likelihood ratio test (eLRT)
```julia
julia> using VCTestMicrobiome
julia> MicrobiomeVCTest(kernelFile = "kernel.csv", covFile = "covariates.csv", responseFile = "y.csv", test = "eLRT")
```
Then the output will be written to "y_baseline_eLRT.out" at the same directory as input files.

* exact restricted likelihood ratio test (eRLRT)
```julia
julia> using VCTestMicrobiome
julia> MicrobiomeVCTest(kernelFile = "kernel.csv", covFile = "covariates.csv", responseFile = "y.csv", test = "eRLRT",yInit=3)
```
* exact score test (eScore)
```julia
julia> using VCTestMicrobiome
julia> MicrobiomeVCTest(kernelFile = "kernel.csv", covFile = "covariates.csv", responseFile = "y.csv", test = "eScore",yInit=3)
```
Note: if the input files are not at the current directory, you should specify the paths correctly.

You can also call MicrobiomeVCTest from command line. For example, to perform eRLRT

```command
$ julia -E 'using VCTestMicrobiome; MicrobiomeVCTest(kernelFile = "kernel.csv", covFile = "covariates.csv", responseFile = "y.csv", test = "eScore",yInit=3)
```

# Option `ZtZ`
If the study is not longitudinal designed or you want to analyze data for one time point, then
```julia
julia> using VCTestMicrobiome
julia> MicrobiomeVCTest(kernelFile = "kernel_baseline.csv", covFile = "covariates_baseline.csv", responseFile = "y_baseline.csv", ZtZ = "none", test = "eRLRT", yInit=3)
```

# Option `pvalueComputing`
Chi squared approximation is recommended (though you don't have to write it out specifically)
```julia
julia> using VCTestMicrobiome
julia> MicrobiomeVCTest(kernelFile = "kernel_baseline.csv", covFile = "covariates_baseline.csv", responseFile = "y_baseline.csv", ZtZ = "none", test = "eRLRT", pvalueComputing = "chi2")
```

Note: Option pvalueComputing is only for eLRT and eRLRT. For eScore, it employs the method of inverting characteristics function to compute the _p_-value.

# Option `KernelLists` and `kernel`
OTUs can be grouped to higher phylogenetic rank such as genus, class, family, phylum _etc_. Therefore, multiple kernel matrix can be calculated based on phylogenetic groups. Function `MMicrobiomeVCTest()` is used for testing association of multiple microbiome clusters.

Specifically, `KernelLists`  indicates a csv file containing the name of the multiple kernel matrix files and `kernel` indicates there are more than one microbiome cluster to be test.
```julia
julia> using VCTestMicrobiome
julia> MMicrobiomeVCTest(nObs=150, KernelLists = "MultiKerName.csv", kernel="multi", covFile = "covariates.csv", responseFile = "y.csv", ZtZ = "intercept", test = "eRLRT")
```
