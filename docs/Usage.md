<style TYPE="text/css">
code.has-jax {font: inherit; font-size: 100%; background: inherit; border: inherit;}
</style>
<script type="text/x-mathjax-config">
MathJax.Hub.Config({
    tex2jax: {
        inlineMath: [['$','$'], ['\\(','\\)']],
        skipTags: ['script', 'noscript', 'style', 'textarea', 'pre'] // removed 'code' entry
    }
});
MathJax.Hub.Queue(function() {
    var all = MathJax.Hub.getAllJax(), i;
    for(i = 0; i < all.length; i += 1) {
        all[i].SourceElement().parentNode.className += ' has-jax';
    }
});
</script>
<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>

# Usage
The core function in package _VCTestMicrobiome.jl_ is MicrobiomeVCTest(), which wraps the three types of exact tests for overall microbiome association test. The usage for MicrobiomeVCTest() is
```julia
MicrobiomeVCTest(Name=Value)
```
Where `Name` is the option name and `Value` is the value of the corresponding option. Note that it's okay if some of the `Value` is not provided since some of the option has a default value.

There are two ways to call _MicrobiomeVCTest()_

*Open up Julia and type
```julia
julia> using VCTestMicrobiome
julia> MicrobiomeVCTest(Name = Value)
```
*From command line, type
```command
$ julia -E 'using VCTestMicrobiome; MicrobiomeVCTest(Name = Value)'
```
# Background
For longitudinal study design, phenotype (outcome variable of interest) and microbiome measurements are collected for a total of $n$ individuals. For each of the individual, they are measured repeatedly at $n_i$ time points. We consider a standard linear mixed model, 

$\bf{y}=\bf{X\beta}+\bf{Zb}+h(\bf{G})+\bf{\epsilon}$

$\bf{h} \sim N(\bf{0}, \sigma_m^2\bf{K}) $

$\bf{b} \sim  \mathcal{N}(0,\bf{Z}\sigma_D^2\bf{Z}')$

$\bf{\epsilon} \sim \mathcal{N}(0,\sigma_e^2\bf{I}) $

where $\bf{y}$ is a vector of continuous phenotype, $\bf{X}$ is covariate matrix, $\bf{\beta}$ is a vector of fixed effects, $h(\bf{G})$ is the random effects contributed by its microbiome profile, $\bf{b}$ is vector of random effects included to control correlation in the repeated measurements, and $\bf{\epsilon}$ is the error term.  It is assumed that $\bf{b}$ and $\bf{h}$ are independent with each other. We have variance-covariance structure 

$\mathrm{Var}(\bf{y})=\bf V = \bf  Z\sigma_D^2 \bf Z'+\sigma_m^2 \bf K +\sigma_e^2\bf {I_n}$

Therefore, $\sigma_D^2$, $\sigma_m^2$ and $\sigma_e^2$ are corresponding variance component parameters from correlation, microbiome and environmental effects.

# Specify input kernel distance matrix file

Option `kernelFile` indicates the file name (with extension) for the input kernel matrix files. 
```julia
MicrobiomeVCTest(kernelFile = "kernel.csv")
```

# Specify input covariates file

Option covFile indicates the file name for the input covariates file. 
```julia
MicrobiomeVCTest(covFile = "covariates.csv")
```
If option `covFile` is not specified, the covariates matrix **X** will be automatically set to a _n_-by-1 matrix with all elements equal to 1, where n is the number of individuals.

#Specify input phenotype file

Option `responseFile` indicates the file name for the input phenotype file. `yInit` indicates the initial column number of reading phenotype file. 
```julia
MicrobiomeVCTest(responseFile = "y.csv", yInit=3)
```


# Specify output file

Option `outFile` indicates the file name for the output file. If the output file name is set to test.out, then use
```julia
MicrobiomeVCTest(outFile = "/PATH/OF/Results.out")
```
If option `outFile` is not specified, the file will be automatically named with phenotype file and type of the exact test. Also, it will be stored at the current directory. 

# Specify type of mixed model 

Option `ZtZ` indicates the type of mixed model. Typically, a random intercept or random intercept and random slope model is most often
used. Now, we simply consider a random intercept model in longitudinal study. If it's a longitudinal microbiome study, then use
```julia
MicrobiomeVCTest(ZtZ= "intercept")
```
If the it's a cross-sectional study, then use
```julia
MicrobiomeVCTest(ZtZ= "none")
```
The default value for option `ZtZ` is intercept.	

# Choose testing scheme

Package **VCTestMicrobiome.jl** utilizes three types of exact tests in  [VarianceComponentTest.jl](https://github.com/Tao-Hu/VarianceComponentTest.jl/blob/master/README.md): exact likelihood ratio test (eLRT), exact restricted likelihood ratio test (eRLRT) and exact score test (eScore). Option `test` indicates which testing scheme you want to perform. The usage is

* MicrobiomeVCTest(test = "eLRT"): perform exact likelihood ratio test
* MicrobiomeVCTest(test = "eRLRT"): perform exact restricted likelihood ratio test
* MicrobiomeVCTest(test = "eScore"): perform exact score test

The default value for option test is eLRT.

# Choose method for computing p-value

Option `pvalueComputing` indicates which method will be used to obtain the _p_-value under null hypothesis. The usage is

MicrobiomeVCTest(pvalueComputing = "MonteCarlo"): use a Monte Carlo method by generating many replicates to obtain the exact null distribution of test statistic and then compute the p-value.
MicrobiomeVCTest(pvalueComputing = "chi2"): use a mixed $\chi^2$ distribution to approximate the null distribution of test statistic and then compute the _p_-value.
The default value for option `pvalueComputing` is _chi2_. The approximation effect of mixed $\chi^2$  distribution has been showed to be good enough, and such approximation will be faster than Monte Carlo method since much less replicates need to be generated. So please use _chi2_ whenever possible.

Note:

Option `pvalueComputing` is only valid for eLRT and eRLRT, since this software employs the method of inverting characteristics function to compute the _p_-value for eScore.
In the approximation method, the exact null distributions of eLRT and eRLRT are approximated by a mixture of the form

$\pi_0\chi_0^2:(1-\pi_0)a\chi_b^2$

where the point mass $\pi_0$, scale parameter $a$, and the degree of freedom $b$ for the $\chi^2$  distribution. First, estimate $\pi_0$ by generating B replicates (the number of replicates B is specified by option `nNullSimPts`), and then estimate $a$ and $b$ using only a small number (300 by default) of replicates.

# Multiple microbiome clusters association test
For testing multiple microbiome clusters, you need to use function MMicrobiomeTest(). Option `KernelLists` need to be specified in this case. Also, you need to set up the number of total observations `nObs`.
```julia
MMicrobiomeVCTest(nObs=150,KernelLists="MultiKerName.csv",kernelFile="",covFile="covariates.csv",responseFile="y.csv",yInit=3,test="eRLRT",ZtZ="intercept")
```

# Choose distance type for computing kernel matrix 

Kernel matrix should be generated in `R`. R package **GUniFrac** provides four types of distance 

* unweighted UniFrac distance: `KernelMatrix(d="d_uw")`
* weighted UniFrac distance: `KernelMatrix(d="d_1")`
* variance adjusted weighted UniFrac distance: `KernelMatrix(d="d_VAW")`
* generalized UniFrac distance: `KernelMatrix(d="d_alpha")`













