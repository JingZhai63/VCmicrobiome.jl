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
The core function in package _VCmicrobiome.jl_ is microvctest(), which wraps the three types of exact tests for overall microbiome association test. The usage for microvctest() is
```julia
microvctest(Name=Value)
```
Where `Name` is the option name and `Value` is the value of the corresponding option. Note that it's okay if some of the `Value` is not provided since some of the option has a default value.

There are two ways to call _microvctest()_

* Open up Julia and type
```julia
julia> using VCmicrobiome
julia> microvctest(Name = Value)
```
* From command line, type
```command
$ julia -E 'using VCmicrobiome; microvctest(Name = Value)'
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
microvctest(kernelFile = "kernel.csv")
```

# Specify input covariates file

Option covFile indicates the file name for the input covariates file.
```julia
microvctest(covFile = "covariates.csv")
```
If option `covFile` is not specified, the covariates matrix **X** will be automatically set to a _n_-by-1 matrix with all elements equal to 1, where n is the number of individuals.

#Specify input phenotype file

Option `responseFile` indicates the file name for the input phenotype file. `yInit` indicates the initial column number of reading phenotype file.
```julia
microvctest(responseFile = "y.csv", yIdx = 3)
```


# Specify output file

Option `outFile` indicates the file name for the output file. If the output file name is set to test.out, then use
```julia
microvctest(outFile = "/PATH/OF/Results.out")
```
If option `outFile` is not specified, the file will be automatically named with phenotype file and type of the exact test. Also, it will be stored at the current directory.

# Specify type of microbiome study

Option `longitudinal` indicates the type of microbiome study. We simply consider a random intercept model in longitudinal study. If it's a longitudinal microbiome study, then use
```julia
microvctest(longitudinal = true)
```
If the it's a cross-sectional study, then use
```julia
microvctest(longitudinal = false)
```
The default value for option `longitudinal` is ``true``.

# Choose testing scheme

Package **VCmicrobiome.jl** utilizes three types of exact tests: exact likelihood ratio test (eLRT), exact restricted likelihood ratio test (eRLRT) and exact score test (eScore). Option `test` indicates which testing scheme you want to perform. The usage is

* microvctest(test = "eLRT"): perform exact likelihood ratio test
* microvctest(test = "eRLRT"): perform exact restricted likelihood ratio test
* microvctest(test = "eScore"): perform exact score test

The default value for option test is eLRT.

# Choose method for computing p-value

Option `pvalueComputing` indicates which method will be used to obtain the _p_-value under null hypothesis. The usage is

microvctest(pvalueComputing = "MonteCarlo"): use a Monte Carlo method by generating many replicates to obtain the exact null distribution of test statistic and then compute the p-value.
microvctest(pvalueComputing = "chi2"): use a mixed $\chi^2$ distribution to approximate the null distribution of test statistic and then compute the _p_-value.
The default value for option `pvalueComputing` is _chi2_. The approximation effect of mixed $\chi^2$  distribution has been showed to be good enough, and such approximation will be faster than Monte Carlo method since much less replicates need to be generated. So please use _chi2_ whenever possible.

Note:

Option `pvalueComputing` is only valid for eLRT and eRLRT, since this software employs the method of inverting characteristics function to compute the _p_-value for eScore.
In the approximation method, the exact null distributions of eLRT and eRLRT are approximated by a mixture of the form

$\pi_0\chi_0^2:(1-\pi_0)a\chi_b^2$

where the point mass $\pi_0$, scale parameter $a$, and the degree of freedom $b$ for the $\chi^2$  distribution. First, estimate $\pi_0$ by generating B replicates (the number of replicates B is specified by option `nNullSimPts`), and then estimate $a$ and $b$ using only a small number (300 by default) of replicates.

# Choose distance type for computing kernel matrix

Kernel matrix should be generated in `Julia` package `PhylogeneticDistance.jl`. It provides four types of distance

* unweighted UniFrac distance: `KernelMatrix(Dtype = "d_UW")`
* weighted UniFrac distance: `KernelMatrix(Dtype = "d_alpha", alpha = 1.0)`
* variance adjusted weighted UniFrac distance: `KernelMatrix(Dtype = "d_VAW")`
* generalized UniFrac distance: `KernelMatrix(Dtype = "d_alpha", alpha = 0.5)`

For more details, please see the documents of  [PhylogeneticDistance.jl](http://phylogeneticdistancejl.readthedocs.io/en/latest/). 
