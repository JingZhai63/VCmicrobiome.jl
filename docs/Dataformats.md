# Data Formats

The type of files involved in this package can be divided into two categories: input files and output file. The summary of different files are listed as

| Name       | Category           | Extension  | Description |
| ------------- |:-------------:|:-----:|-------------|
| Kernel distance matrix file      | Input | .csv | Csv file that contains the information of microbiome count and phylogenetic relationship|
| Kernel name list file      | Input | .csv | Csv file that contains the file names of the multiple kernel distance matrix files |
| Covariates file      | Input      |  .csv | Csv file that contains the values of covariates|
| Phenotype file | Input     |    .csv | Csv file that contains the values of phenotypes |
| Output file | Output     |    .out | Flat file that contains p-value for each group of OTUs under certain testing shceme|

# Kernel Distance Matrix file

UniFrac is a distance metric used for comparing biological communities. More information about the UniFrac distance is [here](https://en.wikipedia.org/wiki/UniFrac). UniFrac distance can be calculated using julia package [_PhylogeneticDistance.jl_](https://github.com/JingZhai63/PhylogeneticDistance.jl). Microbiome count file and phylogenetic tree file are needed. The distance matrix need be transformed to positive-definite kernel matrix.


# Covariates file

Covariates file is a csv file that contains the values of covariates. In covariates file, each line represents one measurement for one subject. The format of covariates file is: the first three fields is Subject ID, Measurement ID and Time Point, the rest fields are the covariates/predictors (e.g. sex, age, weight and etc). The intercept should NOT be included. Fields are separated by a comma. A typical covariates file look like

If a specific covariate value is missing, please use "NA" to indicate that missing value

If no covariates file is provided, the covariates matrix **X** will be automatically set to a _n_-by-1 matrix with all elements equal to 1, where n is the number of individuals.

```csv
"subjectID" "microbiomeID" "age" "sex" "smoker"X001	X001T1	57	0	1X001	X001T2	57	0	1X001	X001T3	57	0	1X002	X002T1	26	1	1X002	X002T2	26	1	1X002	X002T3	26	1	1X003	X003T1	27	0	1X003	X003T2	27	0	1X003	X003T3	27	0	1X004	X004T1	44	0	1X004	X004T2	44	0	1X004	X004T3	44	0	1X005	X005T1	38	0	1X005	X005T2	38	0	1X005	X005T3	38	0	1
```

# Phenotype file

Phenotype file is a csv file that contains the values of phenotypes. In phenotype file, each line represents one measurement at one time point. One subject may have several repeated measurements. The format of phenotype file is the same as covariate file. The fields in phenotype file are

* Subject ID
* Measurement ID
* Time Point
* Phenotypes

```csv
"subjectID" "microbiomeID	" "y1" "y2" "y3" "y4" "y5"X001	X001T1	25.47625505	31.68266284	27.63805982	28.65881349	23.93537075X001	X001T2	24.64079399	31.2348356	26.94997667	29.55297632	26.71868439X001	X001T3	24.93957563	30.80621089	30.21994675	28.9160571	26.61401479X002	X002T1	12.05899015	15.64851403	14.65772794	11.88751962	14.36242423X002	X002T2	14.33583163	14.00079888	16.12400259	11.93163584	15.0027127X002	X002T3	11.45724733	15.7380116	15.12183368	12.66567158	15.15858036X003	X003T1	12.93804869	9.977409278	13.75288042	14.32158349	12.65765104X003	X003T2	10.49035571	12.08767453	14.26815982	12.13905067	13.3291567X003	X003T3	12.42898212	10.59607919	13.9706309	11.73030896	12.47932868X004	X004T1	20.26107286	24.75744212	22.69600831	21.645191	20.18550653X004	X004T2	21.25620395	25.76615493	23.17473835	21.50946595	20.30114608X004	X004T3	20.65251182	26.15018955	22.12612909	20.72239947	20.5354966X005	X005T1	13.94408194	20.73833674	18.9831627	15.65224559	19.46986625X005	X005T2	15.2615871	21.50167569	17.06161367	15.91519641	20.53206818X005	X005T3	16.43634319	19.54563858	18.40833528	16.71131287	19.86249035
```

# Output file

Output file is a flat file that contains _p_-values and other information for each phenotypes under certain testing scheme. In output file, each line represents for one phenotype in phenotypes file(the first line is the header).
