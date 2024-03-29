## rateratio.test: An Exact Rate Ratio Test Assuming Poisson Counts ##

#### Description ####
Performs the uniformly most powerful unbiased test on the ratio of rates
of two Poisson counts with given time (e.g., persons-years) at risk for
each count. This module is a python re-implement of the 
[rateratio.test](https://cran.r-project.org/package=rateratio.test)
R package.

#### Usage ####
```python
from rateratio import test as rateratio_test
rateratio_test(x, n, RR = 1, 
    alternative = "two.sided",
    conf_level = 0.95)
```

#### Arguments ####
```
x               A vector of length 2 with counts for the two rates.
n               A vector of length 2 with time at risk in each rate.
RR              The null rate ratio (two.sided) or the rate ratio on
                boundary between null and alternative.
alternative     A character string specifying the alternative hypothesis,
                must be one of "two.sided" (default), "greater" or
                "less".  You can specify just the initial letter.
conf_level      Confidence level of the returned confidence interval.
                Must be a single number between 0 and 1.
```

#### Details ####
The rateratio.test tests whether the ratio of the first rate (estimated
by x[0]/n[0] over the second rate (estimated by x[1]/n[1]) is either
equal to, less, or greater than RR. The two-sided p-value is defined as
either 1 or twice the minimum of the one-sided p-values.
See Lehmann (1986, p. 152). For full discussion of the p-value and
confidence interval consistency of inferences, see Fay (2010).

#### Value ####
An object of class "htest" containing the following components:
```
p_value         The p-value of the test.  Due to numerical restriction,
                the smallest p-value that can be returned is 2.22E-16.
estimate        A vector with the rate ratio and the two individual rates
null_value      the null rate ratio (two.sided) or the rate ratio on
                boundary between null and alternative.
conf_int        Confidence interval.
alternative     Type of alternative hypothesis.
methods         Description of method.
data_name       Description of data.
```

#### References ####
* Wei, X., Zhang, C., Freddolino, P. L., & Zhang, Y. (2020).
  [Detecting Gene Ontology misannotations using taxon-specific rate ratio
  comparisons](https://doi.org/10.1093/bioinformatics/btaa548). Bioinformatics, 36(16), 4383-4388.
* Fay, M. P. (2010). [Two-sided exact tests and matching confidence
  intervals for discrete data](https://doi.org/10.32614/RJ-2010-008). R Journal, 2(1), 53-58.

#### Examples ####
The following python code
```python
from rateratio import test as rateratio_test
print(rateratio_test((2,9), (17877,16660)))
```
will output
```

	Exact Rate Ratio Test, assuming Poisson counts

data:  c(2, 9) with time of c(17877, 16660), null rate ratio 1
p-value = 0.0501064956141
alternative hypothesis: true rate ratio is not equal to 1
95.0 percent confidence interval:
 0.0217740627591 1.00054910449
sample estimates:
  Rate Ratio       Rate 1       Rate 2
0.207094155743	0.000111875594339	0.000540216086435
```

#### Dependency ####
[scipy](https://www.scipy.org/)
