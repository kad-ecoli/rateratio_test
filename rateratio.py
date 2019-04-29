#!/usr/bin/env python
docstring='''
rateratio.test: An Exact Rate Ratio Test Assuming Poisson Counts

Description
    Performs the uniformly most powerful unbiased test on the ratio of rates
    of two Poisson counts with given time (e.g., persons-years) at risk for
    each count. This module is a python re-implement of the rateratio.test
    R package.

Usage
    from rateratio import test as rateratio_test
    rateratio_test(x, n, RR = 1, 
        alternative = "two.sided",
        conf_level = 0.95)

Arguments
    x               A vector of length 2 with counts for the two rates.
    n               A vector of length 2 with time at risk in each rate.
    RR              The null rate ratio (two.sided) or the rate ratio on
                    boundary between null and alternative.
    alternative     A character string specifying the alternative hypothesis,
                    must be one of "two.sided" (default), "greater" or
                    "less".  You can specify just the initial letter.
    conf_level      Confidence level of the returned confidence interval.
                    Must be a single number between 0 and 1.

Details
    The rateratio.test tests whether the ratio of the first rate (estimated
    by x[0]/n[0] over the second rate (estimated by x[1]/n[1]) is either
    equal to, less, or greater than RR. The two-sided p-value is defined as
    either 1 or twice the minimum of the one-sided p-values.
    See Lehmann (1986, p. 152). For full discussion of the p-value and
    confidence interval consistency of inferences, see Fay (2010).

Value
    An object of class "htest" containing the following components:
    p_value         The p-value of the test. Due to numerical restriction,
                    the smallest p-value that can be returned is 2.22E-16.
    estimate        A vector with the rate ratio and the two individual rates
    null_value      the null rate ratio (two.sided) or the rate ratio on
                    boundary between null and alternative.
    conf_int        Confidence interval.
    alternative     Type of alternative hypothesis.
    methods         Description of method.
    data_name       Description of data.

References
    Fay, M. P. (2010). Two-sided exact tests and matching confidence
        intervals for discrete data. R Journal, 2(1), 53-58.
    Lehmann, E.L. (1986). Testing Statistical Hypotheses (second edition).
        Wadsworth and Brooks/Cole, Pacific Grove, California.

Examples
    from rateratio import test as rateratio_test
    print(rateratio_test((2,9), (17877,16660)))
'''

from scipy.stats import binom
pbinom = binom.cdf
from scipy.stats import beta
qbeta = beta.ppf
from numpy import Inf

class htest:
    def __init__(self, x, n, RR, alternative):
        self.p_value=0
        self.estimate=[] # "Rate Ratio","Rate 1","Rate 2"
        self.null_value=RR # RR "rate ratio"
        self.conf_int=[]
        self.alternative=alternative
        self.method="Exact Rate Ratio Test, assuming Poisson counts"
        self.data_name="data:  c(%d, %d) with time of c(%d, %d), null rate ratio "%(x[0],x[1],n[0],n[1])+str(RR)
        return

    def __repr__(self):
        hypothesis = "not equal to "+str(self.null_value)
        if   self.alternative == "greater":
            hypothesis = "greater than "+str(self.null_value)
        elif self.alternative == "less":
            hypothesis = "less than "+str(self.null_value)
        return '''
\t%s

%s
p-value = %s
alternative hypothesis: true rate ratio is %s
%s percent confidence interval:
 %s %s
sample estimates:
  Rate Ratio       Rate 1       Rate 2
%s\t%s\t%s
'''%(self.method,
    self.data_name,
    str(self.p_value),
    hypothesis,
    str(100*self.conf_int[-1]),
    self.conf_int[0],self.conf_int[1],
    self.estimate[0],self.estimate[1],self.estimate[2],
    )

    def __str__(self):
        return self.__repr__()

# Modify p.L and p.U from binom.test function
def p_L(x,n,alpha):
    if x == 0:
       return 0
    return qbeta(alpha, x, n - x + 1)

def p_U(x,n,alpha):
    if x == n:
       return 1
    return qbeta(1 - alpha, x + 1, n - x)

def test(x, n, RR = 1, alternative = "two.sided", conf_level = 0.95):
    '''
    x           - a vector of length 2 with counts for the two rates
    n           - a vector of length 2 with time at risk in each rate
    RR          - the null rate ratio (two.sided) or the rate ratio on
                  boundary between null and alternative
    alternative - a string specifying the alternative hypothesis, must be
                  one of "two.sided" (default), "greater" or "less".
                  You can specify just the initial letter.
    conf_level  - confidence level of the returned confidence interval.
                  Must be a single number between 0 and 1.
    '''
    # modify checks from prop.test
    k = len(x)
    assert(k == 2),     "x must have a length 2"
    assert(k == len(n)),"'x' and 'n' must have the same length"
    assert(min(n) > 0), "elements of 'n' must be positive"
    assert(min(x) >= 0),"elements of 'x' must be nonnegative"
    assert(RR > 0),     "RR must be greater than 0"
    
    #alternative <- match.arg(alternative)
    if   alternative == "t" or alternative == "two_sided":
        alternative = "two.sided"
    elif alternative == "g":
        alternative = "greater"
    elif alternative == "l":
        alternative = "less"
    assert(alternative in ("two.sided","less","greater")
        ),'''Error in match.arg(alternative) : 
  'arg' should be one of "two.sided", "less", "greater"'''

    assert(conf_level >0 and conf_level < 1
        ),"'conf.level' must be a single number between 0 and 1"

    RVAL=htest(x, n, RR, alternative)

    Y = x[0]
    N = n[0]
    X = x[1]
    M = n[1]
    RVAL.estimate = [ (1.*Y/N)/(1.*X/M), 1.*Y/N, 1.*X/M ]
    pRR =  (1.*N*RR)/(N*RR + M)
    pval_less =  pbinom(Y, X+Y, pRR)
    pval_greater = 1 - pbinom(Y-1, Y+X, pRR)

    if   alternative == "less":
        RVAL.p_value  = pval_less
        RVAL.conf_int = [
            0, (p_U(Y,X+Y,1-conf_level)*M)/(N*(1-p_U(Y,X+Y,1-conf_level) ))]
    elif alternative == "greater":
        RVAL.p_value  = pval_greater
        RVAL.conf_int = [
            (p_L(Y,X+Y,1-conf_level)*M)/(N*(1-p_L(Y,X+Y,1-conf_level))), Inf]
    elif alternative == "two.sided":
        RVAL.p_value  = min( (1, 2*min(pval_less,pval_greater) ) )
        RVAL.conf_int = [
            (p_L(Y,X+Y,(1-conf_level)/2)*M)/(N*(1-p_L(Y,X+Y,(1-conf_level)/2))), 
            (p_U(Y,X+Y,(1-conf_level)/2)*M)/(N*(1-p_U(Y,X+Y,(1-conf_level)/2)))]
    
    RVAL.conf_int.append(conf_level)
    return RVAL

if __name__=="__main__":
    #n = 17877
    #m = 16660
    #print(test((2,9), (n,m)))

    print(docstring)
