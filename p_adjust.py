#!/usr/bin/env python
docstring='''
p.adjust: P-values for Multiple Comparisons

Description
    Given a set of p-values, return q-values (i.e. p-values adjusted).
    This module is a python re-implement of p.adjust R package.

Usage
    from p_adjust import adjust as p_adjust
    p_adjust(p, method, n)

Arguments
    p           numpy array of p-values
    method      correction method. default: fdr
    n           number of comparisons, must be at least len(p); only set
                this (to non-default) when you know what you are doing!
                default: len(p)

Details
    The adjustment methods include the Bonferroni correction ("bonferroni")
    in which the p-values are multiplied by the number of comparisons. Less
    conservative corrections are also included by Benjamini & Hochberg
    (1995) ("BH" or its alias "fdr"). A pass-through option ("none") is also 
    included. The set of methods are contained in the "p_adjust_methods"
    list for the benefit of methods that need to have the method as an option
    and pass it on to "p_adjust".

    The Bonferroni method is designed to give strong control of the
    family-wise error rate.

    The "BH" (aka "fdr") method of Benjamini & Hochberg controls the false 
    discovery rate, the expected proportion of false discoveries amongst
    the rejected hypotheses. The false discovery rate is a less stringent
    condition than the family-wise error rate, so these methods are more
    powerful than the others.

    Note that you can set 'n' larger than 'len(n)' which means the
    unobserved p-values are assumed to be greater than all the observed
    p for 'bonferroni' methods and equal to 1 for the other methods.

Value
    A numpy array of corrected p-values (of the same length as 'p').

References
    Benjamini, Y., and Hochberg, Y. (1995).  Controlling the false
        discovery rate: a practical and powerful approach to multiple
        testing.  Journal of the Royal Statistical Society Series B
        57, 289-300.

Examples
    from p_adjust import p_adjust
    p_adjust([0.03,0.2,0.4], method="fdr")
'''

import numpy as np
from scipy import interpolate

p_adjust_methods=["fdr","bonferroni","none"]

# This function is partly based on Nicolo Fusi's qvalue code at
# https://github.com/nfusi/qvalue/blob/master/qvalue/qvalue.py
def p_adjust(p, method="fdr", n=None):
    '''
    p       - numpy array of p-values
    method  - correction method. default: fdr
    n       - number of comparisons, must be at least len(p); only set
              this (to non-default) when you know what you are doing!
              default: len(p)
    '''
    method=str(method).lower()
    if method=="bh":
        method="fdr"
    assert(method in p_adjust_methods
        ), "ERROR! method must be one of "+','.join(p_adjust_methods)

    if isinstance(p,list):
        p=np.array(p)

    if method=="none":
        return p

    original_shape = p.shape
    p = p.ravel()  # flattens the array in place

    if n>0:
        n *= 1.
    else:
        n = 1.*len(p)

    if method=="bonferroni":
        qv = n * p
        qv[qv>1]=1.
        qv = qv.reshape(original_shape)
        return qv

    # if the number of hypotheses is small, just set pi0 to 1
    if len(p) < 100:
        pi0 = 1.0
    else:
        # evaluate pi0 for different lambdas
        pi0 = []
        lam = np.arange(0, 0.90, 0.01)
        counts = np.array([(p > i).sum() for i in sp.arange(0, 0.9, 0.01)])
        for l in range(len(lam)):
            pi0.append(counts[l]/(m*(1-lam[l])))

        pi0 = np.array(pi0)

        # fit natural cubic spline
        tck = interpolate.splrep(lam, pi0, k=3)
        pi0 = interpolate.splev(lam[-1], tck)

        if pi0 > 1:
            pi0 = 1.0
    if pi0 < 0:
        print("WARNING! pi0 is not between 0 and 1: %f" % pi0)

    p_ordered = np.argsort(p)
    p = p[p_ordered]
    qv = pi0 * n/len(p) * p
    qv[-1] = min(qv[-1], 1.0)

    for i in xrange(len(p)-2, -1, -1):
        qv[i] = min(pi0*n*p[i]/(i+1.0), qv[i+1])

    # reorder qvalues
    qv_temp = qv.copy()
    qv = np.zeros(qv.shape)
    qv[p_ordered] = qv_temp

    # reshape qvalues
    qv = qv.reshape(original_shape)
    return qv

if __name__=="__main__":
    print(docstring)
