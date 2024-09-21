#!/usr/bin/env python
# coding=utf-8
# ==============================================================================
# title           : circ_corrcc_uniform.py
# description     : Adjusted circular correlation coefficient
# author          : Guillaume Dumas, 2023-2024 (guillaume.dumas@ppsp.team)
# notes           : Adapted for the Matlab function "circ_corrcc"
#                   original by Philipp Berens, 2009 (berens@tuebingen.mpg.de)
#                   edited by Marius Zimmermann, 2020-2024 (marz@dtu.dk)
# ==============================================================================

import numpy as np
from scipy.stats import norm, binom


def circ_mean(alpha):
    return np.angle(np.sum(np.exp(1j * alpha)))


def hodges_ajne_test(alpha, alpha_threshold=0.05):
    """Performs Hodges-Ajne test for circular uniformity.
    H0: population is uniformly distributed around the circle.
    Returns the p-value.
    """
    n = len(alpha)
    alpha = np.mod(alpha, 2 * np.pi)
    alpha = np.sort(alpha)
    U = np.array([np.sum((alpha_sorted >= alpha_sorted[i]) & 
                         (alpha_sorted < (alpha_sorted[i] + np.pi))) for i in range(n)])
    m = np.min(U)
    p_val = 2 * binom.cdf(m, n, 0.5)
    if p_val < alpha_threshold:
        print(f"Input signal is not uniformly distributed! p-value = {p_val:.3f}")
    return p_val


def circ_corrcc_uniform(alpha1, alpha2, test_uniform=True, test_uniform_alpha=0.05):
    # Ensure inputs are arrays
    alpha1 = np.asarray(alpha1)
    alpha2 = np.asarray(alpha2)
    
    # Perform Hodges-Ajne test for non-uniformity if requested
    if test_uniform:
        pHodgesAjne_alpha1 = hodges_ajne_test(alpha1, test_uniform_alpha)
        pHodgesAjne_alpha2 = hodges_ajne_test(alpha2, test_uniform_alpha)
    
    # Adjusted circular mean for alpha1 and alpha2
    n = -1 * (circ_mean(np.angle(alpha1 - alpha2)) - circ_mean(np.angle(alpha1 + alpha2))) / 2
    m = circ_mean(np.angle(alpha1 - alpha2)) + n

    x_sin = np.sin(alpha1 - m)
    y_sin = np.sin(alpha2 - n)

    # Calculation of r_minus and r_plus
    r_minus = np.abs(np.sum(np.exp(1j * (alpha1 - alpha2))))
    r_plus = np.abs(np.sum(np.exp(1j * (alpha1 + alpha2))))

    # Numerator and corrected denominator
    num = r_minus - r_plus
    den = 2 * np.sqrt(np.sum(x_sin**2) * np.sum(y_sin**2))
    
    rho = num / den

    # Compute p-value
    n = len(alpha1)
    alpha1_bar = circ_mean(alpha1)
    alpha2_bar = circ_mean(alpha2)
    
    l20 = np.mean(np.sin(alpha1 - alpha1_bar)**2)
    l02 = np.mean(np.sin(alpha2 - alpha2_bar)**2)
    l22 = np.mean(np.sin(alpha1 - alpha1_bar)**2 * np.sin(alpha2 - alpha2_bar)**2)
    
    ts = np.sqrt((n * l20 * l02)/l22) * rho
    pval = 2 * (1 - norm.cdf(abs(ts)))
    
    return rho, pval
