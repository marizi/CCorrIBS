#!/usr/bin/env python
# coding=utf-8
# ==============================================================================
# title           : circ_corrcc_uniform.py
# description     : Adjusted circular correlation coefficient
# author          : Guillaume Dumas, 2023 (guillaume.dumas@ppsp.team)
# notes           : Adapted for the Matlab function "circ_corrcc"
#                   original by Philipp Berens, 2009 (berens@tuebingen.mpg.de)
#                   edited by Marius Zimmermann, 2020 (marz@dtu.dk)
# ==============================================================================

import numpy as np
from scipy.stats import norm

def circ_mean(alpha):
    return np.angle(np.sum(np.exp(1j * alpha)))

def circ_corrcc_uniform(alpha1, alpha2):
    x_sin = np.sin(alpha1 - circ_mean(alpha1))
    y_sin = np.sin(alpha2 - circ_mean(alpha2))
    r_minus = np.abs(np.sum(np.exp((alpha1-alpha2) * 1j)))
    r_plus  = np.abs(np.sum(np.exp((alpha1+alpha2) * 1j)))
    num = (r_minus - r_plus)
    den = 2*np.sqrt(np.sum(x_sin**2) * np.sum(y_sin**2))
    rho = num / den
    # Compute p-value
    n = len(alpha1)
    alpha1_bar = circ_mean(alpha1)
    alpha2_bar = circ_mean(alpha2)
    l20 = np.mean(np.sin(alpha1 - alpha1_bar)**2);
    l02 = np.mean(np.sin(alpha2 - alpha2_bar)**2);
    l22 = np.mean((np.sin(alpha1 - alpha1_bar)**2) * (np.sin(alpha2 - alpha2_bar)**2));
    ts = np.sqrt((n * l20 * l02)/l22) * rho
    pval = 2 * (1 - norm.cdf(abs(ts)))
    return rho, pval