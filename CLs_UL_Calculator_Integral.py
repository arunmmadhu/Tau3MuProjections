# Run with: python3 CLs_UL_Calculator_Integral.py <s_exp> <b_exp>

import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy.core.getlimits")

import numpy as np
from scipy.stats import poisson
from scipy.integrate import quad
from scipy.optimize import brentq
import argparse
import math
from scipy.special import gammaln

def poisson_likelihood(n_obs, mu):
    """
    Generalized Poisson PMF that supports non-integer n_obs.
    """
    return (mu ** n_obs) * np.exp(-mu) / np.exp(gammaln(n_obs + 1))

def posterior_pdf(r, n_obs, b, s0):
    mu = b + r * s0
    return poisson_likelihood(n_obs, mu)

def normalization_constant(n_obs, b, s0):
    integral, _ = quad(posterior_pdf, 0, np.inf, args=(n_obs, b, s0))
    return integral

def normalized_posterior(r, n_obs, b, s0, norm):
    return posterior_pdf(r, n_obs, b, s0) / norm

def tail_integral(r_limit, n_obs, b, s0, norm):
    integral, _ = quad(normalized_posterior, r_limit, np.inf, args=(n_obs, b, s0, norm))
    return integral

def find_upper_limit(b, s0, cl=0.9):
    n_obs = b  # not rounded!
    norm = normalization_constant(n_obs, b, s0)

    def tail_diff(r):
        return tail_integral(r, n_obs, b, s0, norm) - (1 - cl)

    #Ensures that the sign of f(r) flips from positive at 0 to -ve at +ve at r_upper_guess. This loop gives the search range.
    r_upper_guess = 20
    while tail_diff(r_upper_guess) > 0:
        r_upper_guess *= 2

    r_limit = brentq(tail_diff, 0, r_upper_guess, xtol=1e-4)
    s_limit = r_limit * s0
    return r_limit, s_limit, n_obs

def main():
    parser = argparse.ArgumentParser(description="Bayesian 90% CL upper limit with n_obs = b (no rounding)")
    parser.add_argument("s0", type=float, help="Reference signal expectation")
    parser.add_argument("b", type=float, help="Expected background")
    args = parser.parse_args()

    r_lim, s_lim, n_obs = find_upper_limit(args.b, args.s0)

    print(f"\nBayesian Upper Limit at 90% CL (alpha = 0.1):")
    print(f"  n_obs (set equal to b) = {n_obs}")
    print(f"  r < {r_lim:.4f}")
    print(f"  s < {s_lim:.4f} (s0 = {args.s0})\n")

if __name__ == "__main__":
    main()
    