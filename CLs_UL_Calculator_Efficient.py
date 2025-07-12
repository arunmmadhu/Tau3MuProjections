# Run with: python3 CLs_UL_Calculator_Efficient.py <s_exp> <b_exp>

import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy.core.getlimits")

import argparse
import numpy as np
from scipy.stats import poisson
from scipy.optimize import brentq

# Global CLs target and quantile
poisson_ppf_val = 0.5
target_cls = 0.1  # 90% CLs

def calculate_poisson_sum(numbers, lambda_val):
    """
    Calculate the sum of Poisson PMFs for given integer points and mean lambda.
    """
    lambda_val = max(lambda_val, 0.001)  # prevent zero
    return np.sum(poisson.pmf(numbers, lambda_val))

def calculate_cls(s_exp, b_exp):
    """
    Calculate CLs = CL(s+b) / CL(b)
    """
    ppf_val = poisson.ppf(poisson_ppf_val, b_exp)
    ks = np.arange(0, int(ppf_val) + 1)  # discrete range of interest
    CL_sb = calculate_poisson_sum(ks, s_exp + b_exp)
    CL_b = calculate_poisson_sum(ks, b_exp)
    return CL_sb / CL_b

def find_best_s_exp(function, s_exp, b_exp):
    
    def cls_diff(s_calc):
        return function(s_calc, b_exp) - target_cls

    start = 0.0
    end = s_exp + 200.0

    try:
        s_exp_calc = brentq(cls_diff, start, end, xtol=1e-3)
        best_cls = function(s_exp_calc, b_exp)
        signal_strength = s_exp_calc / s_exp
        return [round(s_exp_calc, 2), round(best_cls, 5), round(signal_strength, 3)]
    except ValueError:
        print("No solution found in the specified range.")
        return [None, None, None]

def main():
    parser = argparse.ArgumentParser(description="Input s_exp and b_exp to compute CLs-based upper limit.")
    parser.add_argument("arg1", type=float, help="s_exp (expected signal)")
    parser.add_argument("arg2", type=float, help="b_exp (expected background)")
    args = parser.parse_args()

    s_exp = args.arg1
    b_exp = args.arg2

    print(f"Starting value scan from 0 to {s_exp + 200.0}...")

    s_calc, cls_val, signal_strength = find_best_s_exp(calculate_cls, s_exp, b_exp)

    if s_calc is not None:
        print(f"The best value of s_exp_calc is {s_calc} events")
        print(f"The best value of CLs is {cls_val}")
        print(f"The best value of signal strength / upper limit is {signal_strength}")
    else:
        print("Could not determine an upper limit within the search range.")

if __name__ == "__main__":
    main()
