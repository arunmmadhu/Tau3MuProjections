#Run with python3

import math
import argparse
from scipy.stats import poisson

poisson_ppf_val = 0.5

def integers_less_than(input_number):
    
    if not isinstance(input_number, (int, float)):
        raise ValueError("Input must be a floating-point number.")
    
    # Using a list comprehension to generate the list of integers
    integer_list = [i for i in range(int(input_number)+1)]

    return integer_list
    
def calculate_poisson_sum(numbers, lambda_val):
    
    if lambda_val==0.0:
        lambda_val=0.001
    
    # Check if input is valid
    if all([x < 0 for x in numbers]) or lambda_val < 0.0:
        raise ValueError("Both k and lambda_val must be non-negative.")
    
    # Using map to calculate the Poisson function for each pair of number and lambda
    #poisson_values = map(lambda x: (math.exp(-lambda_val) * (lambda_val ** x)) / math.factorial(x), numbers)
    #poisson_values = map(lambda x: math.exp(x * math.log(lambda_val) - lambda_val - lgamma(x + 1)), numbers)
    poisson_values = poisson.pmf(numbers, lambda_val)
    
    # Using sum to calculate the sum of the Poisson values
    total_sum = sum(poisson_values)
    
    return total_sum
    
def calculate_cls(s_exp, b_exp):
    # Check if input is valid
    CL_s_p_b = calculate_poisson_sum(integers_less_than(poisson.ppf(poisson_ppf_val, b_exp)), s_exp+b_exp)
    CL_b = calculate_poisson_sum(integers_less_than(poisson.ppf(poisson_ppf_val, b_exp)), b_exp)
    CL_s = CL_s_p_b / CL_b
    return CL_s
    
    
def find_domain_for_value(function,s_exp,b_exp):
    range_for_s = 20.0
    start = max(s_exp-range_for_s,0.0)
    end=s_exp+range_for_s
    step=0.01
    
    print(f"Starting value for the s_exp_calc scan is {start}, end is {end} with a step size of {step}.")
    
    """
    Find the domain of a function for which it returns a specific target value.

    Parameters:
    - target_value: The specific value the function should return
    - function: The function to analyze
    - start: The starting point for the search range
    - end: The ending point for the search range
    - step: The step size for the search range

    Returns:
    - List of input values (domain) for which the function returns the target value
    """
    target_value = 0.1 # 90% CLs confidence
    
    domain_value = 20.0
    targ_diff = 100.0
    actual_target = 200.0
    current_value = start
    
    

    while current_value <= end:
        if abs(function(current_value,b_exp)-target_value) < targ_diff:
            domain_value = current_value
            targ_diff = abs(function(current_value,b_exp)-target_value)
            actual_target = function(current_value,b_exp)
        current_value += step

    return [round(domain_value, 2),actual_target,round(domain_value/s_exp, 2)]
    
    
def main():
    parser = argparse.ArgumentParser(description="Input s_exp and b_exp.")
    parser.add_argument("arg1", type=float, help="s_exp")
    parser.add_argument("arg2", type=float, help="b_exp")
    args = parser.parse_args()

    var1, var2, var3 = find_domain_for_value(calculate_cls, args.arg1, args.arg2)
    
    print(f"The best value of s_exp_calc is {var1} events")
    print(f"The best value of CLs is {var2}")
    print(f"The best value of signal strength / upper limit is {var3}")
    
main()    

