from math import factorial
from numpy import log

def polylog_integrand(s, z, t): 
    sgn = -1 if s % 2 == 0 else 1 
    return sgn * z * pow(log(t), s-1) / (factorial(s-1) * (1.0 - t * z))

def polylog_wrapper(func, parameters, x):
    """Wrapper function that evaluates the function and returns updated parameters with results."""
    result = func(parameters["s"], parameters["z"], x)
    updated_params = parameters.copy()  # Copy to avoid modifying original dictionary
    updated_params["x"] = x
    updated_params["func"] = result
    return updated_params