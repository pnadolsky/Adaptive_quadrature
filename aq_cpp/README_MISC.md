# Polylog Test Module

## Overview
The `PolyTest` module implements the Borwein, Borwein & Girgensohn integral representation of the polylogarithm function $ Li_s(z) $. It provides an integrand function, `polylog_integrand`, and a wrapper mechanism for evaluating the function based on the parameters $ s $ and $ z $.

## Files
- **polylog_port.hpp**: Header file defining the function interfaces.
- **polylog_port.cpp**: Implementation of the polylogarithm integrand and wrapper functions.
- **polylog_test.cpp**: Test file to validate the functionality of `polylog_wrapper`.

## Functionality
- **`polylog_integrand`**:
  ```cpp
  double polylog_integrand(int s, double z, double t);
  ```
  - Computes the integrand in the integral representation of the polylogarithm function using the Borwein, Borwein & Girgensohn representation:
$$
f(s,z; t) = \frac{z}{s!} \cdot \frac{\ln^s\left[ 1 \over t \right]}{1-zt} 
$$

- **Wrapper Function: `polylog_wrapper`**:
  ```cpp
  double polylog_wrapper(
      std::unordered_map<std::string, std::variant<int, double, std::string>> parameters, 
      double t);
  ```
  - Extracts `s` and `z` from an unordered map of parameters.
  - Calls `polylog_integrand` to compute the result.
  - Returns the evaluated function value.

## Testing
The `polylog_test.cpp` file provides test cases that:
1. Define a parameter map containing:
   ```cpp
   std::unordered_map<std::string, std::variant<int, double, std::string>> params;
   params["s"] = 3;
   params["z"] = 0.5;
   ```
2. Evaluate `polylog_wrapper` for a sample \( t \) value.
3. Output the computed result.

### Running the Test
To compile and run the test:
```sh
 g++ -o polylog_test polylog_port.cpp polylog_test.cpp -std=c++17 -Wall
 ./polylog_test
```

