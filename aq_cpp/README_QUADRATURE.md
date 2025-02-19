# C++ Quadrature Library
##  Weights Loader

### Overview
The `weights_loader` module is responsible for loading and managing weight matrices for the neural network implementation. It provides functionality to read weights from external files and convert them into a usable format within the C++ program.

### Files
- **weights_loader.hpp**: Header file declaring the `WeightsLoader` class and its methods.
- **weights_loader.cpp**: Implementation file containing the logic for loading and handling weights.
- **weights_loader_test.cpp**: Unit test file that verifies the functionality of the `WeightsLoader` class.

### Features
- Load weights from JSON-formatted files
- Convert data into a format compatible with Eigen or other matrix operations
- Handle file reading errors gracefully

### Dependencies
- C++17 or later
- [Eigen](https://eigen.tuxfamily.org/) for matrix operations
- `json.hpp` (single-header JSON parser) for handling JSON files

### Compile weights_loader_test
~~~
g++ -std=c++17 -o weights_loader_test weights_loader.cpp weights_loader_test.cpp 
~~~
## Quadrature

### Overview
The `Quadrature` module defines an abstract base class for numerical integration methods. It provides a standard interface for different quadrature techniques by enforcing method signatures that must be implemented in derived classes.

### Files
- **quadrature.hpp**: Header file defining the `Quadrature` abstract class, its methods, and protected data members.
    - the header includes these definitions that are used extensively in the code in this library (almost everything includes the header for this class):  
    ```cpp
    using ParamType = std::variant<int, double, std::string>;
    using ParamMap = std::unordered_map<std::string, ParamType>;
    using ParamCollection = std::unordered_map<std::string, std::variant<std::vector<int>, std::vector<double>, std::vector<std::string>>>; 
    ```
    - in addition to these, 
        - ```operator<<``` is overloaded to print ParamType and ParamMaps to ostreams.
        - Hashing and equality structs are provided to assist the functionality of building maps with AdaptiveQuadrature trees in them (see adaptive_gauss_tree.hpp for more details).
- **quadrature.cpp**: Implementation file containing shared functionality for derived quadrature classes.

### Virtual Methods to Overload
Derived classes must implement the following virtual methods:

#### `integrate`
```cpp
virtual double integrate(
    std::function<double(ParamMap, double)> func,
    ParamMap parameters) = 0;
```
- This pure virtual function must be overridden to provide a numerical integration method.
- Takes a function `func` as input, along with a set of parameters.
- Returns the computed integral as a `double`.

#### `transformVariable`
```cpp
virtual double transformVariable(double t) const;
```
- Provides optional variable transformation logic.
- Can be overridden by derived classes for specialized transformations.

### Protected Data Members
These members are available to derived classes:
- `const WeightsLoader& weightsLoader`: Reference to an external weight loader.
- `int order1, order2`: Orders of quadrature methods.
- `double error, result`: Stores integration error and result.
- `std::optional<double> lowerLimit, upperLimit`: Optional limits of integration.
- `std::vector<double> nodes1, weights1`: Nodes and weights for the first quadrature scheme.
- `std::vector<double> nodes2, weights2`: Nodes and weights for the second quadrature scheme.
- `const std::string method`: A compile-time enforced method name.

### Public Getter Methods
The class provides the following getters:

```cpp
std::string get_method() const;
int getOrder1() const;
int getOrder2() const;
double getResult() const;
double getError() const;
std::optional<double> getLowerLimit() const;
std::optional<double> getUpperLimit() const;
```

#### Functionality:
- `get_method()`: Returns the name of the quadrature method.
- `getOrder1()`, `getOrder2()`: Retrieve the quadrature orders.
- `getResult()`: Returns the computed integral result.
- `getError()`: Returns the computed integration error.
- `getLowerLimit()`, `getUpperLimit()`: Retrieve optional integration limits.

### Extending the Quadrature Class
To create a new quadrature method, inherit from `Quadrature` and implement the `integrate` method:

```cpp
class MyQuadratureMethod : public Quadrature {
public:
    double integrate( std::function<double(ParamMap, double)> func, ParamMap parameters) override {
        // Implement integration logic here
    }
};
```


## Legendre Quadrature

### Overview
The `LegendreQuadrature` module implements Gauss-Legendre quadrature, a numerical integration technique that approximates the integral of a function using Legendre polynomials as weight functions. This module inherits from the `Quadrature` base class and provides concrete implementations of the required integration methods.

### Files
- **legendre_quadrature.hpp**: Header file defining the `LegendreQuadrature` class.
- **legendre_quadrature.cpp**: Implementation file containing the logic for integration.
- **legendre_quadrature_test.cpp**: Test file to validate the integration method.

### Inherited Methods from `Quadrature`
The `LegendreQuadrature` class overrides the following methods from `Quadrature`:

#### `integrate`
```cpp
virtual double integrate( std::function<double(ParamMap, double)> func, ParamMap parameters) override;
```
- Implements numerical integration using Gauss-Legendre quadrature.
- Evaluates the function at quadrature nodes and applies weight factors.
- Uses two different quadrature orders for error estimation.

#### `transformVariable`
```cpp
virtual double transformVariable(double t) const override;
```
- Transforms a variable from the standard Legendre domain `[-1,1]` to the given integration interval `[lower, upper]`.
- Used to map quadrature nodes to the actual integration range.

### Class Structure
#### Constructor
```cpp
LegendreQuadrature(const WeightsLoader& loader, int n1, int n2, double lower, double upper);
```
- Calls the `Quadrature` constructor with the method name "Gauss-Legendre".
- Initializes the integration limits and orders.

#### Testing
The `legendre_quadrature_test.cpp` file provides a test case that:
1. Loads quadrature weights from a JSON file (`legendre.json`).
2. Defines function parameters and selects orders for integration.
3. Calls the `integrate` method and prints the computed integral result and error estimate.


###  Compile legendre_quadrature_test
~~~
g++ -std=c++17 -o legendre_quadrature_test weights_loader.cpp polylog_port.cpp quadrature.cpp legendre_quadrature.cpp legendre_quadrature_test.cpp
~~~

## Laguerre Quadrature

### Overview
The `LaguerreQuadrature` module implements Gauss-Laguerre quadrature, a numerical integration technique for approximating integrals over the semi-infinite domain `[0, ∞)`. This method is particularly useful for functions with an exponential decay. The module inherits from the `Quadrature` base class and provides concrete implementations of the required integration methods.

### Files
- **laguerre_quadrature.hpp**: Header file defining the `LaguerreQuadrature` class.
- **laguerre_quadrature.cpp**: Implementation file containing the logic for integration.
- **laguerre_quadrature_test.cpp**: Test file to validate the integration method.

### Inherited Methods from `Quadrature`
The `LaguerreQuadrature` class overrides the following methods from `Quadrature`:

#### `integrate`
```cpp
virtual double integrate( std::function<double(ParamMap, double)> func, ParamMap parameters) override;
```
- Implements numerical integration using Gauss-Laguerre quadrature.
- Evaluates the function at quadrature nodes and applies weight factors.
- Uses two different quadrature orders for error estimation.

#### `transformVariable`
```cpp
virtual double transformVariable(double t) const override;
```
- Maps quadrature nodes directly to the domain `[0, ∞)`, as Laguerre quadrature does not require variable transformation.

### Additional Functionality
The `LaguerreQuadrature` class introduces a weight function option:

#### `laguerre_weight_function`
```cpp
virtual double laguerre_weight_function(double t) const;
```
- Returns the Laguerre weight function `exp(t)`, which can be optionally applied during integration.

### Class Structure
#### Constructor
```cpp
LaguerreQuadrature(const WeightsLoader& loader, int n1, int n2, bool use_weight_function = true);
```
- Calls the `Quadrature` constructor with the method name "Gauss-Laguerre".
- Supports an optional weight function flag.
    - Recall that Laguerre quadrature assumes that you are integrating:  $\int_0^\infty f(t) e^{-t} dt$.
    - use_weight_function = True tells the routine that you are integrating $\int_0^\infty f(t)  dt$, and the extra $e^t$ must be accounted for in the code.

### Testing
The `laguerre_quadrature_test.cpp` file provides test cases that:
1. Load quadrature weights from a JSON file (`laguerre.json`).
2. Define test functions for integration.
3. Call the `integrate` method with and without the weight function.
4. Print the computed integral result and error estimate.


###  Compile laguerre_quadrature_test
~~~
g++ -std=c++17 -o laguerre_quadrature_test weights_loader.cpp quadrature.cpp laguerre_quadrature.cpp laguerre_quadrature_test.cpp
~~~

## Laguerre Singular Quadrature

### Overview
The `LaguerreSingularEndpoint` module extends `LaguerreQuadrature` to handle functions with integrable singularities, such as functions of the form:

$$
 \frac{\log^n(x)}{x^\alpha}, \quad 0 < \alpha < 1
$$

This class introduces a transformation to better approximate integrals with singularities at either endpoint of the domain.  The value of $\alpha$ (or at least and idea of its value ) must be known beforehand, or the integral won't converge. 

### Files
- **laguerre_singular_endpoint.hpp**: Header file defining the `LaguerreSingularEndpoint` class.
- **laguerre_singular_endpoint.cpp**: Implementation file containing the logic for singularity handling.
- **laguerre_singular_test.cpp**: Test file to validate integration of singular functions.

### Inherited Methods from `LaguerreQuadrature`
The `LaguerreSingularEndpoint` class overrides the following methods from `LaguerreQuadrature`:

#### `laguerre_weight_function`
```cpp
virtual double laguerre_weight_function(double t) const override;
```
- Modifies the weight function to account for the singularity.
- Applies a transformation based on the expected exponent `alpha`.
- This replacement is to account for the power law singular behavior when alpha is not zero:
$$
w(t) = \frac{\left( b-a \right)^{1-\alpha}}{1-\alpha} \left( t-\alpha \right)^{\alpha}.
$$

#### `transformVariable`
```cpp
virtual double transformVariable(double t) const override;
```
- Implements a transformation to map nodes to the interval `[lower, upper]`.
- Uses an exponential transformation to adapt integration points to singular behavior.
- The form of the transformation depends on the value of $\alpha$:
$$
x(t) = a \cdot (1-z) + b \cdot z, 
$$
where
$$
z(t) = \exp\left(-{1 \over 1- \alpha}t\right). 
$$

### Class Structure
#### Constructor
```cpp
LaguerreSingularEndpoint(const WeightsLoader& loader, int n1, int n2, double lower, double upper, bool leftIsSingular = true, double alpha = 0);
```
- Calls the `LaguerreQuadrature` constructor with the method name "Laguerre-Singular".
- Sets whether the singularity is at the left or right endpoint.
- Takes `alpha`, the exponent governing the singularity's strength.

### Testing
The `laguerre_singular_test.cpp` file provides test cases that:
1. Load quadrature weights from `laguerre.json`.
2. Define a singular test function such as:
   ```cpp
   double singular_test_function(ParamMap params, double x) {
       return std::log(x) / std::sqrt(x);
   }
   ```
3. Compare results between:
   - Standard `LegendreQuadrature`
   - `LaguerreSingularEndpoint` with a singular left endpoint
   - `LaguerreSingularEndpoint` with a singular right endpoint
4. Output computed integral results and error estimates.



###  Compile laguerre_singular_test
~~~
g++ -std=c++17 -o laguerre_singular_test weights_loader.cpp quadrature.cpp legendre_quadrature.cpp laguerre_quadrature.cpp laguerre_singular_endpoint.cpp laguerre_singular_test.cpp
~~~




##  Polylog Port
### Compile polylog_test
~~~
g++ -std=c++17 -o polylog_test polylog_port.cpp polylog_test.cpp
~~~
