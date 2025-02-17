# C++ Adaptive Quadrature Classes
##  Adaptive Quadrature

### Overview
The `AdaptiveGaussTree` module implements an adaptive quadrature method to numerically integrate functions over a given interval. The adaptive approach refines the integration nodes dynamically, using either Gauss-Legendre or Gauss-Laguerre quadrature based on singularity detection.

### Files
- **adaptive_gauss_tree.hpp**: Header file defining the `AdaptiveGaussTree` class and related functions.
- **adaptive_test.cpp**: Test file demonstrating the usage of the `AdaptiveGaussTree` class.

### Features
- Supports adaptive refinement for accurate numerical integration.
- Utilizes **Gauss-Legendre** and **Gauss-Laguerre** quadrature methods.
- Can construct the tree using direct function parameters or load from a JSON file.
- Saves computed integration results to a JSON file.

### Class Structure
#### **AdaptiveGaussTree**
This class manages the adaptive quadrature process and stores results in a hierarchical tree structure.

##### **Constructors**
1. **From Parameters and an Function to Integrate:**
```cpp
AdaptiveGaussTree(
    std::function<ParamMap, double)> f,
    double lower, double upper, double tol, int minD, int maxD,
    int n1, int n2,
    double alphaA, double alphaB, bool singularA, bool singularB,
    WeightsLoader rl1, WeightsLoader rl2, WeightsLoader ll1, WeightsLoader ll2,
    std::string name="Project", std::string author="Author",  std::string description="project description", 
    std::string reference="references", std::string version="1.0", update_log_message="Initial Train" 
    );
```
- See README_QUADRATURE for information about ParamMap
- Initializes the quadrature tree based on user-defined parameters.
- Uses `WeightsLoader` instances to provide quadrature weights.
- Now includes optional json header fields for project name, author, references, version, and project description

2. **From a JSON File and a Function to Integrate:**
```cpp
AdaptiveGaussTree(
    std::function<double(ParamMap, double)> f,
    WeightsLoader rl1, WeightsLoader rl2, WeightsLoader ll1, WeightsLoader ll2, std::string filename);
```
- Loads a previously saved quadrature tree from a JSON file.

##### **Methods**
- `std::pair<double, double> get_integral_and_error() const;`
  - Returns the total integral and error by traversing the quadrature tree.
- `void save_to_json(std::string filename);`
  - Saves the tree structure, computed integrals, and metadata to a JSON file.
- `void load_from_json(std::string filename);`
  - Loads a quadrature tree from a JSON file.
- `void add_update_log(const std::string& message)`
  - add an entry to the update log
- `void print_update_log()`
  - prints the update log

##### **Tree Node Structure**
Each node in the tree represents an interval of the integration domain and contains:
- `lower, upper`: Integration bounds for this node.
- `depth`: Depth of the node in the tree.
- `tolerance`: Error tolerance at this node.
- `error, result`: Computed integration error and result.
- `method`: Either **Gauss-Legendre** or **Gauss-Laguerre**.
- `left, right`: Pointers to child nodes for further refinement.

### Usage Example
#### **Using the Adaptive Gauss Tree**
```cpp
#include "adaptive_gauss_tree.hpp"
#include <iostream>

// Example function: log(x) / sqrt(x)
double test_function(ParamMap params, double x) {
    return x > 0.0 ? std::log(x) / std::sqrt(x) : 0.0;
}

int main() {
    WeightsLoader legendre_n1("../model_json/legendre_n1.json");
    WeightsLoader legendre_n2("../model_json/legendre_n2.json");
    WeightsLoader laguerre_n1("../model_json/laguerre_n1.json");
    WeightsLoader laguerre_n2("../model_json/laguerre_n2.json");

    AdaptiveGaussTree adaptive_tree(test_function, 0.0, 1.0, 1e-6, 2, 10,
                                    50, 100, 0.5, 0.5, true, false,
                                    legendre_n1, legendre_n2, laguerre_n1, laguerre_n2);
    
    auto [integral, error] = adaptive_tree.get_integral_and_error();
    std::cout << "Computed Integral: " << integral << "\n";
    std::cout << "Estimated Error: " << error << "\n";

    adaptive_tree.save_to_json("adaptive_output.json");
    return 0;
}
```
Note that it should be ok to pass the same WeightsLoader for the same quadrature types (e.g. laguerre_n1, laguerre_n2).
### Running the Test
```sh
 g++ -o adaptive_test adaptive_test.cpp  weights_loader.cpp legendre_quadrature.cpp laguerre_singular_endpoint.cpp laguerre_quadrature.cpp quadrature.cpp  -std=c++17 
 ./adaptive_test
```
##  Batch Integration

WIP