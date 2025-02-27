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
- Defines **QuadCollection**:  
```cpp 
using QuadCollection = std::unordered_map<ParamMap, std::unique_ptr<AdaptiveGaussTree>, ParamMapHash, ParamMapEqual>;
```
- Overloads **operator<<** to print integral and error to ostreams
- **void printQuadCollectionKeys**:  print a list of batch keys.  This will be moved to adaptive_gauss_batch when it is functional.
- **void printQuadCollectionResults**: prints a list of keys and integral results.   This will be moved to adaptive_gauss_batch when it is functional.

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
    ParamMap args={},
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
    WeightsLoader rl1, WeightsLoader rl2, WeightsLoader ll1, WeightsLoader ll2, ParamMap args={},std::string filename);
```
- Loads a previously saved quadrature tree from a JSON file.

##### **Methods**
- `std::pair<double, double> get_integral_and_error() const;`
  - Returns the total integral and error by traversing the quadrature tree.
- `void save_to_json(std::string filename, overwrite = False);`
  - Saves the tree structure, computed integrals, and metadata to a JSON file  (set to True to overwrite file).
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
#include <adaptive_gauss_tree.hpp>
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
g++ -o adaptive_test -Iinclude source/*.cpp test/adaptive_test.cpp  -std=c++17 
 ./adaptive_test
```

# AdaptiveGaussTreeBatch

## Overview
`AdaptiveGaussTreeBatch` is a C++ class designed to manage collections of `AdaptiveGaussTree` objects, which perform adaptive Gaussian quadrature calculations. This class supports parameterized tree reconstruction, merging, and serialization from JSON files.

## Features
- **Batch processing of adaptive Gaussian quadrature trees**
- **Loading and reconstructing trees from serialized JSON files**
- **Merging multiple batches with automatic parameter resolution**
- **Serialization and deserialization of quadrature trees and parameter sets**
- **Sorting and filtering of parameter combinations**

## Dependencies
This class depends on the following external libraries:
- [`nlohmann/json`](https://github.com/nlohmann/json) for JSON serialization
- Custom classes:
  - `WeightsLoader`
  - `AdaptiveGaussTree`

## Usage

### Initialization
You can construct an instance using either direct parameters or by loading from a JSON file.

#### Construct from Parameters
```cpp
AdaptiveGaussTreeBatch batch(
    function,                 // Quadrature function
    0.0, 1.0,                 // Integration bounds
    1e-6,                     // Tolerance
    2, 10,                    // Min and max depth
    5, 10,                    // Orders for Gauss quadrature
    0.0, 0.0,                 // Alpha parameters
    false, false,              // Singular endpoints
    legendre_n1, legendre_n2,  // Weight loaders for Legendre
    laguerre_n1, laguerre_n2,  // Weight loaders for Laguerre
    parameters,                // Parameter collection
    "ExampleProject", "Author", "Description", "Reference", "1.0"
);
```

#### Construct from JSON
```cpp
AdaptiveGaussTreeBatch batch(func, "trees.json");
```

### Merging Two Batches
```cpp
AdaptiveGaussTreeBatch batch1(func, "batch1.json");
AdaptiveGaussTreeBatch batch2(func, "batch2.json");
batch1 += batch2; // Merges batch2 into batch1 and returns the tree (for chaining multiple)
AdaptiveGaussTreeBatch batch3 = batch1+batch2; // Constructs batch3 from deep copies of batch1 and batch2
batch1.merge(batch2) // adds batch2 to the batch1 object (e.g. same as merge, but with no chaining)
```

### Saving to JSON
```cpp
batch.save_to_json("output.json");
```

## Key Methods
- **`void merge(const AdaptiveGaussTreeBatch& other)`**
  - Merges another batch, ensuring unique parameter sets.
- **`void save_to_json(const std::string & filename, bool overwrite = false, bool write_roots=false, bool write_trees =true)`**
  - Saves the batch data into a JSON file.
- **`json parameter_serializer(bool dump_nodes = false)`**
  - Serializes parameter sets and trees.
- **`void printCollection()`**
  - Prints the parameter set with values for the integral and error. 
- **`add_update_log(update_log_message)`**
  - Adds a log message and timestamp to the class object.
- **`const QuadCollection& getCollection() const`**
  - Returns the collection of trees. 


## Data Structure
- `QuadCollection`: Stores mappings of parameter sets to `AdaptiveGaussTree` instances.
- `ParamCollection`: A map storing parameter values of type `int`, `double`, or `string`.
- `results`: Stores all valid parameter combinations.
- `update_log`: Tracks modifications to the batch.

## Error Handling
- Throws `std::runtime_error` if a JSON file fails to load or contains invalid structures.
- Ensures uniqueness in parameter combinations before merging.

## Running the Test
```sh
g++ -o adaptive_gauss_batch_test  -Iinclude source/*.cpp test/adaptive_gauss_batch_test.cpp  -std=c++17 
 ./adaptive_gauss_batch_test
```



