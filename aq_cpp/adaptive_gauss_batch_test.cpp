#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <variant>
#include "polylog_port.hpp"
#include "quadrature.hpp"
#include "adaptive_gauss_tree.hpp"
#include "adaptive_gauss_batch.hpp"
#include <typeinfo>
#include <chrono>

int main() {
    std::function<double(ParamMap, double)> func = polylog_wrapper;
    double lower = 0.0, upper = 1.0, tol = 1e-12;
    int n1 = 100, n2=150, minD = 2,  maxD =20;
    double alphaA =0.0,  alphaB=0.0;
    bool singularA = true, singularB = false;
    WeightsLoader legendre_n1("../model_json/legendre.json");
    WeightsLoader legendre_n2("../model_json/legendre.json");
    WeightsLoader laguerre_n1("../model_json/laguerre.json");
    WeightsLoader laguerre_n2("../model_json/laguerre.json");    
    std::string name="Project",  author="Author",   description="project description"; 
    std::string reference="references",  version="1.0", update_log_message="Initial Train" ;

    // Define input parameter vectors
    std::vector<int> s = {2, 3, 4, 5,6, 7,8,9,10};
    std::vector<double> z = {0.1, 0.2,0.3, 0.4, 0.5,0.6, 0.7, 0.8,0.9, 1.0};
    std::vector<std::string> l ={"one", "two"};

    ParamMap test;
    test["s"] = 2; // s[0]
    test["z"] = 1.0;  // z[9]  
//    test["lab"] = "one"; // some label
    std::cout <<"Some parameter vals from our vectors:\t" << s[8] <<'\t'<<z[9]<<std::endl;
    std::cout <<"Map values:\t" << test["s"] <<'\t'<< test["z"] <<std::endl;
    // Store them in an unordered_map
    ParamCollection params;  // last in first out
//    params["lab"] = l;    
    params["s"] = s;
    params["z"] = z;

    auto start = std::chrono::high_resolution_clock::now();
    AdaptiveGaussTreeBatch batch = AdaptiveGaussTreeBatch(
        func, lower,  upper, tol,  minD,  maxD,  n1,  n2,
         alphaA,  alphaB, singularA,  singularB,
        legendre_n1, legendre_n2, laguerre_n1, laguerre_n2,
        params      
    );
    auto stop = std::chrono::high_resolution_clock::now();
     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);    
    std::cout << "Time to generate in milliseconds: \t" << duration.count() <<" ms"<<std::endl;
    batch.printCollection();
    //  void save_to_json(const std::string & filename, bool overwrite = false, bool write_roots=false, bool write_trees =true)
    batch.save_to_json("test.json", true,  true,  true );
    return 0;
}
