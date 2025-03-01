#include <iostream>
#include <cmath>
#include "adaptive_gauss_tree.hpp"
#include "polylog_port.hpp"

// Example function: log(x) / sqrt(x)
double test_function(ParamMap params, double x) {
    (void) params;
    if (x == 0.0) return 0.0;  // Avoid log(0) undefined case
    return std::log(x) / std::sqrt(x);
}

int main() {
    try {
        // Load quadrature weights
        WeightsLoader legendre_n1("../model_json/legendre.json");
        WeightsLoader legendre_n2("../model_json/legendre.json");
        WeightsLoader laguerre_n1("../model_json/laguerre.json");
        WeightsLoader laguerre_n2("../model_json/laguerre.json");
        
        // Create AdaptiveGaussTree instance from parameters
        AdaptiveGaussTree adaptive_tree(test_function, 0.0, 1.0, 1e-6, 2, 10, 
                                        100, 150,  // Pass explicit n1 and n2
                                        0.5, 0.5, true, false,  
                                        legendre_n1, legendre_n2, laguerre_n1, laguerre_n2);
        
        // Get integral and error
        auto [integral, error] = adaptive_tree.get_integral_and_error();
        std::cout << "Computed Integral: " << integral << "\n";
        std::cout << "Estimated Error: " << error << "\n";
        
        // Save to JSON 
        adaptive_tree.save_to_json("adaptive_output.json");  // exists after first time running this!
        adaptive_tree.save_to_json("adaptive_output.json", true);    // set flag to override  
        std::cout << "Adaptive quadrature tree saved to adaptive_output.json" << std::endl;

        // Load from JSON
        AdaptiveGaussTree loaded_tree(test_function, legendre_n1, legendre_n2, laguerre_n1, laguerre_n2, "adaptive_output.json");
        std::cout << "Adaptive quadrature tree loaded from JSON file." << std::endl;

        auto [integral2, error2] = loaded_tree.get_integral_and_error();
        std::cout << "Loaded Integral: " << integral2 << "\n";
        std::cout << "Loaded Estimated Error: " << error2 << "\n";        

        // Load from JSON generated by python NB

        AdaptiveGaussTree loaded_tree_2(test_function, legendre_n1, legendre_n2, laguerre_n1, laguerre_n2, "../test_dump.json");
        auto [integral3, error3] = loaded_tree_2.get_integral_and_error();
        std::cout << "Loaded Integral (from python output): " << integral3 << "\n";
        std::cout << "Loaded Estimated Error (from python output): " << error3 << "\n";     

        AdaptiveGaussTree loaded_tree_3(test_function, legendre_n1, legendre_n2, laguerre_n1, laguerre_n2, "../test.json");
        auto [integral4, error4] = loaded_tree_3.get_integral_and_error();
        std::cout << "Loaded Integral (from python output): " << integral4 << "\n";
        std::cout << "Loaded Estimated Error (from python output): " << error4 << "\n";             


        // Create AdaptiveGaussTree instance with parameters (e.g. polylog_wrapper; expects 's' and 'z')
        ParamMap poly_args;
        poly_args["s"] = 2;    // dilog
        poly_args["z"] = 1.0;  // z= 1 (e.g. Li_2(1)=zeta(2) )

        AdaptiveGaussTree poly_tree(polylog_wrapper, 0.0, 1.0, 1e-12, 2, 10, // min amnd max depth
            100, 150,  // Pass explicit n1 and n2
            0.0, 0.0, true, false, 
            legendre_n1, legendre_n2, laguerre_n1, laguerre_n2,
            poly_args  // <-------------------   pass optional args at the end !
        );
        auto [integral5, error5] = poly_tree.get_integral_and_error();
        std::cout << "Polylog Integral (e.g. zeta(2)): " << integral5 << "\n";
        std::cout << "Polylog Estimated Error: ~ < 10^-12: " << error5 << "\n";   

        // test copy constructor

        AdaptiveGaussTree poly_tree_2 = poly_tree;
        auto [integral6, error6] = poly_tree_2.get_integral_and_error();
        std::cout << "Polylog (Copy) Integral (e.g. zeta(2)): " << integral6 << "\n";
        std::cout << "Polylog (Copy) Estimated Error: ~ < 10^-12: " << error6 << "\n";           

    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }
    return 0;
}
