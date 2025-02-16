#include <iostream>
#include <cmath>
#include "weights_loader.hpp"
#include "laguerre_quadrature.hpp"

// Example function: exp(-x^2) for Laguerre integration
double gaussian_function(std::unordered_map<std::string, std::variant<int, double, std::string>> params, double t) {
    return std::exp(-t * t);
}

double gaussian_function_with_weight_included(std::unordered_map<std::string, std::variant<int, double, std::string>> params, double t) {
    return std::exp(t-t * t);
}

int main() {
    try {
        // Load Laguerre quadrature weights
        WeightsLoader laguerre_loader("..\\model_json\\laguerre.json");

        // Create Laguerre Quadrature instance (default: use weight function)
        // LaguerreQuadrature laguerre_weighted(laguerre_loader, 10, 20, true); 
        LaguerreQuadrature laguerre_weighted(laguerre_loader, 10, 20);
        // Perform integration with weight function
        double result_weighted = laguerre_weighted.integrate(gaussian_function, {});

        std::cout << "\n=== Laguerre Quadrature with Weight Function ===\n";
        std::cout << "Integral of exp(-x^2) * exp(-x) from 0 to Infinty: " << result_weighted << std::endl;
        std::cout << "Estimated Error: " << laguerre_weighted.getError() << std::endl;

        // Create Laguerre Quadrature instance (without weight function)
        LaguerreQuadrature laguerre_unweighted(laguerre_loader, 10, 20, false);

        // Perform integration without weight function
        double result_unweighted = laguerre_unweighted.integrate(gaussian_function_with_weight_included, {});

        std::cout << "\n=== Laguerre Quadrature without Weight Function ===\n";
        std::cout << "Integral of exp(x-x^2) from 0 to Infinity: " << result_unweighted << std::endl;
        std::cout << "Estimated Error: " << laguerre_unweighted.getError() << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }

    return 0;
}
