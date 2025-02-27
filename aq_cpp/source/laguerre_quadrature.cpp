#include <laguerre_quadrature.hpp>
#include <cmath>

// Constructor: Lower limit = 0, upper limit = infinity, with weight function flag
LaguerreQuadrature::LaguerreQuadrature(const WeightsLoader& loader, int n1, int n2, bool use_weight_function)
    : Quadrature(loader, n1, n2, 0.0, std::nullopt, "Gauss-Laguerre"), use_weight_function(use_weight_function) {}

// Transform variable for Laguerre quadrature (maps nodes directly to [0,∞))
double LaguerreQuadrature::transformVariable(double t) const {
    return t;  // No transformation needed for Laguerre quadrature
}

// Perform integration using two orders for error estimation
double LaguerreQuadrature::integrate( std::function<double(ParamMap, double)> func, ParamMap parameters) {

    double integral1 = 0.0, integral2 = 0.0;

    // Compute integral using order1
    for (size_t i = 0; i < nodes1.size(); ++i) {
        double t = transformVariable(nodes1[i]);
        double weight = use_weight_function ? laguerre_weight_function(t) : 1.0;
        integral1 += weights1[i] * func(parameters, t) * weight;
    }

    // Compute integral using order2
    for (size_t i = 0; i < nodes2.size(); ++i) {
        double t = transformVariable(nodes2[i]);
        double weight = use_weight_function ? laguerre_weight_function(t) : 1.0;
        integral2 += weights2[i] * func(parameters, t) * weight;
    }

    // Store results
    result = integral1;
    error = std::abs(integral1 - integral2);  // Compute error estimation

    return result;
}
