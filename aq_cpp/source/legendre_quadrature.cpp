#include "legendre_quadrature.hpp"
#include <iostream>
#include <cmath>

// Constructor
LegendreQuadrature::LegendreQuadrature(const WeightsLoader& loader, int n1, int n2, double lower, double upper)
    : Quadrature(loader, n1, n2, lower, upper, "Gauss-Legendre") {}

// Transform variable from [-1,1] to [lower,upper]
double LegendreQuadrature::transformVariable(double t) const {
    return (upperLimit.value() - lowerLimit.value()) / 2.0 * t + (upperLimit.value() + lowerLimit.value()) / 2.0;
}

// Perform integration using two orders for error estimation
double LegendreQuadrature::integrate( std::function<double(ParamMap, double)> func, ParamMap parameters) {

    double integral1 = 0.0, integral2 = 0.0;

    double half_length = (upperLimit.value() - lowerLimit.value()) / 2.0;

    // Compute integral using order1
    for (size_t i = 0; i < nodes1.size(); ++i) {
        double t = transformVariable(nodes1[i]);
        integral1 += weights1[i] * func(parameters, t);
    }
    integral1 *= half_length;  // Adjust for interval change

    // Compute integral using order2
    for (size_t i = 0; i < nodes2.size(); ++i) {
        double t = transformVariable(nodes2[i]);
        integral2 += weights2[i] * func(parameters, t);
    }
    integral2 *= half_length;  // Adjust for interval change

    // Store results
    result = integral1;
    error = std::abs(integral1 - integral2);  // Compute error estimation

    return result;
}
