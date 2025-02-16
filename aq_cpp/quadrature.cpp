#include "quadrature.hpp"

// Constructor: Allows infinite limits using std::nullopt
Quadrature::Quadrature(const WeightsLoader& loader, int n1, int n2, std::optional<double> lower, std::optional<double> upper, std::string methodName)
    : weightsLoader(loader), order1(n1), order2(n2), lowerLimit(lower), upperLimit(upper), result(0.0), error(0.0), method(methodName) {

    if (!weightsLoader.hasOrder(order1) || !weightsLoader.hasOrder(order2)) {
        throw std::invalid_argument("Requested quadrature orders not found in WeightsLoader.");
    }

    nodes1 = weightsLoader.getNodes(order1);
    weights1 = weightsLoader.getWeights(order1);

    nodes2 = weightsLoader.getNodes(order2);
    weights2 = weightsLoader.getWeights(order2);
}

// Default transformation: Handles finite limits only
double Quadrature::transformVariable(double t) const {
    if (lowerLimit.has_value() && upperLimit.has_value()) {
        // Standard change of interval: [-1,1] → [lower, upper]
        return (upperLimit.value() - lowerLimit.value()) / 2.0 * t + (upperLimit.value() + lowerLimit.value()) / 2.0;
    }

    throw std::logic_error("transformVariable() must be overridden for infinite limits.");
}
