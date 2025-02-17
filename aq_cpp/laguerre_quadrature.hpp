#ifndef LAGUERRE_QUADRATURE_HPP
#define LAGUERRE_QUADRATURE_HPP

#include "quadrature.hpp"

class LaguerreQuadrature : public Quadrature {
private:
    bool use_weight_function;  // Flag to apply weight function (default: true)

    // Laguerre weight function 
    virtual double laguerre_weight_function(double t) const {
        return std::exp(t);
    }

public:
    // Constructor with optional use_weight_function argument
    LaguerreQuadrature(const WeightsLoader& loader, int n1, int n2, bool use_weight_function = true);

    double integrate( std::function<double(ParamMap, double)> func, ParamMap parameters) override;

    // Laguerre maps [-1,1] to [0,∞) (no transformation needed)
    virtual double transformVariable(double t) const override;
};

#endif // LAGUERRE_QUADRATURE_HPP
