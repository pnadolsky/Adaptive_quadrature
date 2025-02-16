#ifndef LEGENDRE_QUADRATURE_HPP
#define LEGENDRE_QUADRATURE_HPP

#include "quadrature.hpp"

class LegendreQuadrature : public Quadrature {
public:
    LegendreQuadrature(const WeightsLoader& loader, int n1, int n2, double lower, double upper);

    double integrate(
        std::function<double(std::unordered_map<std::string, std::variant<int, double, std::string>>, double)> func,
        std::unordered_map<std::string, std::variant<int, double, std::string>> parameters) override;

    // Transform variable from [-1,1] to [lower,upper]
    double transformVariable(double t) const override;
};

#endif // LEGENDRE_QUADRATURE_HPP
