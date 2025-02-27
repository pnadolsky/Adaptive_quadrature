#ifndef LAGUERRE_SINGULAR_ENDPOINT_HPP
#define LAGUERRE_SINGULAR_ENDPOINT_HPP

#include "laguerre_quadrature.hpp"

class LaguerreSingularEndpoint : public LaguerreQuadrature {
private:
    bool leftIsSingular;
    double alpha;

protected:
    // ✅ Override Laguerre weight function (adjusting for singularity)
    double laguerre_weight_function(double t) const override;

public:
    // ✅ Constructor requires lower & upper and sets `leftIsSingular` and `alpha`
    LaguerreSingularEndpoint(const WeightsLoader& loader, int n1, int n2, double lower, double upper, bool leftIsSingular = true, double alpha = 0);

    // ✅ Override transformVariable based on singularity position
    double transformVariable(double t) const override;
};

#endif // LAGUERRE_SINGULAR_ENDPOINT_HPP
