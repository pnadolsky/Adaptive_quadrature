#include <laguerre_singular_endpoint.hpp>
#include <cmath>

//  Constructor: Passes lower & upper to `Quadrature`, sets singularity and alpha
LaguerreSingularEndpoint::LaguerreSingularEndpoint(const WeightsLoader& loader, int n1, int n2, double lower, double upper, bool leftIsSingular, double alpha)
    : LaguerreQuadrature(loader, n1, n2, true), leftIsSingular(leftIsSingular), alpha(alpha) {

    // Pass lower & upper to base Quadrature class
    this->lowerLimit = lower;
    this->upperLimit = upper;
}

double LaguerreSingularEndpoint::transformVariable(double t) const {
    double a = leftIsSingular ? lowerLimit.value() : upperLimit.value();
    double b = leftIsSingular ? upperLimit.value() : lowerLimit.value();
    double z = std::exp(-t/(1-alpha));
    return a*(1-z) + b * z; 
}

// weight function:
//  w(x) = (t-a)
double LaguerreSingularEndpoint::laguerre_weight_function(double t) const {
    double a = leftIsSingular ? lowerLimit.value() : upperLimit.value();
    double b = leftIsSingular ? upperLimit.value() : lowerLimit.value(); 
    double sign = leftIsSingular ? 1.0 : -1.0;   
    return sign*std::pow(t-a, alpha)*std::pow(b-a,1-alpha)/(1-alpha);  
}
