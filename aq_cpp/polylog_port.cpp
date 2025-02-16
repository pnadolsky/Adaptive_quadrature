#include "polylog_port.hpp"
#include <cmath>
#include <unordered_map>  // Ensure this is included!
#include <variant>

// Function to compute the polylogarithm integrand
double polylog_integrand(int s, double z, double t) {
    int sgn = (s % 2 == 0) ? -1 : 1;
    return sgn * z * std::pow(std::log(t), s - 1) / (std::tgamma(s) * (1.0 - t * z));
}

// Wrapper function to evaluate and return updated parameters
double  polylog_wrapper(
    std::unordered_map<std::string,  std::variant<int, double, std::string>> parameters, 
    double t) {
    int s = std::get<int>(parameters["s"]);
    double z = std::get<double>(parameters["z"]);
    double result =  polylog_integrand(s, z, t);
    return result;
}
