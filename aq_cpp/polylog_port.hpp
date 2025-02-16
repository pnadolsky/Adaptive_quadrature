#ifndef POLYLOG_PORT_H
#define POLYLOG_PORT_H

#include <unordered_map>
#include <string>
#include <unordered_map>  // Ensure this is included!
#include <variant>

double polylog_integrand(int s, double z, double t);
double polylog_wrapper(
    std::unordered_map<std::string,  std::variant<int, double, std::string>> parameters, 
    double t);

#endif // POLYLOG_PORT_H
