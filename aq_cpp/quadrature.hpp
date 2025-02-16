#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP

#include "weights_loader.hpp"
#include <functional>
#include <unordered_map>
#include <variant>
#include <optional>
#include <string>

class Quadrature {
protected:
    const WeightsLoader& weightsLoader;

    int order1, order2;
    double error, result;
    std::optional<double> lowerLimit, upperLimit;
    std::vector<double> nodes1, weights1;
    std::vector<double> nodes2, weights2;

    const std::string method;  // Now enforced at compile time!

public:
    // Constructor requires `method` assignment in derived classes
    Quadrature(const WeightsLoader& loader, int n1, int n2, std::optional<double> lower, std::optional<double> upper, std::string methodName);

    virtual ~Quadrature() = default;

    virtual double integrate(
        std::function<double(std::unordered_map<std::string, std::variant<int, double, std::string>>, double)> func,
        std::unordered_map<std::string, std::variant<int, double, std::string>> parameters) = 0;

    virtual double transformVariable(double t) const;

    // Getters 
    std::string get_method() const { return method; }
    int getOrder1() const { return order1; }
    int getOrder2() const { return order2; }
    double getResult() const { return result; }
    double getError() const { return error; }
    std::optional<double> getLowerLimit() const { return lowerLimit; }
    std::optional<double> getUpperLimit() const { return upperLimit; }
};

#endif // QUADRATURE_HPP
