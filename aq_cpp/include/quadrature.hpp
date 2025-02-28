#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP

#include <weights_loader.hpp>
#include <functional>
#include <unordered_map>
#include <map>
#include <variant>
#include <optional>
#include <string>

using ParamType = std::variant<int, double, std::string>;  // allows for mutable parameter types
using ParamMap = std::unordered_map<std::string, ParamType>;  // ParamMap: single set of values, e.g.:  {s: 1, z: 0.1, label: "A"}
// ParamCollection:  labels and vectors of values , e.g.:  {s: {1,2,3}, z: {0.1,0.2}, label: {"A","B","C"}}
using ParamCollection = std::map<std::string, std::variant<std::vector<int>, std::vector<double>, std::vector<std::string>>>;

std::ostream& operator<<(std::ostream& os, const ParamType& param);    //overload to print ParamType
std::ostream& operator<<(std::ostream& os, const ParamMap& paramMap);  //overload to print ParamMap

struct ParamTypeHash { std::size_t operator()(const ParamType& param) const;}; // Custom hasher for ParamType
struct ParamMapHash { std::size_t operator()(const ParamMap& paramMap) const ;}; // Custom hasher for ParamMap
struct ParamMapEqual { bool operator()(const ParamMap& lhs, const ParamMap& rhs) const ; }; // Custom equality for ParamMap 

class Quadrature {
protected:
    const WeightsLoader& weightsLoader;

    int order1, order2;
    double result, error;
    std::optional<double> lowerLimit, upperLimit;
    const std::string method;  // Now enforced at compile time!    
    std::vector<double> nodes1, weights1;
    std::vector<double> nodes2, weights2;

public:
    // Constructor requires `method` assignment in derived classes
    Quadrature(const WeightsLoader& loader, int n1, int n2, std::optional<double> lower, std::optional<double> upper, std::string methodName);

    virtual ~Quadrature() = default;
    virtual double integrate( std::function<double(ParamMap, double)> func, ParamMap parameters) = 0;
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
