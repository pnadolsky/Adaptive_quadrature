#ifndef WEIGHTS_LOADER_HPP
#define WEIGHTS_LOADER_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <nlohmann/json.hpp>

using json = nlohmann::ordered_json;

class WeightsLoader {
private:
    std::unordered_map<int, std::vector<double>> nodes;
    std::unordered_map<int, std::vector<double>> weights;
    std::string method;  //  e.g. "Legendre" from header
    int n_max;  // New field to store maximum order

public:
    WeightsLoader(){};
    // Constructor that loads a weights JSON file
    WeightsLoader(const std::string& filename);
    // This Constructor loads a single weight set as generated in the adaptive quadrature jsons. Uses the dictonary keys
    // WeightsLoader(..., "legendre_roots_n1", "Legendre","n1")
    WeightsLoader(json js, const std::string& key, const std::string& method, const std::string& n_key); 
    // Getter functions
    std::vector<double> getNodes(int n) const;
    std::vector<double> getWeights(int n) const;
    std::string getMethod() const;
    int getNMax() const;  // New method to retrieve max order
// Add this inside `WeightsLoader` in weights_loader.hpp:
    WeightsLoader(const WeightsLoader&) = default;  // Explicitly use the default copy constructor
    WeightsLoader& operator=(const WeightsLoader&) = default;  // Default assignment operator

    // Function to check if order exists
    bool hasOrder(int n) const;
};

#endif // WEIGHTS_LOADER_HPP
