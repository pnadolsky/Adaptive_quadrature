#ifndef WEIGHTS_LOADER_HPP
#define WEIGHTS_LOADER_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include "json.hpp"

using json = nlohmann::ordered_json;

class WeightsLoader {
private:
    std::unordered_map<int, std::vector<double>> nodes;
    std::unordered_map<int, std::vector<double>> weights;
    std::string method;  //  e.g. "Legendre" from header
    int n_max;  // New field to store maximum order

public:
    // Constructor that loads JSON file
    WeightsLoader(const std::string& filename);

    // Getter functions
    std::vector<double> getNodes(int n) const;
    std::vector<double> getWeights(int n) const;
    std::string getMethod() const;
    int getNMax() const;  // New method to retrieve max order

    // Function to check if order exists
    bool hasOrder(int n) const;
};

#endif // WEIGHTS_LOADER_HPP
