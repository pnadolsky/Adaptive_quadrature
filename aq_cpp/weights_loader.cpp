#include "weights_loader.hpp"

WeightsLoader::WeightsLoader(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Error opening JSON file: " + filename);
    }

    json j;
    file >> j;
    file.close();

    // Load method and n_max
    method = j["method"].get<std::string>();
    n_max = j["n_max"].get<int>();  // Load n_max from JSON

    // Load nodes and weights for each order
    for (const auto& [key, value] : j["n"].items()) {
        int order = std::stoi(key);
        nodes[order] = value["0"].get<std::vector<double>>();
        weights[order] = value["1"].get<std::vector<double>>();
    }
}

std::vector<double> WeightsLoader::getNodes(int n) const {
    if (nodes.find(n) != nodes.end()) {
        return nodes.at(n);
    } else {
        throw std::out_of_range("Order not found in JSON data: " + std::to_string(n));
    }
}

std::vector<double> WeightsLoader::getWeights(int n) const {
    if (weights.find(n) != weights.end()) {
        return weights.at(n);
    } else {
        throw std::out_of_range("Order not found in JSON data: " + std::to_string(n));
    }
}

std::string WeightsLoader::getMethod() const {
    return method;
}

int WeightsLoader::getNMax() const {
    return n_max;
}

bool WeightsLoader::hasOrder(int n) const {
    return nodes.find(n) != nodes.end() && weights.find(n) != weights.end();
}
