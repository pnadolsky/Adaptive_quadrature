#include <iostream>
#include "weights_loader.hpp"

int main() {
    try {
        WeightsLoader loader("..\\model_json\\legendre.json");
//        WeightsLoader loader("..\\model_json\\laguerre.json");
        std::cout << "Method: " << loader.getMethod() << std::endl;
        std::cout << "Max Order (n_max): " << loader.getNMax() << std::endl;

        int order = 10;
        if (loader.hasOrder(order)) {
            std::cout << "Nodes for order " << order << ": ";
            for (double node : loader.getNodes(order)) {
                std::cout << node << " ";
            }
            std::cout << std::endl;

            std::cout << "Weights for order " << order << ": ";
            for (double weight : loader.getWeights(order)) {
                std::cout << weight << " ";
            }
            std::cout << std::endl;
        } else {
            std::cout << "Order " << order << " not found in JSON file." << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}


