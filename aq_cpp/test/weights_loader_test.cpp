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

//   Testing the json functionality.
     
    //WeightsLoader loader("..\\model_json\\legendre.json");
    std::ifstream file("test.json");
    json file_json;
    file >> file_json;
    file.close();

    WeightsLoader wl = WeightsLoader(file_json, "legendre_roots_n1", "Legendre", "n1");
    std::vector<double> nodes = wl.getNodes(100);
    std::vector<double> weights = wl.getWeights(100);
// Print the values
    std::cout <<"from test.json:" << wl.getMethod() <<std::endl;
    std::cout <<"n_max (= one value):" << wl.getNMax() <<std::endl;
    std::cout << "Nodes: (" << std::endl;
    for (double val : nodes) {
            std::cout << val << " ";
    }
    std::cout << std::endl;  
    std::cout << "Weights: (" << std::endl;
    for (double val : weights) {
            std::cout << val << " ";
    }
    std::cout << std::endl;  

    return 0;
}


