#include <iostream>
#include <unordered_map>
#include <variant>
#include "polylog_port.hpp"

int main() {
    try {
        // Define parameters as a map
        std::unordered_map<std::string, std::variant<int, double, std::string>> params;
        params["s"] = 3;    // Integer exponent
        params["z"] = 0.5;  // Double input parameter

        // Test value of t
        double t = 0.9;

        // Call the wrapper function
        double result = polylog_wrapper(params, t);

        // Output the result
        std::cout << "Polylogarithm result: " << result << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }

    return 0;
}
