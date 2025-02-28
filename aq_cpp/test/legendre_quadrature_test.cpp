#include <iostream>
#include <iomanip> // Required for setprecision
#include "weights_loader.hpp"
#include "legendre_quadrature.hpp"
#include "polylog_port.hpp"


int main() {
    try {
        // Load quadrature weights
        WeightsLoader loader("..\\model_json\\legendre.json");

        // Define function parameters
        ParamMap params;
        params["s"] = 2;
        params["z"] = 1.0;

        // Create LegendreQuadrature with two orders (100 and 180 for error estimation)
        LegendreQuadrature legendre(loader, 100, 180, 0.0, 1.0);

        // Perform integration
        //double result = legendre.integrate(polylog_wrapper, params);

        // Output result and error estimation
        std::cout << "Legendre Quadrature Integral result: " << std::setprecision(15) << legendre.getResult() << std::endl;
        std::cout << "Estimated Error: " << legendre.getError() << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }

    return 0;
}
