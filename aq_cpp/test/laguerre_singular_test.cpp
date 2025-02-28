#include <iostream>
#include <cmath>
#include "weights_loader.hpp"
#include "laguerre_quadrature.hpp"
#include "laguerre_singular_endpoint.hpp"
#include "legendre_quadrature.hpp"

// ✅ Updated test function: log(x) / sqrt(x)
double singular_test_function(ParamMap params, double x) {
    (void) params;
    if (x == 0.0) return 0.0;  // Avoid log(0) undefined case
    return std::log(x) / std::sqrt(x);
//    return std::log(x);
}

int main() {
    try {
        // Load Laguerre quadrature weights
        WeightsLoader laguerre_loader("..\\model_json\\laguerre.json");
        WeightsLoader legendre_loader("..\\model_json\\legendre.json");
        double alpha =  0.5;
        // ✅ Standard Legnedre Quadrature
        LegendreQuadrature legendre(legendre_loader, 200, 250, 0.0001, 1.0);
        //double result_standard = legendre.integrate(singular_test_function, {});

        std::cout << "\n=== Standard " << legendre.get_method() << " Quadrature ===\n";
        std::cout << "Integral Result: " << legendre.getResult() << std::endl;
        std::cout << "Estimated Error: " << legendre.getError() << std::endl;

        // ✅ Laguerre Quadrature with Singular Left Endpoint, alpha = 0.5
        LaguerreSingularEndpoint laguerre_singular_left(laguerre_loader, 200, 250, 0.0, 1.0, true, alpha);
        //double result_singular_left = laguerre_singular_left.integrate(singular_test_function, {});

        std::cout << "\n=== Laguerre Quadrature with Singular Left Endpoint (alpha = 0.5) ===\n";
        std::cout << "Integral Result: " << laguerre_singular_left.getResult() << std::endl;
        std::cout << "Estimated Error: " << laguerre_singular_left.getError() << std::endl;

        // ✅ Laguerre Quadrature with Singular Right Endpoint, alpha = 0.5.  reversed limits
        LaguerreSingularEndpoint laguerre_singular_right(laguerre_loader, 200, 250, 1.0, 0.0, false, alpha);
        //double result_singular_right = laguerre_singular_right.integrate(singular_test_function, {});

        std::cout << "\n=== Laguerre Quadrature with Singular Right Endpoint (alpha = 0.5) ===\n";
        std::cout << "Integral Result: " << laguerre_singular_right.getResult() << std::endl;
        std::cout << "Estimated Error: " << laguerre_singular_right.getError() << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }

    return 0;
}
