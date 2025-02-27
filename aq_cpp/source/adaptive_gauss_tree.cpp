#include "adaptive_gauss_tree.hpp"

std::ostream& operator<<(std::ostream& os, const AdaptiveGaussTree& tree) {
    auto [integral, error] = tree.get_integral_and_error();
    os << "( integral: " << integral << ", error: " << error << " )";
    return os;
}