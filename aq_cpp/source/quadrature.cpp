#include <quadrature.hpp>

// Constructor: Allows infinite limits using std::nullopt
Quadrature::Quadrature(const WeightsLoader& loader, int n1, int n2, std::optional<double> lower, std::optional<double> upper, std::string methodName)
    : weightsLoader(loader), order1(n1), order2(n2), result(0.0), error(0.0), lowerLimit(lower), upperLimit(upper), method(methodName) {

    if (!weightsLoader.hasOrder(order1) || !weightsLoader.hasOrder(order2)) {
        throw std::invalid_argument("Requested quadrature orders not found in WeightsLoader.");
    }

    nodes1 = weightsLoader.getNodes(order1);
    weights1 = weightsLoader.getWeights(order1);

    nodes2 = weightsLoader.getNodes(order2);
    weights2 = weightsLoader.getWeights(order2);
}

// Default transformation: Handles finite limits only
double Quadrature::transformVariable(double t) const {
    if (lowerLimit.has_value() && upperLimit.has_value()) {
        // Standard change of interval: [-1,1] → [lower, upper]
        return (upperLimit.value() - lowerLimit.value()) / 2.0 * t + (upperLimit.value() + lowerLimit.value()) / 2.0;
    }

    throw std::logic_error("transformVariable() must be overridden for infinite limits.");
}

//overload to print ParamType
std::ostream& operator<<(std::ostream& os, const ParamType& param) {
    std::visit([&os](auto&& value) {
        os << value;
    }, param);
    return os;
}

//overload to print ParamMap
std::ostream& operator<<(std::ostream& os, const ParamMap& paramMap) {
    os << "{ ";
    for (const auto& [key, value] : paramMap) {
        os << key << ": ";
        std::visit([&os](auto&& val) { os << val; }, value);
        os << ", ";
    }
    os << "}";
    return os;
}

// Custom hasher for ParamType
//struct ParamTypeHash {
    std::size_t ParamTypeHash::operator()(const ParamType& param) const {
        return std::visit([](auto&& arg) -> std::size_t {
            return std::hash<std::decay_t<decltype(arg)>>{}(arg);
        }, param);
    }
//};

// Hash function for ParamMap
//struct ParamMapHash {
    std::size_t ParamMapHash::operator()(const ParamMap& paramMap) const {
        std::size_t hashValue = 0;
        for (const auto& [key, value] : paramMap) {
            std::hash<std::string> hashString;
            std::hash<int> hashInt;
            std::hash<double> hashDouble;

            hashValue ^= hashString(key);

            std::visit([&](auto&& val) {
                using T = std::decay_t<decltype(val)>;
                if constexpr (std::is_same_v<T, int>) {
                    hashValue ^= hashInt(val);
                } else if constexpr (std::is_same_v<T, double>) {
                    hashValue ^= hashDouble(val);
                } else if constexpr (std::is_same_v<T, std::string>) {
                    hashValue ^= hashString(val);
                }
            }, value);
        }
        return hashValue;
    }
//};

// Equality function for ParamMap
//struct ParamMapEqual {
    bool ParamMapEqual::operator()(const ParamMap& lhs, const ParamMap& rhs) const {
        if (lhs.size() != rhs.size()) return false;
        for (const auto& [key, value] : lhs) {
            auto it = rhs.find(key);
            if (it == rhs.end()) return false;
            if (value.index() != it->second.index()) return false;
            bool valuesEqual = std::visit(
                [](const auto& a, const auto& b) -> bool {
                    using A = std::decay_t<decltype(a)>;
                    using B = std::decay_t<decltype(b)>;
                    if constexpr (std::is_same_v<A, B>) {
                        return a == b;
                    } else {
                        return false;
                    }
                },
                value, it->second
            );
            if (!valuesEqual) return false;
        }
        return true;
    }
//};

