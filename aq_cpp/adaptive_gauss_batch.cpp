#include "adaptive_gauss_batch.hpp"

// recursively generate all parameter combinations
void AdaptiveGaussTreeBatch::generate_combinations(
    const std::vector<std::string>& keys,
    const ParamCollection& params,
    std::vector<size_t>& indices,
    std::vector<ParamMap>& results,
    size_t depth
) {
    if (depth == keys.size()) {
        // Construct the current combination
        ParamMap combination;
        for (size_t i = 0; i < keys.size(); ++i) {
            const std::string& key = keys[i];

            // Get the correct type from params
            if (std::holds_alternative<std::vector<int>>(params.at(key))) {
                combination[key] = std::get<std::vector<int>>(params.at(key))[indices[i]];
            } else if (std::holds_alternative<std::vector<double>>(params.at(key))) {
                combination[key] = std::get<std::vector<double>>(params.at(key))[indices[i]];
            } else if (std::holds_alternative<std::vector<std::string>>(params.at(key))) {
                combination[key] = std::get<std::vector<std::string>>(params.at(key))[indices[i]];
            }
        }
        results.push_back(combination);
        return;
    }

    // Iterate through current parameter vector
    const std::string& key = keys[depth];
    size_t size = std::visit([](auto&& vec) { return vec.size(); }, params.at(key));

    for (size_t i = 0; i < size; ++i) {
        indices[depth] = i;
        generate_combinations(keys, params, indices, results, depth + 1);
    }
}

bool AdaptiveGaussTreeBatch::compareParamMaps(const ParamMap& a, const ParamMap& b) {
    for (const auto& key : keys) {
        if (a.count(key) == 0 || b.count(key) == 0) continue; // Skip if key missing

        const ParamType& valA = a.at(key);
        const ParamType& valB = b.at(key);

        if (valA == valB) continue; // If equal, check next key

        return compareVariant(valA, valB);
    }
    return false; // If all keys are equal, maintain order
}

void AdaptiveGaussTreeBatch::sortResults() {
    std::sort(results.begin(), results.end(), 
        [this](const ParamMap& a, const ParamMap& b) {
            return compareParamMaps(a, b);
        }
    );
}

bool AdaptiveGaussTreeBatch::compareVariant(const ParamType& a, const ParamType& b) {
    return std::visit([](const auto& lhs, const auto& rhs) -> bool {
        using L = std::decay_t<decltype(lhs)>;
        using R = std::decay_t<decltype(rhs)>;

        if constexpr (std::is_same_v<L, R>) {
            return lhs < rhs; // Compare if same type
        } else {
            return AdaptiveGaussTreeBatch::getTypeRank<L>() < AdaptiveGaussTreeBatch::getTypeRank<R>(); // Order: int < double < string
        }
    }, a, b);
}

// Define getTypeRank as a static function
template <typename T>
int AdaptiveGaussTreeBatch::getTypeRank() {
    if constexpr (std::is_same_v<T, int>) return 1;
    if constexpr (std::is_same_v<T, double>) return 2;
    if constexpr (std::is_same_v<T, std::string>) return 3;
    return 4; // Default case
}

