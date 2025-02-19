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


