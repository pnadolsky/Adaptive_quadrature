#include <adaptive_gauss_batch.hpp>

// Constructor for loading from a serialized JSON tree
AdaptiveGaussTreeBatch::AdaptiveGaussTreeBatch(
    std::function<double(ParamMap, double)> func,
    std::string filename
): func(func), lower(0.0), upper(1.0), alphaA(0.0), alphaB(0.0) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Error opening file: " + filename);
    }

    json data;
    file >> data;

    // Load basic metadata
    name = data["name"];
    reference = data["reference"];
    description = data["description"];
    author = data["author"];
    version = data["version"];
    tol = data["tol"].get<double>();
    min_depth = data["min_depth"].get<int>();
    max_depth = data["max_depth"].get<int>();
    order1 = data["n1"].get<int>();
    order2 = data["n2"].get<int>();
    a_singular = data["a_singular"];
    b_singular = data["b_singular"];

    // Load weights for quadrature
   
    legendre_n1 = WeightsLoader(data, "legendre_roots_n1", "Legendre", "n1");
    legendre_n2 = WeightsLoader(data, "legendre_roots_n2", "Legendre", "n2");
    laguerre_n1 = WeightsLoader(data, "laguerre_roots_n1", "Laguerre", "n1");
    laguerre_n2 = WeightsLoader(data, "laguerre_roots_n2", "Laguerre", "n2");

    // Load update log

    if (data.contains("update_log") && data["update_log"].is_array()) {
        update_log.clear();
        for (const auto& entry : data["update_log"]) {
            if (entry.contains("timestamp") && entry.contains("message")) {
                std::string timestamp = entry["timestamp"].get<std::string>();
                std::string message = entry["message"].get<std::string>();
                update_log.emplace_back(timestamp, message);
            }
        }
    } else {
        throw std::runtime_error("Invalid format for update_log: Expected an array of objects.");
    }
    
    // Parse the parameters section

    if (data.contains("parameters")) {
        json parameters_json = data["parameters"];
        
        // Extract keys
        std::set<std::string> key_set;
        extract_keys(parameters_json, key_set);
        keys.assign(key_set.begin(), key_set.end());

        // Recursively extract parameter sets and trees

        extract_parameters(parameters_json, parameters);
    }
   
    std::vector<int> s_values = std::get<std::vector<int>>(parameters["s"]);

    std::vector<double> z_values = std::get<std::vector<double>>(parameters["z"]);
  

    // Generate all possible parameter combinations

    std::vector<size_t> indices(keys.size(), 0);

  
    generate_combinations(keys, parameters, indices, results);

    sortResults();




    // Reconstruct trees from stored JSON
    for (const auto& param_map : results) {
 //       std::cout << "reconstructing tree" << std::endl;
//        const json& tree_json = find_tree_json(data["parameters"], param_map);
        json json_head;
        json tree_json = find_tree_json(data["parameters"], param_map);
//        std::cout << "tree dump" << std::endl;
//        std::cout << tree_json.dump(4) << std::endl; 
        // add tree fields to tree_json.  need these:
/*
        name = data["name"];
        reference=data["reference"];
        description=data["description"];
        author=data["author"];
        version=data["version"];       
        tolerance = data["tolerance"];
        min_depth = data["min_depth"];
        max_depth = data["max_depth"];
        order1 = data["n1"];
        order2 = data["n2"];

                            "a": 0.0,
                            "b": 1.0,
                            "depth": 0,
                            "tol": 1e-12,
                            "error": 5.773159728050814e-15,
                            "integral": 1.0009945751278209,
                            "method": "Gauss-Laguerre"        
*/        
        json_head["name"] = name;
        json_head["reference"]= reference;
        json_head["description"]= description;
        json_head["author"]=author;
        json_head["version"]=version;
        json_head["tolerance"]= tol;
        json_head["min_depth"]=min_depth;
        json_head["max_depth"]=max_depth;
        json_head["n1"]=order1;
        json_head["n2"]=order2;
//        std::cout << "assign tree" <<std::endl;
        json_head["tree"] = tree_json;
//        std::cout << "good" <<std::endl;
//        std::cout << "head dump" << std::endl;        
//        std::cout << json_head.dump(4) << std::endl; 
//        std::cout << "completed" << std::endl;
/*
        auto node = std::make_unique<Node>(data["a"], data["b"], data["depth"],
                                           data["tol"],order1, order2, data["method"] == "Gauss-Laguerre");
        node->error = data["error"];
        node->result = data["integral"];
*/

        std::unique_ptr<AdaptiveGaussTree> tree_ptr = std::make_unique<AdaptiveGaussTree>(
            json_head, func, legendre_n1, legendre_n2, laguerre_n1, laguerre_n2 , param_map
        );
        AdaptiveGaussTree & the_tree = *tree_ptr;
//        std::cout << "tree" << std::endl;
//        std::cout << the_tree << std::endl;
        quad_coll[param_map] = std::move(tree_ptr);
    }
}

void AdaptiveGaussTreeBatch::merge(const AdaptiveGaussTreeBatch& other) {
    // Merge quad_coll (deep copy of AdaptiveGaussTree)
    for (const auto& pair : other.quad_coll) {
        const ParamMap& param_key = pair.first;

        // Check if key already exists
        if (quad_coll.find(param_key) != quad_coll.end()) {
            std::cout << "Warning: Duplicate key found in merge(). Keeping existing tree for key: " << param_key << std::endl;
            continue; // Skip or modify this behavior if needed
        }

        // Deep copy and insert new AdaptiveGaussTree
        quad_coll[param_key] = std::make_unique<AdaptiveGaussTree>(*pair.second);
    }

    // Merge results
    results.insert(results.end(), other.results.begin(), other.results.end());

    // Merge update log (keeping chronological order)
    update_log.insert(update_log.end(), other.update_log.begin(), other.update_log.end());

    // Preserve unique keys
    std::set<std::string> key_set(keys.begin(), keys.end());
    key_set.insert(other.keys.begin(), other.keys.end());
    keys.assign(key_set.begin(), key_set.end());

    // Merge parameters (combine values for matching keys)
    for (const auto& [key, value] : other.parameters) {
        if (parameters.find(key) != parameters.end()) {
            // Append new values if they don't already exist
            std::visit([&](auto& vec) {
                using T = std::decay_t<decltype(vec)>;
                auto& existing_vec = std::get<T>(parameters[key]);
                for (const auto& val : vec) {
                    if (std::find(existing_vec.begin(), existing_vec.end(), val) == existing_vec.end()) {
                        existing_vec.push_back(val);
                    }
                }
            }, value);
        } else {
            parameters[key] = value; // Add new parameter entry
        }
    }

    // Log the merge operation
    add_update_log("Merged with another AdaptiveGaussTreeBatch instance.");
}




void AdaptiveGaussTreeBatch::extract_keys(const json& param_json, std::set<std::string>& key_set, int depth) {
    for (auto it = param_json.begin(); it != param_json.end(); ++it) {
        if (it.key() == "tree") continue;  // Ignore "tree"

        // **Only store non-numeric keys at the first depth**
        if (depth % 2 == 0 ) {
            key_set.insert(it.key());
        }

        // Recursively process deeper levels
        if (it->is_object()) {
            extract_keys(*it, key_set, depth + 1);
        }
    }
}

// using ParamCollection = std::map<std::string, std::variant<std::vector<int>, std::vector<double>, std::vector<std::string>>>

void AdaptiveGaussTreeBatch::extract_parameters(
    const json& param_json, 
    ParamCollection& parameters, 
    std::optional<std::string> current_key, 
    bool expecting_value, 
    std::map<std::string, std::set<int>> int_sets, 
    std::map<std::string, std::set<double>> double_sets , 
    std::map<std::string, std::set<std::string>> string_sets ) 
{
    for (auto it = param_json.begin(); it != param_json.end(); ++it) {
        std::string key = it.key();

        // 🚨 Skip "tree" since it contains no parameters
        if (key == "tree") continue;

        //std::cout << "Processing key: " << key << " | Type: " << it->type_name() 
        //          << " | Expecting " << (expecting_value ? "Value" : "Key") << std::endl;

        // If expecting a key (even depth), recurse deeper with key as parameter name
        if (!expecting_value) {
            extract_parameters(*it, parameters, key, true, int_sets, double_sets, string_sets);
        }
        else { // Expecting a VALUE
            if (!current_key) continue;  // Skip if no parameter name is set

            // Determine if the key is a number
            bool is_number = key.find_first_not_of("0123456789.") == std::string::npos;
            bool is_float = key.find('.') != std::string::npos;

            // Ensure the key exists in the correct set for uniqueness tracking
            if (is_number) {
                if (is_float) {
                    double_sets[*current_key].insert(std::stod(key));
                } else {
                    int_sets[*current_key].insert(std::stoi(key));
                }
            } else {
                string_sets[*current_key].insert(key);
            }

            // ✅ Debug: Show what was added
            //std::cout << "✔ Inserted unique value into parameters[" << *current_key << "]: " << key << std::endl;

            // Now switch to expecting a key on next recursion
            extract_parameters(*it, parameters, current_key, false, int_sets, double_sets, string_sets);
        }
    }

    // Convert sets to vectors and store in `parameters`
    for (const auto& [key, values] : int_sets) {
        parameters[key] = std::vector<int>(values.begin(), values.end());
    }
    for (const auto& [key, values] : double_sets) {
        parameters[key] = std::vector<double>(values.begin(), values.end());
    }
    for (const auto& [key, values] : string_sets) {
        parameters[key] = std::vector<std::string>(values.begin(), values.end());
    }
}



const json& AdaptiveGaussTreeBatch::find_tree_json(const json& param_json, const ParamMap& param_map) {
    const json* current = &param_json;
    for (const auto& key : keys) {
        auto it = param_map.find(key);
        if (it == param_map.end()) continue; // Key not found in param_map, move on

        std::string key_value;
        if (std::holds_alternative<int>(it->second)) {
            key_value = std::to_string(std::get<int>(it->second));
        } else if (std::holds_alternative<double>(it->second)) {
            key_value = std::to_string(std::get<double>(it->second));
        } else if (std::holds_alternative<std::string>(it->second)) {
            key_value = std::get<std::string>(it->second);
        }

        if (!current->contains(key)) {
            throw std::runtime_error("❌ Key missing in JSON: " + key);
        }
        if (!(*current)[key].contains(key_value)) {
            throw std::runtime_error("❌ Key value missing in JSON: " + key_value);
        }

        current = &((*current)[key][key_value]);
    }

    // 🚨 Debug before returning "tree"
//    std::cout << "🔍 Final JSON before accessing 'tree':\n" << current->dump(4) << std::endl;

    // ✅ Check if "tree" exists and is valid
    if (!current->contains("tree")) {
        throw std::runtime_error("❌ 'tree' key missing in JSON");
    }
    if ((*current)["tree"].is_null()) {
        throw std::runtime_error("❌ 'tree' key is null in JSON");
    }

    return (*current)["tree"];
}



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

void AdaptiveGaussTreeBatch::save_to_json(const std::string & filename, bool overwrite, bool write_roots, bool write_trees ) {
    json data;

    data["name"] = name;
    data["author"] =author;
    data["version"] =version;        
    data["reference"] = reference;
    data["description"] = description;
    data["tol"] = tol;
    data["min_depth"] = min_depth;
    data["max_depth"] = max_depth;
    data["n1"] = order1;
    data["n2"] = order2;
    data["a_singular"] = a_singular ;
    data["b_singular"] = b_singular;
    data["write_trees"] = write_trees;
    // Serialize update log
    json log_json = json::array();
    for (const auto& entry : update_log) {
        log_json.push_back({{"timestamp", entry.first}, {"message", entry.second}});
    }
    data["update_log"] = log_json;
    // Serialize weights if requested.
    if (write_roots) {
        data["legendre_roots_n1"] = {legendre_n1.getNodes(order1), legendre_n1.getWeights(order1)};
        data["laguerre_roots_n1"] = {laguerre_n1.getNodes(order1), laguerre_n1.getWeights(order1)};
        data["legendre_roots_n2"] = {legendre_n2.getNodes(order2), legendre_n2.getWeights(order2)};
        data["laguerre_roots_n2"] = {laguerre_n2.getNodes(order2), laguerre_n2.getWeights(order2)};
    }         
    data["parameters"] = parameter_serializer(write_trees);
    std::ofstream file(filename);
    file << data.dump(4);        
} 

json AdaptiveGaussTreeBatch::parameter_serializer(bool dump_nodes) {
    json result;

    for (const auto& [param_map, tree_ptr] : quad_coll) {
        json* current = &result;  // Pointer to navigate the JSON structure

        for (const auto& key : keys) {
            auto it = param_map.find(key);
            if (it == param_map.end()) continue; // Skip if key not found

            // Convert variant value to a string key (for hierarchy)
            std::string key_value;
            if (std::holds_alternative<int>(it->second)) {
                key_value = std::to_string(std::get<int>(it->second));
            } else if (std::holds_alternative<double>(it->second)) {
                key_value = std::to_string(std::get<double>(it->second));
            } else if (std::holds_alternative<std::string>(it->second)) {
                key_value = std::get<std::string>(it->second);
            }

            // Ensure parameter label exists in JSON
            current = &((*current)[key]);  // Create or navigate label (e.g., "s", "z")
            current = &((*current)[key_value]);  // Create/navigate value under label
        }

        // Assign the serialized tree to the final position
        (*current)["tree"] = tree_ptr->get_tree_serialized(dump_nodes);
    }

    return result;
}

AdaptiveGaussTreeBatch AdaptiveGaussTreeBatch::operator+(const AdaptiveGaussTreeBatch& other) const {
    AdaptiveGaussTreeBatch result(*this);  // Step 1: Deep copy lhs (`this`)
    result.merge(other);                   // Step 2: Merge rhs (`other`)
    return result;                          // Step 3: Return new merged batch
}

AdaptiveGaussTreeBatch& AdaptiveGaussTreeBatch::operator+=(const AdaptiveGaussTreeBatch& other) {
    this->merge(other);  // Merge `other` into `this`
    return *this;        // Return the updated lhs (`this`)
}

