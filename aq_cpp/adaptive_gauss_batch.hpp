#ifndef ADAPTIVE_GAUSS_BATCH_HPP
#define ADAPTIVE_GAUSS_BATCH_HPP

#include "weights_loader.hpp"
#include "adaptive_gauss_tree.hpp"
#include "json.hpp"
#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <variant>
#include <functional>
#include <algorithm>

using json = nlohmann::ordered_json;
// using ParamCollection = std::map<std::string, std::variant<std::vector<int>, std::vector<double>, std::vector<std::string>>>;
// using  QuadCollection = std::unordered_map<ParamMap, std::unique_ptr<AdaptiveGaussTree>, ParamMapHash, ParamMapEqual> 
class AdaptiveGaussTreeBatch {
private:
    QuadCollection quad_coll;
    std::function<double(ParamMap, double)> func;
    double tol; 
    double lower, upper;
    double alphaA, alphaB;
    int min_depth; int max_depth; int n1;int n2;
    bool a_singular; bool b_singular;
    WeightsLoader legendre_n1; WeightsLoader legendre_n2; WeightsLoader laguerre_n1; WeightsLoader laguerre_n2;
    ParamCollection parameters;
    std::string name; std::string author;  std::string description; 
    std::string reference; std::string version;         
    std::vector<ParamMap> results;
    std::vector<std::pair<std::string, std::string>> update_log;
    std::vector<std::string> keys;

    void generate_combinations(
        const std::vector<std::string>& keys,
        const ParamCollection& params,
        std::vector<size_t>& indices,
        std::vector<ParamMap>& results,
        size_t depth = 0
    );

    void add_update_log(const std::string& message) {
        std::time_t now = std::time(nullptr);
        char timestamp[20];
        std::strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
        update_log.emplace_back(timestamp, message);
    }

    bool compareParamMaps(const ParamMap& a, const ParamMap& b);    // Comparator for sorting based on "keys"
    void sortResults();                                             // Sorting function
    static bool compareVariant(const ParamType& a, const ParamType& b);    // Function to compare std::variant<int, double, std::string>

    // Helper function to rank types for sorting (int < double < string)    
    template <typename T> static int getTypeRank();

public:
    AdaptiveGaussTreeBatch(
        std::function<double(ParamMap, double)> func,
        double lower, double upper,        
        double tol, int min_depth, int max_depth, int n1, int n2,
        double alphaA, double alphaB,
        bool a_singular, bool b_singular,
        WeightsLoader legendre_n1, WeightsLoader legendre_n2, WeightsLoader laguerre_n1, WeightsLoader laguerre_n2,
        ParamCollection parameters,
        std::string name="Project", std::string author="Author",  std::string description="project description", 
        std::string reference="references", std::string version="1.0", std::string update_log_message="Initial Batch Creation"         
    ) : func(func), 
        lower(lower),upper(upper), tol(tol),
        min_depth(min_depth), max_depth(max_depth), n1(n1), n2(n2),
        alphaA(alphaA), alphaB(alphaB),
        a_singular(a_singular),b_singular(b_singular), 
        legendre_n1(legendre_n1),legendre_n2(legendre_n2),laguerre_n1(laguerre_n1),laguerre_n2(laguerre_n2),
        parameters(parameters), 
        name(name), author(author),description(description),reference(reference), version(version)
    {
    
        add_update_log(update_log_message);
        for (const auto& pair : parameters) {
            keys.push_back(pair.first);
        }
        std::vector<size_t> indices(keys.size(), 0);
        // printKeys_internal(); 
        generate_combinations(keys, parameters, indices, results);
        sortResults();
//        printKeys_internal(); 
        for (const auto& combo : results) {
            quad_coll[combo]= std::make_unique<AdaptiveGaussTree>(
                func, lower, upper, tol, min_depth, max_depth, n1, n2,  
                alphaA, alphaB, a_singular, b_singular,
                legendre_n1,  legendre_n2, laguerre_n1, laguerre_n2,
                combo,
                name, author, description, 
                reference, version, update_log_message            
            );                 
        }
    }


    void printCollection() {
           for (const auto& result : results) {
               std::cout << result << *quad_coll[result] <<std::endl;
           }
       };   

    const QuadCollection& getCollection() const { return quad_coll; };

    void save_to_json(const std::string & filename, bool overwrite = false, bool write_roots=false, bool write_trees =true) {
        json data;
    
        data["name"] = name;
        data["author"] =author;
        data["version"] =version;        
        data["reference"] = reference;
        data["description"] = description;
        data["tol"] = tol;
        data["min_depth"] = min_depth;
        data["max_depth"] = max_depth;
        data["n1"] = n1;
        data["n2"] = n2;
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
            data["legendre_roots_n1"] = {legendre_n1.getNodes(n1), legendre_n1.getWeights(n1)};
            data["laguerre_roots_n1"] = {laguerre_n1.getNodes(n1), laguerre_n1.getWeights(n1)};
            data["legendre_roots_n2"] = {legendre_n2.getNodes(n2), legendre_n2.getWeights(n2)};
            data["laguerre_roots_n2"] = {laguerre_n2.getNodes(n2), laguerre_n2.getWeights(n2)};
        }         
        data["parameters"] = parameter_serializer(write_trees);
        std::ofstream file(filename);
        file << data.dump(4);        
    }; 



    json parameter_serializer(bool dump_nodes = false) {
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
    

};



#endif // ADAPTIVE_GAUSS_BATCH_HPP
