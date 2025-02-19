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

using json = nlohmann::ordered_json;
//using ParamCollection = std::unordered_map<std::string, std::variant<std::vector<int>, std::vector<double>, std::vector<std::string>>>;

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
        generate_combinations(keys, parameters, indices, results);
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

    void printKeys() {
        for (const auto& [key, _] : quad_coll) {
            std::cout << key << std::endl;
        }
     };
    
    void printResults() {
           for (const auto& [key, val] : quad_coll) {
               std::cout << key << *val<<std::endl;
           }
       };


/*
    void save_to_json(const std::string& filename) {
        json data;
        data["update_log"] = update_log;
        for (const auto& [key, val_map] : results) {
            for (const auto& [val, tree] : val_map) {
                data["parameters"][key][val] = tree->get_integral_and_error();
            }
        }
        std::ofstream file(filename);
        file << data.dump(4);
    }
*/
/*
    void load_from_json(const std::string& filename, WeightsLoader rl1, WeightsLoader rl2, WeightsLoader ll1, WeightsLoader ll2) {
        std::ifstream file(filename);
        json data;
        file >> data;
        update_log = data["update_log"].get<std::vector<std::pair<std::string, std::string>>>();
        for (auto& [key, val_map] : data["parameters"].items()) {
            for (auto& [val, _] : val_map.items()) {
                ParamMap param_map = {{key, val}};
                results[key][val] = std::make_unique<AdaptiveGaussTree>(func, rl1, rl2, ll1, ll2, filename, param_map);
            }
        }
    }
*/

};



#endif // ADAPTIVE_GAUSS_BATCH_HPP
