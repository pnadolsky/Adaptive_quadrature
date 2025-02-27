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
#include <set>


using json = nlohmann::ordered_json;
// using ParamMap = std::unordered_map<std::string, ParamType>
// using ParamCollection = std::map<std::string, std::variant<std::vector<int>, std::vector<double>, std::vector<std::string>>>;
// using QuadCollection = std::unordered_map<ParamMap, std::unique_ptr<AdaptiveGaussTree>, ParamMapHash, ParamMapEqual> 
class AdaptiveGaussTreeBatch {
private:
    QuadCollection quad_coll;
    std::function<double(ParamMap, double)> func;
    double tol; 
    double lower, upper;
    double alphaA, alphaB;
    int min_depth; int max_depth; int order1;int order2;
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
    void extract_keys(const json& param_json, std::set<std::string>& key_set, int depth = 0) ;
    void extract_parameters(
        const json& param_json, 
        ParamCollection& parameters, 
        std::optional<std::string> current_key = std::nullopt, 
        bool expecting_value = false, 
        std::map<std::string, std::set<int>> int_sets = {}, 
        std::map<std::string, std::set<double>> double_sets = {}, 
        std::map<std::string, std::set<std::string>> string_sets = {})  ;

    const json& find_tree_json(const json& param_json, const ParamMap& param_map);     

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
        min_depth(min_depth), max_depth(max_depth), order1(n1), order2(n2),
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
        // printKeys_internal(); 
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

    AdaptiveGaussTreeBatch(
        std::function<double(ParamMap, double)> func,
        std::string filename
    );    

    void printCollection() {
           for (const auto& result : results) {
               std::cout << result << *quad_coll[result] <<std::endl;
           }
       };   

    const QuadCollection& getCollection() const { return quad_coll; };
    void save_to_json(const std::string & filename, bool overwrite = false, bool write_roots=false, bool write_trees =true); 
    json parameter_serializer(bool dump_nodes = false); 
};



#endif // ADAPTIVE_GAUSS_BATCH_HPP
