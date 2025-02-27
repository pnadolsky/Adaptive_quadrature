#ifndef ADAPTIVE_GAUSS_TREE_HPP
#define ADAPTIVE_GAUSS_TREE_HPP

#include "quadrature.hpp"
#include "legendre_quadrature.hpp"
#include "laguerre_singular_endpoint.hpp"
#include "weights_loader.hpp"
#include "json.hpp"
#include <iostream>
#include <fstream>
#include <filesystem> // For checking file existence
#include <memory>
#include <unordered_map>
#include <variant>
#include <optional>
#include <vector>
#include <utility>
#include <ctime>

using json = nlohmann::ordered_json;

class AdaptiveGaussTree {
    private:
    struct Node {
        double lower, upper;
        int depth;
        double tolerance, error, result;
        int order1, order2;
        bool is_singular;
        std::unique_ptr<Node> left, right;
        
        Node(double lower, double upper, int depth, double tol, int o1, int o2, bool singular)
            : lower(lower), upper(upper), depth(depth), tolerance(tol), order1(o1), order2(o2), is_singular(singular),
              error(0.0), result(0.0), left(nullptr), right(nullptr) {}
    };


    public:
    // Constructor from parameters
    AdaptiveGaussTree(
        std::function<double(ParamMap, double)> f,
        double lower, double upper, double tol, int minD, int maxD,
        int n1, int n2,  // Explicitly pass n1 and n2
        double alphaA, double alphaB, bool singularA, bool singularB,
        WeightsLoader rl1, WeightsLoader rl2, WeightsLoader ll1, WeightsLoader ll2,
        ParamMap args={},  // function argumnets (optional)
        std::string name="Project", std::string author="Author",  std::string description="project description", 
        std::string reference="references", std::string version="1.0", std::string update_log_message="Initial Train" 
    )
        : func(f), tolerance(tol), min_depth(minD), max_depth(maxD),
          order1(n1), order2(n2),  // Explicitly set orders
          alpha_a(alphaA), alpha_b(alphaB), a_singular(singularA), b_singular(singularB),
          roots_legendre_n1(rl1), roots_legendre_n2(rl2),
          roots_laguerre_n1(ll1), roots_laguerre_n2(ll2),
          args(args),
          name(name), author(author), description(description), reference(reference), version(version) {
        
        root = build_tree(lower, upper, 0, tol);
        add_update_log(update_log_message);
    }
        
    // Constructor from JSON file
    AdaptiveGaussTree(
        std::function<double(ParamMap, double)> f,
        WeightsLoader rl1, WeightsLoader rl2, WeightsLoader ll1, WeightsLoader ll2, std::string filename,ParamMap args={})
        : func(f), roots_legendre_n1(rl1), roots_legendre_n2(rl2),
          roots_laguerre_n1(ll1), roots_laguerre_n2(ll2), args(args) {
        load_from_json(filename);
    }

    // Constructor from JSON data
    AdaptiveGaussTree(json jsn,
        std::function<double(ParamMap, double)> f,
        WeightsLoader rl1, WeightsLoader rl2, WeightsLoader ll1, WeightsLoader ll2, ParamMap args={} )
        : func(f), roots_legendre_n1(rl1), roots_legendre_n2(rl2),
          roots_laguerre_n1(ll1), roots_laguerre_n2(ll2), args(args) {
        load_from_json_stream(jsn);
    }

// copy constructor

    AdaptiveGaussTree(const AdaptiveGaussTree& other)
        : func(other.func),
        tolerance(other.tolerance),
        min_depth(other.min_depth), max_depth(other.max_depth),
        order1(other.order1), order2(other.order2),
         a_singular(other.a_singular), b_singular(other.b_singular),
        alpha_a(other.alpha_a), alpha_b(other.alpha_b),
        args(other.args),
        roots_legendre_n1(other.roots_legendre_n1), roots_legendre_n2(other.roots_legendre_n2),
        roots_laguerre_n1(other.roots_laguerre_n1), roots_laguerre_n2(other.roots_laguerre_n2),
        name(other.name), reference(other.reference), description(other.description),
        author(other.author), version(other.version),
        update_log(other.update_log) {

    // Deep copy the tree structure
        root = clone_tree(other.root.get());
}

// Helper function to recursively copy tree nodes
    std::unique_ptr<Node> clone_tree(const Node* node) {
        if (!node) return nullptr;

        auto new_node = std::make_unique<Node>(node->lower, node->upper, node->depth, 
                                            node->tolerance, node->order1, node->order2, 
                                            node->is_singular);
        new_node->result = node->result;
        new_node->error = node->error;

        new_node->left = clone_tree(node->left.get());
        new_node->right = clone_tree(node->right.get());

        return new_node;
    }

    // Method to get the total integral and error
    std::pair<double, double> get_integral_and_error() const {
        return traverse_and_sum(root.get());
    }

    void add_update_log(const std::string& message) {
        // Get current time
        std::time_t now = std::time(nullptr);
        char timestamp[20];
        std::strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
    
        // Append timestamped message
        update_log.emplace_back(timestamp, message);
    }

    void save_to_json(std::string filename, bool overwrite = false, bool dump_log = false) {
        // Check if the file exists
        if (std::filesystem::exists(filename) && !overwrite) {
            std::cerr << "File \"" << filename << "\" exists. Set overwrite = true to overwrite." << std::endl;
            return;
        }
        json data;
        // name(name), author(author), description(description), reference(reference), version(version)        
        data["name"] = name;
        data["reference"] = reference;
        data["description"] = description;
        data["author"] =author;
        data["version"] = version;
        data["tolerance"] = tolerance;
        data["min_depth"] = min_depth;
        data["max_depth"] = max_depth;
        data["n1"] = order1;
        data["n2"] = order2;
        // Serialize update log
        json log_json = json::array();
        for (const auto& entry : update_log) {
            log_json.push_back({{"timestamp", entry.first}, {"message", entry.second}});
        }
        data["update_log"] = log_json;        
        data["tree"] = serialize_tree(root.get());
        std::ofstream file(filename);
        file << data.dump(4);
    }

    void load_from_json(std::string filename) {
        json data;
        std::ifstream file(filename);
        file >> data;
        load_from_json_stream(data);                
    }

    void load_from_json_stream(json data) {
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
        // Deserialize update log
        update_log.clear();
        if (data.contains("update_log")) {
            for (const auto& entry : data["update_log"]) {
                update_log.emplace_back(entry["timestamp"], entry["message"]);
            }
        }
//        std::cout << "deserialize" << std::endl;
        root = deserialize_tree(data["tree"]);
    }

    void print_update_log() const {
        for (const auto& entry : update_log) {
            std::cout << "[" << entry.first << "] " << entry.second << std::endl;
        }
    }
    json get_tree_serialized(bool dump_nodes = false){
        return serialize_tree(root.get(), dump_nodes);
    }
private:


    std::function<double(ParamMap, double)> func;
    double tolerance;
    int min_depth, max_depth;
    int order1, order2;
    bool a_singular, b_singular;
    double alpha_a, alpha_b;
    ParamMap args;
    
    WeightsLoader roots_legendre_n1, roots_laguerre_n1, roots_legendre_n2, roots_laguerre_n2;
    std::unique_ptr<Node> root;

    // json header info
    std::string name;
    std::string reference;
    std::string description;
    std::string author;
    std::string version;
    std::vector<std::pair<std::string, std::string>> update_log;

    std::unique_ptr<Node> build_tree(double lower, double upper, int depth, double tol) {
        bool use_laguerre = (lower == 0 && a_singular) || (upper == 1 && b_singular);
        std::unique_ptr<Quadrature> quadrature;
        if (use_laguerre) {
            quadrature = std::make_unique<LaguerreSingularEndpoint>(roots_laguerre_n1, order1, order2, lower, upper);
        } else {
            quadrature = std::make_unique<LegendreQuadrature>(roots_legendre_n1, order1, order2, lower, upper);
        }
                
        double I2 = quadrature->integrate(func, args);
 //       double I1 = quadrature->integrate(func, {});  
 //       double err = I2-I1  // ChatGPT needs a vacay.
        double err = quadrature->getError();
        
        auto node = std::make_unique<Node>(lower, upper, depth, tol, order1, order2, use_laguerre);
        node->result = I2;
        node->error = err;
        
        if (depth < min_depth || (err >= tol && depth < max_depth)) {
            double mid = (lower + upper) / 2;
            node->left = build_tree(lower, mid, depth + 1, tol / 2);
            node->right = build_tree(mid, upper, depth + 1, tol / 2);
        }
        return node;
    }

    json serialize_tree(Node* node, bool dump_nodes = false) {
        if (!node) return nullptr;
        if (dump_nodes) return {
            {"a", node->lower},
            {"b", node->upper},
            {"depth", node->depth},
            {"tol", node->tolerance},
            {"error", node->error},
            {"integral", node->result},
            {"method", node->is_singular ? "Gauss-Laguerre" : "Gauss-Legendre"}            
        };
        return {
            {"a", node->lower},
            {"b", node->upper},
            {"depth", node->depth},
            {"tol", node->tolerance},
            {"error", node->error},
            {"integral", node->result},
            {"method", node->is_singular ? "Gauss-Laguerre" : "Gauss-Legendre"},
            {"left", serialize_tree(node->left.get())},
            {"right", serialize_tree(node->right.get())}
        };
    }
// Node(double lower, double upper, int depth, double tol, int o1, int o2, bool singular)    
    std::unique_ptr<Node> deserialize_tree(const json& data) {
        if (data.is_null()) return nullptr;
        auto node = std::make_unique<Node>(data["a"], data["b"], data["depth"],
                                           data["tol"],order1, order2, data["method"] == "Gauss-Laguerre");
        node->error = data["error"];
        node->result = data["integral"];
        if (data.contains("left")){
            node->left = deserialize_tree(data["left"]);
            node->right = deserialize_tree(data["right"]);
        } 
        return node;
    }

    std::pair<double, double> traverse_and_sum(const Node* node) const {
        if (!node) return {0.0, 0.0};
        if (!node->left && !node->right) return {node->result, node->error};
        
        auto [left_integral, left_error] = traverse_and_sum(node->left.get());
        auto [right_integral, right_error] = traverse_and_sum(node->right.get());
        
        return {left_integral + right_integral, left_error + right_error};
    }    
    
};

using QuadCollection = std::unordered_map<ParamMap, std::unique_ptr<AdaptiveGaussTree>, ParamMapHash, ParamMapEqual>;


//overload to print AdaptiveGaussTree:  (integral, error)

std::ostream& operator<<(std::ostream& os, const AdaptiveGaussTree& tree);

#endif // ADAPTIVE_GAUSS_TREE_HPP
