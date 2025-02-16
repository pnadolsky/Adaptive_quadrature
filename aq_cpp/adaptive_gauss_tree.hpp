#ifndef ADAPTIVE_GAUSS_TREE_HPP
#define ADAPTIVE_GAUSS_TREE_HPP

#include "quadrature.hpp"
#include "legendre_quadrature.hpp"
#include "laguerre_singular_endpoint.hpp"
#include "weights_loader.hpp"
#include "json.hpp"
#include <iostream>
#include <fstream>
#include <memory>
#include <unordered_map>
#include <variant>
#include <optional>

using json = nlohmann::json;

class AdaptiveGaussTree {
public:
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

private:
    std::function<double(std::unordered_map<std::string, std::variant<int, double, std::string>>, double)> func;
    double tolerance;
    int min_depth, max_depth;
    int order1, order2;
    bool a_singular, b_singular;
    double alpha_a, alpha_b;
    
    WeightsLoader roots_legendre_n1, roots_laguerre_n1, roots_legendre_n2, roots_laguerre_n2;
    std::unique_ptr<Node> root;

public:
    // Constructor from parameters
    AdaptiveGaussTree(
        std::function<double(std::unordered_map<std::string, std::variant<int, double, std::string>>, double)> f,
        double lower, double upper, double tol, int minD, int maxD,
        int n1, int n2,  // Explicitly pass n1 and n2
        double alphaA, double alphaB, bool singularA, bool singularB,
        WeightsLoader rl1, WeightsLoader rl2, WeightsLoader ll1, WeightsLoader ll2)
        : func(f), tolerance(tol), min_depth(minD), max_depth(maxD),
          order1(n1), order2(n2),  // Explicitly set orders
          alpha_a(alphaA), alpha_b(alphaB), a_singular(singularA), b_singular(singularB),
          roots_legendre_n1(rl1), roots_legendre_n2(rl2),
          roots_laguerre_n1(ll1), roots_laguerre_n2(ll2) {
        
        root = build_tree(lower, upper, 0, tol);
    }
    
    
    // Constructor from JSON file
    AdaptiveGaussTree(
        std::function<double(std::unordered_map<std::string, std::variant<int, double, std::string>>, double)> f,
        WeightsLoader rl1, WeightsLoader rl2, WeightsLoader ll1, WeightsLoader ll2, std::string filename)
        : func(f), roots_legendre_n1(rl1), roots_legendre_n2(rl2),
          roots_laguerre_n1(ll1), roots_laguerre_n2(ll2) {
        load_from_json(filename);
    }
    // Method to get the total integral and error
    std::pair<double, double> get_integral_and_error() const {
        return traverse_and_sum(root.get());
    }
private:
    std::unique_ptr<Node> build_tree(double lower, double upper, int depth, double tol) {
        bool use_laguerre = (lower == 0 && a_singular) || (upper == 1 && b_singular);
        std::unique_ptr<Quadrature> quadrature;
        if (use_laguerre) {
            quadrature = std::make_unique<LaguerreSingularEndpoint>(roots_laguerre_n1, order1, order2, lower, upper);
        } else {
            quadrature = std::make_unique<LegendreQuadrature>(roots_legendre_n1, order1, order2, lower, upper);
        }
                
        double I1 = quadrature->integrate(func, {});
        double I2 = quadrature->integrate(func, {});
        double err = std::abs(I2 - I1);
        
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

public:
    void save_to_json(std::string filename) {
        json data;
        data["tolerance"] = tolerance;
        data["min_depth"] = min_depth;
        data["max_depth"] = max_depth;
        data["n1"] = order1;
        data["n2"] = order2;
        data["tree"] = serialize_tree(root.get());
        std::ofstream file(filename);
        file << data.dump(4);
    }

    void load_from_json(std::string filename) {
        std::ifstream file(filename);
        json data;
        file >> data;
        tolerance = data["tolerance"];
        min_depth = data["min_depth"];
        max_depth = data["max_depth"];
        order1 = data["n1"];
        order2 = data["n2"];
        root = deserialize_tree(data["tree"]);
    }

private:
    json serialize_tree(Node* node) {
        if (!node) return nullptr;
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
        node->left = deserialize_tree(data["left"]);
        node->right = deserialize_tree(data["right"]);
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

#endif // ADAPTIVE_GAUSS_TREE_HPP
