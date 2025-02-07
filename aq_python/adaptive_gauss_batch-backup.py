import json
import os
import sys
import numpy as np
import pandas as pd
from itertools import product
from aq_python.adaptive_quadrature import AdaptiveGaussTree

class AdaptiveGaussTreeBatch:
    def __init__(self, func, parameters=None, file=None, tol=1e-6, min_depth=2, max_depth=10, n1=5, n2=10, 
                 a_singular=False, b_singular=False, name="Batch Quadrature", author=None, version="1.0", 
                 reference=None, description=None, update_log_message='Initial Batch Creation'):
        """
        Iterates through parameter sets and creates an AdaptiveGaussTree for each, or loads from a file if specified.
        """
        self.func = func
        self.update_log = []  
        
        if file:
            write_trees = self._load_from_json(file)
            if not write_trees:
                print(f"Warning: Tree structure is not available in {file}.", file=sys.stderr)
            if not self.results:
                print(f"Warning: Results not available in {file}.", file=sys.stderr)
            self.add_update_log(update_log_message)
        else:
            self.parameters = {key: list(value) if isinstance(value, (range, list)) else [value] for key, value in parameters.items()}
            self.param_combinations = list(product(*self.parameters.values()))
            self.param_keys = list(self.parameters.keys())
            
            self.tol = tol
            self.min_depth = min_depth
            self.max_depth = max_depth
            self.n1 = n1
            self.n2 = n2
            self.a_singular = a_singular
            self.b_singular = b_singular
            
            self.name = name
            self.author = author
            self.version = version
            self.reference = reference
            self.description = description
            
            self.add_update_log(update_log_message)
            self.results = self._generate_trees()
    
    def _generate_trees(self):
        """Generates and stores AdaptiveGaussTree results in a hierarchical dictionary."""
        tree_data = {}
        
        for param_set in self.param_combinations:
            param_dict = dict(zip(self.param_keys, param_set))
            nested_tree = tree_data
            
            for key, value in param_dict.items():
                if key not in nested_tree:
                    nested_tree[key] = {}
                if value not in nested_tree[key]:
                    nested_tree[key][value] = {}
                nested_tree = nested_tree[key][value]
            
            quad_tree = AdaptiveGaussTree(
                lambda t: self.func(param_dict, t), 0, 1, tol=self.tol, 
                max_depth=self.max_depth, n1=self.n1, n2=self.n2, a_singular=self.a_singular, min_depth=self.min_depth
            )
            
            nested_tree['tree'] = quad_tree._serialize_tree(quad_tree.root)
        
        return tree_data
    
    def get_values(self):
        """Returns a copy of self.results with 'left' and 'right' removed from each tree node."""
        import copy
        local_results = copy.deepcopy(self.results)
        
        def prune_tree(tree):
            if isinstance(tree, dict):
                if "tree" in tree and isinstance(tree["tree"], dict):
                    pruned_tree = tree["tree"].copy()
                    pruned_tree.pop("left", None)
                    pruned_tree.pop("right", None)
                    tree["tree"] = pruned_tree
                for key in tree:
                    prune_tree(tree[key])
        
        prune_tree(local_results)
        return local_results
    
    def save_to_json(self, filename, prevent_overwrite=True, write_trees=True):
        """Saves the hierarchical parameterized tree structure to a JSON file, preventing overwrite if specified."""
        if prevent_overwrite and os.path.exists(filename):
            print(f"Error: {filename} exists. If you wish to overwrite the file, set prevent_overwrite=False")
            return
        
        data = {
            "name": self.name,
            "author": self.author,
            "version": self.version,
            "reference": self.reference,
            "description": self.description,
            "tol": self.tol,
            "min_depth": self.min_depth,
            "max_depth": self.max_depth,
            "n1": self.n1,
            "n2": self.n2,
            "a_singular": self.a_singular,
            "b_singular": self.b_singular,
            "write_trees": write_trees,
            "update_log": self.update_log,
            "parameters": self.results if write_trees else self.get_values()
        }
        with open(filename, "w") as f:
            json.dump(data, f, indent=4)
    
    def _load_from_json(self, filename) -> bool:
        """Loads AdaptiveGaussTreeBatch from a JSON file."""
        with open(filename, "r") as f:
            data = json.load(f)
        
        self.name = data.get("name", "Batch Quadrature")
        self.author = data.get("author")
        self.version = data.get("version", "1.0")
        self.reference = data.get("reference")
        self.description = data.get("description")
        self.tol = data.get("tol", 1e-6)
        self.min_depth = data.get("min_depth", 2)
        self.max_depth = data.get("max_depth", 10)
        self.n1 = data.get("n1", 5)
        self.n2 = data.get("n2", 10)
        self.a_singular = data.get("a_singular", False)
        self.b_singular = data.get("b_singular", False)
        self.update_log = data.get("update_log", [])
        write_trees = data.get("write_trees", True)
        self.results = data.get("parameters", {})
        return write_trees
    
    def print_trees(self, nested_dict=None, level=0):
        """Recursively prints the hierarchical parameter tree structure."""
        if nested_dict is None:
            nested_dict = self.results
        
        indent = "  " * level
        for key, value in nested_dict.items():
            if isinstance(value, dict):
                print(f"{indent}{key}: ")
                self.print_trees(value, level + 1)
            else:
                print(f"{indent}- {key}: {value}")
    
    def add_update_log(self, message):
        """Adds a timestamped message to the update log."""
        import datetime
        self.update_log.append({
            "timestamp": datetime.datetime.now().isoformat(),
            "message": message
        })
    
    def to_dataframe(self):
        """Converts the pruned results into a pandas DataFrame."""
        pruned_results = self.get_values()
        
        def flatten_results(d, params=None):
            if params is None:
                params = {}
            
            for key, value in d.items():
                if key == "tree" and isinstance(value, dict):
                    yield {**params, **value}
                elif isinstance(value, dict):
                    for sub_key, sub_value in value.items():
                        yield from flatten_results(sub_value, {**params, key: sub_key})
        
        data = list(flatten_results(pruned_results))
        return pd.DataFrame(data)
