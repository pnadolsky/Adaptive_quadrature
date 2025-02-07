import json
import datetime
import numpy as np
from scipy.special import roots_legendre, roots_laguerre
from math import factorial

class AdaptiveGaussTree:
    class Node:
        """Represents a node in the adaptive quadrature tree."""
        def __init__(self, a, b, depth, tol, error=0.0, integral=None, n1=None, n2=None, is_singular=False):
            self.a = a
            self.b = b
            self.depth = depth
            self.tol = tol
            self.error = error
            self.n1 = n1
            self.n2 = n2
            self.is_singular = is_singular
            self.integral = integral
            self.left = None
            self.right = None

    def __init__(self, f=None, a=None, b=None, tol=1e-6, min_depth=2, max_depth=10, n1=5, n2=10, 
                 alpha_a=0, alpha_b=0, a_singular=False, b_singular=False,
                 roots_legendre_n1=None, roots_laguerre_n1=None, roots_legendre_n2=None, roots_laguerre_n2=None,
                 name="Adaptive Quadrature Tree", reference=None, description=None, 
                 author=None, version="1.0", update_log_message='Initial Train', filename=None):
        if filename:
            self._load_from_json(filename)
        else:
            self.f = f
            self.a = a
            self.b = b
            self.tol = tol
            self.min_depth = min_depth
            self.max_depth = max_depth
            
            # Determine n1 and n2 from provided quadrature roots
            self.n1 = len(roots_legendre_n1[0]) if roots_legendre_n1 else n1
            self.n2 = len(roots_legendre_n2[0]) if roots_legendre_n2 else n2
            
            self.alpha_a = alpha_a
            self.alpha_b = alpha_b
            self.a_singular = a_singular
            self.b_singular = b_singular

            # Store quadrature roots for efficiency
            self.legendre_roots_n1 = roots_legendre(self.n1) if roots_legendre_n1 is None else roots_legendre_n1
            self.laguerre_roots_n1 = roots_laguerre(self.n1) if roots_laguerre_n1 is None else roots_laguerre_n1
            self.legendre_roots_n2 = roots_legendre(self.n2) if roots_legendre_n2 is None else roots_legendre_n2
            self.laguerre_roots_n2 = roots_laguerre(self.n2) if roots_laguerre_n2 is None else roots_laguerre_n2

            # Metadata
            self.name = name
            self.reference = reference
            self.description = description
            self.author = author
            self.version = version
            self.update_log = []  
            self.add_update_log(update_log_message)
            
            # Build tree
            self.root = self._build_tree(a, b, tol, depth=0)

    def _build_tree(self, a, b, tol, depth):
        """Recursively builds the quadrature tree."""
        use_laguerre_a = (a == self.a and self.a_singular)
        use_laguerre_b = (b == self.b and self.b_singular)

        if use_laguerre_a:
            integrate = lambda f, a, b, n: self._gauss_laguerre_integrate(f, self.a, b - self.a, self.alpha_a, n)
        elif use_laguerre_b:
            integrate = lambda f, a, b, n: self._gauss_laguerre_integrate(f, self.b, self.b - a, self.alpha_b, n)
        else:
            integrate = self._gauss_legendre_integrate

        I1 = integrate(self.f, a, b, self.n1)
        I2 = integrate(self.f, a, b, self.n2)
        error = abs(I2 - I1)

        node = self.Node(a, b, depth, tol, error, integral=I2, n1=self.n1, n2=self.n2, is_singular=(use_laguerre_a or use_laguerre_b))

        if depth < self.min_depth or (error >= tol and depth < self.max_depth):
            mid = (a + b) / 2
            node.left = self._build_tree(a, mid, tol / 2, depth + 1)
            node.right = self._build_tree(mid, b, tol / 2, depth + 1)

        return node

    def _gauss_legendre_integrate(self, f, a, b, deepset):
        """Performs Gauss-Legendre quadrature."""
        x, w = self.legendre_roots_n2 if deepset  else self.legendre_roots_n1 
        x_mapped = 0.5 * (b - a) * x + 0.5 * (a + b)
        return 0.5 * (b - a) * np.sum(w * f(x_mapped))

    def _gauss_laguerre_integrate(self, f, sing_value, epsilon, alpha, deepset):
        """Performs Gauss-Laguerre quadrature for singularities."""
        y, w = self.laguerre_roots_n2 if deepset  else self.laguerre_roots_n1 
        def transformed_f(y):
            x = sing_value + epsilon * np.exp(-y)
            return f(x) / np.exp(-y * alpha)
        return epsilon**(1 - alpha) * np.sum(w * transformed_f(y))

    def _serialize_tree(self, node):
        """Recursively serializes the tree to a dictionary for JSON output."""
        if node is None:
            return None
        return {
            "a": node.a,
            "b": node.b,
            "depth": node.depth,
            "tol": node.tol,
            "error": node.error,
            "integral": node.integral,
            "method": "Gauss-Laguerre" if node.is_singular else "Gauss-Legendre",
            "left": self._serialize_tree(node.left),
            "right": self._serialize_tree(node.right)
        }

    def _deserialize_tree(self, data):
        """Recursively reconstructs the tree from a dictionary."""
        if data is None:
            return None
        node = self.Node(
            a=data["a"], 
            b=data["b"], 
            depth=data["depth"], 
            tol=data["tol"], 
            error=data["error"], 
            integral=data["integral"],
            is_singular=(data["method"] == "Gauss-Laguerre")
        )
        node.left = self._deserialize_tree(data["left"])
        node.right = self._deserialize_tree(data["right"])
        return node

    
    def integrate(self):
        """Computes the total integral."""
        return self._traverse_and_sum(self.root)

    def _traverse_and_sum(self, node):
        """Recursively sums the integral values from leaf nodes and accumulates total error."""
        if node is None:
            return 0.0, 0.0  # Return both integral and error as 0 for empty nodes

        if node.left is None and node.right is None:
            return node.integral, node.error  # Return integral and error for leaf nodes

        left_integral, left_error = self._traverse_and_sum(node.left)
        right_integral, right_error = self._traverse_and_sum(node.right)

        total_integral = left_integral + right_integral
        total_error = left_error + right_error  # Accumulate errors from all subnodes

        return total_integral, total_error


    def save_to_json(self, filename, dump_roots=False):
        """Saves the quadrature tree to a JSON file."""
        data = {
            "name": self.name,
            "reference": self.reference,
            "description": self.description,
            "author": self.author,
            "version": self.version,
            "update_log": self.update_log,
            "tolerance": self.tol,
            "min_depth": self.min_depth,
            "max_depth": self.max_depth,
            "n1":  self.n1,
            "n2":  self.n2,           
        }
        if dump_roots:
            data["legendre_roots_n1"] = [self.legendre_roots_n1[0].tolist(), self.legendre_roots_n1[1].tolist()]
            data["laguerre_roots_n1"] = [self.laguerre_roots_n1[0].tolist(), self.laguerre_roots_n1[1].tolist()]
            data["legendre_roots_n2"] = [self.legendre_roots_n2[0].tolist(), self.legendre_roots_n2[1].tolist()]
            data["laguerre_roots_n2"] = [self.laguerre_roots_n2[0].tolist(), self.laguerre_roots_n2[1].tolist()]
        
        data["tree"] = self._serialize_tree(self.root)           
        
        with open(filename, "w") as f:
            json.dump(data, f, indent=4)

    def _load_from_json(self, filename):
        """Loads a quadrature tree from a JSON file and reconstructs the tree structure."""
        with open(filename, "r") as f:
            data = json.load(f)

        # Load metadata
        self.name = data.get("name", "Adaptive Quadrature Tree")
        self.reference = data.get("reference")
        self.description = data.get("description")
        self.author = data.get("author")
        self.version = data.get("version", "1.0")
        self.update_log = data.get("update_log", [])

        # Load quadrature parameters
        self.tol = data.get("tolerance", 1e-6)
        self.min_depth = data.get("min_depth", 2)
        self.max_depth = data.get("max_depth", 10)

        self.n1 = data.get("n1", 5)
        self.n2 = data.get("n2", 10)          
        # Load quadrature roots if present
        if "legendre_roots_n1" in data:
            self.legendre_roots_n1 = [np.array(data["legendre_roots_n1"][0]), np.array(data["legendre_roots_n1"][1])]
        else:
            self.legendre_roots_n1 = roots_legendre(self.n1)
        if "laguerre_roots_n1" in data:
            self.laguerre_roots_n1 = [np.array(data["laguerre_roots_n1"][0]), np.array(data["laguerre_roots_n1"][1])]
        else:
            self.laguerre_roots_n1 = roots_legendre(self.n2)           
        if "legendre_roots_n2" in data:
            self.legendre_roots_n2 = [np.array(data["legendre_roots_n2"][0]), np.array(data["legendre_roots_n2"][1])]
        else:
            self.legendre_roots_n2 = roots_legendre(self.n2)
        if "laguerre_roots_n2" in data:
            self.laguerre_roots_n2 = [np.array(data["laguerre_roots_n2"][0]), np.array(data["laguerre_roots_n2"][1])]
        else:
            self.laguerre_roots_n2 = roots_laguerre(self.n2)
            
        # Load integration limits
        self.a = data["tree"]["a"]
        self.b = data["tree"]["b"]

        # Reconstruct the tree structure
        self.root = self._deserialize_tree(data["tree"])


    
    def print_tree(self, node=None, level=0):
        """Prints the tree structure."""
        if node is None:
            node = self.root
        indent = "  " * level
        method = "Gauss-Laguerre" if node.is_singular else "Gauss-Legendre"
        print(f"{indent}Node ({method}): [{node.a}, {node.b}], Depth: {node.depth}, "
              f"Integral: {node.integral:.10f}, Error: {node.error:.5e}")
        if node.left:
            self.print_tree(node.left, level + 1)
        if node.right:
            self.print_tree(node.right, level + 1)

    def add_update_log(self, message):
        """Adds a timestamped message to the update log."""
        self.update_log.append({
            "timestamp": datetime.datetime.now().isoformat(),
            "message": message
        })
