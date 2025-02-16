import json
import scipy.special

def generate_quadrature_json(method="Legendre", n_max=100, filename="quadrature_data.json"):
    data = {
        "method": method,
        "n_max": n_max,
        "n": {}
    }

    for n in range(1, n_max + 1):
        if method == "Legendre":
            abscissas, weights = scipy.special.roots_legendre(n)
        elif method == "Laguerre":
            abscissas, weights = scipy.special.roots_laguerre(n)
        else:
            raise ValueError("Unsupported method. Choose 'Legendre' or 'Laguerre'.")

        data["n"][str(n)] = {
            "0": abscissas.tolist(),
            "1": weights.tolist()
        }

    with open(filename, "w") as f:
        json.dump(data, f, indent=4)

    return filename



