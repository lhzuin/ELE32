from sympy import symbols, Poly, GF, factor, factor_list
import itertools
import numpy as np

def poly_to_vector(p):
    if p.is_zero:
        return [0]
    coeffs = p.all_coeffs()
    # Pad with zeros if necessary to match the degree of the polynomial
    vector = [0] * (p.degree() + 1 - len(coeffs)) + coeffs
    return vector

def vec_to_poly(v, x):
    poly = Poly(np.poly1d(v.astype(int)), x, domain=GF(2))
    return poly
#"""
def poly_to_vector_fixed_len(p, desired_len):
    if p.is_zero:
        return np.array([0]*desired_len)
    coeffs = p.all_coeffs()
    size = max(p.degree() + 1, desired_len)

    # Pad with zeros if necessary to match the degree of the polynomial
    vector = [0] * (size - len(coeffs)) + coeffs
    return np.array(vector)

def vec_to_poly_with_degree(v, x):
    poly = Poly(np.poly1d(v), x, domain=GF(2))
    return poly


def find_irreducible_factors_berlekamp(n):
    x = symbols("x")
    p_x = Poly(x**n + 1, x, domain=GF(2))
    factored_p_x = factor_list(p_x.as_expr(), domain=GF(2))
    #print(factored_p_x)

    # Extract the irreducible factors
    irreducible_factors = []
    for base, exp in factored_p_x[1]:
        for i in range(exp):
            irreducible_factors.append(Poly(base, x, domain=GF(2)))#**exp)
    #print(f"irreducible_factors: {irreducible_factors}")
    return irreducible_factors
def generator_polynomials(n, k):
    irreducible_factors = find_irreducible_factors_berlekamp(n)
    degree = n - k
    gen_polys = []

    for r in range(1, len(irreducible_factors) + 1):
        for combination in itertools.combinations(irreducible_factors, r):
            candidate = Poly(1, symbols("x"), domain=GF(2))

            for factor in combination:
                candidate *= factor
            
            if candidate.degree() == degree:
                v = poly_to_vector(candidate)
                if v not in gen_polys:
                    gen_polys.append(v)

    return gen_polys

def find_best_generator(gen_polys):
    if not gen_polys:
        return []
    min_module = len(gen_polys[0])
    best_g = None
    for g in gen_polys:
        mod = sum(g)
        #print(mod)
        if mod <= min_module:
            min_module = mod
            best_g = g
    #print(best_g)
    return best_g
    
if __name__ == "__main__":
    n = 7
    k = 3
    gen_polys = generator_polynomials(n, k)

    for g in gen_polys:
        print(g)

    g = find_best_generator(gen_polys)    
