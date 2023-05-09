import math
from generator_polys import generator_polynomials, find_best_generator
from find_minimal_distance import get_minimum_distance
from cyclic_encode import CyclicEncode
import numpy as np
from bsc import BSC

def find_possible_ks(n):
    lower_limit = 4*n*0.9/7
    upper_limit = 4*n*1.1/7

    lower_bound = math.ceil(lower_limit)
    upper_bound = math.floor(upper_limit)
    integers = list(range(lower_bound, upper_bound + 1))
    return integers

"""
if __name__ == "__main__":
    found_solution = False
    for n in range(8, 41):
        if found_solution:
            break
        for k in find_possible_ks(n):
            print(f"({n},{k})")
            gen_polys = generator_polynomials(n, k)
            print(f"gen_polys: {gen_polys}")
            if gen_polys:
                g = find_best_generator(gen_polys)
                distance = get_minimum_distance(n, g)
                if distance == 5:
                    print(f"Found Solution! n = {n}, k = {k}")
                    found_solution = True
                    break

"""


#vector_length = 1000000
p = 0.01

if __name__ == "__main__":
    n = 17
    k = 9

    vector_length = 10*k
    info = np.zeros(vector_length)
    gen_polys = generator_polynomials(n, k)
    g = find_best_generator(gen_polys)
    print(f"g: {g}")
    cyclic = CyclicEncode(info, n, g)

    cyclic_code = cyclic.encoder()

    bsc_cyclic = BSC(p, cyclic_code)

    received_code_cyclic = bsc_cyclic.transform_code()

    inferred_info_cyclic = cyclic.decoder(received_code_cyclic)

    num_of_errors_cyclic = np.sum(inferred_info_cyclic)
    print(num_of_errors_cyclic)

