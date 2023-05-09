# Preciso pegar todas as possíveis palavras do código e achar sua distância
import numpy as np
from itertools import product
from cyclic_encode import CyclicEncode
from generator_polys import generator_polynomials, find_best_generator

def get_distance(a, b):
    differences = np.bitwise_xor(a, b)
    return np.sum(differences)



def get_messages(k):
    messages = list(product([0, 1], repeat=k))
    return np.array(messages)

def get_minimum_distance(n, g):
    k = len(g)
    code = get_messages(k).flatten()
    cyclic = CyclicEncode(code, n, g)

    cyclic_code = cyclic.encoder()

    row_sums = cyclic_code.sum(axis=1)

    second_lowest_index = 1  # Index of the second lowest value in a sorted array

    second_lowest = np.partition(row_sums, second_lowest_index)[second_lowest_index]
    return second_lowest

if __name__ == "__main__":
    n = 7
    k = 4
    #code = np.zeros(k)
    code = get_messages(k).flatten()

    gen_polys = generator_polynomials(n, k)
    g = find_best_generator(gen_polys)
    print(f"k = {len(g)}")

    cyclic = CyclicEncode(code, n, g)

    cyclic_code = cyclic.encoder()

    row_sums = cyclic_code.sum(axis=1)
    #sorted_cyclic_code = cyclic_code[np.argsort(row_sums)]

    second_lowest_index = 1  # Index of the second lowest value in a sorted array

    second_lowest = np.partition(row_sums, second_lowest_index)[second_lowest_index]

    print(second_lowest)

