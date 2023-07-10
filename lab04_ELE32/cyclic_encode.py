import numpy as np
from itertools import product
from sympy import symbols, Poly, GF, div
from generator_polys import  vec_to_poly, poly_to_vector_fixed_len
from collections import deque
from generator_polys import generator_polynomials, find_best_generator

class CyclicEncode:
    def __init__(self, codified_code, n, k):
        self.codified_code = codified_code
        self.size = len(self.codified_code)
        self.codified_word_size = n # n
        self.gen_polys = generator_polynomials(n, k)
        self.g = find_best_generator(self.gen_polys)
        self.x = symbols('x')
        self.poly_g = vec_to_poly(np.array(self.g), self.x)
        self.g_degree = len(self.g)-1
        self.word_size = self.codified_word_size - self.g_degree # k
        self.rotate_constant = self.g[:0:-1]
        if self.size % self.codified_word_size != 0:
            raise ValueError(f"Code size must be divisible by {self.codified_word_size}")
        self.s_list = self.get_syndrome_list()  # Cache the results of get_syndrome_list in __init__
        self.n = n
        self.k = k
        self.r = k/n
        self.name = f"Cyclic Encode n={n}, k={k}"
        self.is_binary = True
        self.receive_L = False


    def decoder(self, received_code):
        return np.array([self.decode_word(received_code[i:i+self.codified_word_size]) for i in range(0, len(received_code), self.codified_word_size)]).flatten()

    def decode_word(self, codified_word):
        s = self.get_syndrome(codified_word)
        info_word = self.get_info_word(s, codified_word)
        return info_word

    def get_syndrome(self, codified_word):
        #poly = Poly.from_list(codified_word.tolist(), self.x, modulus=2)
        poly = vec_to_poly(codified_word, self.x)
        quotient, remainder = div(poly, self.poly_g)
        s = poly_to_vector_fixed_len(remainder, self.codified_word_size - self.word_size)
        return np.array(s)

    def get_syndrome_list(self):
        combinations = list(product([0, 1], repeat=(self.codified_word_size-1)))
        combinations = [[1] + list(sublist) for sublist in combinations]
        binary_array = np.array(combinations)
        row_sums = binary_array.sum(axis=1)
        mask = row_sums < 3
        filtered_array = binary_array[mask]

        sorted_binary_array = filtered_array[np.argsort(filtered_array.sum(axis=1))]
        s_list = []
        for i in range(len(sorted_binary_array)):
            s = tuple(self.get_syndrome(sorted_binary_array[i]).tolist())
            if s not in s_list:
                s_list.append(s)
        return s_list

    def get_info_word(self, s, v):
        s = deque(s)
        v = deque(v)
        num_of_rotates = 0
        num_of_changes = 0
        while(sum(s) != 0 and num_of_rotates < self.codified_word_size and num_of_changes < 2):
            if tuple(s) in self.s_list:
                v[0] = int(not v[0])
                s = self.get_syndrome(np.array(v))
                s = deque(s)
                num_of_changes += 1
            else:
                s = self.rotate_s(s)
                v.rotate(1)
                num_of_rotates += 1

        v.rotate(self.codified_word_size - num_of_rotates)
        codified_poly = vec_to_poly(np.array(v), self.x)
        info_poly, remainder = div(codified_poly, self.poly_g)
        info_word = poly_to_vector_fixed_len(info_poly, self.word_size)
        return np.array(info_word)

    def rotate_s(self, s):
        s_deque = deque(s)
        s_deque.rotate(1)
        if s_deque[0] == 1:
            s_deque[0] = 0
            s_deque = deque(map(lambda x, y: (x + y) % 2, s_deque, self.rotate_constant))
        return s_deque
    
    def rotate_v(self, v):
        return np.roll(v, shift=1)