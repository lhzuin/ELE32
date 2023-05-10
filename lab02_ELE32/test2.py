import numpy as np
from itertools import product
from sympy import symbols, Poly, GF, div
from generator_polys import  vec_to_poly, poly_to_vector_fixed_len

class CyclicEncode:
    def __init__(self, code, n, g):
        self.code = code
        self.size = len(self.code)
        self.codified_word_size = n
        self.g = g
        self.x = symbols('x')
        self.poly_g = vec_to_poly(np.array(self.g), self.x)
        self.g_degree = len(g)-1
        self.word_size = self.codified_word_size - self.g_degree
        self.rotate_constant = self.g[:0:-1]
        if self.size % self.word_size != 0:
            raise ValueError(f"Code size must be divisible by {self.word_size}")
        self.s_list = self.get_syndrome_list()  # Cache the results of get_syndrome_list in __init__

    @property
    def divide_code(self):
        matrix = np.array([self.code[i:i+self.word_size] for i in range(0, self.size, self.word_size)])
        return matrix

    def encoder(self):
        grouped_code = self.divide_code
        encoded_code = [self.encode_word(u) for u in grouped_code]
        return np.array(encoded_code, dtype=int)

    def encode_word(self, u):
        poly_u = vec_to_poly(u, self.x)
        result = self.poly_g * poly_u
        result = Poly(result.as_expr(), self.x, domain=GF(2))
        result = poly_to_vector_fixed_len(result, self.codified_word_size)
        return result

    def decoder(self, received_code):
        return np.array([self.decode_word(codified_word) for codified_word in received_code]).flatten()

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
"""
    def get_info_word(self, s, v):
        num_of_rotates = 0
        num_of_changes = 0
        while(sum(s)!= 0 and num_of_rotates < self.codified_word_size and num_of_changes<2):
            if tuple(s) in self.s_list:
                v[0] = int(not v[0])
                s = self.get_syndrome(v)
                num_of_changes += 1
            else:
                s = self.rotate_s(s)
                v = self.rotate_v(v)
                num_of_rotates += 1

        v = np.roll(v, shift=(self.codified_word_size - num_of_rotates))
        codified_poly = vec_to_poly(v, self.x)
        info_poly, remainder = div(codified_poly, self.poly_g)
        info_word = poly_to_vector_fixed_len(info_poly, self.word_size)
        return np.array(info_word)
"""
    def get_info_word(self, s, v):
        num_of_rotates = 0
        v_cycle = cycle(v)  # Use itertools.cycle for rotations

        while sum(s) != 0 and num_of_rotates < 2 * self.codified_word_size:
            if tuple(s) in self.s_list:
                v[0] = int(not v[0])
                s = self.get_syndrome(v)
            else:
                s = self.rotate_s(s)
                next(v_cycle)  # Rotate by calling next() on the cycle
                num_of_rotates += 1

        # Update the value of v using v_cycle before calling np.roll
        v = [next(v_cycle) for _ in range(self.codified_word_size)]

        # Perform the necessary rotations in one call to np.roll
        v = np.roll(v, shift=(self.codified_word_size - num_of_rotates % self.codified_word_size))

        codified_poly = vec_to_poly(v, self.x)
        info_poly, remainder = div(codified_poly, self.poly_g)
        info_word = poly_to_vector_fixed_len(info_poly, self.word_size)

        return np.array(info_word)
    def rotate_s(self, s):
        #print("init")
        rotated = np.roll(s, shift=1)
        if rotated[0] == 1:
            rotated[0] = 0
            #rotated = rotated + self.rotate_constant
            np.add(rotated, self.rotate_constant, out=rotated)
            np.remainder(rotated, 2, out=rotated)
            #rotated = np.remainder(rotated, 2)
        #print("finish")
        return rotated

    def rotate_v(self, v):
        return np.roll(v, shift=1)
