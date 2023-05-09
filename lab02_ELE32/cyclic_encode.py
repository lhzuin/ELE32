import numpy as np
from itertools import product
from sympy import symbols, Poly, GF, div
from generator_polys import poly_to_vector



class CyclicEncode:
    # garante acerto se até um erro for cometido
    def __init__(self, code, n, g):
        self.code = code
        self.size = len(self.code)
        self.codified_word_size = n
        self.g = g
        self.word_size = len(g)
        self.x = symbols('x')
        self.s_list = self.get_syndrome_list()
        if self.size % self.word_size != 0:
            raise ValueError(f"Code size must be divisable by {self.word_size}")
        
    @property
    def divide_code(self):
        matrix = np.array([self.code[i:i+self.word_size] for i in range(0, self.size, self.word_size)])
        return matrix

    def get_dict(self):
        s_dict = {}
        combinations = list(product([0, 1], repeat=(self.codified_word_size-1)))

        combinations = [[1] + list(sublist) for sublist in combinations]
        print(f"combinations: {combinations}")
        #combinations.append([0]*self.codified_word_size)
        binary_array = np.array(combinations)
        row_sums = binary_array.sum(axis=1)
        sorted_binary_array = binary_array[np.argsort(row_sums)]
        s_list = []
        for i in range(len(sorted_binary_array)):
            s = tuple(self.get_syndrome(sorted_binary_array[i]).tolist())
            if s not in s_list:
                s_list.append(s)
        
        return s
        """
        for i in range(len(sorted_binary_array)):
            s = tuple(self.get_syndrome(sorted_binary_array[i]).tolist())
            if s not in list(s_dict.keys()):
                s_dict[s] = tuple(sorted_binary_array[i])
        print(f"s_dict: {s_dict}")
        return s_dict
        """


    def get_syndrome_list(self):
        combinations = list(product([0, 1], repeat=(self.codified_word_size-1)))

        combinations = [[1] + list(sublist) for sublist in combinations]
        #print(f"combinations: {combinations}")
        #combinations.append([0]*self.codified_word_size)
        binary_array = np.array(combinations)
        row_sums = binary_array.sum(axis=1)
        mask = row_sums < 5
        filtered_array = binary_array[mask]
        
        sorted_binary_array = binary_array[np.argsort(filtered_array.sum(axis=1))]
        s_list = []
        for i in range(len(sorted_binary_array)):
            s = tuple(self.get_syndrome(sorted_binary_array[i]).tolist())
            if s not in s_list:
                s_list.append(s)
        print(f"s_list: {s_list}")
        return s_list
    
    def get_syndrome(self, codified_word):
        #s = []
        #for line in received_code:
        poly = Poly.from_list(codified_word.tolist(), self.x, modulus=2)
        #g = Poly.from_list(self.g.to_list(), self.x, modulus=2)
        poly_g = Poly(np.poly1d(self.g), self.x, domain=GF(2))
        quotient, remainder = div(poly, poly_g)
        s = poly_to_vector(remainder)
        s = [0]*(self.codified_word_size - self.word_size - len(s)) + s
        return np.array(s)
            #s.append(poly_to_vector(remainder))
   
        #return np.array(s)

    def encoder(self):
        self.grouped_code = self.divide_code
        encoded_code = []
        for u in self.grouped_code:
            poly_g = Poly(np.poly1d(self.g), self.x, domain=GF(2))
            poly_u = Poly(np.poly1d(u), self.x, domain=GF(2))
            #print(f"u: {poly_u}")
            #print(f"g: {poly_g}")

            result = poly_g * poly_u
            #print(f"result1: {result}")
            result = Poly(result.as_expr(), self.x, domain=GF(2))
            #print(f"result2: {result}")
            result = poly_to_vector(result)

            result = [0]*(self.codified_word_size - len(result)) + result
            #print(f"result3: {result}")
            encoded_code.append(result)

        return np.array(encoded_code, dtype=int)

    def decoder(self, received_code):
        information = []
        for codified_word in received_code:
            s = self.get_syndrome(codified_word)
            info_word = self.get_info_word(s, codified_word)
            information.append(info_word)
        information = np.array(information)
        return information.flatten()

    """
    def get_info_word(self, s, v):
        num_of_rotates = 0
        while(sum(s)!= 0):
            if tuple(s) in self.get_dict().keys():
                v[0] = int(not v[0])
                s = self.get_syndrome(v)
            s = self.rotate_s(s)
            v = self.rotate_v(v)
            num_of_rotates += 1
        
        for i in range(self.codified_word_size - num_of_rotates):
            self.rotate_v(v)
        codified_poly = Poly.from_list(v.tolist(), self.x, modulus=2)
        poly_g = Poly(np.poly1d(self.g), self.x, domain=GF(2))
        info_poly, remainder = div(codified_poly, poly_g)
        info_word = poly_to_vector(info_poly)
        return np.array(info_word)
    """
    def get_info_word(self, s, v):
        num_of_rotates = 0
        while(sum(s)!= 0 and num_of_rotates < self.codified_word_size):
            if tuple(s) in self.s_list:
                v[0] = int(not v[0])
                s = self.get_syndrome(v)
            s = self.rotate_s(s)
            v = self.rotate_v(v)
            num_of_rotates += 1
        
        for i in range(self.codified_word_size - num_of_rotates):
            self.rotate_v(v)
        codified_poly = Poly.from_list(v.tolist(), self.x, modulus=2)
        poly_g = Poly(np.poly1d(self.g), self.x, domain=GF(2))
        info_poly, remainder = div(codified_poly, poly_g)
        info_word = poly_to_vector(info_poly)
        info_word = [0]*(self.word_size - len(info_word)) + info_word
        return np.array(info_word)

    def rotate_s(self, s):
        rotated = np.roll(s, shift=1)
        if rotated[0] == 1:
            rotated[0] = 0
            rotated = rotated + self.g[:-1]
            rotated = np.remainder(rotated, 2)
        return rotated

    def rotate_v(self, v):
        return np.roll(v, shift=1)

"""
if __name__ == '__main__':
    information = np.random.randint(2, size=vector_length)
    print(f"information: {information}")
    hamming = Hamming(information)
    hamming_code = hamming.encoder()
    bsc = BSC(p, hamming_code)
    received_code = bsc.transform_code()
    print(f"received_code: {received_code}")
    inferred_code = hamming.decoder(received_code)
    print(f"inferred_code: {inferred_code}")
    inferred_info = inferred_code.flatten()
    print(f"inferred_info: {inferred_info}")
    differences = np.bitwise_xor(np.array(inferred_info, dtype=np.int64), information)
    num_of_errors = np.sum(differences)
    print(f"num_of_errors: {num_of_errors}")
"""




    
