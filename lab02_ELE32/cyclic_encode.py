import numpy as np
from itertools import product
from sympy import symbols, Poly, GF, div
from generator_polys import poly_to_vector, vec_to_poly, poly_to_vector_fixed_len



class CyclicEncode:
    # garante acerto se até um erro for cometido
    def __init__(self, code, n, g):
        self.code = code
        self.size = len(self.code)
        self.codified_word_size = n
        self.g = g
        self.x = symbols('x')
        self.poly_g = vec_to_poly(self.g, self.x)
        self.g_degree = len(g)-1
        self.word_size = self.codified_word_size - self.g_degree
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
        #print(f"combinations: {combinations}")
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
        mask = row_sums < 3
        filtered_array = binary_array[mask]
        
        sorted_binary_array = filtered_array[np.argsort(filtered_array.sum(axis=1))]
        #print(sorted_binary_array)
        #print(sorted_binary_array)
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
        #print(f"poly: {poly}")
        #g = Poly.from_list(self.g.to_list(), self.x, modulus=2)
        #poly_g = Poly(np.poly1d(self.g), self.x, domain=GF(2))
        #print(f"g: {self.poly_g}")
        quotient, remainder = div(poly, self.poly_g)
        #print(f"remainder: {remainder}")
        #s = poly_to_vector(remainder)
        #s = [0]*(self.codified_word_size - self.word_size - len(s)) + s
        s = poly_to_vector_fixed_len(remainder, self.codified_word_size - self.word_size)
        #print(f"s: {s}")
        return np.array(s)
            #s.append(poly_to_vector(remainder))
   
        #return np.array(s)

    def encoder(self):
        self.grouped_code = self.divide_code
        encoded_code = []
        for u in self.grouped_code:
            #poly_g = Poly(np.poly1d(self.g), self.x, domain=GF(2))
            #poly_u = Poly(np.poly1d(u), self.x, domain=GF(2))
            poly_u = vec_to_poly(u, self.x)
            #print(f"u: {poly_u}")
            #print(f"g: {poly_g}")

            result = self.poly_g * poly_u
            #print(f"result1: {result}")
            result = Poly(result.as_expr(), self.x, domain=GF(2))
            #print(f"result2: {result}")
            #result = poly_to_vector(result)

            #result = [0]*(self.codified_word_size - len(result)) + result
            result = poly_to_vector_fixed_len(result, self.codified_word_size)
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
        #print(f"original_v: {v}")
        while(sum(s)!= 0 and num_of_rotates < 2*self.codified_word_size):
            if tuple(s) in self.s_list:
                #print(1)
                v[0] = int(not v[0])
                #print(f"new_v: {v}")
                s = self.get_syndrome(v)
                #print(f"new_s: {s}")
            else:
                #print(2)
                s = self.rotate_s(s)
                v = self.rotate_v(v)
                #print(f"rotated_v: {v}")
                #print(f"rotated_s: {s}")
                num_of_rotates += 1
        #print(f"num_of_rotates: {num_of_rotates}")
        num_of_rotates = num_of_rotates%self.codified_word_size
        for i in range(self.codified_word_size - num_of_rotates):
            #print(f"rotated_v2: {v}")
            v = self.rotate_v(v)
        #codified_poly = Poly.from_list(v.tolist(), self.x, modulus=2)
        v = v.astype(int)
        codified_poly = vec_to_poly(v, self.x)

        info_poly, remainder = div(codified_poly, self.poly_g)
        #info_word = poly_to_vector(info_poly)
        #info_word = [0]*(self.word_size - len(info_word)) + info_word
        info_word = poly_to_vector_fixed_len(info_poly, self.word_size)
        return np.array(info_word)

    def rotate_s(self, s):
        rotated = np.roll(s, shift=1)
        if rotated[0] == 1:
            rotated[0] = 0
            #rotated = rotated + self.g[:-1] # rever índice
            rotated = rotated + self.g[:0:-1]
            rotated = np.remainder(rotated, 2)
        #print(f"rotated_s: {rotated}")
        return rotated.astype(int)

    def rotate_v(self, v):
        return np.roll(v, shift=1).astype(int)
    
    def reset(self, code):
        self.code = code
        self.size = len(self.code)
        if self.size % self.word_size != 0:
            raise ValueError(f"Code size must be divisable by {self.word_size}")


if __name__ == '__main__':
    n = 17
    k = 9
    #vector_length = find_vector_length()
    vector_length = 1*k
    info = np.zeros(vector_length)
    g = np.array([1, 0, 0, 1, 1, 1, 0, 0, 1])
    print(f"g: {g}")
    cyclic = CyclicEncode(info, n, g)

    cyclic_code = cyclic.encoder()
    print(f"s: {cyclic.get_syndrome(np.array([1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))}")

    #bsc_cyclic = BSC(p, cyclic_code)
    """
    received_code_cyclic = bsc_cyclic.transform_code()
    received_code_cyclic = np.array([[1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0]])
    inferred_info_cyclic = cyclic.decoder(received_code_cyclic)

    num_of_errors_cyclic = np.sum(inferred_info_cyclic)
    print(num_of_errors_cyclic)
    """





    
