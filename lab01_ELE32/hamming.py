import numpy as np
import random
from itertools import product
import matplotlib.pyplot as plt
import math


vector_length = 1000000
p = 0.5


class BSC:
    def __init__(self, p, code):
        self.p = p
        self.code = code
    
    def generate_randomness(self, x):
        x = x%2
        random_num = random.random()
        if random_num < self.p:
             return not x
        else:
             return x
     
    def transform_code(self):
        result = np.vectorize(self.generate_randomness)(self.code)
        return result
        """
        
        new_code = np.empty([])
        for word in self.code
        for bit in self.code:
            
                new_code.append(not bit)
            else: 
                new_code.append(bit)
        """
class Hamming:
    def __init__(self, code):
        self.code = code
        self.size = len(self.code)
        if self.size % 4 != 0:
            raise ValueError("Code size must be divisable by 4")
        self.transformation = np.array([[1, 0, 0, 0, 1, 1, 1], [0, 1, 0, 0, 1, 0, 1], [0, 0, 1, 0, 1, 1, 0], [0, 0, 0, 1, 0, 1, 1]])
        self.verification = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 0], [0, 1, 1], [1, 0, 0], [0, 1, 0], [0, 0, 1]])

    
    @property
    def divide_code(self):
        matrix = np.array([self.code[i:i+4] for i in range(0, self.size, 4)])
        return matrix
    
    @property
    def generate_hamming(self):
        result = np.matmul(self.grouped_code, self.transformation)
        return np.remainder(result, 2)

    def get_dict(self):
        s_dict = {}
        combinations = list(product([0, 1], repeat=7))
        binary_array = np.array(combinations)
        row_sums = binary_array.sum(axis=1)
        sorted_binary_array = binary_array[np.argsort(row_sums)]

        for i in range(len(sorted_binary_array)):
            s = tuple(self.get_syndrome(sorted_binary_array[i]).tolist())
            if s not in list(s_dict.keys()):
                s_dict[s] = sorted_binary_array[i]
        return s_dict

    def get_syndrome(self, received_code):
        s = np.matmul(received_code, self.verification)
        return np.remainder(s, 2)

    def encoder(self):
        self.grouped_code = self.divide_code
        self.hamming_code = self.generate_hamming
        return self.hamming_code

    def decoder(self, received_code):
        s_array = self.get_syndrome(received_code)
        s_dict = self.get_dict()
        for i in range(len(s_array)):
            s = tuple(s_array[i].tolist())
            e = s_dict[s]
            received_code[i] = np.remainder(received_code[i] + e, 2)
        
        information = received_code[:, :4]
        return information

        


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


if __name__ == '__main__':
    p_list = [0.5]
    last_num = 0.5
    for i in range(30):
        power = pow(10,round(i/3)+1)
        last_num = math.floor(last_num*power/2)/power
        p_list.append(last_num)
    print(p_list)

    prob_list = []
    for p in p_list:
        information = np.random.randint(2, size=vector_length)
        hamming = Hamming(information)
        hamming_code = hamming.encoder()
        bsc = BSC(p, hamming_code)
        received_code = bsc.transform_code()
        inferred_code = hamming.decoder(received_code)
        inferred_info = inferred_code.flatten()
        differences = np.bitwise_xor(np.array(inferred_info, dtype=np.int64), information)
        num_of_errors = np.sum(differences)
        error_prob = num_of_errors/vector_length
        prob_list.append(error_prob)
        #print(f"num_of_errors: {num_of_errors}")
    
    print(prob_list)
    plt.figure()
    plt.plot(p_list, prob_list)
    plt.gca().invert_xaxis()
    plt.xscale('log',  basex=10)
    plt.yscale('log',  basey=10)
    plt.savefig("mygraph.png")
    #plt.show()



    
