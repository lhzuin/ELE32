import numpy as np
from itertools import product



class CustomEncode:
    # garante acerto se até um erro for cometido
    def __init__(self, code):
        self.code = code
        self.size = len(self.code)
        self.word_size = 5
        self.codified_word_size = 10
        if self.size % self.word_size != 0:
            raise ValueError(f"Code size must be divisable by {self.word_size}")
        self.transformation = np.array([[1, 0, 0, 0, 0, 1, 1, 0, 0, 0], [0, 1, 0, 0, 0, 0, 1, 1, 0, 0], [0, 0, 1, 0, 0, 0, 0, 1, 1, 0], [0, 0, 0, 1, 0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 1, 1, 0, 0, 0, 1]])
        self.verification = np.array([[1, 1, 0, 0, 0], [0, 1, 1, 0, 0], [0, 0, 1, 1, 0], [0, 0, 0, 1, 1], [1, 0, 0, 0, 1], [1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]])

    
    @property
    def divide_code(self):
        matrix = np.array([self.code[i:i+self.word_size] for i in range(0, self.size, self.word_size)])
        return matrix
    
    @property
    def generate_transf(self):
        result = np.matmul(self.grouped_code, self.transformation)
        return np.remainder(result, 2)

    def get_dict(self):
        s_dict = {}
        combinations = list(product([0, 1], repeat=self.codified_word_size))
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
        self.hamming_code = self.generate_transf
        return self.hamming_code

    def decoder(self, received_code):
        s_array = self.get_syndrome(received_code)
        s_dict = self.get_dict()
        for i in range(len(s_array)):
            s = tuple(s_array[i].tolist())
            e = s_dict[s]
            received_code[i] = np.remainder(received_code[i] + e, 2)
        
        information = received_code[:, :self.word_size]
        return information.flatten()

        


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




    
