import numpy as np
from itertools import product



class Hamming:
    def __init__(self, codified_code):
        self.codified_code = codified_code
        self.size = len(self.codified_code)
        self.k = 4
        self.n = 7
        self.r = self.k/self.n
        self.name = "Hamming"
        self.is_binary = True
        self.receive_L = False
        if self.size % 7 != 0:
            raise ValueError("Code size must be divisable by 7")
        self.transformation = np.array([[1, 0, 0, 0, 1, 1, 1], [0, 1, 0, 0, 1, 0, 1], [0, 0, 1, 0, 1, 1, 0], [0, 0, 0, 1, 0, 1, 1]])
        self.verification = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 0], [0, 1, 1], [1, 0, 0], [0, 1, 0], [0, 0, 1]])

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


    def decoder(self, received_code):
        if received_code.ndim == 1:
            received_code = np.array([received_code[i:i+7] for i in range(0, len(received_code), 7)])

        s_array = self.get_syndrome(received_code)
        s_dict = self.get_dict()
        for i in range(len(s_array)):
            s = tuple(s_array[i].tolist())
            e = s_dict[s]
            received_code[i] = np.remainder(received_code[i] + e, 2)
        
        information = received_code[:, :4]
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




    
