import numpy as np

class BPSK:
    def __init__(self, N0, Eb, code):
        self.N0 = N0
        self.code = code
        self.Eb = Eb
        self.size = len(code)
        self.get_values()

    def get_values(self):
        ones = np.ones((1, self.size))*np.sqrt(self.Eb)
        self.values = ones - 2*self.code
    
    def transform_code(self, is_binary = True):
        n = np.random.normal(0, self.N0/2, self.size)

        if is_binary:
            return convert_to_binary(n+self.code)

        return n+self.code
    


def convert_to_binary(code):
        positive_mask = code >= 0
        negative_mask = code < 0

        # Use the mask to set positive values to 0 and negative ones to 1
        code[positive_mask] = 0
        code[negative_mask] = 1

        return code