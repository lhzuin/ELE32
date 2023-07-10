import numpy as np

class BPSK:
    def __init__(self, N0, Eb, codified_code):
        self.N0 = N0
        self.codified_code = codified_code
        self.Eb = Eb
        self.size = len(codified_code)
        self.get_values()

    def get_values(self):
        ones = np.ones((1, self.size))
        self.values = (ones - 2*self.codified_code)*np.sqrt(self.Eb) # transforma zeros e uns em -1 e +1
    
    def transform_code(self, is_binary = True):
        n = np.random.normal(0, self.N0/2, self.size)
        transformed_code = n+self.values
        if is_binary:
            return convert_to_binary(transformed_code)

        return transformed_code.flatten()
    


def convert_to_binary(codified_code):
        positive_mask = codified_code >= 0
        negative_mask = codified_code < 0

        # Use the mask to set positive values to 0 and negative ones to 1
        codified_code[positive_mask] = 0
        codified_code[negative_mask] = 1

        return codified_code.flatten()