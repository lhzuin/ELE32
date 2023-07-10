import numpy as np


class NoEncode:
    def __init__(self, codified_code):
        self.codified_code = np.array(codified_code)
        self.size = len(self.codified_code)
        self.k = 1
        self.n = 1
        self.r = 1
        self.name = "No Encode"
        self.is_binary = True
        self.receive_L = False
    
    def encoder(self):
        return self.code
    
    def decoder(self, received_code):
        return received_code

