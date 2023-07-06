import numpy as np


class NoEncode:
    def __init__(self, code):
        self.code = np.array(code)
        self.size = len(self.code)
        self.k = 1
        self.name = "No Encode"
        self.is_binary = True
    
    def encoder(self):
        return self.code
    
    def decoder(self, received_code):
        return received_code

