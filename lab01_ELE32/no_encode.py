import numpy as np
import math


class NoEncode:
    def __init__(self, code):
        self.code = np.array(code)
        self.size = len(self.code)
    
    def encoder(self):
        return self.code
    
    def decoder(self, received_code):
        return received_code

