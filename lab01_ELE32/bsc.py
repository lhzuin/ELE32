import numpy as np
import random

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