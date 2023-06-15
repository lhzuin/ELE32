import numpy as np
import random

class BSC:
    def __init__(self, p, code):
        self.p = p
        self.code = code
    
    def generate_randomness(self, x):
        x = x%2
        #print(6)
        random_num = random.random()
        #print(7)
        if random_num < self.p:
             #print(8)
             return int(not x)
        else:
             #print(9)
             return int(x)
    """
    def transform_code(self):
        #print(5)
        result = np.vectorize(self.generate_randomness)(self.code)
        return result
    """
    def transform_code(self):
        x = self.code % 2
        random_num = np.random.rand(*x.shape)
        mask = random_num < self.p
        result = np.where(mask, 1-x, x)
        return result