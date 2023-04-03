import numpy as np
import random
from itertools import product
import matplotlib.pyplot as plt
import math


class NoEncode:
    def __init__(self, code):
        self.code = np.array(code)
        self.size = len(self.code)
        #self.transformation = np.array([[1, 0, 0, 0, 1, 1, 1], [0, 1, 0, 0, 1, 0, 1], [0, 0, 1, 0, 1, 1, 0], [0, 0, 0, 1, 0, 1, 1]])
        #self.verification = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 0], [0, 1, 1], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
