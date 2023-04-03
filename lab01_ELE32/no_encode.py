import numpy as np
import math


class NoEncode:
    def __init__(self, code):
        self.code = np.array(code)
        self.size = len(self.code)

