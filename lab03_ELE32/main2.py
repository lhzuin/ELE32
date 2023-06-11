import numpy as np
from scipy.sparse import csr_matrix, lil_matrix
from random import randint
import pandas as pd



N = 5
M = 4

sparse_matrix = np.array([[1, 1, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 1, 1]])
sparse_matrix = csr_matrix(sparse_matrix)
codified_word = np.array([0, 0, 1, 0, 0])
inferred_word = codified_word
#print(inferred_word)
error_list = np.array([0]*N)
#print("decoding word")
for i in range(10):
    print(f"idx: {i}")
    #equation_line = codified_word @ sparse_matrix

    # Descubro quais equações não foram satisfeitas
    equation_line = inferred_word @ sparse_matrix
    equation_line = np.remainder(equation_line, 2)
    #equation_matrix = np.repeat(equation_line, N, axis=0)

    # Repito a suposta informação, replicando-a em N linhas
    equation_matrix = np.tile(equation_line, (N, 1))
    #print(equation_matrix)
    sparse_equation_matrix = csr_matrix(equation_matrix)

    # Encontro número de equações erradas que cada nó participa
    result = sparse_equation_matrix.multiply(sparse_matrix)

    error_list = result.sum(axis=1).flatten()

    result_array = np.zeros(N)

    max_value = np.amax(error_list)

    if max_value == 0:
        print("*** COND 1 ***")
        break

    # Get the indices of elements in error_list equal to max_value
    indices = np.where(error_list == max_value)[1]
    print(indices)

    result_array[indices] = 1

    if sum(result_array) == 0:
        print("*** COND 2 ***")
        break

    inferred_word = inferred_word + result_array

    inferred_word = np.remainder(inferred_word, 2)
    print(inferred_word)

print(inferred_word)