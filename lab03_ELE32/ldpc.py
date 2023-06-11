from scipy.sparse import csr_matrix
import numpy as np
from random import randint
import pandas as pd


def find_M(N, dc, dv):
    M = N*dv/dc
    print(M)
    M = int(M)
    print(M)
    return M

def create_graph(N, dc, dv):
    M = find_M(N, dc, dv)
    dc_list = [0]*M
    dv_list = [0]*N

    min_dv = 0
    min_dc = 0

    num_cols = M
    num_rows = N

    zero_data = [[0]*num_cols for _ in range(num_rows)]

    # Convert the 2D array to a sparse matrix
    sparse_matrix = csr_matrix(zero_data)

    while(min_dc < dc and min_dv < dv):

        temp_dc = [i for i in range(M) if dc_list[i] == min_dc]
        temp_dv = [i for i in range(N) if dv_list[i] == min_dv]

        random_c = temp_dc[randint(0, len(temp_dc)-1)]
        random_v = temp_dv[randint(0, len(temp_dv)-1)]

        if sparse_matrix[random_v, random_c] != 1:
            sparse_matrix[random_v, random_c] = 1
            dc_list[random_c] += 1
            dv_list[random_v] += 1

        min_dc = min(dc_list)
        min_dv = min(dv_list)
    
    return sparse_matrix
    


def generate_csv(sparse_matrix, N):
    dense_matrix = sparse_matrix.toarray()

    # Create a pandas DataFrame
    df = pd.DataFrame(dense_matrix)

    # Save DataFrame to CSV
    df.to_csv(f"sparse_matrix_N_{N}.csv", index=False)





    # Crio coluna por coluna de verificação de paridade


    #cada linha tem dv uns e cada coluna tem dc uns


def find_dc_dv(T):
    for dv in range(1, 100):
        dc = dv/(1-T)
        if abs(dc - round(dc)) < 10**(-4):
            return dc, dv
        


class LDPC:
    def __init__(self, code, dc, dv, N):
        self.code = code
        self.dc = dc
        self.dv = dv
        self.N = N
        self.M = find_M(N, dc, dv)
        self.size = len(self.code)
        self.word_size = self.M
        self.codified_word_size = N
        self.sparse_matrix = create_graph(N, dc, dv)
    
    @property
    def divide_code(self):
        matrix = np.array([self.code[i:i+self.word_size] for i in range(0, self.size, self.word_size)])
        return matrix

    def encoder(self):
        grouped_code = self.divide_code
        encoded_code = [self.encode_word(u) for u in grouped_code]
        return np.array(encoded_code, dtype=int)

    """
    def encode_word(self, u):
        poly_u = vec_to_poly(u, self.x)
        result = self.poly_g * poly_u
        result = Poly(result.as_expr(), self.x, domain=GF(2))
        result = poly_to_vector_fixed_len(result, self.codified_word_size)
        return result
    """
    def decoder(self, received_code):
        return np.array([self.decode_word(received_code[i:i+self.codified_word_size]) for i in range(0, len(received_code), self.codified_word_size)]).flatten()

    def decode_word(self, codified_word):

        inferred_word = codified_word
        error_list = np.array([0]*self.N)
        for i in range(4):
            # Descubro quais equações não foram satisfeitas
            equation_line = inferred_word @ self.sparse_matrix
            equation_line = np.remainder(equation_line, 2)

            # Repito a suposta informação, replicando-a em N linhas
            equation_matrix = np.tile(equation_line, (self.N, 1))

            sparse_equation_matrix = csr_matrix(equation_matrix)

            # Encontro número de equações erradas que cada nó participa
            result = sparse_equation_matrix.multiply(self.sparse_matrix)

            error_list = result.sum(axis=1).flatten()

            result_array = np.zeros(self.N)

            max_value = np.amax(error_list)

            if max_value == 0:
                #print("*** COND 1 ***")
                break

            # Get the indices of elements in error_list equal to max_value
            indices = np.where(error_list == max_value)[1]

            result_array[indices] = 1

            if sum(result_array) == 0:
                #print("*** COND 2 ***")
                break

            inferred_word = inferred_word + result_array

            inferred_word = np.remainder(inferred_word, 2)

        return inferred_word