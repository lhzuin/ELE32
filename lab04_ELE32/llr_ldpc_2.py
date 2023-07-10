import numpy as np
#import cupy as np
from random import randint
import pandas as pd
from numba import njit



def find_L(r, N0):
    sigma = N0/2
    L = 2*r/(sigma**2)
    return L



def find_M(N, dc, dv):
    M = N*dv/dc
    #print(M)
    M = int(M)
    #print(M)
    return M

def create_graph(N, dc, dv):
    M = find_M(N, dc, dv)
    dc_list = np.zeros(M, dtype=int)
    dv_list = np.zeros(N, dtype=int)

    min_dv = 0
    min_dc = 0

    num_cols = M
    num_rows = N

    array_matrix = np.zeros((num_rows, num_cols))

    while(min_dc < dc and min_dv < dv):

        temp_dc = np.where(dc_list == min_dc)[0]
        temp_dv = np.where(dv_list == min_dv)[0]

        random_c = np.random.choice(temp_dc)
        random_v = np.random.choice(temp_dv)

        if array_matrix[random_v, random_c] != 1:
            array_matrix[random_v, random_c] = 1
            dc_list[random_c] += 1
            dv_list[random_v] += 1

        min_dc = np.min(dc_list)
        min_dv = np.min(dv_list)
    
    return array_matrix 

    


def generate_csv(sparse_matrix, N):
    dense_matrix = sparse_matrix

    # Create a pandas DataFrame
    df = pd.DataFrame(dense_matrix)

    # Save DataFrame to CSV
    df.to_csv(f"sparse_matrix_N_{N}.csv", index=False, header=False)





    # Crio coluna por coluna de verificação de paridade


    #cada linha tem dv uns e cada coluna tem dc uns


def find_dc_dv(T):
    for dv in range(1, 100):
        dc = dv/(1-T)
        if abs(dc - round(dc)) < 10**(-4):
            return dc, dv
        

@njit
def v_iteration(N, M, L_array, c2v, sparse_matrix):
    sum_lines = np.sum(c2v, axis=1)
    temp = sum_lines + L_array
    temp = temp.reshape(-1, 1)
    v2c = temp - c2v
    v2c = np.multiply(v2c, sparse_matrix)

    return v2c

@njit
def c_iteration(N, M, v2c, c2v):
    is_valid_code = True

    for j in range(M):
        column = v2c[:, j]
        non_zero_elems = column[np.nonzero(column)]
        product = np.sign(non_zero_elems).prod()
        if product < 0:
            is_valid_code = False

        for i in range(N):
            if column[i] != 0:
                min_abs_val = np.min(np.abs(column[(np.arange(N) != i) & (column != 0)]))
                c2v[i, j] = np.sign(product*column[i]) * min_abs_val

    return is_valid_code

class LlrLDPC:
    def __init__(self, codified_code, dc, dv, N, sparse_matrix_name = None):
        self.codified_code = codified_code
        self.dc = dc
        self.dv = dv
        self.N = N
        self.M = find_M(N, dc, dv)
        self.size = len(self.codified_code)
        self.word_size = self.M
        self.codified_word_size = N
        if sparse_matrix_name is None:
            self.sparse_matrix = create_graph(N, dc, dv)
        else:
            self.sparse_matrix = np.genfromtxt(sparse_matrix_name, delimiter=',')
        #print(self.sparse_matrix)
        self.v2c = np.zeros((self.N, self.M))
        self.c2v = np.zeros((self.N, self.M))
        self.r = (self.dc-self.dv)/self.dc # taxa
        self.k = int(self.N * self.r)
        self.name = f'LDPC LLR N={N}'
        self.is_binary = False
        self.receive_L = True

    def decoder(self, received_code):
        return np.array([self.decode_word(received_code[i:i+self.codified_word_size]) for i in range(0, len(received_code), self.codified_word_size)]).flatten()

    def decode_word(self, L_initial):
        self.v2c.fill(0)# np.zeros((self.N, self.M))
        self.c2v.fill(0) #np.zeros((self.N, self.M))
        #inferred_word = codified_word
        max_num_of_iterations = 5#20
        #print(L_initial)

        for _ in range(max_num_of_iterations):
            self.v_iteration(L_initial)
            if self.c_iteration():
                break
        sum_lines = np.sum(self.c2v, axis=1)
        inferred_word = np.transpose(sum_lines + L_initial)
        # Create a boolean mask for positive and negative values, including zero
        positive_mask = inferred_word >= 0
        negative_mask = inferred_word < 0

        # Use the mask to set positive values to 0 and negative ones to 1
        inferred_word[positive_mask] = 0
        inferred_word[negative_mask] = 1

        return inferred_word

    
    def v_iteration(self, L_array):
        self.v2c = v_iteration(self.N, self.M, L_array, self.c2v, self.sparse_matrix)
        """
        sum_lines = np.sum(self.c2v, axis=1) # Array com N elementos
        #print(sum_lines)
        temp = np.transpose(sum_lines + L_array) # corresponde a estimativas para decisão (positivo ou negativo)
        temp = temp.reshape(-1, 1)
        self.v2c = np.tile(temp, (1, self.M)) - self.c2v
        self.v2c = np.multiply(self.v2c, self.sparse_matrix)
        #print(self.v2c)
        """

    #@njit
    def c_iteration(self):
        #v_iteration(self.sum_lines, self.c2v, self.v2c, L_array, self.sparse_matrix)
        return c_iteration(self.N, self.M, self.v2c, self.c2v)
        is_valid_code = True

        for j in range(self.M):
            column = self.v2c[:, j]
            non_zero_elems = column[column != 0]
            product = np.prod(np.sign(non_zero_elems))
            
            if product < 0:
                is_valid_code = False

            for i in range(self.N):
                if column[i] != 0:
                    column_temp = column.copy()
                    column_temp[i] = 0
                    min_abs_val = np.min(np.abs(column_temp[column_temp != 0]))
                    self.c2v[i, j] = np.sign(product*column[i]) * min_abs_val
                else:
                    self.c2v[i, j] = 0

        return is_valid_code
        

    def get_error_prob(self, received_code):
        inferred_info = self.decoder(received_code)
        num_of_errors = np.sum(inferred_info)
        error_prob = num_of_errors/self.size
        return error_prob