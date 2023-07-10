from scipy.sparse import csr_matrix
import numpy as np
from random import randint
import pandas as pd


def find_M(N, dc, dv):
    M = N*dv/dc
    #print(M)
    M = int(M)
    #print(M)
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
        


class BitFlippingLDPC:
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
            self.sparse_matrix = create_graph(N, dc, dv)  # Assuming this function returns a csr_matrix
        else:
            dense_matrix = np.genfromtxt(sparse_matrix_name, delimiter=',')
            self.sparse_matrix = csr_matrix(dense_matrix)
        self.r = (self.dc-self.dv)/self.dc # taxa
        self.k = int(self.N * self.r)
        self.name = f'LDPC Bit Flipping N={N}'
        self.is_binary = True
        self.receive_L = False

    
    def decoder(self, received_code):
        #print(received_code.shape)
        return np.array([self.decode_word(received_code[i:i+self.codified_word_size]) for i in range(0, len(received_code), self.codified_word_size)]).flatten()

    def decode_word(self, codified_word):
        #print(codified_word.shape)
        inferred_word = codified_word
        error_list = np.zeros(self.N)

        for _ in range(20):
            # Descubro quais equações não foram satisfeitas
            equation_line = inferred_word @ self.sparse_matrix
            equation_line %= 2

            # Repito a suposta informação, replicando-a em N linhas
            equation_matrix = np.tile(equation_line, (self.N, 1))
            sparse_equation_matrix = csr_matrix(equation_matrix)

            # Encontro número de equações erradas que cada nó participa
            result = sparse_equation_matrix.multiply(self.sparse_matrix)

            error_list = result.sum(axis=1).flatten()

            if np.amax(error_list) == 0:
                break

            # Get the indices of elements in error_list equal to max_value
            indices = np.where(error_list == np.amax(error_list))[1]

            inferred_word[indices] += 1
            inferred_word %= 2

        return inferred_word
    
    def get_error_prob(self, received_code):
        inferred_info = self.decoder(received_code)
        num_of_errors = np.sum(inferred_info)
        error_prob = num_of_errors/self.size
        return error_prob