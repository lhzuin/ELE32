import numpy as np
from random import randint
import pandas as pd



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
        self.v2c = np.zeros((self.N, self.M))
        self.c2v = np.zeros((self.N, self.M))
        self.k = (self.dc-self.dv)/self.dc
        self.name = f'LDPC N={N}'
        self.is_binary = False
    
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
    """
    def decoder(self, received_code):
        return np.array([self.decode_word(received_code[i:i+self.codified_word_size]) for i in range(0, len(received_code), self.codified_word_size)]).flatten()

    def decode_word(self, codified_word):
        inferred_word = codified_word
        max_num_of_iterations = 20

        for _ in range(max_num_of_iterations):
            self.v_iteration(codified_word)
            if self.c_iteration():
                break
        sum_lines = np.sum(self.c2v, axis=1)
        inferred_word = np.transpose(sum_lines + codified_word)
        # Create a boolean mask for positive and negative values, including zero
        positive_mask = inferred_word >= 0
        negative_mask = inferred_word < 0

        # Use the mask to set positive values to 0 and negative ones to 1
        inferred_word[positive_mask] = 0
        inferred_word[negative_mask] = 1

        return inferred_word

    
    def v_iteration(self, L_array):
        sum_lines = np.sum(self.c2v, axis=1)
        temp = np.transpose(sum_lines + L_array) # corresponde a estimativas para decisão (positivo ou negativo)
        temp = temp.reshape(-1, 1)
        #print(temp)
        self.v2c = np.tile(temp, (1, self.M)) - self.c2v

    def c_iteration(self):
        N, M = self.v2c.shape
        self.c2v = np.zeros((N, M))
        is_valid_code = True

        for j in range(M):
            # Select column j
            column = self.v2c[:, j]
            
            # Compute the product of non-zero elements for each line
            non_zero_elems = column[column != 0]
            product = np.prod(non_zero_elems)
            
            if product < 0:
                is_valid_code = False
            
            for i in range(N):
                # Compute minimum absolute non-zero value in column excluding i-th row
                min_abs_val = np.min(np.abs(column[np.where((np.arange(N) != i) & (column != 0))]))
                
                # Compute sign based on the multiplication of all non null elements of column j that are not in line i
                if column[i] != 0:
                    self.c2v[i, j] = np.sign(product/column[i]) * min_abs_val
                else:
                    self.c2v[i, j] = 0


        return is_valid_code

    def get_error_prob(self, received_code):
        inferred_info = self.decoder(received_code)
        num_of_errors = np.sum(inferred_info)
        error_prob = num_of_errors/self.size
        return error_prob