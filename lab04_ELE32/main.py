from hamming import Hamming
import numpy as np
from bpsk import BPSK
from llr_ldpc_2 import LlrLDPC, find_L, find_dc_dv, generate_csv
from bit_flipping_ldpc import BitFlippingLDPC
from no_encode import NoEncode
from chebyshev import find_vector_length2
import matplotlib.pyplot as plt
from cyclic_encode import CyclicEncode

min_p = 0.00001
error_margin = 0.05
error_freq = 0.05

if __name__ == "__main__":
    Ei_N0_db_list = [x*0.5 for x in range(11)]
    Ei_N0_list = [10**(x/10) for x in Ei_N0_db_list]
    Ei =1
    N0_list = [Ei/(x) for x in Ei_N0_list]

    print(N0_list)
    N_list_ldpc = [1001]

    n_cyclic = 17
    k_cyclic = 9
    N_list = N_list_ldpc + [n_cyclic]
    print(N_list)

    T = 4/7
    dc, dv = find_dc_dv(T)
    vector_length = find_vector_length2(min_p, error_margin, error_freq, N_list)

    
    transmited_code = np.array([0]*vector_length)
    ldpc_list = []
    for i in range(len(N_list_ldpc)):
        ldpc_list.append(LlrLDPC(transmited_code, dc, dv, N_list_ldpc[i], f"sparse_matrix_N_{N_list_ldpc[i]}.csv"))
        ldpc_list.append(BitFlippingLDPC(transmited_code, dc, dv, N_list_ldpc[i], f"sparse_matrix_N_{N_list_ldpc[i]}.csv"))
        #generate_csv(ldpc_list[i].sparse_matrix, N_list_ldpc[i])
    
    #transmited_code = np.zeros(vector_length*N)
    test_list = ldpc_list + [Hamming(transmited_code), CyclicEncode(transmited_code, n_cyclic, k_cyclic), NoEncode(transmited_code)]
    print(test_list)

    prob_list = [[] for _ in range(len(test_list))]
    for j in range(len(N0_list)):
        N0 = N0_list[j]
        Eb_list = [Ei*test_list[i].r for i in range(len(test_list))] 
        print(f"Ei/N0 = {Ei_N0_db_list[j]}")
        bpsk_list = [BPSK(N0, Eb_list[i], transmited_code) for i in range(len(test_list))]
        received_code_list = [bpsk_list[i].transform_code(test_list[i].is_binary) for i in range(len(test_list))]
        L_list_of_arrays = [find_L(received_code, N0) for received_code in received_code_list]
        print(L_list_of_arrays[0])
        inferred_info = []
        num_of_errors = []
        error_prob = []

        for i in range(len(test_list)):
            print(f"i: {i}")
            decode_input = L_list_of_arrays[i] if test_list[i].receive_L else received_code_list[i]
            inferred_info.append(test_list[i].decoder(decode_input)) # len(codified_code)*r
            num_of_errors.append(np.sum(inferred_info[i])) # qtd de bits errados
            #print(inferred_info[i])
 
            error_prob.append(num_of_errors[i]/len(inferred_info[i])) 
            
            prob_list[i].append(error_prob[i])
        
        #prob_list[4][j] = p
        
    #print(prob_list)
    print(prob_list[0])
    plt.figure()
    for i in range(len(test_list)):
        plt.plot(Ei_N0_db_list, prob_list[i], label=test_list[i].name)
        plt.scatter(Ei_N0_db_list, prob_list[i])
    #plt.gca().invert_xaxis()
    #plt.xscale('log')
    plt.yscale('log')
    plt.title('Incidence of errors for different p values')
    plt.xlabel('Ei/N0')
    plt.ylabel("Pb")
    plt.legend()
    plt.savefig("lab4.png")

    with open('lab4.txt', 'w') as f:
        for sublist in prob_list:
            # Join each integer in the sublist to a string, separated by a space, and write it to the file
            f.write(' '.join(map(str, sublist)) + '\n')



