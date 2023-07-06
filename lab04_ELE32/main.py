import math
from hamming import Hamming
import numpy as np
from bpsk import BPSK
from ldpc import LDPC, find_L
from no_encode import NoEncode
from chebyshev import find_vector_length2
import matplotlib.pyplot as plt

from ldpc import find_dc_dv, generate_csv, LDPC

min_p = 0.00001
error_margin = 0.05
error_freq = 0.05

if __name__ == "__main__":
    Ei_N0_list = [x*0.5 for x in range(11)]
    Ei =1
    N0_list = [Ei/(10**(x/10)) for x in Ei_N0_list]

    print(N0_list)
    N_list = [98]

    T = 4/7
    dc, dv = find_dc_dv(T)
    vector_length = find_vector_length2(min_p, error_margin, error_freq, N_list)
    information_length = vector_length*T

    #sparse_matrix_list = []
    
    information = np.array([0]*int(information_length))
    transmited_code = np.array([0]*vector_length)
    ldpc_list = []
    for i in range(len(N_list)):
        ldpc_list.append(LDPC(transmited_code, dc, dv, N_list[i]))
        generate_csv(ldpc_list[i].sparse_matrix, N_list[i])
    
    #transmited_code = np.zeros(vector_length*N)
    test_list = ldpc_list + [Hamming(information), NoEncode(information)]
    print(test_list)

    prob_list = [[], [], []]
    for j in range(len(N0_list)):
        N0 = N0_list[j]
        Eb_list = [test_list[i].k for i in range(len(test_list))]
        #Eb = test_list[j].k
        print(f"Ei/N0 = {Ei_N0_list[j]}")
        bpsk_list = [BPSK(N0, Eb_list[i], transmited_code) for i in range(len(test_list))]
        received_code_list = [bpsk_list[i].transform_code(test_list[i].is_binary) for i in range(len(test_list))]
        L_list_of_arrays = [find_L(received_code, N0) for received_code in received_code_list]
        inferred_info = []
        num_of_errors = []
        error_prob = []

        for i in range(len(test_list)):
            print(f"i: {i}")
            inferred_info.append(test_list[i].decoder(received_code_list[i]))
            num_of_errors.append(np.sum(inferred_info[i]))
            #print(inferred_info[i])
 
            error_prob.append(num_of_errors[i]/len(inferred_info[i])) 
            
            prob_list[i].append(error_prob[i])
        
        #prob_list[4][j] = p
        
    #print(prob_list)

    plt.figure()
    for i in range(len(test_list)):
        plt.plot(Ei_N0_list, prob_list[i], label=test_list[i].name)
        plt.scatter(Ei_N0_list, prob_list[i])
    plt.gca().invert_xaxis()
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Incidence of errors for different p values')
    plt.xlabel('Ei/N0')
    plt.ylabel("Pb")
    plt.legend()
    plt.savefig("lab4.png")


