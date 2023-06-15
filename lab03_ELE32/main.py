import math
from hamming import Hamming
import numpy as np
from bsc import BSC
from no_encode import NoEncode
from chebyshev import find_vector_length2
import matplotlib.pyplot as plt

from ldpc import find_dc_dv, generate_csv, LDPC

min_p = 0.00001
error_margin = 0.05
error_freq = 0.05

if __name__ == "__main__":
    p_list = [0.5]
    last_num = 0.5
    for i in range(14):
        power = pow(10,round(i/3)+1)
        last_num = math.floor(last_num*power/2)/power
        p_list.append(last_num)
    print(p_list)
    N_list = [98, 203, 497]

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

    prob_list = [[], [], [], [], []]
    for j in range(len(p_list)):
        p = p_list[j]
        print(f"p = {p}")
        bsc = BSC(p, transmited_code)
        received_code = bsc.transform_code()
        inferred_info = []
        num_of_errors = []
        error_prob = []

        for i in range(len(test_list)):
            print(f"i: {i}")
            inferred_info.append(test_list[i].decoder(received_code))
            num_of_errors.append(np.sum(inferred_info[i]))
            #print(inferred_info[i])
 
            error_prob.append(num_of_errors[i]/len(inferred_info[i]))
            
            prob_list[i].append(error_prob[i])
        
        prob_list[4][j] = p
        
    #print(prob_list)

    plt.figure()
    plt.plot(p_list, prob_list[3], label="Hamming")
    plt.scatter(p_list, prob_list[3])
    plt.plot(p_list, prob_list[4], label="No Encode")
    plt.scatter(p_list, prob_list[4])
    plt.plot(p_list, prob_list[0], label=f"LDPC N={N_list[0]}")
    plt.scatter(p_list, prob_list[0])
    plt.plot(p_list, prob_list[1], label=f"LDPC N={N_list[1]}")
    plt.scatter(p_list, prob_list[1])
    plt.plot(p_list, prob_list[2], label=f"LDPC N={N_list[2]}")
    plt.scatter(p_list, prob_list[2])
    plt.gca().invert_xaxis()
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Incidence of errors for different p values')
    plt.xlabel('p')
    plt.ylabel("Pb")
    plt.legend()
    plt.savefig("lab3.png")


