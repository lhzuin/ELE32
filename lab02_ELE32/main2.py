import numpy as np
import matplotlib.pyplot as plt
import math
from bsc import BSC
import math
from generator_polys import generator_polynomials, find_best_generator
from find_minimal_distance import get_minimum_distance
from cyclic_encode import CyclicEncode
from hamming import Hamming
from no_encode import NoEncode
import numpy as np
from bsc import BSC
from chebyshev import find_vector_length

vector_length = 1000000
min_p = 0.000001
error_margin = 0.05
error_freq = 0.05


if __name__ == '__main__':
    p_list = [0.5]
    last_num = 0.5
    for i in range(12):
        power = pow(10,round(i/3)+1)
        last_num = math.floor(last_num*power/2)/power
        p_list.append(last_num)
    print(p_list)

    prob_list_hamming = []
    prob_list_no_encode = []
    prob_list_cyclic = []

    k_hamming = 4
    n_cyclic = 17
    k_cyclic = 9

    gen_polys = generator_polynomials(n_cyclic, k_cyclic)
    g = find_best_generator(gen_polys)
    print(f"g: {g}")

    vector_length = find_vector_length(min_p, error_margin, error_freq, [k_hamming, k_cyclic])
    info = np.zeros(vector_length)

    hamming = Hamming(info)
    no_encode = NoEncode(info)
    cyclic = CyclicEncode(info, n_cyclic, g)

    hamming_code = hamming.encoder()
    no_encode_code = no_encode.encoder()
    cyclic_code = cyclic.encoder()

    for p in p_list:
        print(f"p = {p}")

        bsc_hamming = BSC(p, hamming_code)
        bsc_no_encode = BSC(p, no_encode_code)
        bsc_cyclic = BSC(p, cyclic_code)
 
        received_code_hamming = bsc_hamming.transform_code()
        received_code_no_encode = bsc_no_encode.transform_code()
        received_code_cyclic = bsc_cyclic.transform_code()
      
        inferred_info_hamming = hamming.decoder(received_code_hamming)
        inferred_info_no_encode = no_encode.decoder(received_code_no_encode)
        inferred_info_cyclic = cyclic.decoder(received_code_cyclic)

        num_of_errors_hamming = np.sum(inferred_info_hamming)
        num_of_errors_no_encode = np.sum(inferred_info_no_encode)
        num_of_errors_cyclic = np.sum(inferred_info_cyclic)
 
        error_prob_hamming = num_of_errors_hamming/vector_length
        error_prob_no_encode = num_of_errors_no_encode/vector_length
        error_prob_cyclic = num_of_errors_cyclic/vector_length

        prob_list_hamming.append(error_prob_hamming)
        prob_list_no_encode.append(error_prob_no_encode)
        prob_list_cyclic.append(error_prob_cyclic)
        
    for i in range(len(p_list)):
        print(f"p: {p_list[i]} -> hamming: {prob_list_hamming[i]}")
    for i in range(len(p_list)):
        print(f"p: {p_list[i]} -> cyclic: {prob_list_cyclic[i]}")
    print(prob_list_hamming)
    print(prob_list_cyclic)
    plt.figure()
    plt.plot(p_list, prob_list_hamming, label="Hamming")
    plt.scatter(p_list, prob_list_hamming)
    plt.plot(p_list, prob_list_no_encode, label="No Encode")
    plt.scatter(p_list, prob_list_no_encode)
    plt.plot(p_list, prob_list_cyclic, label="Cyclic (n=17, k=9)")
    plt.scatter(p_list, prob_list_cyclic)
    plt.gca().invert_xaxis()
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Incidence of errors for different p values')
    plt.xlabel('p')
    plt.ylabel("Pb")
    plt.legend()
    plt.savefig("lab2.png")