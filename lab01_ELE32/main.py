import numpy as np
import matplotlib.pyplot as plt
import math
from bsc import BSC
from hamming import Hamming
from no_encode import NoEncode

vector_length = 1000000

if __name__ == '__main__':
    p_list = [0.5]
    last_num = 0.5
    for i in range(15):
        power = pow(10,round(i/3)+1)
        last_num = math.floor(last_num*power/2)/power
        p_list.append(last_num)
    print(p_list)

    prob_list_hamming = []
    prob_list_no_encode = []
    for p in p_list:
        information = np.random.randint(2, size=vector_length)
        hamming = Hamming(information)
        no_encode = NoEncode(information)
        hamming_code = hamming.encoder()
        no_encode_code = no_encode.code
        bsc_hamming = BSC(p, hamming_code)
        bsc_no_encode = BSC(p, no_encode)
        received_code_hamming = bsc_hamming.transform_code()
        received_code_no_encode = bsc_no_encode.transform_code()
        inferred_code_hamming = hamming.decoder(received_code_hamming)
        inferred_info_hamming = inferred_code_hamming.flatten()
        inferred_info_no_encode = received_code_no_encode
        differences_hamming = np.bitwise_xor(np.array(inferred_info_hamming, dtype=np.int64), information)
        differences_no_encode = np.bitwise_xor(np.array(inferred_info_no_encode, dtype=np.int64), information)
        num_of_errors_hamming = np.sum(differences_hamming)
        num_of_errors_no_encode = np.sum(differences_no_encode)
        error_prob_hamming = num_of_errors_hamming/vector_length
        error_prob_no_encode = num_of_errors_no_encode/vector_length
        prob_list_hamming.append(error_prob_hamming)
        prob_list_no_encode.append(error_prob_no_encode)
    
    print(prob_list_hamming)
    plt.figure()
    plt.plot(p_list, prob_list_hamming)
    plt.plot(p_list, prob_list_no_encode)
    plt.gca().invert_xaxis()
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig("mygraph.png")