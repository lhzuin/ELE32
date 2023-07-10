from math import lcm, ceil
from functools import reduce

def find_vector_length(min_p, error_margin, error_freq, k_list):
    var = 1
    mmc = reduce(lcm, k_list)
    #return mmc#10000008#1000*mmc
    return 980


def find_vector_length2(min_p, error_margin, error_freq, k_list):
    num = (1-min_p)/(min_p*error_freq*error_margin**2) 
    num = min(num, 10**7)
    num = min(num, 10**6)
    #print(num)
    mmc = reduce(lcm, k_list)
    #print(mmc)
    num = ceil(num/mmc)*mmc
    #num = 98000
    print(f"vector length: {num}")
    return num



if __name__ == "__main__":
    min_p = 0.00001
    error_margin = 0.075
    error_freq = 0.075
    k_hamming = 4
    n_cyclic = 17
    k_cyclic = 9
    print(find_vector_length2(min_p, error_margin, error_freq, [k_hamming, k_cyclic]))