from math import lcm, ceil
from functools import reduce

def find_vector_length(min_p, error_margin, error_freq, k_list):
    var = 1
    mmc = reduce(lcm, k_list)
    return 500*mmc


def find_vector_length2(min_p, error_margin, error_freq, k_list):
    num = (1-min_p)/(min_p*error_freq*error_margin**2) 
    print(num)
    mmc = reduce(lcm, k_list)
    print(mmc)
    num = ceil(num/mmc)*mmc

    return num



if __name__ == "__main__":
    min_p = 0.000001
    error_margin = 0.075
    error_freq = 0.075
    k_hamming = 4
    n_cyclic = 17
    k_cyclic = 9
    print(find_vector_length2(min_p, error_margin, error_freq, [k_hamming, k_cyclic]))