
import os
import ctypes
from ctypes import CDLL, RTLD_GLOBAL, c_double, c_size_t
import numpy as np
import time


if __name__ == '__main__':
    lib_path = '/home/daniele/dynasty/dynasty.so'
    try:
        lib = CDLL(lib_path, mode=RTLD_GLOBAL)
    except:
        print('Cannot load library from {}.'.format(lib_path))

    R = c_double(np.log(0.07))
    E = c_double(2.1)
    B = c_double(0.17)
    D = c_double(0.42)
    G = c_double(0.09)
    H = c_double(0.1)
    Q = c_double(0.4)
    ttran = c_double(1e4)
    tstop = c_double(1e5)
    nev = c_size_t(256)
    
    start = time.time()
    lib.integrate(R, E, B, D, G, H, Q, ttran, tstop, nev)
    stop = time.time()

    print('Elapsed time = {:.3f} s'.format(stop - start))
    
