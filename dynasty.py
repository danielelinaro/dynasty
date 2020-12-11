
import os
import time
import ctypes
from ctypes import CDLL, RTLD_GLOBAL, c_double, c_size_t, byref
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    lib_path = '/home/daniele/dynasty/dynasty.so'
    try:
        lib = CDLL(lib_path, mode=RTLD_GLOBAL)
    except:
        print('Cannot load library from {}.'.format(lib_path))

    n_eq = 3
    n_pars = 7

    y0 = (c_double * n_eq)()
    y0[0] = 0.8
    y0[1] = 0.1
    y0[2] = 0.1

    R = 0.1
    E = 2.3
    B = 0.17
    D = 0.42
    G = 0.09
    H = 0.1
    Q = 0.4

    pars = (c_double * n_pars)()
    pars[0] = R
    pars[1] = E
    pars[2] = B
    pars[3] = D
    pars[4] = G
    pars[5] = H
    pars[6] = Q

    t_tran = c_double(1e4)
    t_stop = c_double(1e5)
    n_ev = c_size_t(1024)

    atol = (c_double * n_eq)()
    for i in range(n_eq):
        atol[i] = 1e-8
    rtol = c_double(1e-10)

    sol = (c_double * n_ev.value * n_eq)()
    
    start = time.time()
    n_ev = lib.integrate(pars, y0, t_tran, t_stop, n_ev, atol, byref(rtol), sol)
    stop = time.time()

    print('Elapsed time = {:.3f} s'.format(stop - start))

    data = np.array(sol)
    data = np.reshape(data, (data.shape[1], data.shape[0]), 'C')
    data = data[:n_ev,:]

    fig,ax = plt.subplots(1, 1, figsize=(4, 4))
    ax.plot(data[:,2], data[:,1], 'ko', markersize=3, markerfacecolor='w', markeredgewidth=0.5)
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    ax.set_xlabel(r'$\mathrm{x}_3$')
    ax.set_ylabel(r'$\mathrm{x}_2$')
    ax.grid('on')
    fig.tight_layout()
    plt.show()
    
