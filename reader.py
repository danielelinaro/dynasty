
import numpy as np

filename = 'dynasty.dat'

with open(filename,'rb') as fid:
    for binary_line in fid:
        line = binary_line.decode('ascii').rstrip()
        if 'NEQ:' in line:
            neq = int(line.split(' ')[1])
        elif 'NPARS:' in line:
            npars = int(line.split(' ')[1])
        elif 'NEV:' in line:
            nev = int(line.split(' ')[1])
        elif 'NINT:' in line:
            nint = int(line.split(' ')[1])
        elif 'Binary:' in line:
            break

    print('neq   =', neq)
    print('npars =', npars)
    print('nev   =', nev)
    print('nint  =', nint)

    cnt = 1

    for i in range(nint):
        if cnt != 0:
            num = np.fromfile(fid, np.int64, count=1)

        if num != 0:
            import ipdb
            ipdb.set_trace()

        pars = np.fromfile(fid, np.float64, count=npars)
        #print(pars)

        for j in range(nev):
            cnt = np.fromfile(fid, np.int64, count=1)
            if cnt == 0:
                break
            x = np.fromfile(fid, np.float64, count=neq+1)
            print('[{:03d}]  {:8.2f}  {:7.5f}  {:7.5f}  {:7.5f}' \
                  .format(j+1, x[0], x[1], x[2], x[3]))

                          
