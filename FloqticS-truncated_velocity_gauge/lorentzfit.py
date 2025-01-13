
import numpy as np

np.set_printoptions(floatmode="fixed")


broadening_lorentz = 0.06                   #in eV

def lrnzfunc(xgrid, sg, freq, intens):
    ygrid = sg * intens / 2.0 / np.pi / ((xgrid - freq) ** 2 + sg ** 2 / 4.0)
    return ygrid

def optimize_code():
    sgm = broadening_lorentz
    start = 0.0
    end = 40.0
    num = 10000

    data_filename = 'spectrum.txt'
    
    freq = thisarray[thisarray[:, 3] > 0.06, 3]
    intens = thisarray[thisarray[:, 3] > 0.06, 4]

    xgrid = np.linspace(start, end, num)
    lrnzplt = np.zeros((num,))
    for freq_val, intens_val in zip(freq, intens):
        lrnzplt += lrnzfunc(xgrid, sgm, freq_val, intens_val)

    lorenz = np.vstack((xgrid, lrnzplt)).T

    output_filename = 'spec_broadened.txt'
    print(output_filename)
    np.savetxt(output_filename, lorenz, fmt='%9.5f')

optimize_code() 
        
