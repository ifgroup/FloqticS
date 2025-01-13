
import numpy as np

np.set_printoptions(linewidth=np.inf)
np.set_printoptions(suppress=True)
np.set_printoptions(floatmode="fixed")


def lrnzfunc(xgrid, sg, freq, intens):
    ygrid = sg * intens / 2.0 / np.pi / ((xgrid - freq) ** 2 + sg ** 2 / 4.0)
    return ygrid

def optimize_code():
    sgm = 0.04
    start = 0.0
    end = 40.0
    num = 10000

    for ka in range(0, 4):
        data_filename = 'res' + str(ka) + '.txt'
        print(data_filename)
        
        try:
            thisarray = np.loadtxt(data_filename)
        except:
            print('cannot read'+'data_filename')
            continue

        freq = thisarray[thisarray[:, 3] > 0.06, 3]
        intens = thisarray[thisarray[:, 3] > 0.06, 4]

        intens = intens*25/500

        xgrid = np.linspace(start, end, num)
        lrnzplt = np.zeros((num,))
        for freq_val, intens_val in zip(freq, intens):
            lrnzplt += lrnzfunc(xgrid, sgm, freq_val, intens_val)

        lorenz = np.vstack((xgrid, lrnzplt)).T

        output_filename = 'gauss' + str(ka) + '.txt'
        print(output_filename)
        np.savetxt(output_filename, lorenz, fmt='%9.5f')

optimize_code() 
        
