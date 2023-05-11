
import numpy as np

np.set_printoptions(linewidth=np.inf)
np.set_printoptions(suppress=True)
np.set_printoptions(floatmode="fixed")

sgm=0.06
start=0.0
end=80.0
num=40000

filename='spec.txt'
with open(filename) as afile:
    thislist = []
    for line in afile:
        if "*" in line:
            continue
        data = np.array([float(i) for i in line.split()])
        thislist.append(data)

thisarray=np.array(thislist)

freq=thisarray[thisarray[:,3] > 0.03 ,3]
intens=thisarray[thisarray[:,4] > 0.03 ,4]

final=np.vstack((freq,intens)).T


def lrnzfunc(sg,fq,it):
    xgrid=np.linspace(start,end,num)
    ygrid=sg*it/2.0/np.pi/((xgrid-fq)**2+sg**2/4.0)
    return ygrid

lrnzplt=np.zeros((num,))
for val1,val2 in final:
    lrnzplt+=lrnzfunc(sgm,val1,val2)


xgr=np.linspace(start,end,num)
lorenz=np.vstack((xgr,lrnzplt/100.0)).T

 
filename='absorptionspectrum'+str(ka)+'.txt'
np.savetxt(filename,lorenz,fmt='%15.5f')
    
        
