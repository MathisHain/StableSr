# Inverse model of stable Sr data

import numpy as np
import matplotlib.pyplot as plt

X = np.genfromtxt('RawData.txt', dtype=float, delimiter='\t')
X = X[1:,:]

dX = np.zeros((47,2))

for row in range(len(X)-1):
    dX[row,0] = (X[row+1,0]+X[row,0])/2
    dX[row,1] = X[row+1,3]-X[row,3]


p = np.poly1d(np.polyfit(X[:,0],X[:,3], 9))
pd = np.polyder(p) 
xp = np.linspace(0, 35, 100)
error = p(X[:,0])-X[:,3]


if (0):
    plt.subplot(311)
    plt.plot(X[:,0],X[:,3],xp, p(xp), '--')

    plt.subplot(312)
    plt.plot(dX[:,0],dX[:,1],xp, pd(xp), '--')

    plt.subplot(313)
    plt.hist(error)
    plt.show()


########### MODEL 
Fin = 38e9
Sr0 = Fin*2.5e6
dsw0 = p(0)
din = 0.18
eps = -0.15

Model = np.zeros((35000,4))
Model[0,0] = 0
Model[0,1] = Sr0
Model[0,2] = dsw0
Model[0,3] = None

for kyr in range(1,35000):
    Myr = kyr/1000.0
    Model[kyr,0] = Myr
    Model[kyr,3] = (pd(Myr)/1000*Model[kyr-1,1] + (din-Model[kyr-1,2])*Fin*1000)/eps
    Model[kyr,1] = Model[kyr-1,1] - Fin*1000 + Model[kyr,3]
    Model[kyr,2] = (Model[kyr-1,2]*Model[kyr-1,1] - din*Fin*1000 + (Model[kyr-1,2]+eps)*Model[kyr,3])/Model[kyr,1]

plt.figure(2)
plt.subplot(311)
plt.plot(X[:,0],X[:,3],xp, p(xp), '--',Model[:,0],Model[:,2],'-')
plt.grid()
plt.ylabel("d88/86Sr")
plt.subplot(312)
plt.plot(Model[:,0],Model[:,1],'-')
plt.ylabel("[Sr]")
plt.grid()
plt.subplot(313)
plt.plot(Model[:,0],Fin-Model[:,3]/1000,'-')
plt.ylabel("In-OUT [mol/yr]")
plt.xlabel("MyrBP")
plt.grid()

plt.savefig("SrInverse_V0.png")
plt.show()
print(Model)    
