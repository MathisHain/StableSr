# Script to run simple Sr forward model ensembles
# Want to learn more? mhain@ucsc.edu
 
import numpy as np
import copy
import matplotlib.pyplot as mpl

class state():
    def __init__(self,Fin,Imb):    
        self.t = 0
        self.eps = -0.22
        self.Fin = Fin
        self.Fout = Fin - Imb
        self.Sr = 88e-6 * 1.4e21
        self.din = 0.345
        self.d = self.din -self.eps

class model():
    def __init__(self,Fin,Imb):
        self.now = state(Fin,Imb)
        self.nxt = state(Fin,Imb)
        self.dt = 1000 # one step is 1000 years
        self.out = output(self.now)
    
    def Set(self,Fin,Imb):
        self.now.Fin = Fin
        self.now.Fout = Fin - Imb
    
    def timestep(self):
        self.nxt.t = self.now.t + self.dt
        self.nxt.Sr = self.now.Sr + (self.now.Fin-self.now.Fout)*self.dt
        self.nxt.d =  (self.now.Sr*self.now.d + (self.now.Fin*self.now.din-self.now.Fout*(self.now.d+self.now.eps))*self.dt)/self.nxt.Sr
        self.nxt.Fin = self.now.Fin
        self.nxt.Fout = self.now.Fout
        self.now = copy.deepcopy(self.nxt)
        
    def reportoutput(self):
        self.out.t = np.append(self.out.t,self.now.t)
        self.out.Sr = np.append(self.out.Sr,self.now.Sr)
        self.out.d = np.append(self.out.d,self.now.d)
        self.out.Fin = np.append(self.out.Fin,self.now.Fin)
        self.out.Fout = np.append(self.out.Fout,self.now.Fout)
        
    def RunModel(self,steps):
        reportinterval = 25
        for step in range(0,int(steps/reportinterval)):
            for i in range(0,reportinterval):
                self.timestep()
                
            self.reportoutput()
        
class output():
    def __init__(self, now):   
        self.t = np.array(round(now.t,2))
        self.Sr = np.array(now.Sr)
        self.d = np.array(now.d)
        self.Fin = np.array(now.Fin)
        self.Fout = np.array(now.Fout)
        
    def Show(self):
        mpl.figure(1)
        mpl.subplot(311)
        mpl.plot(self.t, self.Sr)
        mpl.grid()
        mpl.ylabel("Sr")
        
        mpl.subplot(312)
        mpl.plot(self.t, self.d)
        mpl.grid()
        mpl.ylabel("d")
        
        mpl.subplot(313)
        mpl.plot(self.t, self.Fin)
        mpl.plot(self.t, self.Fout)
        mpl.grid()
        mpl.ylabel("F")
        
        mpl.show()

class ensemble():
    def __init__(self,N,Sr0,Fin,Imb,id):
        self.run = []
        self.Fin = Fin
        self.Imb = Imb
        self.Sr0 = Sr0
        self.id = id
        self.N = N
        
        if id == "ex1":
            print("Running ex1")
            for n in range(0,self.N):
                self.run.append(model((n/5)*Fin,0))
                self.run[n].RunModel(100)
                self.run[n].Set((n/5)*Fin,Imb)
                self.run[n].RunModel(2000)
                self.run[n].Set((n/5)*Fin,0)
                self.run[n].RunModel(7900)
                print("Done with run {}.".format(n))
                
        if id == "ex2":
            print("Running ex2")
            for n in range(0,self.N):
                nn = n-3
                self.run.append(model(Fin,0))
                self.run[n].RunModel(100)
                self.run[n].Set(Fin,(nn/3*Fin/2))
                self.run[n].RunModel(2000)
                self.run[n].Set(Fin,0)
                self.run[n].RunModel(7900)
                print("Done with run {}.".format(n))
                
        if id == "ex3":
            print("Running ex3")
            for n in range(0,self.N):
                nn = n-3
                self.run.append(model(Fin,0))
                self.run[n].RunModel(100)
                self.run[n].Set(Fin+(nn/3*Fin/2),(nn/3*Fin/2))
                self.run[n].RunModel(2000)
                self.run[n].Set(Fin,0)
                self.run[n].RunModel(7900)
                print("Done with run {}.".format(n))
        
    
    def ShowAll(self,Fig,title,filename):
        
        mpl.figure(Fig)
        mpl.suptitle(title)
        mpl.subplot(312)
        mpl.grid()
        mpl.ylabel("Sr inventory [Pmol]")
        for n in range(0,self.N):
            mpl.plot(self.run[n].out.t/1e6, self.run[n].out.Sr/1e15,'-k')
    
        mpl.subplot(313)
        mpl.grid()
        mpl.ylabel("$δ^{88/86}Sr$")
        mpl.xlabel("Time [Myrs]")
        for n in range(0,self.N):
            mpl.plot(self.run[n].out.t/1e6, self.run[n].out.d,'k')
    
        mpl.subplot(311)
        mpl.grid()
        mpl.ylabel(r'$F_{IN}, F_{OUT}  [Gmol/yr]$')
        for n in range(0,self.N):
            mpl.plot(self.run[n].out.t/1e6, self.run[n].out.Fin/1e9,'-b')
            mpl.plot(self.run[n].out.t/1e6, self.run[n].out.Fout/1e9,'--r')
        
        mpl.savefig(filename, format='ps')
        
        #mpl.show()      

############################################################################### 
ex1 = ensemble(6,None,38e9,-0.1*38e9,"ex1")
ex1.ShowAll(11,"Constant imbalance, Variable throughflow","FigSuppInf_constImbalance.ps")

ex2 = ensemble(7,None,38e9,None,"ex2")
ex2.ShowAll(12,"Constant $F_{IN}$, Variable $F_{OUT}$","FigSuppInf_constIN.ps")

ex3 = ensemble(7,None,38e9,None,"ex3")
ex3.ShowAll(13,"Constant $F_{OUT}$, Variable $F_{IN}$","FigSuppInf_constOUT.ps") 