# Inverse model of stable Sr data

import numpy as np

class StableStrontium():
    
    def __init__(self,Fin=38e9, Sr0=88e-6 * 1.4e21, eps=-0.18, polyN=7, Nruns=30):
        
        self.Nruns = Nruns
        self.Data = self.LoadDataSc(0)
        col = 1
        (p,pd,xp,error) = self.FitPoly(self.Data, col, polyN)
        dmean = np.mean(p(xp))
        print("The mean d88Sr of the polynomial fit at sample times: ", dmean)

        dsw0 = p(0)
        print("The d88Sr of the polynomial fit at 0 Myr: ", dsw0)

        din = dmean + eps
        print("The assumed d88Sr Sr input required for secular balance: ",din)
        
        SingleModel = self.IntegrateConstFin(pd,Fin,Sr0,dsw0,din,eps)
        
        self.MCModel = np.zeros((351,7,Nruns+1))
        self.MCModelDT = np.zeros((351,7,Nruns+1))
        self.MCddiff = np.zeros(Nruns+1)
        self.MCtrend = np.zeros(Nruns+1)
        (self.MCModel[:,:,0],self.MCModelDT[:,:,0]) = self.IntegrateConstFin(pd,Fin,Sr0,dsw0,din,eps)
        self.MCddiff[0] = din - (dmean + eps)
        
        for run in range(1,Nruns+1):
            while True:
                dinE = din + 0.01*np.random.randn()
                epsE = eps + 0.005*np.random.randn()
                self.MCddiff[run] = dinE - (dmean + epsE)
                (self.MCModel[:,:,run],self.MCModelDT[:,:,run]) = self.IntegrateConstFin(pd,Fin,Sr0,dsw0,dinE,epsE)
                #Model = np.stack((Model,NewModel),axis=2)
                if (np.min(self.MCModel[:,1,run])>0):
                    break  
            print("Done with run {0}".format(run))
            
    
    def PublishedFigure(self):
        Model=self.MCModel
        ModelDT=self.MCModelDT
        Data=self.Data
        ddiff = self.MCddiff
        col = 1
        
        import matplotlib.pyplot as plt
        plt.figure(1,figsize=(10,8))
    
        plt.subplot(321)
        #for run in range(1,Nruns+1):
        #    if (ddiff[run]>0):
        #        plt.plot(Model[:,0,run],Model[:,4,run],'-r')
        #    else:
        #        plt.plot(Model[:,0,run],Model[:,4,run],'-b')
        plt.plot(Model[:,0,0],Model[:,4,0],'-k',linewidth=3.0)
        plt.plot(Data[:,0],Data[:,col],'+k')
        plt.ylim((0.295,0.405)) 
        plt.ylabel("seawater d88/86Sr [permil]")
        plt.xlim((-2,37))
        plt.grid()
    
    
        plt.subplot(322)
        for run in range(1,self.Nruns+1):
            if (ddiff[run]>0):
                plt.plot(ModelDT[:,0,run],ModelDT[:,4,run],'-r')
            else:
                plt.plot(ModelDT[:,0,run],ModelDT[:,4,run],'-b')
        plt.plot(ModelDT[:,0,0],ModelDT[:,4,0],'-k',linewidth=3.0)
        plt.plot(Data[:,0],Data[:,col],'+k')
        plt.ylim((0.295,0.405))
        plt.xlim((-2,37))
        plt.grid()
    
    
        plt.subplot(323)
        for run in range(1,self.Nruns+1):
            if (ddiff[run]>0):
                plt.plot(Model[:,0,run],-1e6*1e-15*(Model[:,2,run]-Model[:,3,run]),'-r')
            else:
                plt.plot(Model[:,0,run],-1e6*1e-15*(Model[:,2,run]-Model[:,3,run]),'-b')
        p = np.polyfit(Model[:,0,0],Model[:,1,0], 1)
        plt.plot(Model[:,0,0],-1e6*1e-15*(Model[:,2,0]-Model[:,3,0]) ,'-k',linewidth=3.0)
        plt.ylim((-1.2e1,1.2e1))
        plt.xlim((-2,37))
        plt.ylabel("Sr in-out imbalance [Pmol/Myr]")
        plt.grid()
    
    
        plt.subplot(324)
        for run in range(1,self.Nruns+1):
            if (ddiff[run]>0):
                plt.plot(ModelDT[:,0,run],-1000*1e-12*(ModelDT[:,2,run]-ModelDT[:,3,run]),'-r')
            else:
                plt.plot(ModelDT[:,0,run],-1000*1e-12*(ModelDT[:,2,run]-ModelDT[:,3,run]),'-b')
        p = np.polyfit(ModelDT[:,0,0],ModelDT[:,1,0], 1)
        plt.plot(ModelDT[:,0,0],-1000*1e-12*(ModelDT[:,2,0]-ModelDT[:,3,0]) ,'-k',linewidth=3.0)
        plt.ylim((-1.2e1,1.2e1))
        plt.xlim((-2,37))
        plt.grid()
    
        L03 = self.LoadLear03()
    
        plt.subplot(325)
        for run in range(1,self.Nruns+1):
            if (ddiff[run]>0):
                plt.plot(Model[:,0,run],Model[:,1,run]/1.4e21*1e6,'-r')
            else:
                plt.plot(Model[:,0,run],Model[:,1,run]/1.4e21*1e6,'-b')
        plt.plot(Model[:,0,0],Model[:,1,0]/1.4e21*1e6 ,'-k',linewidth=3.0)
        plt.ylim((-0.1e2,2e2))
        plt.xlim((-2,37))
        plt.xlabel("Age [MyrBP]")
        plt.ylabel("seawater Sr concentration [µmol/kg]")
        plt.grid()
    
        plt.twinx()
        plt.plot(L03[:,1],L03[:,2:9],'+k')
        plt.ylim((0,3))
        plt.xlim((-2,37))
    
        plt.subplot(326)
        for run in range(1,self.Nruns+1):
            if (ddiff[run]>0):
                plt.plot(ModelDT[:,0,run],ModelDT[:,1,run]/1.4e21*1e6,'-r')
            else:
                plt.plot(ModelDT[:,0,run],ModelDT[:,1,run]/1.4e21*1e6,'-b')
        plt.plot(ModelDT[:,0,0],ModelDT[:,1,0]/1.4e21*1e6 ,'-k',linewidth=3.0)
        plt.ylim((-0.1e2,2e2))
        plt.xlim((-2,37))
        plt.xlabel("Age [MyrBP]")
        plt.grid()
    
        plt.twinx()
        plt.plot(L03[:,1],L03[:,2:9],'+k')
        plt.ylim((0,3))
        plt.xlim((-2,37))
        
        plt.savefig('Published_Paytan_etal_2021.pdf')
    
        plt.show()
        

    def SingleFigure(self):
        Model=self.Model
        Data=self.Data
        
        import matplotlib.pyplot as plt
        plt.figure(2)
        plt.subplot(311)
        plt.plot(Data[:,0],Data[:,3],xp, p(xp), '--',Model[:,0],Model[:,4],'-')
        plt.grid()
        plt.ylabel("d88/86Sr")
        plt.subplot(312)
        plt.plot(Model[:,0],Model[:,1],'-')
        plt.ylabel("[Sr]")
        plt.grid()
        plt.subplot(313)
        plt.plot(Model[:,0],Fin-Model[:,2],'-')
        plt.ylabel("In-OUT [mol/yr]")
        plt.xlabel("MyrBP")
        plt.grid()
        plt.savefig("Figure_Sc5_eps018_N30.pdf")

    def LoadLear03(self):
        Data = np.genfromtxt('Lear2003_SrOverCa_RAW.txt', dtype=float, delimiter='\t')
        Data = Data[1:,:]
    
        DataD = np.copy(Data)
        DataD[:,2:9] *= 0
        idx = np.where(Data == 522)
        DataD[idx,2:9] = Data[idx,2:9] + (3.000 - 1.000) * 0.101
        idx = np.where(Data == 523)
        DataD[idx,2:9] = Data[idx,2:9] + (3.100 - 1.000) * 0.101
        idx = np.where(Data == 525)
        DataD[idx,2:9] = Data[idx,2:9] + (1.500 - 1.000) * 0.101
        idx = np.where(Data == 573)
        DataD[idx,2:9] = Data[idx,2:9] + (3.650 - 1.000) * 0.101
        idx = np.where(Data == 608)
        DataD[idx,2:9] = Data[idx,2:9] + (3.470 - 1.000) * 0.101
        idx = np.where(Data == 689)
        DataD[idx,2:9] = Data[idx,2:9] + (1.365 - 1.000) * 0.101
        idx = np.where(Data == 690)
        DataD[idx,2:9] = Data[idx,2:9] + (2.055 - 1.000) * 0.101
        idx = np.where(Data == 926)
        DataD[idx,2:9] = Data[idx,2:9] + (3.525 - 1.000) * 0.101
        idx = np.where(Data == 1052)
        DataD[idx,2:9] = Data[idx,2:9] + (0.750 - 1.000) * 0.101

        DataDS = np.copy(DataD)
        DataDS[:,2] += 0.29 
        DataDS[:,3] += 0.11
        DataDS[:,4] += 0.42 
        DataDS[:,5] += 0.28 
        DataDS[:,6] += 0.15 
        DataDS[:,8] += 0.08
        DataDS[:,9] += 0.16
    
        DataDSCa = np.copy(DataDS)
        DataDSCa[:,2] *=  (1.0 + 0.7* DataDSCa[:,1]/50.0)
        DataDSCa[:,3] *=  (1.0 + 0.7* DataDSCa[:,1]/50.0)
        DataDSCa[:,4] *=  (1.0 + 0.7* DataDSCa[:,1]/50.0)
        DataDSCa[:,5] *=  (1.0 + 0.7* DataDSCa[:,1]/50.0)
        DataDSCa[:,6] *=  (1.0 + 0.7* DataDSCa[:,1]/50.0)
        DataDSCa[:,7] *=  (1.0 + 0.7* DataDSCa[:,1]/50.0)
        DataDSCa[:,8] *=  (1.0 + 0.7* DataDSCa[:,1]/50.0)
        DataDSCa[:,9] *=  (1.0 + 0.7* DataDSCa[:,1]/50.0)
        if (False):
            plt.figure(99)
            plt.subplot(311)
            plt.plot(Data[:,1],Data[:,2:9],'.')
            plt.xlim((0,40))
            plt.grid()
            plt.subplot(312)
            plt.plot(DataD[:,1],DataD[:,2:9],'.')
            plt.xlim((0,40))
            plt.grid()
            plt.subplot(313)
            plt.plot(DataDS[:,1],DataDS[:,2:9],'.')
            plt.xlim((0,40))
            plt.grid()
            plt.show()
    
        return DataDS

    def LoadData(self,skiptoprows):
        Data = np.genfromtxt('RawData.txt', dtype=float, delimiter='\t')
        Data = Data[skiptoprows:,:]
        Data[:,3] = Data[:,1]+0.538
        return Data

    def LoadDataSc(self,skiptoprows):
        #Data = np.genfromtxt('RawData.txt', dtype=float, delimiter='\t')
        #Data = np.genfromtxt('Sc1.txt', dtype=float, delimiter='\t')
        Data = np.genfromtxt('Sc5.txt', dtype=float, delimiter='\t') # This is the scenario we picked for the 07/2020 submission
        Data = Data[skiptoprows:,:]
        return Data

    def DataDerivative(self,Data):
        dX = np.zeros((Data.shape(0),2))
        for row in range(len(X)-1):
            dX[row,0] = (X[row+1,0]+X[row,0])/2
            dX[row,1] = X[row+1,3]-X[row,3]
        return dX

    def FitPoly(self,Data, col,n):
        p = np.poly1d(np.polyfit(Data[:,0],Data[:,col], n))
        pd = np.polyder(p) 
        xp = np.linspace(0, 35, 100)
        error = p(Data[:,0])-Data[:,col]
        return (p,pd,xp,error)

    def PlotPolyFit(self,Data,p,pd,xp):
        plt.subplot(311)
        plt.plot(Data[:,0],Data[:,3],xp, p(xp), '--')

        plt.subplot(312)
        plt.plot(dX[:,0],dX[:,1],xp, pd(xp), '--')

        plt.subplot(313)
        plt.hist(error)
        plt.show()

    def IntegrateConstFin(self,pd,Fin,Sr0,dsw0,din,eps):

        dt = 1000   # timestep
    
        dout0 = dsw0+eps
        Fout0 = (pd(0)/1e6*Sr0 + (din-dsw0)*Fin)/eps
        state0 = np.array([0,Sr0,Fout0,Fin,dsw0,dout0,din])
        Model = np.zeros((351,7))
        Model = state0
        Fout = np.zeros(35001)
        Fout[0] = Fout0
    

        for kyr in range(1,35001):
            state1 = state0*0
            Myr = kyr/1000.0
            state1[0] = Myr
            state1[1] = state0[1] + state0[2]*dt - state0[3]*dt
            state1[4] = (state0[1]*state0[4] + state0[2]*dt*state0[5] - state0[3]*dt*state0[6])/state1[1] 
            state1[2] = (pd(Myr)/1e6*state1[1] + (din-state1[4])*Fin)/eps
            state1[3] = Fin
            state1[5] = state1[4]+eps
            state1[6] = din
        
            if (kyr%100 == 0):
                Model = np.vstack([Model,state1])
            
            Fout[kyr] = state1[2]
            state0 = state1
    
        del state0,state1
    
        
        
        FoutDT = Fout - np.mean(Fout) + Fin
        state0 = np.array([35,Sr0,FoutDT[-1],Fin,0.3924,0.3924+eps,din])
        ModelDT = state0

        for kyr in range(1,35001):
            state1 = state0*0
            Myr = 35 - kyr/1000.0

            state1[0] = Myr
            state1[1] = state0[1] - state0[2]*dt + state0[3]*dt
            state1[4] = (state0[1]*state0[4] - state0[2]*dt*state0[5] + state0[3]*dt*state0[6])/state1[1] 
            state1[2] = FoutDT[-kyr]
    
            state1[3] = Fin
            state1[5] = state1[4]+eps
            state1[6] = din
    
            if (kyr%100 == 0):
                #print(kyr, Myr)
                ModelDT = np.vstack([ModelDT,state1])
        
            state0 = state1
        
        ModelDT
    
        del state0,state1

        #print(ModelDT[0:3,:])
        #print(ModelDT[-3:-0,:])
    
        return Model, ModelDT
        
        
    

########### MODEL #####################
if __name__ == "__main__":
    
    
    import InverseModel_Vfinal as ModelLib
    
    ModelObject = ModelLib.StableStrontium()
    ModelObject.PublishedFigure()
    
    if (False):
        # this is a legacy code fragment that serves as an illustration of how the old V0/V1 interface was used
        # leaving this in here for additional tracability and debuging 
        
        #Data = LoadData(0)
        col = 3
        Data = LoadDataSc(0)
        col = 1

        (p,pd,xp,error) = FitPoly(Data, col, 7)

        dmean = np.mean(Data[:,col])
        dmean = np.mean(p(xp))
        print(dmean)

        Fin = 38e9
        Sr0 = 88e-6 * 1.4e21 #Fin*2.5e6
        dsw0 = p(0)
        print("d initial",dsw0)
        print(p(35))
        eps = -0.18 #-0.22 #-0.167
        din = dmean + eps #0.183
        print("dmean: ",dmean)
        print("din: ",din)

        if (True):
    
            Nruns = 3
            Model = np.zeros((351,7,Nruns+1))
            ModelDT = np.zeros((351,7,Nruns+1))
            ddiff = np.zeros(Nruns+1)
            trend = np.zeros(Nruns+1)
            (Model[:,:,0],ModelDT[:,:,0]) = IntegrateConstFin(pd,Fin,Sr0,dsw0,din,eps)
            ddiff[0] = din - (dmean + eps)

            for run in range(1,Nruns+1):
                while True:
                    dinE = din + 0.01*np.random.randn()
                    epsE = eps + 0.005*np.random.randn()
                    ddiff[run] = dinE - (dmean + epsE)
                    (Model[:,:,run],ModelDT[:,:,run]) = IntegrateConstFin(pd,Fin,Sr0,dsw0,dinE,epsE)
                    #Model = np.stack((Model,NewModel),axis=2)
                    if (np.min(Model[:,1,run])>0):
                        break  
                print("Done with run {0}".format(run))

    
            plt.figure(1,figsize=(10,8))
    
            plt.subplot(321)
            #for run in range(1,Nruns+1):
            #    if (ddiff[run]>0):
            #        plt.plot(Model[:,0,run],Model[:,4,run],'-r')
            #    else:
            #        plt.plot(Model[:,0,run],Model[:,4,run],'-b')
            plt.plot(Model[:,0,0],Model[:,4,0],'-k',linewidth=3.0)
            plt.plot(Data[:,0],Data[:,col],'+k')
            plt.ylim((0.295,0.405)) 
            plt.ylabel("seawater d88/86Sr [permil]")
            plt.xlim((-2,37))
            plt.grid()
    
    
            plt.subplot(322)
            for run in range(1,Nruns+1):
                if (ddiff[run]>0):
                    plt.plot(ModelDT[:,0,run],ModelDT[:,4,run],'-r')
                else:
                    plt.plot(ModelDT[:,0,run],ModelDT[:,4,run],'-b')
            plt.plot(ModelDT[:,0,0],ModelDT[:,4,0],'-k',linewidth=3.0)
            plt.plot(Data[:,0],Data[:,col],'+k')
            plt.ylim((0.295,0.405))
            plt.xlim((-2,37))
            plt.grid()
    
    
            plt.subplot(323)
            for run in range(1,Nruns+1):
                if (ddiff[run]>0):
                    plt.plot(Model[:,0,run],-1e6*1e-15*(Model[:,2,run]-Model[:,3,run]),'-r')
                else:
                    plt.plot(Model[:,0,run],-1e6*1e-15*(Model[:,2,run]-Model[:,3,run]),'-b')
            p = np.polyfit(Model[:,0,0],Model[:,1,0], 1)
            plt.plot(Model[:,0,0],-1e6*1e-15*(Model[:,2,0]-Model[:,3,0]) ,'-k',linewidth=3.0)
            plt.ylim((-1.2e1,1.2e1))
            plt.xlim((-2,37))
            plt.ylabel("Sr in-out imbalance [Pmol/Myr]")
            plt.grid()
    
    
            plt.subplot(324)
            for run in range(1,Nruns+1):
                if (ddiff[run]>0):
                    plt.plot(ModelDT[:,0,run],-1000*1e-12*(ModelDT[:,2,run]-ModelDT[:,3,run]),'-r')
                else:
                    plt.plot(ModelDT[:,0,run],-1000*1e-12*(ModelDT[:,2,run]-ModelDT[:,3,run]),'-b')
            p = np.polyfit(ModelDT[:,0,0],ModelDT[:,1,0], 1)
            plt.plot(ModelDT[:,0,0],-1000*1e-12*(ModelDT[:,2,0]-ModelDT[:,3,0]) ,'-k',linewidth=3.0)
            plt.ylim((-1.2e1,1.2e1))
            plt.xlim((-2,37))
            plt.grid()
    
            L03 = LoadLear03()
    
            plt.subplot(325)
            for run in range(1,Nruns+1):
                if (ddiff[run]>0):
                    plt.plot(Model[:,0,run],Model[:,1,run]/1.4e21*1e6,'-r')
                else:
                    plt.plot(Model[:,0,run],Model[:,1,run]/1.4e21*1e6,'-b')
            plt.plot(Model[:,0,0],Model[:,1,0]/1.4e21*1e6 ,'-k',linewidth=3.0)
            plt.ylim((-0.1e2,2e2))
            plt.xlim((-2,37))
            plt.xlabel("Age [MyrBP]")
            plt.ylabel("seawater Sr concentration [µmol/kg]")
            plt.grid()
    
            plt.twinx()
            plt.plot(L03[:,1],L03[:,2:9],'+k')
            plt.ylim((0,3))
            plt.xlim((-2,37))
    
            plt.subplot(326)
            for run in range(1,Nruns+1):
                if (ddiff[run]>0):
                    plt.plot(ModelDT[:,0,run],ModelDT[:,1,run]/1.4e21*1e6,'-r')
                else:
                    plt.plot(ModelDT[:,0,run],ModelDT[:,1,run]/1.4e21*1e6,'-b')
            plt.plot(ModelDT[:,0,0],ModelDT[:,1,0]/1.4e21*1e6 ,'-k',linewidth=3.0)
            plt.ylim((-0.1e2,2e2))
            plt.xlim((-2,37))
            plt.xlabel("Age [MyrBP]")
            plt.grid()
    
            plt.twinx()
            plt.plot(L03[:,1],L03[:,2:9],'+k')
            plt.ylim((0,3))
            plt.xlim((-2,37))
    
            #plt.show()
    

            plt.savefig('VectorOUT.pdf')
    
            plt.show()
    
            if (False):
                plt.figure(3)
                plt.subplot(411)
                plt.plot(ddiff,trend ,'+')
                plt.grid()
    
                plt.subplot(412)
                plt.plot(Data[:,0],Data[:,3],'+k')
                for run in range(1,Nruns+1):
                    if (ddiff[run]>0):
                        plt.plot(Model[:,0,run],Model[:,4,run],'-r')
                    else:
                        plt.plot(Model[:,0,run],Model[:,4,run],'-b')
                plt.grid()
    
                plt.subplot(413)
                plt.plot(Data[:,0],Data[:,3],'+k')
                for run in range(1,Nruns+1):
                    if (ddiff[run]>0):
                        plt.plot(ModelDT[:,0,run],ModelDT[:,4,run],'-r')
                    else:
                        plt.plot(ModelDT[:,0,run],ModelDT[:,4,run],'-b')
                plt.grid()
    
                plt.subplot(414)
                for run in range(1,Nruns+1):
                    if (ddiff[run]>0):
                        plt.plot(Model[:,0,run],Model[:,2,run],'-r')
                    else:
                        plt.plot(Model[:,0,run],Model[:,2,run],'-b')
                plt.grid()
                plt.show()

        if (False):
            Model = IntegrateConstFin(pd,Fin,Sr0,dsw0,din,eps)
    
            plt.figure(2)
            plt.subplot(311)
            plt.plot(Data[:,0],Data[:,3],xp, p(xp), '--',Model[:,0],Model[:,4],'-')
            plt.grid()
            plt.ylabel("d88/86Sr")
            plt.subplot(312)
            plt.plot(Model[:,0],Model[:,1],'-')
            plt.ylabel("[Sr]")
            plt.grid()
            plt.subplot(313)
            plt.plot(Model[:,0],Fin-Model[:,2],'-')
            plt.ylabel("In-OUT [mol/yr]")
            plt.xlabel("MyrBP")
            plt.grid()

            plt.savefig("SrInverse_V0.png")
            plt.savefig("Figure_Sc5_eps018_N30.png")
            plt.savefig("Figure_Sc5_eps018_N30.ps")
            plt.show()
            print(Model)    



