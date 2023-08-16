from simul import Simul
import matplotlib.pyplot as plt
from animatesimul import AnimateSimul
import numpy as np
import pickle


def b_VdW():
    
    np.random.seed(10)
    
    sample_time = 0.01
    Ni = np.linspace(10,190,10)
    Li = np.logspace(0,1,11)
    
    frames = 500
    
    t = np.zeros((frames,len(Ni),len(Li)))
    p = np.zeros((frames,len(Ni),len(Li)))
    e = np.zeros((frames,len(Ni),len(Li)))
    pi = np.zeros((len(Ni),len(Li)))
    ei = np.zeros((len(Ni),len(Li)))
    p_prom = np.array([[]])
    V_N = np.array([[]])
    e_prom = np.array([[]])
    
    for n in range(len(Ni)):
        for l in range(len(Li)):

            simulation = Simul(sample_time, sigma=0.01, L = Li[l] , N = Ni[n])  #  sigma particle radius  
        
            for i in range(frames): 
                        
                t[i,n,l] = i*sample_time
                Pressure, Energy = simulation.md_step()
                p[i,n,l] = Pressure
            
            V_N = np.append(V_N,(Li[l]**2)/Ni[n])
            pi[n,l] = np.mean(p[:,n,l])
            ei[n,l] = Energy/Ni[n]
            p_prom = np.append(p_prom,pi[n,l])
            e_prom = np.append(e_prom,ei[n,l])
            
        print(";)")
    
    sort_index = np.argsort(V_N)
    V_N = np.sort(V_N)
    e_prom = e_prom[sort_index]
    p_prom = p_prom[sort_index]
    b = V_N - e_prom/p_prom
    
    Va = Li**2
    Na = np.transpose(np.array([Ni]))
    
    Data = {'V': Va, 'N': Ni, 'V/N': V_N, 'E': e_prom, 'P': p_prom, 'bi': b, 'ei': ei, 'pi': pi, 't': t, 'p': p}
    
    with open('Data.pickle', 'wb') as file:  # this sends all the state to a file
        pickle.dump(Data,file)

#b_VdW()

with open('Data.pickle','rb') as file:
    Data = pickle.load(file)
    
V_N = Data['V/N']
P = Data['P']
pi = Data['pi']
V = Data['V']
N = Data['N']
b = Data['bi']
ei = Data['ei']
E = Data['E']
print("The average b-value for VdW equation is:",np.mean(b))

grid = plt.GridSpec(2, 2, wspace=0.3, hspace=0.5)

plt.subplot(grid[0, :])
plt.plot(V_N,P,color='k',linewidth=0.5,label="P") 
plt.xlabel('V/N') 
plt.ylabel('Pressure')
plt.ylim(0,20)
plt.grid()

plt.subplot(grid[1, 0])
plt.plot(V,pi[0,:],linewidth=0.5,label = "N = " + str(N[0])) 
plt.plot(V,pi[1,:],linewidth=0.5,label = "N = " + str(N[1])) 
plt.plot(V,pi[3,:],linewidth=0.5,label = "N = " + str(N[3])) 
plt.plot(V,pi[5,:],linewidth=0.5,label = "N = " + str(N[5])) 
plt.plot(V,pi[8,:],linewidth=0.5,label = "N = " + str(N[8])) 
plt.ylim(0,100)
plt.xlabel('Volumen') 
plt.ylabel('Pressure') 
plt.title('P vs. V')
plt.legend()
plt.grid()

plt.subplot(grid[1,1])
plt.plot(N,pi[:,0],linewidth=0.5,label = "V = " + str(round(V[0]))) 
plt.plot(N,pi[:,2],linewidth=0.5,label = "V = " + str(round(V[2]))) 
plt.plot(N,pi[:,4],linewidth=0.5,label = "V = " + str(round(V[4]))) 
plt.plot(N,pi[:,7],linewidth=0.5,label = "V = " + str(round(V[7]))) 
plt.plot(N,pi[:,10],linewidth=0.5,label = "V = " + str(round(V[10]))) 
plt.ylim(0,100)
plt.xlabel('No. Particles') 
plt.ylabel('Pressure') 
plt.title('P vs. N') 
plt.legend()
plt.grid()

plt.show()

#print(Data)