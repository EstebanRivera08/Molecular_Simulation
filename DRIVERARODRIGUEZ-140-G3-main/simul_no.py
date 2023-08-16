import numpy as np
import math

class Simul:
    """ 
    This is the prototype of the simulation code
    It moves the particles with at _velocity, using a vector notation: numpy should be used.
    """

    
    def __init__(self, sample_time, sigma, Lx, Ly , N, m):
        print("Simul init")
        
        self.sigma = sigma  # particle radius
        self._sample_time = sample_time
        self.mass = m
        self.boite = np.array([[Lx,Ly]])
        self.N = N 
        
        Nx = round(np.sqrt(N*Lx/Ly))
        Ny = round(Nx*Ly/Lx)
       
        self.position = np.zeros((Nx*Ny,2))
        
        posx = np.linspace(2*self.sigma, Lx - 2*self.sigma,Nx)
        posy = np.transpose(np.linspace(2*self.sigma, Ly - 2*self.sigma, Ny))
        
        #[np.newaxis]
        n = 0
        for i in range(Nx) :
            for j in range(Ny) :
                self.position[n,:] = np.array([[posx[i],posy[j]]])
                n += 1
 
        
       #self.position = np.zeros([20,2]) + np.array([[0.5,0.5]])2*
        #self.position = np.array([[0.25, 0.25], [0.75, 0.75]])
        #self._velocity = np.array([[1 ,1], [-1.5,-1.5]])
        
        #self.position = np.array([[0.25, 0.25], [0.75, 0.25], [0.25, 0.75], [0.75, 0.75]])  # starting positions
        
        self._i, self._j = np.triu_indices(self.position.shape[0], k=1)  # all pairs of indices between particles
        
        self._velocity = np.random.normal(size = self.position.shape)*2 #random velocities

        

    def _wall_time(self):
    # Looking for vertical and horizontal intersceptions with axes? (x=0;x=Lx;y=0;y=Ly); Lx=Ly=1
    # r = r0 + v*t = 
        tc = np.where(self._velocity > 0, (self.boite - self.sigma-self.position) / self._velocity, (0 + self.sigma - self.position) / self._velocity)
        
        i,j = np.where(tc < self._sample_time)
        print(i,j)
        
        imin, jmin = np.unravel_index(tc.argmin(), tc.shape)
        
        return tc[i,j], i, j
        
    def _pair_time(self):
        
        rij = self.position[self._i]-self.position[self._j]  # set of all 6 separation vectors        
        vij = self._velocity[self._i]-self._velocity[self._j]
        
        a = np.sum(vij*vij,1)
        b = 2*np.sum(rij*vij,1)
        c = np.sum(rij*rij,1) - (2*self.sigma)**2
        
        tcol = (-b - np.sqrt(b**2-4*a*c))/(2*a)
        #print("t,b",tcol,b)
        
       #tcol[tcol < 0] = self._sample_time + 1
        #print("2da Mod",tcol)
        tcol[c < 0] = self._sample_time*0.5  #???
        tcol[b > 0] = self._sample_time + 1

        index = np.where(tcol < self._sample_time)[0]

        return self._i[index], self._j[index], b[index], vij[index], index
        
    def md_step(self):
        # Dans cette méthode, il n'y pas la gestion temporelle attendue, à savoir que l'on puisse observer plusieurs évènements de type choc contre les parois et/ou chocs entre les particules dans un intervalle de temps sample_time, de sorte que l'on peut avoir l'illusion d'un fonctionnement correct lorsque sample_time est petit devant les temps des évènements physiques caractéristiques de ce système. 
        Enertot=np.sum(np.sum(self._velocity**2))
        print('Energie totale',  Enertot)
        pressure = 0
        ic, jc, b, vij, index = self._pair_time()
        tc, imin, jmin = self._wall_time()
        print('Simul::md_step')
        
        if imin.size != 0 :
            self._velocity[imin,jmin] *= -1
            sumV = sum(abs(self._velocity[imin,jmin]))
            pressure += 2*self.mass*sumV/(2*self.boite[0])/self._sample_time
            print("P",pressure)
            
        elif index.size != 0 :
            print('ic=',ic)
            print('jc=',jc)
            # on s'aperçoit que la gestion des vitesses peut concerner un nombre de particules >2, ce qui n'est pas gérable avec lapproche considérée, c'est à ce moment là que la conservation d'énergie n'est plus respectée (d'où le ralentissement des particules observées)...
            K11 = np.sum(self._velocity[ic]**2+self._velocity[jc]**2)
            K12 = np.sum(self._velocity[ic]+self._velocity[jc])
            print("Conservation 1: ",K11,K12)
            
            ni = (self.position[jc]-self.position[ic])/np.linalg.norm(self.position[jc]-self.position[ic]) 
            holi = np.transpose(np.sum(ni*vij,1)[np.newaxis])

            self._velocity[ic] -= ni*holi
            self._velocity[jc] += ni*holi
            
            K21 = np.sum(self._velocity[ic]**2+self._velocity[jc]**2)
            K22 = np.sum(self._velocity[ic]+self._velocity[jc])
            #print("Conservation 2: ",K21,K22)
            
            print("Comparation E: ",K11/K21,K11-K21,(K11-K21)/K11)
            
            

        self.position = self.position + self._sample_time * self._velocity
            
        return pressure
        

    def __str__(self):   # this is used to print the position and velocity of the particles
        p = np.array2string(self.position)
        v = np.array2string(self._velocity)
        return 'pos= '+p+'\n'+'vel= '+v+'\n'

    # Pour remettre le progamme dans une forme qui permette l'exécution fiable de cette simulation (quel que soit le pas de temps choisi et avec la certitude de la conservation d'énergie) il faudrait le réformer en profondeur, nous choisissons de corriger celui de Suzana Giraldo Betancur. 