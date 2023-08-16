import numpy as np
import math


class Simul:
    """ 
    This is the prototype of the simulation code
    It moves the particles with at _velocity, using a vector notation: numpy should be used.
    """
    def __init__(self, sample_time, sigma, L, N):
        
        self.sigma = sigma  # particle radius
        
        Nx = round(np.sqrt(N*L/L))
        Ny = round(Nx*L/L)
        self.position = np.zeros((Nx*Ny,2))
        
        posx = np.linspace(2*self.sigma, L - 2*self.sigma,Nx)
        posy = np.transpose(np.linspace(2*self.sigma, L - 2*self.sigma, Ny))
        
        #[np.newaxis]
        n = 0
        for i in range(Nx) :
            for j in range(Ny) :
                self.position[n,:] = np.array([[posx[i],posy[j]]])
                n += 1
                
        self._velocity = np.random.normal(size=self.position.shape)  # random velocities
        self._i, self._j = np.triu_indices(self.position.shape[0], k=1)  # all pairs of indices between particles
        self._sample_time=sample_time
        self.L=L
        self.boite = np.array([[L,L]])
        self.p_ =  np.array([])
        self.N = Nx*Ny

    def _wall_time(self):
        
        tc = np.where(self._velocity >0, ( self.L- self.sigma - self.position)/self._velocity , (0 + self.sigma - self.position)/self._velocity )
        disk, direction = np.unravel_index(tc.argmin(), tc.shape)
        
        return tc, disk, direction

    def _pair_time(self):
        
        i, j = np.triu_indices(self.position.shape[0], k=1)
        rij = self.position[i]-self.position[j]
        vij = self._velocity[i] - self._velocity[j]
        
        a = np.sum(vij*vij,1)
        b = 2*np.sum(rij*vij,1)
        c = np.sum(rij*rij,1) - 4*self.sigma*self.sigma
        d = b*b - 4*a*c
        
        tcp = np.where((d > 0) & (b<0) & (c>0), (-b - np.sqrt(d))/(2*a) , 1000 )
        

        tcp_p = tcp[tcp.argmin()]
        ic = i[tcp.argmin()]
        jc = j[tcp.argmin()]
        
        vi = self._velocity[ic]
        vj = self._velocity[jc]
        dif_v = vi - vj
        
        return tcp_p, ic, jc, dif_v
        
    def md_step(self):
        d
        tc,disk,direction = self._wall_time()
        tcp_p, ic, jc, dif_v = self._pair_time()        
        tmin = min(tc[disk, direction], tcp_p)
        
        pressure = 0
        
        now = 0
        while now + tmin < self._sample_time:
            now+=tmin
            if tc[disk, direction] < tcp_p:
                
                self.position += tc[disk, direction]*self._velocity
                self._velocity[disk, direction] = -self._velocity[disk, direction]             
                pressure += 2*1*abs(self._velocity[disk, direction])*(1/(4*self.L*self._sample_time))
                          
            else:
                self.position += tcp_p* self._velocity
                u = (self.position[jc] - self.position[ic])/np.linalg.norm(self.position[jc] - self.position[ic])
                
                self._velocity[ic] -= u*np.sum(u*dif_v)
                self._velocity[jc] += u*np.sum(u*dif_v)   
                
            tc,disk,direction= self._wall_time()
            tcp_p, ic, jc, dif_v = self._pair_time()
            tmin = min(tc[disk, direction], tcp_p)   
               
        self.position = self.position + (self._sample_time-now) * self._velocity
        
        Energy = 1/2*1*np.sum(self._velocity**2)
        
        return pressure, Energy

    def __str__(self):   # this is used to print the position and velocity of the particles
        print("Step")
        p = np.array2string(self.position)
        v = np.array2string(self._velocity)
 
        return 