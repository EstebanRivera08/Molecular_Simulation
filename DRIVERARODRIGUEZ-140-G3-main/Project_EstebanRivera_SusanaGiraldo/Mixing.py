import numpy as np
import math


class Mix:
    """ 
    This is the prototype of the simulation code
    It moves the particles with at _velocity, using a vector notation: numpy should be used.
    """
    def __init__(self, sample_time, sigma1, sigma2, Lx, Ly, N):
        
        self.sigma1 = sigma1  # particle radius
        self.sigma2 = sigma2
        
        Nx= round(np.sqrt(N*Ly/(Lx/2)))
        Ny = round(Nx*Ly/(Lx/2))
        
        self.position1 = np.zeros((Nx*Ny,2))
        self.position2 = np.zeros((Nx*Ny,2))
        
        posx1 = np.linspace(0, Lx/2, Nx+2)
        posx2 = np.linspace(Lx/2, Lx, Nx+2)
        posy = np.linspace(0, Ly, Ny+2)
        
        #[np.newaxis]
        n = 0
        for i in range(Nx) :
            for j in range(Ny) :
                self.position1[n,:] = np.array([[posx1[i+1],posy[j+1]]])
                self.position2[n,:] = np.array([[posx2[i+1],posy[j+1]]])
                n += 1

        self._velocity1 = np.random.normal(size=self.position1.shape)  # random velocities
        self._velocity2 = np.random.normal(size=self.position2.shape)  # random velocities
        
        self._sample_time = sample_time
        self.boite = np.array([[Lx,Ly]])

    def _wall_time(self):
        
        tc1 = np.where(self._velocity1 >0, ( self.boite - self.sigma1 - self.position1)/self._velocity1 , (0 + self.sigma1 - self.position1)/ self._velocity1 )
        
        tc2 = np.where(self._velocity2 >0, ( self.boite - self.sigma2 - self.position2)/self._velocity2 , (0 + self.sigma2 - self.position2)/ self._velocity2 )
        

        disk1, direction1 = np.unravel_index(tc1.argmin(), tc1.shape)
        disk2, direction2 = np.unravel_index(tc2.argmin(), tc2.shape)
        
        if tc1[disk1,direction1] < tc2[disk2,direction2] :
            tc = tc1
            Int = 1
        
        else :
            tc = tc2
            Int = 2
            
        disk, direction = np.unravel_index(tc.argmin(), tc.shape)
        
        return tc, disk, direction, Int

    def _pair_time(self):
        
        i, j = np.triu_indices(self.position2.shape[0], k=1)
        
        rij11 = self.position1[i] - self.position1[j]
        vij11 = self._velocity1[i] - self._velocity1[j]
        
        rij21 = self.position2[i] - self.position1[j]
        vij21 = self._velocity2[i] - self._velocity1[j]
        
        rij12 = self.position1[i] - self.position2[j]
        vij12 = self._velocity1[i] - self._velocity2[j]
        
        rij22 = self.position2[i] - self.position2[j]
        vij22 = self._velocity2[i] - self._velocity2[j]
        
        a11 = np.sum(vij11 * vij11,1)
        b11 = 2*np.sum(rij11 * vij11,1)
        c11 = np.sum(rij11 * rij11,1) - 4*self.sigma1*self.sigma1
        d11 = b11*b11 - 4*a11*c11
        
        a12 = np.sum(vij12 * vij12,1)
        b12 = 2*np.sum(rij12 * vij12,1)
        c12 = np.sum(rij12 * rij12,1) - (self.sigma2+self.sigma1)**2
        d12 = b12*b12 - 4*a12*c12
        
        a21 = np.sum(vij21 * vij21,1)
        b21 = 2*np.sum(rij21 * vij21,1)
        c21 = np.sum(rij21 * rij21,1) - (self.sigma2+self.sigma1)**2
        d21 = b21*b21 - 4*a21*c21
        
        a22 = np.sum(vij22 * vij22,1)
        b22 = 2*np.sum(rij22 * vij22,1)
        c22 = np.sum(rij22 * rij22,1) - 4*self.sigma2*self.sigma2
        d22 = b22*b22 - 4*a22*c22
        
        #print(d12,b12,c12)
        
        tcp11 = np.where((d11 > 0) & (b11 < 0) & (c11 > 0), (-b11 - np.sqrt(d11))/(2*a11) , 1000)
        tcp12 = np.where((d12 > 0) & (b12 < 0) & (c12 > 0), (-b12 - np.sqrt(d12))/(2*a12) , 1000 )
        tcp21 = np.where((d21 > 0) & (b21 < 0) & (c21 > 0), (-b21 - np.sqrt(d21))/(2*a21) , 1000 )
        tcp22 = np.where((d22 > 0) & (b22 < 0) & (c22 > 0), (-b22 - np.sqrt(d22))/(2*a22) , 1000 )
        
        
        tcp_p11 = tcp11[tcp11.argmin()]
        
        tcp_p12 = tcp12[tcp12.argmin()]
        
        tcp_p21 = tcp21[tcp21.argmin()]
        
        #tcp22[np.isnan(tcp22)] = 1 + tcp22[tcp22.argmax()]
        #tcp22[tcp22 < 0] = 1 + tcp22[tcp22.argmax()]
        tcp_p22 = tcp22[tcp22.argmin()]
        
        #print(tcp_p11, tcp_p12, tcp_p21, tcp_p22)
        tcp_p = min(tcp_p11, tcp_p12,tcp_p21, tcp_p22)
        
        if tcp_p == tcp_p11 :  
            Interaction = 11
            tcp = tcp11
            ic = i[tcp.argmin()]
            jc = j[tcp.argmin()]

            vi = self._velocity1[ic]
            vj = self._velocity1[jc]
            dif_v = vi - vj
            
        elif tcp_p == tcp_p12 :
            Interaction = 12
            tcp = tcp12
            ic = i[tcp.argmin()]
            jc = j[tcp.argmin()]

            vi = self._velocity1[ic]
            vj = self._velocity2[jc]
            dif_v = vi - vj
            
        elif tcp_p == tcp_p21 :        
            Interaction = 21
            tcp = tcp21
            ic = i[tcp.argmin()]
            jc = j[tcp.argmin()]

            vi = self._velocity2[ic]
            vj = self._velocity1[jc]
            dif_v = vi - vj
        
        else :           
            Interaction = 22
            tcp = tcp22
            ic = i[tcp.argmin()]
            jc = j[tcp.argmin()]

            vi = self._velocity2[ic]
            vj = self._velocity2[jc]
            dif_v = vi - vj
        
        
        return tcp_p, ic, jc, dif_v, Interaction
        
    def md_step(self):
        
        #print('Simul::md_step')
        
        tc,disk,direction, Int = self._wall_time()
        tcp_p, ic, jc, dif_v, Interaction = self._pair_time()        
        tmin = min(tc[disk, direction], tcp_p)
        
        pressure = 0
        
        now = 0
        while now + tmin < self._sample_time:
            now+=tmin
            
            if tc[disk, direction] < tcp_p:
                #print("Choque pared")
                if Int == 1 :
                    #print(disk, direction)
                    self.position1 += tc[disk, direction]*self._velocity1
                    self._velocity1[disk, direction] = -self._velocity1[disk, direction]   
                    #print("Azul")
                else :
                    self.position2 += tc[disk, direction]*self._velocity2
                    self._velocity2[disk, direction] = -self._velocity2[disk, direction]
                    #print("Roja")
            else:
                #print("Particulas")
                if Interaction == 11 :
                    self.position1 += tcp_p* self._velocity1
                    
                    u = (self.position1[jc] - self.position1[ic])/np.linalg.norm(self.position1[jc] - self.position1[ic])

                    self._velocity1[ic] -= u*np.sum(u*dif_v)
                    self._velocity1[jc] += u*np.sum(u*dif_v)
                    #print("Azul-Azul")
                
                elif Interaction == 12 :
                    self.position1 += tcp_p* self._velocity1
                    self.position2 += tcp_p* self._velocity2
                    
                    u = (self.position2[jc] - self.position1[ic])/np.linalg.norm(self.position2[jc] - self.position1[ic])

                    self._velocity1[ic] -= u*np.sum(u*dif_v)
                    self._velocity2[jc] += u*np.sum(u*dif_v)
                    #print("Azul-Rojo")
                
                elif Interaction == 21 :
                    self.position1 += tcp_p* self._velocity1
                    self.position2 += tcp_p* self._velocity2
                    
                    u = (self.position1[jc] - self.position2[ic])/np.linalg.norm(self.position1[jc] - self.position2[ic])

                    self._velocity2[ic] -= u*np.sum(u*dif_v)
                    self._velocity1[jc] += u*np.sum(u*dif_v)
                    #print("Rojo-Azul")
                    
                else :
                    self.position2 += tcp_p* self._velocity2
                    
                    u = (self.position2[jc] - self.position2[ic])/np.linalg.norm(self.position2[jc] - self.position2[ic])

                    self._velocity2[ic] -= u*np.sum(u*dif_v)
                    self._velocity2[jc] += u*np.sum(u*dif_v)
                   # print("Rojo-Rojo")
                    
            tc, disk, direction, Int = self._wall_time()
            tcp_p, ic, jc, dif_v, Interaction = self._pair_time()
            tmin = min(tc[disk, direction], tcp_p)   
               
        self.position1 += (self._sample_time-now) * self._velocity1
        self.position2 += (self._sample_time-now) * self._velocity2
        
        return pressure

    def __str__(self):   # this is used to print the position and velocity of the particles

        return 'Se logró :D'