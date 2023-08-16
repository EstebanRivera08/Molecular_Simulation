import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.collections import PatchCollection, EllipseCollection
from matplotlib.patches import Ellipse, Rectangle, Circle
import numpy as np


class circulos:
    def __init__(self, radii, position):
        
        self.angles = np.linspace(0, 2*np.pi, 100)
        self.y = np.sin(angles)[np.newaxis] + position[:,0]
        self.x = np.cos(angles)[np.newaxis] + position[:,1]
    
    def 
        

class Simul:
    def __init__(self, sample_time, c , a, b, Positions):
        
        self.position = Positions
        #self.x = Positions[:,0]
        #self.y = Positions[:,1]
        self.sample_time = sample_time
        self.c = c
        self.sigma = 0
        self.boite = np.array([[a,b]])

    def md_step(self):
        
        self.sigma += self.c*self.sample_time
        

class AnimateSimul:
    def __init__(self, simulation):

        self.simulation = simulation
        self.fig, self.ax = plt.subplots(figsize=(5, 5))  # initialise  graphics
        
        #self.circles = Ellipse(simulation.position, 2*simulation.sigma,  2*simulation.sigma)  # circles at pos
        #self.ax.add_patch(self.circles)
        
        #patches = [Circle((x_, y_), s_) for x_, y_, s_ in np.broadcast(simulation.x, simulation.y, simulation.sample_time*simulation.c)]
        #collection = PatchCollection(patches, color = 'k', fc = 'none')
        
        self.circles = EllipseCollection(widths=2*simulation.sigma, heights=2*simulation.sigma, angles=0, units='x',
                                         offsets=simulation.position, transOffset=self.ax.transData)  # circles at pos
        
        self.ax.add_collection(self.circles)
        #self.ax.add_collection(collection)

        rect = Rectangle((0, 0), simulation.boite[0,0], simulation.boite[0,1], ec='black', facecolor='none')   #   enclosing box
        self.ax.add_patch(rect)

        self.ax.set_xlim(left=-100, right=simulation.boite[0,0] + 100)
        self.ax.set_ylim(bottom=-100, top=simulation.boite[0,1] + 100)
        
    def _anim_step(self, m):  # m is the number of calls tshat have occurred to this function
        
        p = self.simulation.md_step()  # perform simulation step
        print('anim_step m = ', m)
        #self.circles.set_width(self.simulation.sigma) 
        #self.circles.set_height(self.simulation.sigma) 
        self.circles = EllipseCollection(widths=2*self.simulation.sigma, heights=2*self.simulation.sigma, angles=0, units='x',
                                         offsets=self.simulation.position, transOffset=self.ax.transData)  # circles at pos
        
        self.ax.add_collection(self.circles)
        #patches = [Circle((x_, y_), s_) for x_, y_, s_ in np.broadcast(self.simulation.x, self.simulation.y, self.simulation.sample_time*self.simulation.c)]
        #collection = PatchCollection(patches, color = 'k', fc = 'none')
        #self.ax.add_collection(collection)

        
    def go(self, nframes):
        
        self._ani = animation.FuncAnimation(self.fig, func=self._anim_step, frames=nframes,
                                      repeat=False, interval=20)  # run animation
        plt.show()

        
def virtual(a,b,x0,y0,order) :  #base; height; source
    
    aux_x = lambda x0, y0 : np.array([[x0 - (2*x0), y0],[x0 + 2*(a-x0), y0]])
    aux_y = lambda x0, y0 : np.array([[x0, y0 - (2*y0)],[x0, y0 + 2*(b-y0)]])
    
    virtual_x = aux_x(x0,y0)
    virtual_y = aux_y(x0,y0)
    point1 = 0
    point2 = len(virtual_x)

    if order > 1 :
        for i in range(order) :
            for j in range(point1,point2) :
                virtual_x = np.concatenate((virtual_x,aux_x(virtual_x[j,0],virtual_x[j,1])), axis = 0)
                virtual_x = np.concatenate((virtual_x,aux_x(virtual_y[j,0],virtual_y[j,1])), axis = 0)
                virtual_y = np.concatenate((virtual_y,aux_y(virtual_y[j,0],virtual_y[j,1])), axis = 0)
                virtual_y = np.concatenate((virtual_y,aux_y(virtual_x[j,0],virtual_x[j,1])), axis = 0)
                virtual_x = np.delete(virtual_x,np.where((virtual_x == [[x0,y0]]).all(axis = 1)), axis = 0)
                virtual_y = np.delete(virtual_y,np.where((virtual_y == [[x0,y0]]).all(axis = 1)), axis = 0)

            point1 = point2
            point2 = len(virtual_x)
        
    dib = np.concatenate((virtual_x[point1:point2,:],virtual_y[point1:point2,:]), axis = 0)    
    return dib

def main():
    
    a = 20
    b = 30
    x0 = 15
    y0 = 20
    order = 6
    c = 2400*100*10e-6 # 2400m/s - > cm/micro s
    sample_time = 100
    Positions = virtual(a,b,x0,y0,order)
    
    simulation = Simul(sample_time , c , a, b, Positions)  #  sigma particle radius
    animate = AnimateSimul(simulation)
    animate.go(nframes= 500)
    print(simulation)  #  print last configuration to screen
    
if __name__ == '__main__':
    main()
    
    
