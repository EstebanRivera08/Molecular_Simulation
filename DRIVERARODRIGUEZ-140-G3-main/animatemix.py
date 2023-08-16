import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.collections import EllipseCollection
from matplotlib.patches import Rectangle
from Mixing import Mix


class AnimateMix:
    def __init__(self, simulation):
        print('AnimateSimul')
        self.simulation = simulation
        
        self.fig, self.ax = plt.subplots(figsize=(10, 5))  # initialise  graphics
        
        self.circles1 = EllipseCollection(widths=2*simulation.sigma1, heights=2*simulation.sigma1, angles=0, units='x',
                                         offsets=simulation.position1, transOffset = self.ax.transData)  # circles at pos
        
        self.circles2 = EllipseCollection(widths=2*simulation.sigma2, heights=2*simulation.sigma2, angles=0, units='x',
                                         offsets=simulation.position2, color = 'r', transOffset=self.ax.transData)  # circles at pos
        
        self.ax.add_collection(self.circles1)
        self.ax.add_collection(self.circles2)
        
        rect = Rectangle((0, 0), simulation.boite[0,0], simulation.boite[0,1], ec='black', facecolor='none')   #   enclosing box
        
        self.ax.set_xlim(left=-.1, right=simulation.boite[0,0] + .1)
        self.ax.set_ylim(bottom=-.1, top=simulation.boite[0,1] + .1)
        self.ax.add_patch(rect)

    def _anim_step(self, m):  # m is the number of calls that have occurred to this function
        print('anim_step m = ', m)
        print('time',(m+1)*0.01)
        
        self.simulation.md_step()  # perform simulation step
        
#  calculation a window containing the particles and reset axes
      #  rmin = np.amin(self.simulation.position) - self.simulation.sigma/2. - .1
       # rmax = np.amax(self.simulation.position) + self.simulation.sigma/2. + .1
      #  rmin = min(-.1, rmin)
      #  rmax = max(1.1, rmax)
      #  self.ax.set_xlim(left=rmin, right=rmax)
      #  self.ax.set_ylim(bottom=rmin, top=rmax)
# signal graphics update
        self.circles1.set_offsets(self.simulation.position1)  #Azules
        self.circles2.set_offsets(self.simulation.position2)  #Rojas

    def go(self, nframes):
        print('go')
        self._ani = animation.FuncAnimation(self.fig, func=self._anim_step, frames=nframes,
                                      repeat=False, interval=20)  # run animation
        plt.show()
