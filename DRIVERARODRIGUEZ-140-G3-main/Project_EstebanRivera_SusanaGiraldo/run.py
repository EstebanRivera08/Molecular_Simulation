from simul import Simul
from animatesimul import AnimateSimul
import numpy as np

def main():
    np.random.seed(10)
    simulation = Simul(sample_time=0.01, sigma=0.01, L = 1, N = 100)  #  sigma particle radius
    print(simulation.__doc__)  # print the documentation from the class

    animate = AnimateSimul(simulation)
    animate.go(nframes= 500)
    print(simulation)  #  print last configuration to screen


if __name__ == '__main__':
    main()
    
    
    
