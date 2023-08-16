from animatemix import AnimateMix
from Mixing import Mix
import numpy as np

def main():
    np.random.seed(10)
    simulation = Mix(sample_time=0.01, sigma1 = 0.02, sigma2 = 0.02, Lx = 2, Ly = 1, N = 20)  #  sigma particle radius
    print(simulation.__doc__)  # print the documentation from the class

    animate = AnimateMix(simulation)
    animate.go(nframes= 500)
    print(simulation)  #  print last configuration to screen


if __name__ == '__main__':
    main()
    
    
    
