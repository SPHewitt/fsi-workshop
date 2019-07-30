from matplotlib import pyplot as plt
import re
import numpy as np


if __name__ == '__main__':
    print("\n ---- Plotting forces ----\n")

    # File name
    filename='history/0/force.dat'


    data=np.loadtxt(filename)
    time = data[:,0]
    forceX = data[:,1]
    forceY = data[:,2]
    
    # Plotting routines
    # Uncomment the three lines below and
    # Try changing some of the parameters 
    
    #plt.title('My Title')
    #plt.xlim(8,10)
    #plt.ylim(-200,200)
    plt.xlabel('Time (s)')
    plt.ylabel('Force, X (N)')
    plt.grid(True)
    plt.plot(time,forceX,'-')
    plt.savefig('forceX.pdf',bbox_inches='tight')
    plt.close()

    #plt.title('My Title')
    #plt.xlim(8,10)
    #plt.ylim(-200,200)
    plt.xlabel('Time (s)')
    plt.ylabel('Force, Y (N)')
    plt.grid(True)
    plt.plot(time,forceY,'-')
    plt.savefig('forceY.pdf',bbox_inches='tight')
    plt.close()


