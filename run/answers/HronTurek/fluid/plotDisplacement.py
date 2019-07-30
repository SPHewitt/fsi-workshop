from matplotlib import pyplot as plt
import numpy as np


if __name__ == '__main__':
    print("\n ---- Plotting Displacements ----\n")

    # File name
    filename='history/0/point_286.dat'

    # Read in data file
    data=np.loadtxt(filename,skiprows=1)

    # Seperate data
    time=data[:,0]
    dispX = data[:,1]
    dispY = data[:,2]

    # Plotting routines
    # Uncomment the three lines below and
    # Try changing some of the parameters 
    
    #plt.title('My Title')
    #plt.xlim(8,10)
    #plt.ylim(-0.03,0.01)
    plt.xlabel('Time (s)')
    plt.ylabel('Displacment, X (m)')
    plt.grid(True)
    plt.plot(time,dispX,'-')
    plt.savefig('dispX.pdf',bbox_inches='tight')
    plt.close()

    #plt.title('My Title')
    #plt.xlim(8,10)
    #plt.ylim(-0.08,0.08)
    plt.xlabel('Time (s)')
    plt.ylabel('Displacement, Y (m)')
    plt.grid(True)
    plt.plot(time,dispY,'-')
    plt.savefig('dispY.pdf',bbox_inches='tight')
    plt.close()


