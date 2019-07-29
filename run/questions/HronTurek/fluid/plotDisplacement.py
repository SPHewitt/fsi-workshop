from matplotlib import pyplot as plt
import numpy as np


if __name__ == '__main__':
    print("\n ---- Plotting Displacements ----\n")

    # File name
    filename='history/0/point_487.dat'

    data=np.loadtxt(filename,skiprows=1)

    time=data[:,0]
    dispX = data[:,1]
    dispY = data[:,2]


    plt.xlabel('Time (s)')
    plt.ylabel('Displacment, X (m)')
    plt.grid(True)
    plt.xlim(2,5)
    plt.plot(time,dispX,'-')
    plt.savefig('dispX.pdf',bbox_inches='tight')
    plt.close()

    plt.xlabel('Time (s)')
    plt.ylabel('Displacement, Y (m)')
    plt.grid(True)
    plt.xlim(2,5)
    plt.plot(time,dispY,'-')
    plt.savefig('dispY.pdf',bbox_inches='tight')
    plt.close()


