from pylab import *
from scipy.integrate import solve_ivp
from matplotlib import cm
import matplotlib.animation as animation
from progress.bar import Bar
from time import sleep
import time
from matplotlib import animation

"""Variables for different parameters in basilar membrane"""
H = 0.67  # cm
l = 3.5  # cm
rho = 1  # g/cm^3
c = 1.43*10**5  # cm/s
m = 0.143  # g/cm^2
A = 0.2  # Amplitude in pascal

"""
500Hz=0.2s, 2000Hz=0.05s, 5000Hz=0.02s, 10000Hz=0.01 ;Precision in x
500Hz=0.002s, 2000Hz=0.0005s, 5000Hz=0.0002s, 10000Hz=0.0001 ;Precision in a amplitude = 0.1
Setting arrays of x and t, defining steps
"""
f = list(map(int, input("List of frequencies in [Hz]: ").split()))
N = input("Number of oscillators [N]:  ")
Time = input("Time duration in seconds [s]:  ")
Time = float(Time)
x = linspace(0.0, l, N)
dx = x[1] - x[0]
dt = dx/(c*2)
M = int(Time/dt)
t = linspace(0, Time, M)
k = 4e9*exp(-4*x)

start = time.time()


class Basilarmembranen:

    def __init__(self, x, t, c, k, dt, dx):

        self.x = x
        self.t = t
        self.c = c
        self.k = k
        self.dt = dt
        self.dx = dx

    def maxInRows(self, f, N):

        x = self.x
        t = self.t
        c = self.c
        dt = self.dt
        dx = self.dx

        num1 = len(t)
        num2 = len(x)
        n = zeros((num2, num1))
        u = zeros((num2, num1))

        # Defining constant for easy interpreting
        C = c * dt/dx
        alpha = (2*rho)/(H)
        beta = dt**2/m
        print ("--------------------------------------------------------------------------------")
        print ('Length of M  = ', num1)
        print ('Length of N  = ', num2)
        print ('C = ', C)
        print ('beta =', beta)
        print ('alpha =', alpha)
        print ('dx =', dx)
        print ('dt =', dt)
        print ("Time: ", Time)
        """
        widgets = ['Progress: ', Percentage(), ' ', Bar(marker='0',left='[',right=']'),
                   ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options
                   """
        pbar = Bar('Loading', fill='#', suffix='%(percent)d%%')
        pbar.start()

        # For loop for N X M matrix
        for j in range(num1-1):
            # Initial force, applied from staples.
            u[0, j] = A*sin(t[j]*2*pi*f)
            pbar.next(j+1)
            time.sleep(0.1)
            for i in range(num2-1):
                k[i] = 4e9*exp(-4*x[i])
                # my statement consisting of n, amplitude of membrane wave
                n[i, j+1] = beta*(u[i, j] - k[i]*n[i, j]) + 2*n[i, j] - n[i, j-1]
                # my statement consisting of u, pressure difference
                u[i, j+1] = C**2*(u[i+1, j] - 2*u[i, j] + u[i-1, j]) + 2*u[i, j]\
                    - u[i, j-1] - alpha * (n[i, j+1] - 2*n[i, j] + n[i, j-1])
        pbar.finish()
        print ("--------------------------------------------------------------------------------")

        maxInRows = amax(n, axis=1)
        normalized = maxInRows / maxInRows.max()

        return normalized

    def Displacement(self, f, N):

        x = self.x
        t = self.t
        c = self.c
        dt = self.dt
        dx = self.dx

        num1 = len(t)
        num2 = len(x)
        n = zeros((num2, num1))
        u = zeros((num2, num1))

        # Defining constant for easy interpreting
        C = c * dt/dx
        alpha = (2*rho)/(H)
        beta = dt**2/m
        print ("--------------------------------------------------------------------------------")
        print ('Length of M  = ', num1)
        print ('Length of N  = ', num2)
        print ('C = ', C)
        print ('beta =', beta)
        print ('alpha =', alpha)
        print ('dx =', dx)
        print ('dt =', dt)
        print ("Time: ", Time)
        """
        bar = progressbar.ProgressBar(num1-1,
                                      widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
        """
        """
        widgets = ['Progress: ', Percentage(), ' ', Bar(marker='2',left='[',right=']'),
                   ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options
                   """
        pbar = Bar('Loading', fill='#', suffix='%(percent)d%%')
        pbar.start()

        # For loop for N X M matrix
        for j in range(num1-1):
            # Initial force, applied from staples.
            u[0, j] = A*sin(t[j]*2*pi*f)
            pbar.next(j+1)
            time.sleep(0.1)
            for i in range(num2-1):
                k[i] = 4e9*exp(-4*x[i])
                # my statement consisting of n, amplitude of membrane wave
                n[i, j+1] = beta*(u[i, j] - k[i]*n[i, j]) + 2*n[i, j] - n[i, j-1]
                # my statement consisting of u, pressure difference
                u[i, j+1] = C**2*(u[i+1, j] - 2*u[i, j] + u[i-1, j]) + 2*u[i, j]\
                    - u[i, j-1] - alpha * (n[i, j+1] - 2*n[i, j] + n[i, j-1])
        pbar.finish()
        print ("--------------------------------------------------------------------------------")
        return n[:, M-1]


def figure_maxInRows():
    #legend(loc='upper left')
    figure(1)
    for i in range(len(f)):
        print ("[Resonance] INITIATING CALCULATIONS FOR FREQUENCY:", f[i], "Hz")
        plot(x, solutions.maxInRows(f[i], N), linewidth=1, label='f=%.f Hz' % f[i])
    legend(loc='best')
    title('Resonance in basilar membrane')
    xlabel('Position along membrane [cm]')
    ylabel('Normalized amplitude displacement')
    tight_layout()
    grid()


def figure_Displacement():
    #legend(loc='upper left')
    for i in range(len(f)):
        figure(i+2)
        print ("[Displacement] INITIATING CALCULATIONS FOR FREQUENCY:", f[i], "Hz")
        plot(x, solutions.Displacement(f[i], N), linewidth=1, label='f=%.f Hz' % f[i])
        legend(loc='best')
        title('Basilar membrane displacement')
        xlabel('Position along membrane [cm]')
        ylabel('Amplitude basilar membrane [cm]')
        tight_layout()
        grid()


solutions = Basilarmembranen(x, t, c, k, dt, dx)

figure_maxInRows()
figure_Displacement()

show()
end = time.time()
print ("Code time:", end - start)
