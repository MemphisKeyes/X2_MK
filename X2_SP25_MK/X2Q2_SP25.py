# region imports
from scipy.integrate import solve_ivp
from math import sin
import math
import numpy as np
from matplotlib import pyplot as plt
# endregion

#region class definitions
class circuit():
    def __init__(self, R=20, L=20, C=0.05, A=20, w=20, p=0):
        '''
        Initializes the RLC circuit with given parameters.
        :param R: Resistance in Ohms
        :param L: Inductance in Henrys
        :param C: Capacitance in Farads
        :param A: Amplitude of input voltage in Volts
        :param w: Frequency of input voltage in rad/s
        :param p: Phase of input voltage in radians
        '''
        #region attributes
        self.R = R
        self.L = L
        self.C = C
        self.A = A
        self.w = w
        self.p = p
        self.t = None
        self.i1 = None
        self.i2 = None
        self.vc = None
        self.vin = None
        #endregion

    #region methods
    def ode_system(self, t, X):
        """
        this is the odeSystem callback I'm using for solve_ivp().
        :param X: the current values of the state variables
        :param t: the current time
        :return: list of derivatives of state variables
        """
        i1 = X[0]
        i2 = X[1]
        v_t = self.A * math.sin(self.w * t + self.p)
        di1dt = (v_t - self.R * (i1 - i2) - (1 / self.C) * i2) / self.L
        di2dt = (i1 - i2) / (self.R * self.C)
        return [di1dt, di2dt]

    def simulate(self, t=10, pts=500):
        """
        For simulating transient behavior of circuit.
        :param: time over which to carry out the simulation in seconds
        :param pts: number of points in the simulation
        :return: nothing, just store I
        """
        time_array = np.linspace(0, t, pts)
        sol = solve_ivp(self.ode_system, [0, t], [0, 0], t_eval=time_array)
        self.t = sol.t
        self.i1 = sol.y[0]
        self.i2 = sol.y[1]
        self.vc = (1 / self.C) * self.i2
        self.vin = self.A * np.sin(self.w * self.t + self.p)

    def doPlot(self, ax=None):
        """
        Re-written on 4/21/2022 to adapt to plotting on GUI if ax is not None
        :param args: contains ((R, list of time values, and results of solve_ivp))
        :param ax:
        :return:
        """
        if ax == None:
            fig, ax = plt.subplots()
            QTPlotting = False  # actually, we are just using CLI and showing the plot
        else:
            QTPlotting = True

        #JES MISSING CODE for making the plot
        ax.plot(self.t, self.i1, 'k-', label='i1(t)')
        ax.plot(self.t, self.i2, 'k--', label='i2(t)')
        ax.set_xlabel('t (s)')
        ax.set_ylabel('i1, i2 (A)', color='k')
        ax.tick_params(axis='y', labelcolor='k')
        ax.grid(True)

        ax2 = ax.twinx()
        ax2.plot(self.t, self.vin, 'k:', label='v(t)')
        ax2.set_ylabel('v(t) (V)', color='k')
        ax2.tick_params(axis='y', labelcolor='k')

        ax.legend(loc='upper left')
        ax2.legend(loc='upper right')

        if not QTPlotting:
            plt.tight_layout()
            plt.show()
    #endregion

#endregion

# region function definitions

def main():
    """
    For solving problem 2 on exam.
    :return:
    """
    goAgain = True
    Circuit = circuit(R=20, L=20, C=0.05, A=20, w=20, p=0)  # create a circuit object with default values
    while goAgain:
        #JES MISSING CODE for soliciting user input.
        try:
            R = float(input("Enter resistance R (Ohms): "))
            L = float(input("Enter inductance L (Henrys): "))
            C = float(input("Enter capacitance C (Farads): "))
            A = float(input("Enter amplitude A of v(t): "))
            w = float(input("Enter frequency w (rad/s): "))
            p = float(input("Enter phase p (radians): "))
            Circuit = circuit(R, L, C, A, w, p)
        except:
            print("Invalid input. Using previous values.")

        Circuit.simulate(t=10, pts=500)
        Circuit.doPlot()

        #JES MISSING CODE for soliciting user input.
        again = input("Simulate again with different values? (y/n): ").lower()
        goAgain = (again == 'y')
    pass
# endregion

# region function calls
if __name__ ==  "__main__":
    main()
# endregion
