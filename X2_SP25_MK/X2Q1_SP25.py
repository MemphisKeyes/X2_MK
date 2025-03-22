# region imports
from scipy.integrate import solve_ivp, quad
import numpy as np
import matplotlib.pyplot as plt
# endregion

# region function definitions
def S(x):
    """
    Computes the Fresnel integral S(x) = ∫₀ˣ sin(t²) dt using quad.
    """
    s = quad(lambda t: np.sin(t**2), 0, x)  # the solution for S(x) found using quad
    return s[0]

def Exact(x):
    """
    Computes the exact solution: y = 1 / (2.5 - S(x)) + 0.01 * x²
    """
    return 1 / (2.5 - S(x)) + 0.01 * x**2  # exact solution for i.v.p. at x

def ODE_System(x, y):
    """
    Defines the differential equation y' = (y - 0.01x²)² * sin(x²) + 0.02x
    """
    Y = y[0]  # rename first state variable into convenient name
    Ydot = (Y - 0.01 * x**2)**2 * np.sin(x**2) + 0.02 * x  # calculate derivatives of state variable(s)
    return [Ydot]

def PlotResults(*args):
    """
    Plots the exact and numerical solutions according to problem formatting instructions.
    """
    xRange_Num, y_Num, xRange_Xct, y_Xct = args  # unpack args
    plt.plot(xRange_Xct, y_Xct, 'k-', label='Exact')  # solid line for exact solution
    plt.plot(xRange_Num, y_Num, 'k^', label='Numerical')  # triangles for numerical solution

    plt.xlabel('x')
    plt.ylabel('y')
    plt.xticks(np.arange(0, 6.1, 1.0))
    plt.yticks(np.arange(0, 1.1, 0.2))
    plt.tick_params(axis='both', direction='in', top=True, right=True)
    plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: f'{val:.1f}'))
    plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: f'{val:.1f}'))
    plt.xlim([0, 6])
    plt.ylim([0, 1])
    plt.title(r"IVP: $y' = (y - 0.01x^2)^2 \sin(x^2) + 0.02x,\ y(0)=0.4$")
    plt.legend()
    plt.grid(True)
    plt.show()

def main():
    """
    Solves the IVP: y' = (y - 0.01x²)²*sin(x²) + 0.02x with y(0) = 0.4 using solve_ivp.
    Plots the exact and numerical solutions.
    """
    xRange = np.arange(0, 5.1, 0.2)  # create a numpy array for the x range to evaluate numerical solution (h=0.2)
    xRange_xct = np.linspace(0,5,500)  # create a numpy array for the x range for the exact solution
    Y0 = [0.4]  # create initial conditions
    sln = solve_ivp(ODE_System, [0,5], Y0, t_eval=xRange)  # numerically solve i.v.p. with default RK45 method
    xctSln = np.array([Exact(x) for x in xRange_xct])  # produce array of y values for exact solution
    PlotResults(xRange, sln.y[0], xRange_xct, xctSln)  # call the plotting function to produce the required plot
    pass
# end region

# region function calls
if __name__ ==  "__main__":
    main()
# end region
