import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import root_scalar
import matplotlib
matplotlib.use('Qt5Agg')

inf = 5e-8

raiop = 10
# refraction index 
def n(r):
    return  1 + .05 * (r>raiop) * np.exp(-(r - raiop) ) + (r<=raiop) * .05
# atenuation constant 
def a(r):
   return   0.5 * np.exp(-(r - raiop))
#
# dtheta/dr is derived from Euler Lagrange equation of 
# Lagrangian L = n(r) * ds, where ds = sqrt(1+(dtheta/dr)**2) * dr
#

def dtheta_dr(theta, r, h):
    return (-1/r) / np.sqrt(max(inf, (n(r) * r / h)**2 - 1))

# Parameters
r0 = 15.5

alpha =np.pi/2 -  np.pi / 6

# Angles of incidence in radians
rinit_values = np.linspace(r0, 2*r0, 1000)

# Plot the results in polar coordinates with background color representing n(r)
plt.figure(figsize=(10, 10))
ax = plt.subplot(111, polar=True)

# Create polar heatmap of n(r)
r_values_heatmap = np.linspace(0, 2*r0, 100)
theta_values_heatmap = np.linspace(0, 2 * np.pi, len(r_values_heatmap))
R, th = np.meshgrid(r_values_heatmap, theta_values_heatmap)
z = n(R) / 2
ax.pcolormesh(th, R, z, cmap='YlGnBu')
# PLotting the figure
ax.plot(th, R, color='b', ls='none') 

for r0 in rinit_values:  # Angle of incidence in radians

    # Calculate from contour conditions
    h = n(r0) * r0 * np.cos(alpha)

    # Calculate rclosest using a solver
    def equation_to_solve(rclosest):
        return n(rclosest) * rclosest - h


    rclosest_guess = h
    solution = root_scalar(equation_to_solve, x0=rclosest_guess)
    rclosest = solution.root
    residue = equation_to_solve(rclosest)

    #  print(f"Integration param: {h}\nrclosest: {rclosest}\n Residue: {residue}" )
    
    if abs(residue) < 0.001 and rclosest > raiop:
        # Solve for theta using odeint


        # varing r values (descending and ascending)

        r_values = np.linspace(r0, rclosest+inf, 1000)
        r_values2 = np.linspace(rclosest+inf, r0, 1000)

        theta_initial = 0.0

        # Adjust tolerance parameters for the solver
        rtol = 1e-6
        atol = 1e-8

        def combined_equations(y, r, h):
            theta, absorbance = y
            dsdt = (n(r) * r / h)**2 - 1
            dthetadr = (-1/r) / np.sqrt(max(inf, dsdt))
            dadt = -a(r) * np.sqrt(1 + (r * dthetadr)**2)
            return [dthetadr, dadt]

        def combined_equations2(y, r, h):
            theta, absorbance = y
            dsdt = (n(r) * r / h)**2 - 1
            dthetadr = (1/r) / np.sqrt(max(inf, dsdt))
            dadt = a(r) * np.sqrt(1 + (r * dthetadr)**2)
            return [dthetadr, dadt]

        #integrate to obtain angles and absorbances, descending r
        solutions = odeint(combined_equations, [theta_initial, 0.0], r_values, args=(h,), rtol=rtol, atol=atol)
     
        theta_values, absorbance_values = solutions.T

        # Convert to polar coordinates
        theta_values = theta_values.flatten()
        absorbance_values = absorbance_values.flatten()
       
        # integrate to obtain angles and absorbances, ascending r
        solutions = odeint(combined_equations2, [theta_values[-1], absorbance_values[-1]], r_values2, args=(h,), rtol=rtol, atol=atol)
        
        theta_values2, absorbance_values2 = solutions.T

        # Convert to polar coordinates
        theta_values2 = theta_values2.flatten()
        absorbance_values2 = absorbance_values2.flatten()
       
 
        totalabs = np.exp(-absorbance_values2[-1])

        print(f"ABsorbance %0.3f" % (absorbance_values2[-1]))

        ax.plot(np.concatenate((theta_values,theta_values2)), np.concatenate((r_values, r_values2)), color=(1-totalabs,0.5,totalabs))


ax.set_title('Path of Light Ray Through Atmosphere')
ax.legend()
ax.grid(True)
plt.show()

