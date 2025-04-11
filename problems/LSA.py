import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
import ternary
import matplotlib.colors as mcolors
import matplotlib.cm as cm

# Python code to perform linear stability analysis for 3 phase cahn hilliard equation

# Define free energy functions and gradients

def f(phi1, phi2, params):
    chi12, chi13, chi23 = params['chi12'], params['chi13'], params['chi23']
    N1, N2, N3 = params['N1'], params['N2'], params['N3']

    # phi1 = max(1e-8, min(phi1, 1 - 1e-8))
    # phi2 = max(1e-8, min(phi2, 1 - 1e-8))
    phi3 = 1 - phi1 - phi2

    f = phi1 * math.log(phi1) / N1 + phi2 * math.log(phi2) / N2 + (phi3) * math.log(phi3) / N3 + chi12 * phi1 * phi2 + chi13 * phi1 * phi3 + chi23 * phi2 * phi3
    return f

def df(phi1, phi2, params):
    chi12, chi13, chi23 = params['chi12'], params['chi13'], params['chi23']
    N1, N2, N3 = params['N1'], params['N2'], params['N3']
    phi3 = 1 - phi1 - phi2

    df1 = math.log(phi1) / N1 + 1.0 / N1 - math.log(phi3) / N3 - 1.0 / N3 + chi12 * phi2 - chi13 * phi1 + chi13 * phi3 - chi23 * phi2
    df2 = math.log(phi2) / N2 + 1.0 / N2 - math.log(phi3) / N3 - 1.0 / N3 + chi12 * phi1 - chi13 * phi1 - chi23 * phi2 + chi23 * phi3 
    
    return [df1, df2]

def d2f(phi1, phi2, params):
    chi12, chi13, chi23 = params['chi12'], params['chi13'], params['chi23']
    N1, N2, N3 = params['N1'], params['N2'], params['N3']
    phi3 = 1 - phi1 - phi2

    d2f11 = 1.0 / (N1 * phi1) + 1.0 / (N3 * phi3) - 2.0 * chi13
    d2f12 = 1.0 / (N3 * phi3) + chi12 - chi13 - chi23
    d2f22 = 1.0 / (N2 * phi2) + 1.0 / (N3 * phi3) - 2.0 * chi23

    return [[d2f11, d2f12], [d2f12, d2f22]]


# Compute eigen values
def compute_eigen_values(k, phi1, phi2, chemistry_params, simulation_params):
    d2f_values = d2f(phi1, phi2, chemistry_params)
    Cn = simulation_params['Cn']

    a11 = d2f_values[0][0]
    a12 = d2f_values[0][1]
    a22 = d2f_values[1][1]

    x = k * math.pi
    B = a11 + a22 - math.sqrt((a11 - a22)**2 + 4 * a12**2)
    C = a11 + a22 + math.sqrt((a11 - a22)**2 + 4 * a12**2)
    b = Cn**2

    lambda1 = - 0.5 * x**2 * (B + 2 * b * x**2)
    lambda2 = - 0.5 * x**2 * (C + 2 * b * x**2)

    return [lambda1, lambda2]

# Chemistry parameters
chemistry_params = {
    'chi12': 1.0,
    'chi13': 0.3,
    'chi23': 0.3,
    'N1': 25.0,
    'N2': 25.0,
    'N3': 1.0
}

# Simulation parameters
simulation_params = {
    'Cn': 1e-5,
    'Ly': 1.0,
    'Nelemy': 150
}

# Compute lambda values
k_min = 1
h = simulation_params['Ly'] / simulation_params['Nelemy']
k_max = math.pi / h

N = 50

k_vec = np.linspace(k_min, k_max, N)
lambda_vec = np.zeros((2,N))

for i in range(N):
    lambda_vec[:,i] = compute_eigen_values(k_vec[i], 0.9, 0.10 - 1e-5, chemistry_params, simulation_params)

#print(lambda_vec[0,:])


# Plot the eigen values
plt.figure()
plt.plot(k_vec, lambda_vec[0,:], label='lambda1')
plt.plot(k_vec, lambda_vec[1,:], label='lambda2')
# Draw a horizontal line at y=0
plt.axhline(y=0, color='k', linestyle='--')
plt.legend()
plt.xlabel('$k$')
plt.ylabel('$\lambda$')
plt.title('Eigen values vs k')
plt.show()
plt.close()


# Create a ternery contour plot for lambda1

# Define the grid
n = 50
phi1_vec = np.linspace(0.01, 0.99, n)
phi2_vec = np.linspace(0.01, 0.99, n)

# Store ternary data
ternary_data = []
lambda_values = []

for phi1 in phi1_vec:
    for phi2 in phi2_vec:
        phi3 = 1 - phi1 - phi2
        if phi3 > 0:  # Ensure valid ternary condition
            lambda1 = compute_eigen_values(k_min, phi1, phi2, chemistry_params, simulation_params)[0]
            ternary_data.append((phi1, phi2, phi3))
            lambda_values.append((lambda1))

# Convert to numpy arrays
ternary_data = np.array(ternary_data)
lambda_values = np.array(lambda_values)

print(lambda_values)

# Create a colormap and normalize
cmap = cm.viridis  # Choose a colormap
norm = mcolors.Normalize(vmin=np.min(lambda_values), vmax=np.max(lambda_values))
sm = cm.ScalarMappable(norm=norm, cmap=cmap)  # ScalarMappable for colorbar

# Initialize ternary plot
figure, tax = ternary.figure(scale=1.0)
tax.boundary(linewidth=2.0)
tax.gridlines(color="gray", multiple=0.1)

# Scatter plot with color mapping
tax.scatter(ternary_data, c=lambda_values, cmap=cmap, s=15)

# Labels
tax.left_axis_label("$\phi_1$", fontsize=12)
tax.right_axis_label("$\phi_2$", fontsize=12)
tax.bottom_axis_label("$\phi_3$", fontsize=12)

# Add colorbar using the ScalarMappable
cbar = plt.colorbar(sm, ax=tax.ax)
cbar.set_label("$\lambda_1$", fontsize=12)

# Title
plt.title("Ternary Contour Plot of $\lambda_1$")
plt.show()





