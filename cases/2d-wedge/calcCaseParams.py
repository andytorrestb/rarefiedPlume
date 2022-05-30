import math
import numpy as np

# Math constants
pi = math.pi
p = 3 # precision

# Bird Breakdown Parameter
B = 0.0025

# Gas Properties
gamma = 5.0/3.0
R = 208 # ??
d = 4.17e-10 # [m], molecular diameter

# Sonic Properties (@ throat of nozzle)
Rt = 1.0 # radius [m]
Mt = 1.0 # Mach number
Tt = 1000 # Temperature [k]

# =========== Solve for sonic conditions ============
# Numerical Density [molecules / m ^3]
nt = 0.5 * math.sqrt(gamma / pi) * (Mt / (B * d * d * Rt))

# Magnitude of Radial Velocity  r_hat [m/s]
vt = round(Mt * math.sqrt(gamma * R * Tt), p)

# Mean free path at inflow 
mfpt = round((5 * Mt) / (2 * d * d * vt * nt) * math.sqrt((3 * gamma * R * Tt) / pi), 3)

print("                                 Sonic conditions")
print()
print('velocity =', vt, 'm/s, numDensity =', '{:0.3e}'.format(nt), 'molecules/m^3, mfp =', mfpt, 'm')
print()

# ============ Sove for inflow conditions ============
def TTt_Ma(Ma, ga = 1.67):
    return (ga + 1) / (2 + (ga - 1) * Ma ** 2)

def rrt_Ma(Ma, ga = 1.67):
    rrt = (1 / TTt_Ma(Ma, ga)) ** ((ga + 1) / (ga - 1))
    return np.sqrt(np.sqrt(rrt) / Ma) 

def nnt_Ma(Ma, ga = 1.67):
    return TTt_Ma(Ma, ga) ** (1 / (ga + 1))

def a(T, ga = 1.67, R = 208):
    return np.sqrt(ga * R * T)
r_i = 1.05
r_o = 20
M = Mt

# Iterate to find an acceptable M value
while rrt_Ma(M) - r_i < 0:
    # print(rrt_Ma(M) - r_i)
    M = M + 1.0e-3

# print(rrt_Ma(M) - r_i)
# print(M)

T = TTt_Ma(M) * Tt
v = M * a(T)
n = nnt_Ma(M) * nt
mfp = round((5 * M) / (2 * d * d * v * n) * math.sqrt((3 * gamma * R * T) / pi), 3)

print("                                 Inflow Conditions")
print()
print('velocity =', v, 'm/s, numDensity =', '{:0.3e}'.format(n), 'molecules/m^3, mfp =', mfp, 'm')
print()

# ============ Sove for mesh resolution settings ============
dx = mfp / 10
# dx = 0.08 / 1
x_domain = r_o - r_i
n_x_cells = int(x_domain / dx)

print("                                 Mesh Resolution Settings")
print('Required cells for x-direction =', n_x_cells)
print('Required cells for y-direction =', 4 * n_x_cells)

# ============ Sove for molecule resolution settings =========
