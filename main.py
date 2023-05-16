# This is the supplementary information file provided along with submission of the
# paper titled "A method to measure effective flexural and transverse shear modulus of
# composite structures with large aspect ratio" to Experimental Mechanics
# This file contains the python source code to perform BVLS method on a set of data
# from 3-point bending test with different span length on unidirectional CF reinforced
# epoxy beam.
# article title: A method to measure effective flexural and transverse shear modulus
# of composite structures with large aspect ratio
# journal name: Experimental Mechanics
# author names: Yang Shen and David J. Branscomb
# affiliation: Highland Composites
# e-mail address of the corresponding author: yshen398@gatech.edu
import numpy as np
from scipy.optimize import lsq_linear
import matplotlib as mpl
import matplotlib.pyplot as plt

# cordpreg data
L = np.array([203.2,
152.4,
101.6,
88.9,
76.2,
63.5,
50.8,
38.1,
25.4,
12.7,
7.62,
]) * 1e-3
force_disp = np.array([8.73E+02,
2.09E+03,
6.83E+03,
1.04E+04,
1.54E+04,
2.65E+04,
5.27E+04,
1.12E+05,
2.95E+05,
9.75E+05,
1.37E+06,
])
disp_force = 1 / force_disp
M = np.array([L**3/48, L/4]).transpose()
r = 1.12e-3
I = np.pi * r**4 / 4
A = np.pi * r**2
v = 0.04
k = 6*(1+v)**2 / (7 + 12 * v + 4 * v**2)
boundE = [[1 / (200e9*I), 1 / (20e9*I)], [1 / (50e9*I), 1 / (20e9*I)], [1 /
(200e9*I), 1 / (20e9*I)]]
boundG = [[1 / (k*10e9*A), 1 / (k*1e9*A)], [1 / (k*10e9*A), 1 / (k*1e9*A)], [1 /
(k*20e9*A), 1 / (k*10e9*A)]]
Ef_list = []


G_list = []
for i in range(len(boundE)):
    res = lsq_linear(M, disp_force,
                    bounds=([boundE[i][0], boundG[i][0]], [boundE[i][1], boundG[i]
[1]]),
                    method="bvls", tol=1e-10, verbose=0)
    Ef = 1 / (res.x[0]*I)
    G = 1/ (k*res.x[1]*A)
    print(res.x)
    print("Ef: ", Ef, "G: ", G)
    Ef_list.append(Ef)
    G_list.append(G)
mpl.rcParams['font.size'] = '18'
fig, ax = plt.subplots()
# ax.plot(L, disp_force, marker='o')
ax.scatter(L, disp_force, marker='o', label='Experimental data', color='black')
ax.set_ylim([-0.00025, 0.002])
ax.set_xlabel("Span length (m)")
ax.set_ylabel("Deflection to force ratio (m/N)")
linestyle = ['dashed', 'dotted', (0, (3, 1, 1, 1))]
for i in range(len(Ef_list)):
    x = np.linspace(0, 210e-3, num=100)
    y = x**3 / (48*Ef_list[i]*I) + (6/5.0) * x / (4*G_list[i]*A)
    ax.plot(x, y, linestyle=linestyle[i],
    label='BVLS: Ef={:.2f}, G={:.2f} with Ef[{:.2f}, {:.2f}] and G[{:.2f}, {:.2f}]'.format(Ef_list[i]*1e-9, G_list[i]*1e-9, 1e-9 / (boundE[i][1]*I), 1e-9 /
(boundE[i][0]*I), 1e-9 / (k*boundG[i][1]*A), 1e-9 / (k*boundG[i][0]*A)),
    color='gray')
ax.legend()
plt.show()