#!/usr/bin/env python
# coding: utf-8

# In[15]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def dAdt(A, t, beta, gamma, N, alpha, phi, kappa, mu, rho, omega, sigma):
    S = A[0]
    E = A[1]
    I = A[2]
    Q_R = A[3]
    Q_D = A[4]
    H_R = A[5]
    H_D = A[6]
    R = A[7]
    D = A[8]
    V = A[9]
    return[
        -beta * (2/np.pi*np.arctan((-t+t0)/k)+1) * S * I - alpha * S,
        beta * (2/np.pi*np.arctan((-t+t0)/k)+1) * S * I - gamma * E,
        gamma * E - sigma * I,
        sigma * (1-phi) * (1-kappa) * I - mu * Q_R,
        sigma * (1-phi) * kappa * I - rho * Q_D,
        sigma * phi * (1-kappa) * I - omega * H_R,
        sigma * phi * kappa * I - rho * H_D,
        mu * Q_R + omega * H_R,
        rho * Q_D + rho * H_D,
        alpha * S
    ]

times = np.arange(0, 34, 1)
gamma = np.log(2)/5
N = 1332898000
alpha = 0.0000003915059061
phi = 0.15
kappa = 0.01338343
mu = 1/5
rho = 1/15
omega = 1/10
sigma = np.log(2)/2
t0 = 10
k = 4
beta = 4e-10

i = 24479
e = 56694
a = 168235
qr = (a-i)*(1-phi)*(1-kappa)
qd = (a-i)*(1-phi)*kappa
hr = (a-i)*phi*(1-kappa)
hd = (a-i)*phi*kappa
r = 10394352
d = 154010
v = 3350260
s = N-a-r-d-v-e

S0, E0, I0, QR0, QD0, HR0, HD0, R0, D0, V0 = s, e, i, qr, qd, hr, hd, r, d, v
sol = odeint(dAdt, y0=[S0, E0, I0, QR0, QD0, HR0, HD0, R0, D0, V0], t=times, args=(beta, gamma, N, alpha, phi, kappa, mu, rho, omega, sigma))

S = sol.T[0]
E = sol.T[1]
I = sol.T[2]
QR = sol.T[3]
QD = sol.T[4]
HR = sol.T[5]
HD = sol.T[6]
A = sol.T[2]+sol.T[3]+sol.T[4]+sol.T[5]+sol.T[6]
R = sol.T[7]
D = sol.T[8]
V = sol.T[9]

plt.ylim(0, 35000000)
plt.plot(times, A, color='red', label='A modelled')
plt.plot(times, R, color='blue', label='R modelled')
plt.plot(times, D, color='green', label='D modelled')
plt.legend(loc='upper left')
plt.xlabel('t / Day')
plt.ylabel('Number of people')
plt.title('Number of active/recovered/deaceased cases')

s1 = S[33]
e1 = E[33]
i1 = I[33]
qr1 = QR[33]
qd1 = QD[33]
hr1 = HR[33]
hd1 = HD[33]
r1 = R[33]
d1 = D[33]
v1 = V[33]

times = np.arange(34, 100, 1)
beta = 9.3e-9

S1, E1, I1, QR1, QD1, HR1, HD1, R1, D1, V1 = s1, e1, i1, qr1, qd1, hr1, hd1, r1, d1, v1
sol = odeint(dAdt, y0=[S1, E1, I1, QR1, QD1, HR1, HD1, R1, D1, V1], t=times, args=(beta, gamma, N, alpha, phi, kappa, mu, rho, omega, sigma))

S1 = sol.T[0]
E1 = sol.T[1]
I1 = sol.T[2]
QR1 = sol.T[3]
QD1 = sol.T[4]
HR1 = sol.T[5]
HD1 = sol.T[6]
A1 = sol.T[2]+sol.T[3]+sol.T[4]+sol.T[5]+sol.T[6]
R1 = sol.T[7]
D1 = sol.T[8]
V1 = sol.T[9]

plt.ylim(0, 35000000)
plt.plot(times, A1, color='red')
plt.plot(times, R1, color='blue')
plt.plot(times, D1, color='green')

s2 = S1[65]
e2 = E1[65]
i2 = I1[65]
qr2 = QR1[65]
qd2 = QD1[65]
hr2 = HR1[65]
hd2 = HD1[65]
r2 = R1[65]
d2 = D1[65]
v2 = V1[65]

times = np.arange(99, 145, 1)
beta = 5.5e-9

S2, E2, I2, QR2, QD2, HR2, HD2, R2, D2, V2 = s2, e2, i2, qr2, qd2, hr2, hd2, r2, d2, v2
sol = odeint(dAdt, y0=[S2, E2, I2, QR2, QD2, HR2, HD2, R2, D2, V2], t=times, args=(beta, gamma, N, alpha, phi, kappa, mu, rho, omega, sigma))

S2 = sol.T[0]
E2 = sol.T[1]
I2 = sol.T[2]
QR2 = sol.T[3]
QD2 = sol.T[4]
HR2 = sol.T[5]
HD2 = sol.T[6]
A2 = sol.T[2]+sol.T[3]+sol.T[4]+sol.T[5]+sol.T[6]
R2 = sol.T[7]
D2 = sol.T[8]
V2 = sol.T[9]

plt.ylim(0, 35000000)
plt.plot(times, A2, color='red')
plt.plot(times, R2, color='blue')
plt.plot(times, D2, color='green')

s3 = S2[44]
e3 = E2[44]
i3 = I2[44]
qr3 = QR2[44]
qd3 = QD2[44]
hr3 = HR2[44]
hd3 = HD2[44]
r3 = R2[44]
d3 = D2[44]
v3 = V2[44]

times = np.arange(144, 171, 1)
beta = 1.73e-8

S3, E3, I3, QR3, QD3, HR3, HD3, R3, D3, V3 = s3, e3, i3, qr3, qd3, hr3, hd3, r3, d3, v3
sol = odeint(dAdt, y0=[S3, E3, I3, QR3, QD3, HR3, HD3, R3, D3, V3], t=times, args=(beta, gamma, N, alpha, phi, kappa, mu, rho, omega, sigma))

S3 = sol.T[0]
E3 = sol.T[1]
I3 = sol.T[2]
QR3 = sol.T[3]
QD3 = sol.T[4]
HR3 = sol.T[5]
HD3 = sol.T[6]
H3 = sol.T[5]+sol.T[6]
A3 = sol.T[2]+sol.T[3]+sol.T[4]+sol.T[5]+sol.T[6]
R3 = sol.T[7]
D3 = sol.T[8]
V3 = sol.T[9]

plt.ylim(0, 35000000)
plt.plot(times, A3, color='red')
plt.plot(times, R3, color='blue')
plt.plot(times, D3, color='green')

import pandas as pd
import matplotlib.pyplot as plt
df = pd.read_excel('/Users/caoliwen/Desktop/Optimisation research/India covid cases copy.xlsx', sheet_name='India')
plt.twinx()
plt.ylim(0, 35000000)
plt.plot(df['Day'], df['Active Cases'], color='pink', label='A actual')
plt.plot(df['Day'], df['Cured/Discharged'], color='skyblue', label='R actual')
plt.plot(df['Day'], df['Death'], color='limegreen', label='D actual')
plt.legend(loc='upper center')
plt.axis('off')

plt.tight_layout()
plt.show()


# In[ ]:




