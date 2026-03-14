import matplotlib.pyplot as plt
import numpy as np
M0 = 1.4

#теоретическое число Маха для площадей
def G(gamma, A, A_crit, M_last):
    p = (gamma + 1) / (gamma - 1)
    degree = 2 * (gamma - 1) / (gamma + 1)
    k = p * ((A) * (M_last) / (A_crit))**degree - 2 / (gamma - 1)
    return np.sqrt(k)

def S(diameter):
    return 0.25 * np.pi * diameter**2

M = [M0]
gamma = 1.4
diameter_out = 2.25
diameter_crit = 1.645
A = S(diameter_out)
A_crit = S(diameter_crit)
count = 100
for i in range(1, count, 1):
    M.append(G(gamma, A, A_crit, M[-1]))
M_res = np.array(M)
print(M_res[99])

#для числа Маха для первого сопла через давления -> p1 / p2
def H(gamma, p1, p2, M_last):
    A = (gamma - 1) / (2 * gamma)
    B = (gamma + 1) / (2 * gamma)
    degree1 = - gamma + 1
    degree2 = -gamma
    D = 2 / (gamma + 1)
    F = (gamma - 1) / 2
    return np.sqrt(A + B * ((p1 / p2)**degree1) * (D * (1 / (M_last**2) + F))**degree2)
M = [M0]
gamma = 1.4
#торможение
#атмосферное
delta_p = 1
p1 = 276 * 0.015 + delta_p
#ресивер
p2 = 7.8 + delta_p
count = 100
for i in range(count):
    M.append(H(gamma, p1, p2, M[-1]))
M_res = np.array(M)
print(M_res[99])


