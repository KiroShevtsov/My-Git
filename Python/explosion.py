import numpy as np
import matplotlib.pyplot as plt

#parameters
Te = 368.0
R = 0.05
eta0 = 0.1
E = 1.6263e5
R0 = 8.3144126
Z = 4 * 1e17
rho = 1670
Q = 1.1972e6
alpha = 5.0
c = 1682.0
lambda_l = 0.29
N = 100
t_max = 1e6
t_start = 0.0
max_iter = 300000

#array for time and radiuses
r = np.linspace(0, R, N)
#step for r
h = R / (N - 1)
dt = 0.5

#data for plot
save_every = 100
rows = int(max_iter / save_every)

r_layer = np.zeros(N)

terms = np.zeros((rows, N))
row_idx = 0
times = []

#max temp
explosion_idx = 0
flag_fixed = False

#sweep method for solving tridiagonal system
def thompson(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray):
    n = len(d)
    alpha = np.zeros(n)
    beta = np.zeros(n)
    alpha[0] = - (c[0] / b[0])
    beta[0] = d[0] / b[0]

    #main loop
    for iter in range(1, n):
        alpha[iter] = - (c[iter]) / (a[iter] * alpha[iter - 1] + b[iter])
        beta[iter] = (d[iter] - beta[iter - 1] * a[iter]) / (a[iter] * alpha[iter - 1] + b[iter])
    
    #solving relative unknown x
    x = np.zeros(n)
    x[n - 1]= beta[n - 1]
    start = n - 2
    for index in range(start, -1, -1):
        x[index] = alpha[index] * x[index + 1] + beta[index]
    return x
#rhs for equation
def kinematic(T, eta):
    magnitude = (-E / (R0 * T))
    return (eta + eta0) * (1 - eta) * Z * np.exp(magnitude)

#runge for kinematic differential equation
def runge4(T_current, eta_current, dt):
    k1 = dt * kinematic(T_current, eta_current)
    k2 = dt * kinematic(T_current, eta_current + 0.5 * k1)
    k3 = dt * kinematic(T_current, eta_current + 0.5 * k2)
    k4 = dt * kinematic(T_current, eta_current + k3)
    res = eta_current + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    #eta < 1
    return np.clip(res, 0.0, 1.0)

Temperatures = np.full(N, Te)
eta_current = np.zeros(N)

print("Начинается расчет взрыва...")
for ext in range(1, max_iter + 1):

    current_time = ext * dt
    times.append(float(dt))
    #source
    eta_next = runge4(Temperatures, eta_current, dt)
    source = rho * Q * (eta_next - eta_current) / dt
    #solve equation of heat
    a = np.zeros(N)
    b = np.zeros(N)
    c_mtx = np.zeros(N)
    d = np.zeros(N)
    a[0] = 0.0
    b[0] = rho * c + (6 * lambda_l * dt) / h**2
    c_mtx[0]= (-6 * lambda_l * dt) / h**2
    d[0] = (rho * c) * Temperatures[0] + source[0] * dt

    for i in range(1, N - 1):
        r_right = r[i] + 0.5 * h
        r_left = r[i] - 0.5 * h
        current_r = r[i]

        #added radiuses for plor
        r_layer[i] = current_r

        #diagonals
        a[i] = -lambda_l * dt / (current_r**2 * h**2) * r_left**2
        b[i] = rho * c + lambda_l * dt / (current_r**2 * h**2) * (r_right**2 + r_left**2)
        c_mtx[i] = -lambda_l * dt / (current_r**2 * h**2) * r_right**2
        d[i] = rho * c * Temperatures[i] + source[i] * dt
    #last coefficients
    a[N - 1] = -lambda_l / h
    b[N - 1] = alpha + lambda_l / h
    c_mtx[N - 1] = 0.0
    d[N - 1] = alpha * Te

    #for plot
    r_layer[N - 1] = R

    #get new array temperatures
    Temperatures_new = thompson(a, b, c_mtx, d)
    Temperatures = Temperatures_new.copy()
    eta_current = eta_next.copy()

    #for plot
    if ext % save_every == 0:
        terms[row_idx, :] = Temperatures

        progress = (ext / max_iter) * 100
        # \r возвращает курсор в начало строки, end='' не дает перенести строку
        center_T = Temperatures[i]
        print(f"\rПрогресс: {progress:5.1f}%, Шаг: {ext}/{max_iter}, t: {current_time:8.1f} с, T в центре: {center_T:6.1f} K", end='', flush=True)
        if center_T > 450.0 and not flag_fixed:
            explosion_idx = row_idx
            flag_fixed = True
        row_idx += 1
print("Рассчитано")
#surface
X, Y = np.meshgrid(np.linspace(0, R, N), np.linspace(0, max_iter * dt, len(terms)))

ax = plt.figure().add_subplot(projection='3d')
surf = ax.plot_surface(X, Y, terms, cmap='plasma')
plt.colorbar(surf, pad=0.1)

ax.set_xlabel('r, м')
ax.set_ylabel('t, с')
plt.show()

#T(t)
time_axis = np.linspace(0, max_iter * dt, len(terms))

plt.figure(figsize=(10, 6))
plt.plot(time_axis, terms[:, 0], color='blue')

plt.xlabel('Время (t), с')
plt.ylabel('Температура в центре T(t, 0), К')
plt.title('Динамика температуры в центре шара')
plt.grid(True, which='both')
plt.show()

#T_max(r)
plt.figure(figsize=(10, 6))
plt.plot(r, terms[explosion_idx, :], color='blue', label='начало взрыва')

plt.xlabel('Радиус r, м')
plt.ylabel('Температура T, К')

explosion_time = explosion_idx * save_every * dt

plt.title(f'Профиль температуры T(r) в момент взрыва (t = {explosion_time:.1f} с)')
plt.grid(True, which='both', alpha=0.5)
plt.legend()

plt.show()