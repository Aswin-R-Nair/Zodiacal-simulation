import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

r = np.array([[5,0, 10, np.pi/2]])
rdot = np.array([[0, 0.5, 0, 0.35]])

y = np.vstack([r, rdot]).astype(float)


def f(_, yi):
    y1 = np.array([])
    y2 = yi[1]

    for i in range(0, len(yi[0]), 2):
        rddot = (-1/(yi[0, i]**2)) * (yi[0, i:i+2] / yi[0, i])
        y1 = np.hstack((y1, rddot))

    return np.vstack((y2, y1))

def rk4_step(fyx, h, t_i, y_i):

    k1 = h * fyx(t_i, y_i)
    k2 = h * fyx((t_i + (h / 2)), (y_i + (k1 / 2)))  # RK4 algorithm
    k3 = h * fyx((t_i + (h / 2)), (y_i + (k2 / 2)))
    k4 = h * fyx((t_i + h), (y_i + k3))
    k = (k1 + (2 * k2) + (2 * k3) + k4) / 6

    y_i = y_i + k

    return y_i

def rk4_int(fun, yv, n, t_ini, t_end):
    dt = (t_end - t_ini) / n

    yi = yv
    ti = t_ini

    for i in range(n):
        yi = rk4_step(fun, dt, ti, yi)
        ti += dt

    return yi

print(rk4_int(f, y, 100, 0, 100))




