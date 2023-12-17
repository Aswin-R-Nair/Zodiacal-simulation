import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from trialObjects import Body


# Defining the r vector modulus function
def rmod(rv):
    return np.sqrt(rv[0]**2 + rv[1]**2 + rv[2]**2)


# Defining the (aceeleration) diff equation f vector
def f(_, y):
    y1 = y[0]
    y2 = y[1]

    return np.array([y2, (-1/rmod(y1)**3)*y1])  # Vector form of gravity


# Defining Runge-Kutta 4th order function for a single timestep
def rk4_step(fyx, h, t_i, y_i):

    k1 = h * fyx(t_i, y_i)
    k2 = h * fyx((t_i + (h / 2)), (y_i + (k1 / 2)))  # RK4 algorithm
    k3 = h * fyx((t_i + (h / 2)), (y_i + (k2 / 2)))
    k4 = h * fyx((t_i + h), (y_i + k3))
    k = (k1 + (2 * k2) + (2 * k3) + k4) / 6

    y_i = y_i + k

    return y_i


# Defining RK4 over the whole integration range
def rk4_int(func, bs, n, t_init, t_end):

    dt = (t_end - t_init) / n   # Defining timestep
    tp = [t_init]  # Initiating timesteps list

    t_i = t_init
    for i in range(n):
        for b in bs:
            b.y = rk4_step(func, dt, t_i, b.y)   # In each timestep, it loops through each body and updates y values
            b.r_points.append(b.y.tolist()[0])  # Store the updated y values (the r coordinates)
        t_i += dt                # Update timestep
        tp.append(t_i)

    return tp



# Initiating bodies (r_naught(xyz), rdot_naught(xyz))
m1 = Body([3, 3, 3], [0.3, -0.3, -0.3])
m2 = Body([10, 0, 0], [0, 0.3, 0])
m3 = Body([6, 0, 0], [0, 0.35, 0])

bodies = [m1]

# Performing rk4 with some initial values (function, bodies, n_steps, t_start, t_end)
t = rk4_int(f, bodies, 200, 0, 2000)

# Setting up plot figure and axes
fig = plt.figure()
axes = plt.axes(projection="3d")


tl = 20  # Length of tail for animated bodies


# Defining function for drawing each frame
def animate(i):
    axes.clear()

    x_coordinate = [b.r_points[i][0] for b in bodies]
    y_coordinate = [b.r_points[i][1] for b in bodies]  # Fetch the scatter points for all bodies for this frame
    z_coordinate = [b.r_points[i][2] for b in bodies]

    axes.scatter3D(x_coordinate, y_coordinate, z_coordinate, s=10, color="red")  # Plot the bodies
    axes.scatter3D(0,0,0, s=50, color="yellow")   # Stationary non interacting origin

    axes.set_xlim(-8, 6)
    axes.set_ylim(-6, 6)
    axes.set_zlim(-8, 8)


# Animating the simulation
anime = FuncAnimation(fig, animate, interval=2000/60, repeat=False)

plt.tight_layout()
plt.show()


