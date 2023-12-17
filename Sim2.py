import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter


# Initializing bodies
class Body:
    def __init__(self, r, rdot):
        self.pos = r
        self.vel = rdot

# Giving inital values for each body
sun = Body(np.array([0, 0, 0]), 0)
mars = Body(np.array([1.639, 0, 0]), np.array([0, -4.72, 0]))     # Body(postiion(xyz), velocity(xyz))
jupiter = Body(np.array([5.367, 0, 0]), np.array([0, -2.66, 0]))


dt = 0.1                 # Timestep

k = 39.478                   # proportionality constant for sun system
kd = 0.000012672           # proportionality constant for dust-mars system
kdj = 0.0377            # proportionality constant for dust-jupiter system
kb = 6327.8


# Initializing dust parameters
generatedust = 1              # set True if you want dust
n = 20                        # number of dust particles


ml = []
sdl = []
ni = 0

# Main simulation fuction
def av_distance_calc(n_iter, rd, on, lr):

    global ni

    kr = kb * (rd**2)  # Defining radiation pressure constant with given radius


    rl_external = [] # initiaing final position list
    dust = []


    for i in range(n_iter):

        rl = sun.pos # initiaing position matrix with sun position


        if i <= n and generatedust:
            # dustv = np.random.uniform(-0.5, 0.5, (1, 3))  # Initializing random velocity(min, max, (matrix))
            dustv = np.array([0.03382911, 0.48789013, -0.23769815])
            dustv_normalised = 1.05 * dustv/np.linalg.norm(dustv)  # Normalizing velocity (number is total velocity)
            dust.append(Body(mars.pos, mars.vel + dustv_normalised))  # Initializing dust body with (position, velocity)

        # Updating mars position using Euler-cromer algorithm
        rv = mars.pos - sun.pos  # Mars-sun vector
        mars.vel = mars.vel - k * dt * rv / (np.linalg.norm(rv) ** 3)  # v_n+1 = v_n - a*dt (a is g-force)
        mars.pos = mars.pos + (mars.vel * dt)  # x_n+1 = x_n + v_n+1*dt

        # Updating jupiter position using Euler-cromer
        rvj = jupiter.pos - sun.pos  # Jupiter-sun vector
        jupiter.vel = jupiter.vel - k * dt * rvj / (np.linalg.norm(rvj) ** 3)  # v_n+1 = v_n - a*dt (a is g-force)
        jupiter.pos = jupiter.pos + (jupiter.vel * dt)  # x_n+1 = x_n + v_n+1*dt


        # Updating dust position using Euler-cromer on each particle
        for o, p in enumerate(dust):
            d_mars_r = p.pos - mars.pos
            d_sun_r = p.pos - sun.pos  # Obtaining dust-body vectors
            d_jup_r = p.pos - jupiter.pos

            if np.linalg.norm(d_sun_r) > 6:  # Killing off dust that gets too far
                dust.pop(o)
                continue

            p.vel = p.vel - (((k - kr) * dt * d_sun_r / (np.linalg.norm(d_sun_r)) ** 3) + (
                    kd * dt * d_mars_r / (np.linalg.norm(d_mars_r) ** 3)) + (  # velocity update of dust
                                     kdj * dt * d_jup_r / (
                                         np.linalg.norm(d_jup_r) ** 3)))  # (Forces from sun, mars, jupiter)

            p.pos = p.pos + (p.vel * dt)  # position update
            rl = np.vstack((rl, p.pos))  # Adding updated dust positions

        # print(f"{on}/{lr}, {i}/{n_iter}")  # Prints out progress of simulation
        if i == n_iter - 1:
            rl_external = rl   # Stores the final position of particles at the end of the simulation

    print(rl_external)

    d = np.array([np.linalg.norm(rl_external[0] - rl_external[i]) for i in range(1, len(rl_external))])
    # print(d)
    ml.append(np.mean(d))    # Calculating mean orbital distance
    sdl.append(np.std(d))    # Calcualting standard deviation from mean



frames = 700

# rds = np.linspace(10, 100, 10)   # Creating range of radii (start, end, step)
rds = [1, 12, 23, 34]

for on, r in enumerate(rds):
    av_distance_calc(frames, r, on, len(rds))  # Main function looped through the radii



# Setting up plots for mean and standard deviation
fig, (ax1, ax2) = plt.subplots(2)

ax1.plot(rds, ml)
ax1.set_ylabel("Mean orbital radius (AU)")
# ax1.set_xscale("log", base=2)
# ax1.xaxis.set_major_formatter(ScalarFormatter())

ax2.plot(rds, sdl)
ax2.set_ylabel("Standard deviation (AU)")
ax2.set_xlabel(r"Radius of dust particle ($\mu$m)")
# ax2.set_xscale("log")
# ax2.xaxis.set_major_formatter(ScalarFormatter())

plt.show()








