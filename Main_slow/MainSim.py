import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import rcParams
rcParams['animation.ffmpeg_path'] = r'C:\Program Files\ImageMagick-7.0.11-Q16-HDRI\ffmpeg.exe'

# Initializing bodies
class Body:
    def __init__(self, r, rdot):
        self.pos = r
        self.vel = rdot

# Giving inital values for each body
sun = Body(np.array([0, 0, 0]), 0)
mars = Body(np.array([1.639, 0, 0]), np.array([0, -4.72, 0]))     # Body(postiion(xyz), velocity(xyz))
jupiter = Body(np.array([5.367, 0, 0]), np.array([0, -2.66, 0]))
dust = []

dt = 0.01                 # Timestep

k = 39.478                   # proportionality constant for sun system
kd = 0.000012672           # proportionality constant for dust-mars system
kdj = 0.0377            # proportionality constant for dust-jupiter system

rho = 3580.986                          # denisty of particle
q_pr = 2                                # radiation coefficient (qpr = 2 is assuming perfect reflectivity)
betaR = 5.7 * (10**(-4)) * (q_pr/rho)   # beta*R (where, beta = F_r/F_g, R = radius of particle)

R = 50 * (10**(-6))                    # radius of particle
beta = betaR / R

# Setting up plotting module
fig = plt.figure()
axes = plt.axes(projection="3d")

minlim = -2                         # Initial zoom
maxlim = 2

b1, = axes.plot3D([], [], [], "*m")  # Sun
b2, = axes.plot3D([], [], [], "ro")  # Mars
b3, = axes.plot3D([], [], [], "oy") # Jupiter
b4, = axes.plot3D([], [], [], ".b")    # Dust

tl = axes.text2D(0.5, 0.9, s="", transform=plt.gcf().transFigure)


axes.set_xlim(minlim, maxlim)
axes.set_ylim(minlim, maxlim)
axes.set_zlim(minlim, maxlim)

elev, azim = 90, 0

# Initializing dust parameters
generatedust = 1              # set True if you want dust
n = 1000                         # number of dust particles
randomvelocity = 1             # set True if dust should have random values, False if you want custom value


# Main animation function
def animate(i):

    global elev, azim

    if i <= n and generatedust:
        if randomvelocity:
            dustv = np.random.uniform(-1, 1, (1, 3))          # Initializing random velocity(min, max, (matrix))
        else:
            dustv = np.array([0.03382911, 0.48789013, -0.23769815])                       # Initializing custom velocity(xyz)

        dustv_normalised = 1.05*dustv/np.linalg.norm(dustv)             # Normalizing velocity (number is total velocity)
        dust.append(Body(mars.pos, mars.vel + dustv_normalised))   # Initializing dust body with (position, velocity)


    # Updating mars position using Euler-cromer algorithm
    rv = mars.pos - sun.pos                                         # Mars-sun vector
    mars.vel = mars.vel - k * dt * rv / (np.linalg.norm(rv) ** 3)   # v_n+1 = v_n - a*dt (a is g-force)
    mars.pos = mars.pos + (mars.vel * dt)                           # x_n+1 = x_n + v_n+1*dt

    # Updating jupiter position using Euler-cromer
    rvj = jupiter.pos - sun.pos                                             #Jupiter-sun vector
    jupiter.vel = jupiter.vel - k * dt * rvj / (np.linalg.norm(rvj) ** 3)   # v_n+1 = v_n - a*dt (a is g-force)
    jupiter.pos = jupiter.pos + (jupiter.vel * dt)                          # x_n+1 = x_n + v_n+1*dt


    # Adding updated positions of the bodies
    rl = np.vstack((sun.pos, mars.pos, jupiter.pos))

    # Updating dust position using Euler-cromer on each particle
    for o, p in enumerate(dust):
        d_mars_r = p.pos - mars.pos
        d_sun_r = p.pos - sun.pos              # Obtaining dust-body vectors
        d_jup_r = p.pos - jupiter.pos

        if np.linalg.norm(d_sun_r) > 6:        # Killing off dust that gets too far
            dust.pop(o)
            continue

        p.vel = p.vel - ((k * (1-beta) * dt * d_sun_r / (np.linalg.norm(d_sun_r)) ** 3) + (
                    kd * dt * d_mars_r / (np.linalg.norm(d_mars_r) ** 3)) +(          # velocity update of dust
                kdj * dt * d_jup_r / (np.linalg.norm(d_jup_r) ** 3)))                 # (Forces from sun, mars, jupiter)


        p.pos = p.pos + (p.vel * dt)         # position update
        rl = np.vstack((rl, p.pos))          # Adding updated dust positions

    b1.set_data_3d(rl[0])
    b2.set_data_3d(rl[1])
    b3.set_data_3d(rl[2])                    # Sending positions to be plotted
    b4.set_data_3d(rl[3:].T)
    tl.set_text(f"{(i*0.01):.02f} years")


    elev -= 0.09                                                         # rotate axes
    azim += 0.18
    axes.view_init(elev, azim)

    return b1, b2, b3, b4



funky = FuncAnimation(fig, animate, interval=1000/60, repeat=False, frames=2000)


# funky.save("Simulation2.mp4", writer="ffmpeg", fps=60, dpi=300)                   # Saving the animation
# print("Done")


plt.tight_layout()                                                                # Plotting the animation
plt.show()








