import Sim2_cmod
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import rcParams
import math
rcParams['animation.ffmpeg_path'] = r'C:\Program Files\ImageMagick-7.0.11-Q16-HDRI\ffmpeg.exe'

# Initializing bodies
class Body:
    def __init__(self, r, rdot):
        self.pos = Sim2_cmod.make_vec(r["x"], r["y"], r["z"])
        self.vel = Sim2_cmod.make_vec(rdot["x"], rdot["y"], rdot["z"])


def dict_mapper(arr):
    a = np.array(["x", "y", "z"])
    li = []
    for i in arr:
        li.append(np.vectorize(i.get)(a))
    li = np.vstack(li)
    return li


# Giving inital values for each body
sun = Body({'x': 0.0, 'y': 0.0, 'z': 0.0}, {'x': 0.0, 'y': 0.0, 'z': 0.0})
mars = Body({'x': 1.639, 'y': 0.0, 'z': 0.0}, {'x': 0.0, 'y': -4.72, 'z': 0.0})    # Body(postiion(xyz), velocity(xyz))
jupiter = Body({'x': 5.367, 'y': 0.0, 'z': 0.0}, {'x': 0.0, 'y': -2.66, 'z': 0.0})
halley = Body({'x': 33.5586, 'y': 0.0, 'z': -10.9276}, {'x': 0.0, 'y': 0.1916, 'z': 0.0})

dust = []

dt = 0.01                 # Timestep

k = 39.478                   # proportionality constant for sun system (G*M_sun)
kd = 0.000012672           # proportionality constant for dust-mars system (G*M_mars)
kdj = 0.0377            # proportionality constant for dust-jupiter system  (G*M_jupiter)

rho = 3580.986                          # denisty of particle
q_pr = 2                                # radiation coefficient (qpr = 2 is assuming perfect reflectivity)
betaR = 5.7 * (10**(-4)) * (q_pr/rho)   # beta*R (where, beta = F_r/F_g, R = radius of particle)

R = 5 * (10**(-6))                    # radius of particle
beta = betaR / R


# Setting up plotting module
fig = plt.figure()
axes = plt.axes(projection="3d")

minlim = -35                         # Initial zoom
maxlim = 35

b1, = axes.plot3D([], [], [], "*m")  # Sun
b2, = axes.plot3D([], [], [], "ro")  # Mars
b3, = axes.plot3D([], [], [], "oy") # Jupiter
b4, = axes.plot3D([], [], [], ".b")    # Dust
b5, = axes.plot3D([], [], [], "og")   # Halley
b6, = axes.plot3D([], [], [])        #halley orbit

tl = axes.text2D(0.5, 0.9, s="", transform=plt.gcf().transFigure)


axes.set_xlim(minlim, maxlim)
axes.set_ylim(minlim, maxlim)
axes.set_zlim(minlim, maxlim)

elev, azim = 90, 0

# Initializing dust parameters
generatedust = 1               # set True if you want dust
n = 800                         # number of dust particles
randomvelocity = 1             # set True if dust should have random values, False if you want custom value

orbitline = []

# Main animation function
def animate(i):
    rl = []

    global elev, azim

    r_sun = Sim2_cmod.mag(Sim2_cmod.sub_v(sun.pos, halley.pos))
    
    if len(dust) <= n and generatedust and r_sun<=5:
        
        dustv = np.random.uniform(-0.5, 0.5, (1, 3))
        dustv = Sim2_cmod.make_vec(dustv[0,0], dustv[0,1], dustv[0,2])
        vmax = 1.05
        dustv_normalized = Sim2_cmod.smul_v((vmax/Sim2_cmod.mag(dustv)), dustv)
        dust.append(Body(halley.pos, Sim2_cmod.add_v(halley.vel, dustv_normalized)))          # Initializing random velocity(min, max, (matrix))
        
                
        
        


    # Updating mars position using Euler-cromer algorithm
    rv = Sim2_cmod.sub_v(mars.pos, sun.pos)
    mars.vel = Sim2_cmod.sub_v(mars.vel, Sim2_cmod.smul_v(k*dt/(Sim2_cmod.mag(rv)**3), rv))
    mars.pos = Sim2_cmod.add_v(mars.pos, Sim2_cmod.smul_v(dt, mars.vel))


    # Updating jupiter position using Euler-cromer
    rvj = Sim2_cmod.sub_v(jupiter.pos, sun.pos)
    jupiter.vel = Sim2_cmod.sub_v(jupiter.vel, Sim2_cmod.smul_v(k*dt/(Sim2_cmod.mag(rvj)**3), rvj))
    jupiter.pos = Sim2_cmod.add_v(jupiter.pos, Sim2_cmod.smul_v(dt, jupiter.vel))


    # Updating halley position using Euler-cromer
    rvh = Sim2_cmod.sub_v(halley.pos, sun.pos)
    halley.vel = Sim2_cmod.sub_v(halley.vel, Sim2_cmod.smul_v(k*dt/(Sim2_cmod.mag(rvh)**3), rvh))
    halley.pos = Sim2_cmod.add_v(halley.pos, Sim2_cmod.smul_v(dt, halley.vel))


    # Adding updated positions of the bodies
    rl.extend((sun.pos, mars.pos, jupiter.pos, halley.pos))
    if Sim2_cmod.mag(Sim2_cmod.sub_v(sun.pos, halley.pos)) <= 2:
        if not i%10:
            orbitline.append(halley.pos)
    else:
        if not i%100:
            orbitline.append(halley.pos)

    # Updating dust position using Euler-cromer on each particle
    for o,p in enumerate(dust):
        d_sun_r = Sim2_cmod.sub_v(p.pos, sun.pos)
        d_mars_r = Sim2_cmod.sub_v(p.pos, mars.pos)
        d_jupiter_r = Sim2_cmod.sub_v(p.pos, jupiter.pos)

        a = Sim2_cmod.add_v(Sim2_cmod.smul_v(kd*dt/(Sim2_cmod.mag(d_mars_r)**3), d_mars_r), Sim2_cmod.smul_v(kdj*dt/(Sim2_cmod.mag(d_jupiter_r)**3), d_jupiter_r))
        b = Sim2_cmod.add_v(a, Sim2_cmod.smul_v((k*(1-beta))*dt/(Sim2_cmod.mag(d_sun_r)**3), d_sun_r))
        p.vel = Sim2_cmod.sub_v(p.vel, b)
        p.pos = Sim2_cmod.add_v(p.pos, Sim2_cmod.smul_v(dt, p.vel))

        rl.append(p.pos)       # Adding updated dust positions
    
    rl = dict_mapper(rl)
    oll = dict_mapper(orbitline)
    # print(oll)

    b1.set_data_3d(rl[0])
    b2.set_data_3d(rl[1])
    b3.set_data_3d(rl[2])
    b5.set_data_3d(rl[3])                    # Sending positions to be plotted
    b4.set_data_3d(rl[4:].T)
    b6.set_data_3d(oll.T)
    tl.set_text(f"{(i*dt):.02f} years")

    # elev -= 0.09                                                         # rotate axes
    # azim += 0.18
    # axes.view_init(elev, azim)

    return b1, b2, b3, b4, b5, b6



funky = FuncAnimation(fig, animate, interval=10, repeat=False, frames=12000)


# funky.save("Simulation2.mp4", writer="ffmpeg", fps=60, dpi=300)                   # Saving the animation
# print("Done")


plt.tight_layout()                                                                # Plotting the animation
plt.show()















