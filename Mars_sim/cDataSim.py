import numpy as np
import Sim2_cmod
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# Initializing bodies
class Body:
    def __init__(self, r, rdot):
        self.pos = Sim2_cmod.make_vec(r["x"], r["y"], r["z"])
        self.vel = Sim2_cmod.make_vec(rdot["x"], rdot["y"], rdot["z"])
    
    
# Giving inital values for each body
sun = Body({'x': 0.0, 'y': 0.0, 'z': 0.0}, {'x': 0.0, 'y': 0.0, 'z': 0.0})
mars = Body({'x': 1.639, 'y': 0.0, 'z': 0.0}, {'x': 0.0, 'y': -4.72, 'z': 0.0})    # Body(postiion(xyz), velocity(xyz))
jupiter = Body({'x': 5.367, 'y': 0.0, 'z': 0.0}, {'x': 0.0, 'y': -2.66, 'z': 0.0})


dt = 0.01                 # Timestep

k = 39.478                   # proportionality constant for sun system
kd = 0.000012672           # proportionality constant for dust-mars system
kdj = 0.0377            # proportionality constant for dust-jupiter system

rho = 3580.986                          # denisty of particle
q_pr = 2                                # radiation coefficient (qpr = 2 is assuming perfect reflectivity)
betaR = 5.7 * (10**(-4)) * (q_pr/rho)   # beta*R (where, beta = F_r/F_g, R = radius of particle)



# Initializing dust parameters
generatedust = 1              # set True if you want dust
n = 800                  # number of dust particles 


ml = []
sdl = []


# Main simulation fuction
def av_distance_calc(n_iter, rd):
    sun = Body({'x': 0.0, 'y': 0.0, 'z': 0.0}, {'x': 0.0, 'y': 0.0, 'z': 0.0})
    mars = Body({'x': 1.639, 'y': 0.0, 'z': 0.0}, {'x': 0.0, 'y': -4.72, 'z': 0.0})    # Body(postiion(xyz), velocity(xyz))
    jupiter = Body({'x': 5.367, 'y': 0.0, 'z': 0.0}, {'x': 0.0, 'y': -2.66, 'z': 0.0})


    R = rd * (10**(-6))                    # radius of particle
    beta = betaR / R

    
    rl_external = [] # initiaing final position list
    dust = []
    
    

    for i in range(n_iter):

        rl = []

        if i <= n and generatedust:
            dustv = np.random.uniform(-0.5, 0.5, (1, 3))  # Initializing random velocity(min, max, (matrix))
            dustv = Sim2_cmod.make_vec(dustv[0,0], dustv[0,1], dustv[0,2])
            dustv_normalized = Sim2_cmod.smul_v((1.05/Sim2_cmod.mag(dustv)), dustv)
            dust.append(Body(mars.pos, Sim2_cmod.add_v(mars.vel, dustv_normalized)))

        rv = Sim2_cmod.sub_v(mars.pos, sun.pos)
        mars.vel = Sim2_cmod.sub_v(mars.vel, Sim2_cmod.smul_v(k*dt/(Sim2_cmod.mag(rv)**3), rv))
        mars.pos = Sim2_cmod.add_v(mars.pos, Sim2_cmod.smul_v(dt, mars.vel))

        rvj = Sim2_cmod.sub_v(jupiter.pos, sun.pos)
        jupiter.vel = Sim2_cmod.sub_v(jupiter.vel, Sim2_cmod.smul_v(k*dt/(Sim2_cmod.mag(rvj)**3), rvj))
        jupiter.pos = Sim2_cmod.add_v(jupiter.pos, Sim2_cmod.smul_v(dt, jupiter.vel))
        
        rl.extend((sun.pos, mars.pos))
        
        for o,p in enumerate(dust):
            d_sun_r = Sim2_cmod.sub_v(p.pos, sun.pos)
            d_mars_r = Sim2_cmod.sub_v(p.pos, mars.pos)
            d_jupiter_r = Sim2_cmod.sub_v(p.pos, jupiter.pos)

            a = Sim2_cmod.add_v(Sim2_cmod.smul_v(kd*dt/(Sim2_cmod.mag(d_mars_r)**3), d_mars_r), Sim2_cmod.smul_v(kdj*dt/(Sim2_cmod.mag(d_jupiter_r)**3), d_jupiter_r))
            b = Sim2_cmod.add_v(a, Sim2_cmod.smul_v((k*(1-beta))*dt/(Sim2_cmod.mag(d_sun_r)**3), d_sun_r))
            p.vel = Sim2_cmod.sub_v(p.vel, b)
            p.pos = Sim2_cmod.add_v(p.pos, Sim2_cmod.smul_v(dt, p.vel))

            rl.append(p.pos)

            

        if i == n_iter - 1:
            rl_external = rl
        
    
    
    d = np.array([Sim2_cmod.mag(Sim2_cmod.sub_v(rl_external[0], rl_external[int(i)])) for i in range(2, len(rl_external))])
    ml.append(np.mean(d))    # Calculating mean orbital distance
    sdl.append(np.std(d)) 
    



frames = 1000


rds = np.linspace(0.1, 100, 500)   # Creating range of radii (start, end, step)

print("Running...")
for on, r in enumerate(rds):
    av_distance_calc(frames, r)  # Main function looped through the radii



# Setting up plots for mean and standard deviation
fig, (ax1, ax2) = plt.subplots(2)

showlog = 1

ax1.plot(rds, ml)
ax1.set_ylabel("Mean orbital radius (AU)")
ax1.grid()


ax2.plot(rds, sdl)
ax2.set_ylabel("Standard deviation (AU)")
ax2.set_xlabel(r"Radius of dust particle ($\mu$m)")
ax2.grid()


if showlog:
    ax1.set_xscale("log")
    ax1.grid(which='minor', color='grey', axis="y", linestyle='--', alpha=0.4)
    ax1.xaxis.set_major_formatter(ScalarFormatter())
    ax1.set_yscale("log")
    ax1.yaxis.set_major_formatter(ScalarFormatter())
    ax2.set_xscale("log")
    # ax2.xaxis.set_major_formatter(ScalarFormatter())


plt.show()
        
        









 


