import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import mpl_toolkits.mplot3d.axes3d as p3

fig= plt.figure()
ax = p3.Axes3D(fig)
ln, = plt.plot([], [], [], 'b.', alpha=0.5)
center, = plt.plot([], [], [], 'ro')
eliv, azim = 0, 0

sunpos = np.array([0.5, 0.5, 0.5])
marspos = np.array([0.2, 0.2, 0.2])
marsvel = np.array([0.02, -0.02, -0.02])
k = 0.001
kdust = 0.00000068
dustgenerate = True
dust=[]


class Dust :
    def __init__(sub, pos, vel):
        sub.pos = pos
        sub.vel = vel

    def update(sub):
        sunvec = sunpos - sub.pos
        marsvec = marspos - sub.pos
        sub.vel += (sunvec*1000*kdust/(np.linalg.norm(sunvec)**3) + marsvec*kdust/(np.linalg.norm(marsvec)**3))
        sub.pos = sub.pos + sub.vel

    def position(sub):
        return sub.pos


def update(frame):

    global sunpos, marspos, marsvel, eliv, azim, dust, dustgenerate     #|      fetch variables to local scope

    if dustgenerate :                                                   #|      generate new dust particles
        initvel = np.random.rand(1,3) - np.array([0.5, 0.5, 0.5])       #|
        initvel = 0.005*initvel/np.linalg.norm(initvel)               #|
        dust.append(Dust(marspos, marsvel+initvel))                     #|

   # eliv += 0.1                                                         #|      rotate axes
   # azim -= 0.15                                                        #|
   # ax.view_init(eliv, azim)                                        #|

    rvec = sunpos - marspos                                             #|      euler-cromer intg on mars-sun sys
    marsvel = marsvel + rvec*k/(np.linalg.norm(rvec)**3)                #|
    marspos = marspos + marsvel                                         #|

    arr = np.stack((sunpos, marspos))                                   #|      prepare array and euler-cromer on dust-sun-mars sys
    center.set_data_3d(arr.T)                                           #|
    for particle in dust :                                              #|
        particle.update()                                               #|
        arr = np.concatenate((arr, particle.position()))                #|
    
    ln.set_data_3d(arr.T)                                               #|      send data for animation
    return ln,center                                                    #|      


ani = FuncAnimation(fig, update, interval=1000/30)
plt.show()