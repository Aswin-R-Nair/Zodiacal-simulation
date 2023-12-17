This was part of the summer internship in 2022 which aimed to study if the dust particles that produce Zodiacal light on earth, originates from mars.
The code simulates dust particles being ejected from mars at some random velocity (besides the escape velocity of mars) and moving in a heliocentric orbit. 
The dust motion is governed by the gravity of the sun, mars and jupiter and radiation pressure from the sun.
There is also a simulation that investigates the motion of dust as ejected by a comet, Halleys comet in this case. 

The code has several versions. The pure python code in the Main_slow directory. As the name suggests the code is not optimised and accrue too much memory overhead causing it to be slow. This only has the dust originating from mars.

The more optimised versions run using Cython, C wrapped in python.
The mars_sim directory contains simulation for dust originating from mars and the root directory contains simulation for dust originating from Halleys comet, which run using Cython.

In these directories, the "MainSim", "cMS2" or "cMainSim" contains the visual simulation and "Sim" ,"cS2" or "cDataSim" contains code that returns a plot of orbital parameters.
