import numpy as np


class Body:
    def __init__(self, r_0, rdot_0):
        # self.mass = mass
        self.r = np.array(r_0)
        self.rdot = np.array(rdot_0)
        self.r_points = [self.r.tolist()]
        self.y = np.array([self.r, self.rdot])

for i in range(10):
    if not i%2:
        continue
    print(i)

