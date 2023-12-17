import numpy as np
cimport numpy as np

cdef struct Vec3:
    double x, y, z


cdef Vec3 vec3(double x, double y, double z):
    cdef Vec3 v
    v.x = x
    v.y = y
    v.z = z
    return v


cdef Vec3 add(Vec3 u, Vec3 v):
    return vec3(u.x + v.x, u.y + v.y, u.z + v.z)

cdef Vec3 sub(Vec3 u, Vec3 v):
    return vec3(u.x-v.x, u.y-v.y, u.z-v.z)

cdef Vec3 smul(double n, Vec3 v):
    return vec3(v.x*n, v.y*n, v.z*n)

cdef Vec3 sdiv(double n, Vec3 v):
    return vec3(v.x/n, v.y/n, v.z/n)


from libc.math cimport sqrt
cdef magnitude(Vec3 v):
    cdef double n
    n = sqrt(v.x*v.x + v.y*v.y + v.z*v.z)
    return n



cdef append(Vec3 v):
    cdef Vec3 vec_list[810]
    vec_list[0] = v
    return vec_list

cdef append_sp(Vec3 v, Vec3 u, Vec3 w):
    cdef Vec3 v_list[810]
    v_list[0] = v
    v_list[1] = u
    v_list[2] = w
    return v_list



def make_vec(double a, double b, double c):
    return vec3(a, b, c)


def add_v(Vec3 a, Vec3 b):
    return add(a, b)

def sub_v(Vec3 a, Vec3 b):
    return sub(a, b)

def smul_v(double i, Vec3 a):
    return smul(i, a)

def sdiv_v(double i, Vec3 a):
    return sdiv(i, a)


def mag(Vec3 x):
    return magnitude(x)

    
def give(Vec3 a):
    return append(a)

def stack(Vec3 a, Vec3 b, Vec3 c):
    return append_sp(a, b, c)
