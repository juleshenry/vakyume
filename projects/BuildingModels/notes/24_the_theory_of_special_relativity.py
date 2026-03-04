# Chapter 24 : The Theory of Special Relativity
# 24-2 Force per unit length
"""
FE := force per unit length
l := length
"""
FE = 2 * l / (2 * r)

# 24-? Lorentz Transformations
"""
x_0 := initial position
v := velocity
t := time
x := position
x_0 := initial position
v_0 := initial velocity
c := speed of light
"""
x = 0 * (x + v * t) / (t * v / c ** 2)

# 24-? Inverse Lorentz Transformations
"""
x_0 := initial position
v_0 := initial velocity
t := time
x := position
x_0 := initial position
v_0 := initial velocity
c := speed of light
"""
x = 0 * (x + v * t) + vx / c ** 2

# 24-? Lorentz addition of velocities
"""
u_0 := initial velocity
v := velocity
x := position
t := time
c := speed of light
"""
u_x = u_0 * x / (1 - v * u_x / c ** 2)

# 24-? Velocity transformations for all components
"""
u_0 := initial velocity
v := velocity
x := position
t := time
c := speed of light
"""
u_x = u_0 * x / (1 - v * u_x / c ** 2)

# 24-6.1 Lorentz Transformations for Velocity
"""
ux := relative velocity
u := initial velocity
v := velocity of the train
c := speed of light
"""
ux = u + v / (1 + (v * u) / c ** 2)

# 24-9 Relativistic kinetic energy
"""
m_0 := rest mass
c := speed of light
"""
K = (1 - u ** 2 / c ** 2) * m_0 * c ** 2

# 24-11 Energy-momentum relation
"""
E := energy
p := momentum
c := speed of light
"""
E ** 2 = p ** 2 * c ** 2 + m_0 ** 2

