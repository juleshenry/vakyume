# Chapter 4 : Describing motion in multiple dimensions
# 4-1 Average velocity vector
"""
t := time
x_1 := initial x-coordinate
y_1 := initial y-coordinate
x_2 := final x-coordinate
y_2 := final y-coordinate
"""
v = (x_2 - x_1) / t + (y_2 - y_1) / t

# 4-2 Instantaneous velocity vector
"""
v := velocity
r := position
t := time
"""
v = d r / dt

# 4-3 Instantaneous velocity vector
"""
v := velocity
x := x-component
y := y-component
t := time
"""
v = vx * x + vy * y

# 4-4 Acceleration vector
"""
a := acceleration
x := x-component
y := y-component
t := time
"""
a = ax * x + ay * y

# 4-5 Velocity vector with constant acceleration
"""
v := velocity
v_0 := initial velocity
a := acceleration
t := time
"""
v = v_0 + a * t

# 4-6 Position vector with constant acceleration
"""
r := position
r_0 := initial position
v_0 := initial velocity
a := acceleration
t := time
"""
r = r_0 + v_0 * t + 0.5 * a * t ** 2

