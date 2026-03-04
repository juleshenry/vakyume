# Chapter 3 : Describing motion in one dimension
# 3-1 Position with constant velocity
"""
x := position
x_0 := initial position
v_x := velocity
t := time
"""
x = x_0 + v_x * t

# 3-2 Velocity with constant acceleration
"""
v := velocity
v_0x := initial velocity
ax := acceleration
t := time
"""
v = v_0x + ax * t

# 3-3 Position with constant acceleration
"""
x := position
x_0 := initial position
v_0x := initial velocity
ax := acceleration
t := time
"""
x = x_0 + v_0x * t + 0.5 * ax * t**2

# 3-4 Instantaneous velocity
"""
x := position
t := time
v_0 := initial velocity
a := acceleration
"""
v = v_0 + a * t

# 3-6 Position of object A in the x-reference frame
"""
x := position
x_0 := initial position
v_0 := initial velocity
a := acceleration
t := time
"""
x = v_0 * t + x_0

# 3-14 Relative acceleration
"""
a := relative acceleration
a_A := acceleration of the passenger
"""
a = a_A

# 3-15 Velocity of the passenger
"""
v := velocity of the passenger
v_B := velocity of the boat
v_A := velocity of the passenger
t := time
"""
v = (v_B + v_A) * t
