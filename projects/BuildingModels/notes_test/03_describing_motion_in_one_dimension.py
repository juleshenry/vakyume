# Chapter 3 : Describing motion in one dimension
# 3-1 Position with constant velocity
"""
x := position
x_0 := initial position
v_x := velocity
t := time
"""
x = x_0 + v_x * t

# 3.1 Velocity with constant speed
"""
v := velocity
x := distance
t := time
"""
v = (x) / t

# 3-2 Velocity with constant acceleration
"""
v := velocity
v_0x := initial velocity
a_x := acceleration
t := time
"""
v = v_0x + a_x * t

# 3-3 Position with constant acceleration
"""
x := position
x_0 := initial position
v_0x := initial velocity
a_x := acceleration
t := time
"""
x = x_0 + v_0x * t + 0.5 * a_x * t * * 2

# 3-4 Instantaneous velocity as time derivative of position
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
v_0 := initial velocity
t := time
v_B := velocity of the train
xA(t) := description of Alice's motion
xB(t) := position of Brice's origin
"""
x = v_B * t + xA

# 3-8 Velocity as measured by Igor (alternative form)
"""
v_A := velocity of object A
v_B := velocity of the train
t := time
xA(t) := description of Alice's motion
"""
v_A = v_B + vA

# 3-10 Velocity as measured by Igor (final form)
"""
v_A := velocity of object A
v_B := velocity of the train
v_A(t) := speed of Alice relative to the train
t := time
"""
v_A = v_B + vA

# 3-14 Relative acceleration
"""
a := relative acceleration
a_A := acceleration of passenger
"""
a = a_A

