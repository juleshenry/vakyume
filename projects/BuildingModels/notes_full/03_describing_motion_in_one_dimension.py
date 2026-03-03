# Chapter 3 : Describing motion in one dimension
# Motion with constant speed (Example)
"""
x := position
t := time
v_0 := initial velocity
dx/dt := acceleration
"""
x = x_0 + v_0 * t + 0.5 * dx / dt * t * * 2
x_0 : = position at time zero, which is given as 0.5 m in this example
v_0 : = initial velocity of the grey ball, not specified in text but required to solve equations
dx / dt : = acceleration of the grey ball, also not provided in text and needs specification

# Motion with constant speed - Equation (3.4)
"""
x := position of the object
v_0 := initial velocity of the object in m/s
t := time in s
c := slope, which is zero for motion at a constant speed and thus can be ignored
m := mass of the object in kg
g := acceleration due to gravity in m/s^2 (constant)
"""
x = x_0 + v_0 * t

# Chapter 3-Motion with Constant Acceleration: Position as a Function of Time for an Object Released from Rest at t=0 s
"""
x := position
v_0 := initial velocity (m/s)
a_g := acceleration due to gravity (m/s^2, negative if downwards)
t := time (s)
x_0 := initial position (m)
"""
x = x_0 + v_0 * t - 0.5 * a_g * t * * 2

# Chapter 3-2 Motion with Constant Acceleration
"""
v := velocity
v_0 := initial velocity
a := acceleration (negative due to gravity)
x_0 := initial position
t := time
"""
dx = v * dt + 0.5 * a * dt * * 2

# Constant acceleration motion equations for one dimension
"""
x := position
a_0 := initial acceleration
C := constant term in velocity function
t := time
m := mass of cricket
g := gravitational acceleration
j := jumping rate constant
k := spring constant (if applicable)
v_0 := initial speed
A_0 := initial angular acceleration (if applicable)
"""
x = x_0 + 1 / 2 * a_0 * t * * 2
v = v_0 * cos + C_x
y = y_0 - 1 / 2 * g * t * * 2 (if jumping upwards or downwards from the ground level, otherwise use kinematic equations for projectile motion with air resistance if applicable)
a_y = a_g + j * t + k * x * * 3 (if spring is involved in horizontal movement of cricket)

# Chapter 3-4 Relative Motion
"""
vA := velocity_of_passenger_on_boat
xA(t) := position_of_passenger_on_boat
vB := velocity_of_boat
aA := acceleration_of_passenger_as_measured_by_Igor
dtvA(t) := change_in_velocity_of_passenger_over_time_from_Igor's_perspective
dadtA(t) := a↔A(t)
"""
x = x_0 + v * t + 0.5 * a * t * * 2
v = v_0 + a * t
xA : = x_B - x_B(0) + x_A(0) - x_B(0)
vA = v_B + v_A - v_B
v = v_0'

# Rob's position with constant velocity
"""
x_r := rob's initial position
v_r := rob's speed
t := time since the start of the problem
a_r := acceleration due to gravity
"""
0 = x_r + v_r * t 
x_v : = velociraptor's initial position
v_v : = velociraptor's speed when Rob passes by
a_v : = velociraptor's acceleration
t_r : = time taken for the velociraptor to catch up with Rob from t = 0
x_v = x_v + v_v * t_r 
v_v : = initial velocity of velociraptor (unknown)
t_v : = time taken for the velociraptor to reach Rob's position from rest at t = 0, which equals x_r / a_v

# Chapter 3-7 Sample Problems and Solutions - Motion with Constant Acceleration
"""
v_0 := initial velocity in m/s
a_R := acceleration of Rob's car R in m/s^2
t_reaction := reaction time in s
x_0R := initial position of Rob's car R in meters
xV := final velocity of velociraptor in m/s
"""
v = v_0 + a_R * t - (3.7) Velocity as a function of time for Rob's car with constant acceleration and reaction time
t : = time after the start signal is given to both cars, including reaction time
x : = position
a_V : = acceleration of velociraptor in m / s * * 2
v = xV - vR 
t_total : = total reaction + travel times for both cars, including the initial delay of robots to start moving in s
x_0V : = final position of velociraptor when it catches up with Rob in meters
vR = v_0 - a_R * t 
t_reaction : = robots initial delay before starting to move, s
xR : = position of Rob's car R at any time t in meters
vV : = final velocity of velociraptor when it catches up with Rob in m / s
xV_0 : = initial position of velociraptor relative to the starting point of Rob, in meters
t_total : = robots reaction + travel times for both cars until they meet at a common location, s

