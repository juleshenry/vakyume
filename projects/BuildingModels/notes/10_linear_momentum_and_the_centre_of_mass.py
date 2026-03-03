# Chapter 10 : Linear momentum and the centre of mass
# Chapter 10-2 Momentum in One Dimension for a Particle with Constant Mass
"""
m := mass
v := velocity
a := acceleration due to gravity (constant)
t := time
"""
v = v_initial + g * t
p : = momentum
dt : = dp / dt
dp = m * dv = m * a * dt = m * g * dt = g * dt
Integrating both sides with respect to time gives us velocity as a function of time. Since initial v = 0 (particle at rest):
v = m * g * t
p = mv = constant * m * g * t
dp / dt = - mg k̂

# Chapter 10: Linear Momentum and the Centre of Mass - Impulse Equation (10.3)
"""
J := impulse
F := force applied during time interval
t_avg := average force over time interval
t_start := initial time of contact
t_end := final time of contact
"""
J = F * t_avg

# Chapter 10-2 Change in Momentum and Impulse
"""
p := momentum
m := mass
v_i := initial velocity
v_f := final velocity (zero for both)
t := time of collision
Jnet := net impulse received by a particle during the period $t$
ϱF netdt := force times change in time, scalar quantity
"""
ϱF netdt = - (1 / 2) * m * v_i
Jnet = - m * v_i
v_i = $v$
↭F = ma = m$v
F = (9.26 m / s) = 740 N
F = (9.26 m / s) = 1853 N

# Chapter 10-5 Conservation of Momentum for a System with Two Particles
"""
p_i := momentum of the i-th particle
F_ext := external force acting on the system
ϱP := total momentum of the system = p_1 + p_2
t := time
"""
d(ϱP) / dt = F_ext

# Linear momentum and the centre of mass - Momentum conservation in a collision
"""
p_train := momentum of train before collision
p_car := initial momentum of car (which is zero)
p_total := total momentum = p_train + p_car
v_initial := speed of train before collision
m_train := mass of one train car
N_cars := number of cars in the train
v_final := final velocity after collision for both train and attached car
"""
p_total = m_train * v_initial + 0 (since p_car is zero)
v_final = N_cars / (N_cars + 1) * v_initial
m_train * v_initial = (N_cars + 1) * v_final

# Inelastic Collisions - Example 10-5: Two objects collide and stick together
"""
m_1 := mass of object 1
v_1i := initial velocity of object 1 before collision
m_2 := mass of object 2 (initially at rest)
v_f := final velocity after collision
r_coll := distance between centres of the two objects during collision
"""
(9.3) m_1 * v_1i = (m_1 + m_2) * v_f 
0 = (1 / 2) * (m_1 + m_2) * v_f * * 2 - KE_initial

# Elastic Collisions - Example Problems
"""
m_p := mass of proton
m_N := mass of nucleus
v_i := initial velocity of proton
v_f := final velocity of proton after collision
v'_N := final velocity of nucleus after elastic collision (unknown)
t := time at which the collision occurs (not needed for calculation, assumed to be instantaneous)
r := distance between colliding particles before and during collision (also not needed)
"""
m_p * v_i = m_p * 0 + m_N * v'_N 
m_N * v'_N = m_N * (v_f - v_i) 
v'_N = v_i 
0.5 * m_p * (0) * * 2 + 0.5 * m_N * v_f * * 2 = 0.5 * m_p * v_i * * 2 
0.5 * m_N * (v'_N) * * 2 + 0.5 * m_p * (0) * * 2 = 0.5 * m_p * v_i * * 2 
(v'_N) * * 2 = v_i * * 2 
v'_N = v_i

# Chapter 10-2 Elastic Collisions and Conservation Laws
"""
P := total momentum before collision
Prime_P := total momentum after collision
vM := velocity of mass M before collision
vM' := velocity of mass M after collision
vm := velocity of smaller mass m before collision
vm' := velocity of smaller mass m after collision
m := mass of the smaller object (mass m)
M := mass of larger object (mass M)
t := time duration or, in this case, not applicable as it is an instantaneous event
v_rel := relative speed between masses before and after collision
"""
Prime_P = Prime_V ↔
PM + pm = PM' + pm'
M + m = M' + m
m(vm' - vm) = m(v'→
↭ P = Prime_P
MvM + mv = Mv'→
Prime_E : = total mechanical energy before collision
Prime_E' : = total mechanical energy after collision
KE : = kinetic energy of a single object
PE : = potential energy of a single object (not needed here as it is constant)
v : = velocity of an individual mass
mV : = momentum of smaller mass m before the collision
MV : = momentum of larger mass M before the collision
Prime_E = Prime_KE' + Prime_KE''
1 / 2M(vM') * * 2 + 1 / 2mv * * 2 = 1 / 2M(v→
m - v) * * 2 = m(v→
↭ Prime_E = Prime_KE' + Prime_KE''
P = P ↔
MvM + mv = Mv→
K : = kinetic energy, K_1 : = kinetic energy of mass M, K_2 : = kinetic energy of smaller object m
PE_1 : = potential energy of larger object M (not needed here as it is constant)

# Linear momentum and centre of mass - Equation (3.1)
"""
p_cm := linear momentum of the system's centre of mass
m_total := total mass of the system
v_cm := velocity of the centre of mass
"""
p = m * v

# Chapter 10: Linear Momentum and Centre of Mass - Velocity Components Equations
"""
v_cmx := velocity of centre of mass in frame X direction
m1 := mass of block 1
m2 := mass of block 2 (identical to m1)
v1 := initial velocity of block 1 along the x axis
vCMx := center of mass velocity component in frame X direction
"""
P_total = m_tot * v_cmx

# Chapter-10: Centre of Mass and Linear Momentum
"""
F_ext := external force acting on a particle i
m_i := mass of the ith point particle in the system
ϱa_i := acceleration of the ith point particle relative to an inertial frame
ϱr_i := position vector of the ith point particle relative to some origin O5. 
ϱv_i := velocity vector of the ith point particle6.
"""
F_ext = m_i * a_i

# Chapter 10-3 Momentum and Centre of Mass Equations
"""
p := momentum
r := mass
v := velocity
t := time
a := acceleration
M_total := sum(ri * ri*v) for all i in system particles
zCM := (1/M_total)*sum(zi*zi*zi*dz) for all zi in bowl disks
"""
p = m * v
r : = mass of particle
t : = time
a : = acceleration acting on particle
M_total : = sum(ri * ri * v) over entire system particles
zCM : = (1 / M_total) * sum(zi * zi * dz * * 3) / h for all zi in bowl disks, h is infinitesimally small thickness of disk
p = m * v
r : = mass of particle i
t : = time
a : = acceleration acting on particle i
M_total : = sum(ri * ri * v) over entire system particles
zCM : = (1 / M_total) * sum(zi * (dz * * 3)) / (h * * 3) for all zi in bowl disks, h is infinitesimally small thickness of disk

# Important Equations Momentum and Impulse for a point particle
"""
p := momentum
d_t p = F net dt
F_net := force_net
r_cm := position_of_center_mass
v := velocity of the center of mass
a_cm := acceleration of the center of mass
M := total_mass
m := particle_mass
"""
p = m * v
d_t p = F net dt
F_net = d_t (m * r_cm)
r_cm = 0
M : = total_mass
a_cm : = acceleration of the center of mass
d_t Mϱv_cm = F net d_t r_cm
F_ext : = external_force
r_cm = 0
x_i : = position of the i - th particle in terms of x_cm
v_i : = velocity of the i - th particle in terms of v_cm
a_i : = acceleration of the i - th particle in terms of a_cm
M_i : = mass of the i - th particle
d_t M_i * v_i = 0
x_r : = position of centre of mass in terms of x
v_r : = velocity of centre of mass in terms of v
a_r : = acceleration of the centre of mass in terms of a
M_r : = total_mass of the system
d_t M_r * v_r = 0

# Problem 10-1 Bullet and Ballistic Pendulum Speed Calculation
"""
m := bullet mass
M := pendulum bob mass
h := height reached by pendulum after collision
L := length of the strings holding the pendulum
v_bullet := initial velocity of the bullet (unknown)
v_pendulum := final velocity of the pendulum and embedded bullet right after impact (unknown)
g := acceleration due to gravity
"""
m * v_bullet = (M + m) * v_pendulum * sin(θ) 
v_pendulum = sqrt(2gh) 
m * v_bullet = (M + m) * g * h 
v_bullet = ((M + m) * g * h) / m

# Chapter 10 - Linear Momentum and Centre of Mass
"""
E_cm := total mechanical energy in centre of mass frame
k := spring constant
x := compression distance
vM_cm := velocity of masses before collision
m := smaller block's mass
M := larger block's mass
g := acceleration due to gravity
h := initial height difference between blocks
"""
kx * * 2 = E_cm
x = sqrt((m + M) / (m + M)) * k * h / g
vM_cm = sqrt(g * h / ((m + M) - m))

# Chapter 10-4 Centre of Mass Y Position Calculation
"""
y_cm := y position of centre of mass
M := total mass of wire
R := radius of semi circle
d_m := dm = (mass per unit length) * ds
ds := differential arc length element on the wire
ϑ := angle corresponding to a point on the wire with respect to the x-axis
"""
y_cm = 1 = R2 dϑ = - R2 cos ϑ dϑ = - 2R * * 2 / 2 = - R * * 2 = M * (π / 2)

