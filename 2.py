#2-1 Reynold's number
'''
Re = rho * D * v / mu
rho := density, lb/ft^3
D := pipe inside diam, ft
v := vel. ft/s
mu := viscosity, lb/ft*s
'''

#2-2 Maxwell-Boltzmann 
'''
lambda = (pi * delta**2 * psi * 2**.5)
lambda := average mean free path , in
delta := mol. diam , in
psi:= mol. density molecules/in^3
'''

#2-3 Knudsen's number
'''
kn = lambda / D
D:= inside diameter, in
lambda:=avg. mean free path, in
'''
#2-4 Internal viscosity
'''
beta = mu * vel_grad
mu:=coefficient of viscosity
'''
#2-5 Hagen-Poiseuille
'''
q = pi*D^4*delta_P / (128 * L)
q:=volumetric flow cm^3/s
D:= pipe diam.,cm
delta_P := upstream-downstream pressure, dyne/cm^3
L:=length, cm
mu = coef. of visco., poise
'''
#2-6 Average Molecular velocity, cm/s
'''
mu = 0.35 * rho * lambda * v_a
mu :=viscosity, poise
rho:= density, g/cm^3
lambda:= mean free path, cm
'''
#2-7
'''
v_a = ((8 * k * T)/(pi * m))**.5
k:=boltz
T:= abs temp
m:= mass of a molecule
 = 
'''

#2-8 Critical point viscosity
'''
mu_c = (7.7 * (M **.5 ) * P_c**(2/3)) / T_c**(1/6)

M = mol. weight
T_c = critical temp, K
P_c = critical pressure, atm
'''


#2-10 Suction pressure
'''
Suc_Pres = oper_press - delta_P
delta_P := pressure loss
'''

#2-11 , Darcy-Fanning isothermal flow
'''
h_r = f * L * v**2 / (D*2*g_c)

f:= Moody friction
L:=length_pipe, ft
v:= velocity, ft/s
D:= inside diameter, ft
g_c:= dimensional constant, 32.2 lb * ft / lb * s
'''

#2-12 Pressure drop
'''
delta_P = 4.31 * rho * f * L * v**2 / (2 * d * g)

delta_P = 2.15 * rho * f * L * q**2 / (d ** 5)

rho:= density, lb/ft^3
d:= pipe inside diameter, in
q:= vol. flow rate, ft^3/min
'''

#2-14
'''
Incompressibility not always valid assumption:
100 micron to 1 torr calculation.  
Rule of thumb: it holds for velocities less than 1/3 sonic velocity

v_s =  (k*g_c * R / M * T)**.5
v_s := sonic_velocity
k:=ratio of specific heat at constant temp to the specific heat at constant volume
'''

#2-15 Turbulent flow smooth pipe, Blausius equation

'''
f =  0.316 / Re**(.25)
given Re < 2e5

'''

#2-16, 2-17 Laminar Flow

'''
f = 64 / Re
delta_P = 0.0345* mu * L * v / d**2
delta_P = 0.105 * mu * L * q / d**4

'''

#2-18,2-19 Noncircular ducts
'''
D_eq = 4 * R_ll
R_ll := cross-sectional area of duct / wetted perimeter
R_ll = w * h / (2 *(w+h))
Re = 4 * R_ll * rho * v / mu = (2 * w * h * rho * v)/ ((w+h) * mu)
'''

#2-20
'''
L ~= sum(pipe) + sum(equiv. length)
L:= laminar flow
'''

#2-22
'''
Q = S_p * P_s
Q:= through_put, sucking pressure P, S_p = dV / Dt
'''

#2-25 Conductance, Reciprocal of resistance, expressed in ft^3/min
'''
C = Q / (P_1 - P_2)
pressure loss
'''

#2-26 Poiseuille's eqn for isothermal flow
'''
q* P` =  pi * D**4 / (128 * mu * L) * P` * (P_upstream - P_downstream)

'''

#2-28 Laminar conductance
'''
C = pi * D**4 / (128 * mu * L) * P`
'''

#2-29 S, pumping speed
'''
S_1**-1 = S_2**-1 + 1/C
'''

#2-31 , General pump formula for overall speed
'''
S = (S_p * C) / (S_p + C)
'''

#2-32
'''
1 / C_series = geometric_sum (C...)
'''

#2-33
'''
1 / C_paralell = arithmetic_sum (C...)
'''

#2-34, Transitional Flow, 1 > Kn > 0.01
'''
C = C_1 * (D**4 / (mu * L)) * P` + C_2 * (D**3 / L) 
'''

#2-35, relating laminar and transitional flow conductance
'''
C_T =  C_L * F_p
F_P:= correction factor for Poiseuille's eqn from Figure 2-11
'''

#2-36, Conductance of any vacuum system component
'''
C =  C_0 * F_t
C_0:=conductance thin walled aperture
F_t:=transmission prob. for component
'''

#2-37
'''
C = 38.3 * (T * A * F_t / M )**.5
F_t:= 1, for an aperture
'''
