# Chapter 2 : Fluid Flow in Vacuum Lines
import math

# 2-1 Reynold's number
"""
rho := density, lb/ft^3
D := pipe inside diam, ft
v := vel. ft/s
mu := viscosity, lb/ft*s
"""
Re = rho * D * v / mu
# 2-2 Maxwell-Boltzmann
"""
lambd := average mean free path , in
delta := mol. diam , in
psi:= mol. density molecules/in^3
"""
lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5

# 2-3 Knudsen's number
"""
D:= inside diameter, in
lambd:=avg. mean free path, in
"""
kn = lambd / D
# 2-4 Internal viscosity
"""
mu:=coefficient of viscosity
"""
_beta = mu * vel_grad
# 2-5 Hagen-Poiseuille
"""
q:=volumetric flow cm^3/s
D:= pipe diam.,cm
delta_P := upstream-downstream pressure, dyne/cm^3
L:=length, cm
mu:= coef. of visco., poise
"""
q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
# 2-6 Average Molecular velocity, cm/s
"""
mu :=viscosity, poise
rho:= density, g/cm^3
lambd:= mean free path, cm
"""
mu = 0.35 * rho * lambd * v_a
# 2-7 Maxwell-Boltzmann Velocity Distribution
"""
k:=boltz
T:= abs temp
m:= mass of a molecule
"""
v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
# 2-8 Critical point viscosity
"""
M:= mol. weight
T_c:= critical temp, K
P_c:= critical pressure, atm
"""
mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)


def eqn_2_8(M: float, P_c: float, T_c: float):
    mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
    return mu_c


# 2-10 Suction pressure
"""
delta_P := pressure loss
"""
Suc_Pres = oper_press - delta_P
# 2-11 , Darcy-Fanning isothermal flow
"""
f:= Moody friction
L:=length_pipe, ft
v:= velocity, ft/s
D:= inside diameter, ft
g_c:= dimensional constant, 32.2 lb * ft / lb * s
"""
h_r = f * L * v ** 2 / (D * 2 * g_c)

# 2-12 Pressure drop
"""
rho:= density, lb/ft^3
d:= pipe inside diameter, in
q:= vol. flow rate, ft^3/min
"""
delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
# 2-13 Pressure drop
"""
rho:= density, lb/ft^3
d:= pipe inside diameter, in
q:= vol. flow rate, ft^3/min
"""
delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
# 2-14 Room Temperature Pressure Drop
"""
Incompressibility not always valid assumption:
100 micron to 1 torr calculation.  
Rule of thumb: it holds for velocities less than 1/3 sonic velocity

v_s := sonic_velocity
k:=ratio of specific heat at constant temp to the specific heat at constant volume
"""
v_s = (k * g_c * R / M * T) ** 0.5

# 2-15 Turbulent flow smooth pipe, Blausius equation
"""
given Re < 2e5
"""
f = 0.316 / Re ** (0.25)
# 2-16, 2-17 Laminar Flow
"""
f = 64 / Re
delta_P = 0.0345* mu * L * v / d**2
delta_P = 0.105 * mu * L * q / d**4
"""

# 2-18a Noncircular ducts equivalent diameter
D_eq = 4 * R_ll
"""
R_ll := cross-sectional area of duct / wetted perimeter
"""
# 2-18b Noncircular ducts Hydraulic radius
R_ll = w * h / (2 * (w + h))
# 2-19a Noncircular Reynold's Number Hydraulic Radius
Re = 4 * R_ll * rho * v / mu
# 2-19b Noncircular Reynold's Number Rectangle
Re = (2 * w * h * rho * v) / ((w + h) * mu)
# 2-20 total Loss
"""
L:= laminar flow
"""
L = sum_pipe + sum_equivalent_length
# 2-22 Throughput Molecular Flow
"""
Q:= through_put, sucking pressure P
S_p:= dV / Dt
"""
Q_throughput = S_p * P_s
# 2-25 Conductance, Reciprocal of resistance, expressed in ft^3/min
"""
pressure loss
"""
C = Q_throughput / (P_1 - P_2)
# 2-26 Poiseuille's eqn for isothermal flow
q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)


# 2-28 Laminar conductance
C = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p


# 2-29 S, pumping speed
S_1 ** -1 = S_2 ** -1 + 1 / C


# 2-31 , General pump formula for overall speed
S_pump_speed = (S_p * C) / (S_p + C)


# 2-32 overall conductance series
1 / C_series = geometric_sum_C


# 2-33 overall conductance parallel

1 / C_paralell = arithmetic_sum_C


# 2-34, Transitional Flow, 1 > Kn > 0.01

C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)


# 2-35, relating laminar and transitional flow conductance
"""
F_P:= correction factor for Poiseuille's eqn from Figure 2-11
"""
C_T = C_L * F_p
# 2-36, Conductance of any vacuum system component
"""
C_0:=conductance thin walled aperture
F_t:=transmission prob. for component
"""
C = C_0 * F_t

# 2-37 Conductance of Cylindrical Line of Cross Sectional Area A
"""
F_t:= 1, for an aperture
"""
C = 38.3 * (T * A * F_t / M) ** 0.5
