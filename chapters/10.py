# 10-1
"""
Rotational speed of liquid ring pumps
sig_R := rotor tip speed ft/s
D_r := rotor Diameter
w := rotational speed
"""
sig_R = 0.00436 * D_r * w
# 10-2
PS = - V * dP / dt + Q
# 10-3
Q = 9.25 * N * T
# 10-4
t = V/S * ln ((SP_1 - Q) / (SP_2 - Q))
# 10-5
t = V/S * ln (P_1 / P_2)
# 10-6
S_a = V/t * ln (P_1 / P_2)
# 10-7
t = V(ln(P_1 / P_2) / S_1_2 + ln(P_2 / P_3) / S_2_3 + ...)
# 10-8 Effective Sealant Temperature
delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
# 10-9
T_c = T_s + delta_T
# 10-10
"""
bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
"""
# 10-11
"""
single-stage water sealed pumps compressing dry noncondensable loads
T_c = T_s + 10
"""
# 10-12
"""
two-stage water sealed pumps compressing dry noncondensable loads
T_c = T_s + 5
"""
# 10-13
"""
Saturated Load single-stage
T_c = T_s + 25
"""
# 10-14
"""
Saturated Load two-stage
T_c = T_s + 12
"""
# 10-15
"""
S = S_Th * (P - p_s) / P
"""
# 10-16
"""
Standard air capacitance
S_Th = S_0 * (P / (P - p_0)) ** 0.6
"""
# 10-17
"""
Dry air capacitance
S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
"""
# 10-18
"""
Pump condensation in Inlet

T_i := inlet  temperature of load

S = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
"""
# 10-19
"""
S = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
"""
# 10-20
"""
S_0 = S ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
"""
# 10-21
"""
Discharge pressure
P_prime := pseudo suction pressure
P_d := actual pump discharge pressure
P_prime = P / P_d * 760

"""
