# Chapter 5: Process Applications I
# 5-1 Phase Equilibria
"""
K_i := volatility
y_i := component concentration, vapor
x_i := component concentration, liquid
"""
K_i = y_i / x_i
# 5-2a Ease of Separation
alpha_1_2 = K_1 / K_2

# 5-2b Ease of Separation
K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)


# 5-3 Raoult's Law
"""
p_i := component partial pressure
x_i := liquid component mole fraction
P_0_i := pure component vapor pressure at equilibrium temperature
"""
p_i = x_i * P_0_i

# 5-4 Dalton's Law + Raoult's Law

y_i * P = x_i * P_0_i


# 5-5, relative volatlity in ideal liquid solution

alpha_12 = P_0_1 / P_0_2


# 5-6, Liquid Phase Nonideality I

p_i = x_i * gamma_i * P_0_i

# 5-7, Liquid Phase Nonideality II

y_i * P = x_i * gamma_i * P_0_i


# 5-8, Liquid Phase Nonideality III

alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)


# 5-9, Distillation

L_0 / V_1 = L_0 / (L_0 + D)

# 5-10a Equilibrium Distilliation Top
L_0 / V_1 = (L_0 / D) / (L_0 / D + 1)
# 5-10b Equilibrium Distilliation Reflux
L_0 / V_1 = R / (R + 1)
# 5-10c Equilibrium Distilliation Reflux
(L_0 / D) / (L_0 / D + 1) = R / (R + 1) 

# 5-11 Equilibrium Distilliation Bottoms Takeoff
L_N / V_0 = (V_0 + B) / V_0

# 5-12 Distillation Column
N_t = N_ES / Eff ** T

# 5-13 Required Packing Height
H_p = N_ES * HETP

# 5-14 Langmuir-Knudsen
W_E = 0.0583 * P_0 * (M / T) ** 0.5

# 5-15 separation efficiency
a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4

# 5-16 Henry's Law
p_i = x_i * H_i

# 5-17 Henry's Law 
"""
 math.log is "natural log, ln"
"""
log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
