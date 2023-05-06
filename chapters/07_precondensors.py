# Chapter 7 : Precondensors
# 7-1 Dalton's Law
y_i = p_i / P
# 7-2 Equilbirium Vapor Pressure
p_i = x_i * P_i_0
# 7-3 Partial Pressure Condensable Vapor
p_i = x_i * epsilon_i * P_i_0
# 7-4 Partial Pressure Noncondensibles
p_nc = P - p_c
n_i / n_nc = p_i / p_nc = p_i / (p - P_c)
# 7-5 Molar Flow Rate
N_i = N_nc * (p_i) / (P - P_c)
# 7-6 Noncondensable Mass Flow Rate
W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
# 7-7 Liquid Phase Nonideality
W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
# 7-8 Spray and Cascade Condensor Required Cooling Liquid
L_c = Q / (c_p * del_T)
# 7-9 Cooling Liquid Volumetric Flow Rate
L_c = Q / (c_p * del_T * rho * 8.02)
# 7-10 Cooling Water Volumetric Flow Rate
L_c_P = Q / (500 * del_T)
# 7-11 Required Condensing Zone Volume
V_c = Q / (U_v * del_T_LM)
# 7-12 Surface Condensor Heat Transfer
Q = U * A * del_T
# 7-14 Temperature Difference
A = Q / (U * del_T_LM) = (Q / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
# 7-15 Heat Transfer Difference
1 / U = SIGMA(R)
# 7-16 Heat Transfer with Tube Wall Resistance
1 / U_0 = (
    1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
)
# 7-17 Heat Transfer Film Difference
R_0 = R_nc + 1 / h_c
# 7-18 Noncondensable Film Resistances
1 / U_0 = (
    R_nc
    + 1 / h_c
    + R_fo
    + (x_w * D_0) / (k_w * D_LM)
    + R_fi * D_0 / D_i
    + D_0 / (h_i * D_i)
)
