# Chapter 6: Process Applications II
# 6-1 Refrigeration Capacity Mass Flow Heat Coefficient
w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
# 6-2 Refrigeration Capacity Q_r
w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
# 6-4 Vapor Load Steam Jet
w_v = 12000 * Q_v / delta_h_v
# 6-6 Makeup Water Flow Required to Replace Evaporative Losses
f_m = 24 * Q_r / delta_h_v = 0.0266 * Q_r
# 6-7 Average Vapor Load to Steam Jet in T
m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
# 6-8 Average vapor load to booster during time interval delta_t, lb/hr
w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
# 6-9 Vacuum Filter Volume Filtrate Poiseuille
dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
# 6-10 Filtrate Volum Neglecting Filter-Medium Resistance
dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
# 6-11a Vacuum Drying Total Residence Time A
t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T) 
# 6-11b Vacuum Drying Total Residence Time B
t_R= (delta_h_v * m_b * (m_i - m_l)) / (A_d * h_d * (T_w - T_s))
