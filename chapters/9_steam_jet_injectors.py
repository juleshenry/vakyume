# Chapter 9 : Steam Jet Injectors
# 9-1 Ejector Operating Principle
"""
w_s := motive steam flow rate, lb/hr
v:= velocity
A:= cross sectional area, ft^2
rhos_s := motive steam density, lb/ft^3
"""
w_s = v * A * rho_s
# 9-2 Motive Steam Flow
"""
!flow coefficient of 0.97
d_n := nozzle throat diameter
P_m := motive steam pressure at point 1, psia
rhos_s := motive steam density at point 1, lb/ft^3
"""
w_s = 865.8 * d_n ** 2 * (P_m * rho_s) ** 0.5
# 9-3 Evacuation Hogger Estimation
"""
Order of magnitue accuracy estimate of time required for point specific ejector to evacuate a system
t_e := time required to evacuate system, minutes
P_s := design suction pressure of the ejector, torr
V := free volume of process system, ft^3
w_j := ejector capacity, 70 deg_F basis, lb/hr
"""
t_e = (2.3 - 0.003 * P_s) * V / w_j
# 9-4 Single-stage jets
"""
w_s:= motive steam requirement
r := pounds of steam required to compress 1 lb air from ejector suction pressure P_s to discharge pressure P_d
SC := size correction factor
"""
w_s = AEL * r * SC
# 9-5 Hogging jet
"""
w_h:= motive steam hogging
r_h:=pounds of 100-psig stream required per cubic foot
V:= process system free volume, ft^3
t_h := time permitted for evatuation, hr
"""
w_h = r_h * V / t_h
