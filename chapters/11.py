import math

# 11-1 Low Pressure Evacuation Time for Rotary-Piston/Rotary-Vane Pumps
"""
Q_0 := throughput of gas flow due to system outgassing
"""
PS = -V * dP / dT + Q + Q_0
# 11-2
"""
Q_0 := throughput of gas flow due to system outgassing, must be constant for 11-2
t = V/S * math.log( (SP_1 - (Q + Q_0))/ (SP_2 - (Q + Q_0)))
"""

# 11-3 System Factor
"""
t:= actual evacuation time
t_c:= calculated evacuation time using Eq 10-4
F_s:= system factor, based on operating experience

t = t_c ( F_s )
"""

# 11-4 Traps and Condensors, Quantity of vapor that can be pumped
"""
p_v := partial pressure of vapor at pump suction, torr
p_g := pressure of permanent gas at pump suction, torr
p_s := pump suction pressure, sum of partial pressure of vapor and partial pressure of permanent gas, torr
P_0_V := saturation pressure of vapor at pump operating temperature, torr
P_D := pump discharge pressure, torr
p_v / (p_v + p_g) = p_v / p_s <= P_0_v / P_D
"""

# 11-5
"""
p_v_max := maximum allowable partial pressure p_v_max of the process vapor at the pump suction
p_v_max = P_0_v * p_g / (P_D - P_0_v)
"""

# 11-6 Oil-sealed pump with a gas ballast
"""
P_0_v := saturation vapor pressure of a condensable vapor
S_B := maximum permissible gas ballast flow rate, ft^3/min
S_D := free air displacement of the vacuum pump, ft^3/min
p_b := partial pressure of vapor in the ballast gas, e.g. partial pressure of water vapor in ATM, torr
p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
"""
