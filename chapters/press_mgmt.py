# Chapter 3: Pressure Measurement
# 3-1 Absolute Pressure
Abs_Pressure = BarometricPressure - Vacuum
# 3-2 Liquid manometer
P = G / (G_C * rho * H)
# 3-3 McLeod Gauge
"""
H1,H2 distances from the fluid levels to a fixed point, top of the left capillary
"""
P_P - P = H_2 - H_1
# 3-4 Boyle's Law

P * V = KAPPA

# 3-5 McLeod Boyle's
P_P = P * (V / V_P)

# 3-6 Combined  McLeod Boyle's
P = V_P * (H_2 - H_1) / (V - V_P)
# 3-8 McLeod Bulb Volume
"""A_C cross sectional area"""
V_P = A_C * H_2
# 3-9 Torr Pressure in  Bulb McLeod
P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
# 3-11 Zeroed McLeod
P = A_C / V * (H_2) ** 2
# 3-12 Zeroed McLeod, Kappa Gauge Constant
"""
KAPPA := A_C / V, THE `GAUGE CONSTANT`
"""
P = KAPPA_1 * H_2 ** 2
# 3-13 Left Capillary Compressed Volume
P = KAPPA_2 * (H_2 - H_1)
# 3-15 Practical McLeod Minimum Volume
"""
V_PMIN := `PRACTICAL MIN, 1982`
"""
V_PMIN = 3.141592653589793 / 4
# 3-16 Maximum Compression Ratio
"""
~255000
"""
V_div_V_P_MAX = 200000 / (3.141592653589793 / 4)
# 3-17 Minimum Detectable McLeod Millimeter Change
"""
3.9*10**-6
"""
P_MIN = (3.141592653589793 / 4) / (200000)
