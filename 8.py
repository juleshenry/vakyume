# 8-1 Steam Jets
"""
Note: Need to inflation adjust
NS:= number ejector stages
NC:= number of condensors
SCON:=steam consumption based on 100-psig motive steam, lb/hr
Installation_Cost = 16000 * (NS + 2 * NC) (SCON / 1000) ** 0.35
"""
# 8-2 Liquid Ring Pumps
"""
hp:= horse power of pump
Installed_Costs = 33000 * (hp / 10) **0.5
"""

# 8-3 Rotary Piston or Two-Stage Rotary Vane Pumps
"""
Installed_Costs = 38000 * (hp / 10) **0.45

"""
# 8-4 Rotary Lobe Blowers
"""
Installed costs = 26000 * (hp / 10)** 0.4
"""

# 8-5
"""
E:= thermal efficiency
E = theoretical_adiabatic_horsepower / actual_brake_horsepower
"""

# 8-6
"""
deg_R:=absolute temperature
M:=molecular weight
R:=gas constant, 1544 ft*lb_f / (lb*mol) * deg_R
T:= absolute temperature, deg_R
P:= absolute pressure, torr
Adiabatic_hp = k / (k - 1) * (w * R *T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k-1)/k) - 1)
"""

# 8-7
"""
Adiabatic_hp = (w / 20) * ((P_2/ P_1) ** 0.286 - 1)
"""

# 8-8
"""
adiabatic_power_watts = f / 12 * ( (P_2 / P_1) **0.286 - 1)
"""

# 8-9
"""
E_j:=ejector thermal efficiency
e:=electrical cost, cents per kWh
s:=steam cost, dollar per 1000 lb
E_m:=mechanical pump thermal efficiency
r = 2.93 ( E_j * e) / (E_m * s)
"""
