# Chapter 4: Air Leakage
# 4-1 Strong Vacuum Acceptable Leakage
"""
1 <= P <= 10 #torr; 
Leaks result from metal porosities and cracks along weld lines from the following eqn
"""
W = 0.026 * P * 0.34 * V ** 0.6
# 4-2 Medium Vacuum Acceptable Leakage
"""
10 <= P <= 100 #torr
Leaks result from metal porosities and cracks along weld lines from the following eqn
"""
W = 0.032 * P * 0.26 * V ** 0.6
"""
# 4-3 Rough Vacuum Acceptable Leakage
100 <= P <= 760 #torr;
Leaks result from metal porosities and cracks along weld lines from the following eqns.
"""
W = 0.106 * V ** 0.6
# 4-4 High Vacuum Leakage around static or rotary seals
"""
D:= sealed diameter, in
theta:= specific leakage rate, lh * in/hr
!:= small_omega <= 10 lb/hr, the max acceptable air leakage rate for a component
Estimate acceptable air leakage resulting from leakage around static or 
rotary seals, valves, access ports, and other features for process. 
1 <= P <= 10 #torr; 
"""
small_omega = 3.141592653589793 * D * theta * P ** 0.34
# 4-5 Medium Vacuum Leakage around static or rotary seals
"""
Estimate acceptable air leakage resulting from leakage around static or 
rotary seals, valves, access ports, and other features for process. 
10 <= P <= 100 #torr
"""
small_omega = 1.2 * 3.141592653589793 * D * theta * P ** 0.26
# 4-6 Rough Vacuum Leakage around static or rotary seals
""" 
Estimate acceptable air leakage resulting from leakage around static or 
rotary seals, valves, access ports, and other features for process. 
100 <= P <= 760 #torr
"""
small_omega = 3.98 * 3.141592653589793 * D * theta
# 4-7 Total Acceptable Leakage Rate
W_T = W + sum_individual_leak_rates
# 4-8 Rough Vacuum Acceptable Air Leakage Rates
""" 
0.1 <= P < = 1 #torr,
"""
W = 0.026 * P ** 0.64 * V ** 0.6
# 4-9 Acceptable Air Leakage Rates
"""
0.1 <= P < = 1 #torr; 
"""
W = 3.141592653589793 * D * theta * P ** 0.64
# 4-10 Air Leakage over time
leakage = 0.0059 * V * del_P / t * 530 / T  # lb/hr
