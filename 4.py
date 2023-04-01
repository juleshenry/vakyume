# 4-1
"""
!:=Estimate acceptable air leakage resulting from metal porosities and cracks along weld lines from the following eqns.
1 <= P <= 10 #torr; W = 0.026 * P*.34 * V **.6
"""
# 4-2
"""
!:=Estimate acceptable air leakage resulting from metal porosities and cracks along weld lines from the following eqns.
10 <= P <= 100 #torr; W = 0.032 * P*.26 * V **.6
"""
# 4-3
"""
!:=Estimate acceptable air leakage resulting from metal porosities and cracks along weld lines from the following eqns.
100 <= P <= 760 #torr; W = 0.106 * V **.6
"""
# 4-4
"""
D:= sealed diameter, in
theta:= specific leakage rate, lh * in/hr
!:=Estimate acceptable air leakage resulting from leakage around static or rotary seals, valves, access ports, and other features for process. 
!:= Small_omega <= 10 lb/hr, the max acceptable air leakage rate for a component
1 <= P <= 10 #torr; Small_omega = pi * D * theta * P **0.34
"""
# 4-5
"""
!:=Estimate acceptable air leakage resulting from leakage around static or rotary seals, valves, access ports, and other features for process. 
!:= Small_omega <= 10 lb/hr, the max acceptable air leakage rate for a component
10 <= P <= 100 #torr; Small_omega = 1.2* pi * D * theta * P **0.26
"""
# 4-6
"""
!:=Estimate acceptable air leakage resulting from leakage around static or rotary seals, valves, access ports, and other features for process. 
!:= Small_omega <= 10 lb/hr, the max acceptable air leakage rate for a component
100 <= P <= 760 #torr; Small_omega = 3.98 * pi * D * theta
"""
# 4-7
"""
W_T = W + Sigma (Small_omega)
"""
# 4-8
"""
Acceptable Air Leakage Rates
0.1 <= P < = 1 #torr, W = 0.026 * P ** 0.64 * V **0.6
"""
# 4-9
"""
Acceptable Air Leakage Rates
0.1 <= P < = 1 #torr; W = pi * D * theta * P ** 0.64 
"""
# 4-10
"""
leakage = 0.0059 * V * del_P / t * 530 / T#lb/hr
"""
# 4-11
"""
W = f( phi , aleph , sigma, del_P, SF ) 
"""
