from sympy import Symbol, solve, sympify, symbols

eqn_str = "V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0))) - t"
print(f"Testing equation: {eqn_str}")

try:
    # Sympy's sympify usually handles 'ln' as 'log' (natural log)
    # But let's check what it does.
    # We should define symbols explicitly to avoid conflicts like Q (assumption key)
    t, V, S_vol_pump_speed, SP_1, SP_2, Q_external_gas_throughput, Q_0, Q = symbols('t V S_vol_pump_speed SP_1 SP_2 Q_external_gas_throughput Q_0 Q')
    
    # We can pass local symbols to sympify
    expr = sympify(eqn_str, locals={'t': t, 'V': V, 'S_vol_pump_speed': S_vol_pump_speed, 
                                   'SP_1': SP_1, 'SP_2': SP_2, 
                                   'Q_external_gas_throughput': Q_external_gas_throughput, 
                                   'Q_0': Q_0, 'Q': Q})
    print(f"Sympified expression: {expr}")
    
    for var_name in ["t", "V", "Q"]:
        print(f"\nSolving for {var_name}...")
        solns = solve(expr, Symbol(var_name))
        print(f"Solutions: {solns}")
except Exception as e:
    import traceback
    traceback.print_exc()
