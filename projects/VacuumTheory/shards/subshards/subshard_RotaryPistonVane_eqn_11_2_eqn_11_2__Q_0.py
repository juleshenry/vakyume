import math
from math import *
import numpy as np
import sympy
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
class UnsolvedException(Exception): pass
from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
import numpy as np




def eqn_11_2__Q_0(Q: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float, **kwargs):
    # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
    # Error during Sympy solve: Sympy solve failed
    def func(Q_0):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        return eval("(V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + x))/ (SP_2 - (Q + x)))) - (t)".replace('x', str(Q_0)))
    raise UnsolvedException("Pending LLM/Manual Repair")