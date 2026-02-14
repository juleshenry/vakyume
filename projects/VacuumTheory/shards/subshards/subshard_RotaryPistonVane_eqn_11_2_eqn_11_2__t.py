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




def eqn_11_2__t(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, **kwargs):
    # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
    # Error during Sympy solve: Sympy solve failed
    def func(t):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        return eval("(V / S_vol_pump_speed * ln( (SP_1 - (Q_exxernal_gas_xhroughpux + Q_0))/ (SP_2 - (Q + Q_0)))) - (x)".replace('x', str(t)))
    raise UnsolvedException("Pending LLM/Manual Repair")