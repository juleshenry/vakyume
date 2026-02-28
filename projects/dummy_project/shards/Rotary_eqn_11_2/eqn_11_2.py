from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from kwasak import kwasak_static

class Rotary:
    @kwasak_static
    def eqn_11_2(Q=None, Q_0=None, Q_external_gas_throughput=None, SP_1=None, SP_2=None, S_vol_pump_speed=None, V=None, t=None, **kwargs):
        return
