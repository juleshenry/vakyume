from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_11_2__Q_cap import eqn_11_2__Q
from .eqn_11_2__Q_cap_0 import eqn_11_2__Q_0
from .eqn_11_2__Q_cap_external_gas_throughput import eqn_11_2__Q_external_gas_throughput
from .eqn_11_2__S_capP_cap_1 import eqn_11_2__SP_1
from .eqn_11_2__S_capP_cap_2 import eqn_11_2__SP_2
from .eqn_11_2__S_cap_vol_pump_speed import eqn_11_2__S_vol_pump_speed
from .eqn_11_2__V_cap import eqn_11_2__V
from .eqn_11_2__t import eqn_11_2__t

class Rotary:
    eqn_11_2__Q = eqn_11_2__Q
    eqn_11_2__Q_0 = eqn_11_2__Q_0
    eqn_11_2__Q_external_gas_throughput = eqn_11_2__Q_external_gas_throughput
    eqn_11_2__SP_1 = eqn_11_2__SP_1
    eqn_11_2__SP_2 = eqn_11_2__SP_2
    eqn_11_2__S_vol_pump_speed = eqn_11_2__S_vol_pump_speed
    eqn_11_2__V = eqn_11_2__V
    eqn_11_2__t = eqn_11_2__t

    @kwasak
    def eqn_11_2(self, Q=None, Q_0=None, Q_external_gas_throughput=None, SP_1=None, SP_2=None, S_vol_pump_speed=None, V=None, t=None):
        return
