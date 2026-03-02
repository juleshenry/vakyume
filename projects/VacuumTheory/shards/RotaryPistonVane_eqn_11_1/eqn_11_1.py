from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_11_1__PS_cap import eqn_11_1__PS
from .eqn_11_1__Q_0_cap import eqn_11_1__Q_0
from .eqn_11_1__Q_external_gas_throughput_cap import eqn_11_1__Q_external_gas_throughput
from .eqn_11_1__V_cap import eqn_11_1__V
from .eqn_11_1__dP_cap import eqn_11_1__dP
from .eqn_11_1__dT_cap import eqn_11_1__dT

class RotaryPistonVane:
    eqn_11_1__PS = eqn_11_1__PS
    eqn_11_1__Q_0 = eqn_11_1__Q_0
    eqn_11_1__Q_external_gas_throughput = eqn_11_1__Q_external_gas_throughput
    eqn_11_1__V = eqn_11_1__V
    eqn_11_1__dP = eqn_11_1__dP
    eqn_11_1__dT = eqn_11_1__dT

    @kwasak_static
    def eqn_11_1(self, PS=None, Q_0=None, Q_external_gas_throughput=None, V=None, dP=None, dT=None):
        return
