from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_2_31__C_cap import eqn_2_31__C
from .eqn_2_31__S_p_cap import eqn_2_31__S_p
from .eqn_2_31__S_pump_speed_cap import eqn_2_31__S_pump_speed


class FluidFlowVacuumLines:
    eqn_2_31__C = eqn_2_31__C
    eqn_2_31__S_p = eqn_2_31__S_p
    eqn_2_31__S_pump_speed = eqn_2_31__S_pump_speed

    @kwasak
    def eqn_2_31(self, C=None, S_p=None, S_pump_speed=None):
        return
