from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_7_11__Q_cap_condensor_heat_duty import eqn_7_11__Q_condensor_heat_duty
from .eqn_7_11__U_cap_v import eqn_7_11__U_v
from .eqn_7_11__V_cap_c import eqn_7_11__V_c
from .eqn_7_11__del_T_cap_L_capM_cap import eqn_7_11__del_T_LM


class Precondensors:
    eqn_7_11__Q_condensor_heat_duty = eqn_7_11__Q_condensor_heat_duty
    eqn_7_11__U_v = eqn_7_11__U_v
    eqn_7_11__V_c = eqn_7_11__V_c
    eqn_7_11__del_T_LM = eqn_7_11__del_T_LM

    @kwasak
    def eqn_7_11(self, Q_condensor_heat_duty=None, U_v=None, V_c=None, del_T_LM=None):
        return
