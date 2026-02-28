from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_7_11__Q_condensor_heat_duty_cap import eqn_7_11__Q_condensor_heat_duty
from .eqn_7_11__U_v_cap import eqn_7_11__U_v
from .eqn_7_11__V_c_cap import eqn_7_11__V_c
from .eqn_7_11__del_T_LM_cap import eqn_7_11__del_T_LM

class Precondensors:
    eqn_7_11__Q_condensor_heat_duty = staticmethod(eqn_7_11__Q_condensor_heat_duty)
    eqn_7_11__U_v = staticmethod(eqn_7_11__U_v)
    eqn_7_11__V_c = staticmethod(eqn_7_11__V_c)
    eqn_7_11__del_T_LM = staticmethod(eqn_7_11__del_T_LM)

    @kwasak_static
    def eqn_7_11(Q_condensor_heat_duty=None, U_v=None, V_c=None, del_T_LM=None, **kwargs):
        return
