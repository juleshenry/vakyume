from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_7_14a__A import eqn_7_14a__A
from .eqn_7_14a__Q_condensor_heat_duty import eqn_7_14a__Q_condensor_heat_duty
from .eqn_7_14a__U import eqn_7_14a__U
from .eqn_7_14a__del_T_LM import eqn_7_14a__del_T_LM

class Precondensors:
    eqn_7_14a__A = eqn_7_14a__A
    eqn_7_14a__Q_condensor_heat_duty = eqn_7_14a__Q_condensor_heat_duty
    eqn_7_14a__U = eqn_7_14a__U
    eqn_7_14a__del_T_LM = eqn_7_14a__del_T_LM

    @kwasak
    def eqn_7_14a(self, A=None, Q_condensor_heat_duty=None, U=None, del_T_LM=None):
        return
