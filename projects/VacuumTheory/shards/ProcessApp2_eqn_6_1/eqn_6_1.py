from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_6_1__T_1_cap import eqn_6_1__T_1
from .eqn_6_1__T_2_cap import eqn_6_1__T_2
from .eqn_6_1__T_R_cap import eqn_6_1__T_R
from .eqn_6_1__c_p import eqn_6_1__c_p
from .eqn_6_1__del_h_v import eqn_6_1__del_h_v
from .eqn_6_1__w_1 import eqn_6_1__w_1
from .eqn_6_1__w_2 import eqn_6_1__w_2
from .eqn_6_1__w_v import eqn_6_1__w_v

class ProcessApp2:
    eqn_6_1__T_1 = staticmethod(eqn_6_1__T_1)
    eqn_6_1__T_2 = staticmethod(eqn_6_1__T_2)
    eqn_6_1__T_R = staticmethod(eqn_6_1__T_R)
    eqn_6_1__c_p = staticmethod(eqn_6_1__c_p)
    eqn_6_1__del_h_v = staticmethod(eqn_6_1__del_h_v)
    eqn_6_1__w_1 = staticmethod(eqn_6_1__w_1)
    eqn_6_1__w_2 = staticmethod(eqn_6_1__w_2)
    eqn_6_1__w_v = staticmethod(eqn_6_1__w_v)

    @kwasak_static
    def eqn_6_1(T_1=None, T_2=None, T_R=None, c_p=None, del_h_v=None, w_1=None, w_2=None, w_v=None, **kwargs):
        return
