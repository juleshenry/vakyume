from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_6_8__C_1_cap import eqn_6_8__C_1
from .eqn_6_8__C_2_cap import eqn_6_8__C_2
from .eqn_6_8__T_1_cap import eqn_6_8__T_1
from .eqn_6_8__T_2_cap import eqn_6_8__T_2
from .eqn_6_8__c_p import eqn_6_8__c_p
from .eqn_6_8__delta_h_c import eqn_6_8__delta_h_c
from .eqn_6_8__delta_h_v import eqn_6_8__delta_h_v
from .eqn_6_8__delta_t import eqn_6_8__delta_t
from .eqn_6_8__m_b import eqn_6_8__m_b
from .eqn_6_8__w_v import eqn_6_8__w_v

class ProcessApp2:
    eqn_6_8__C_1 = staticmethod(eqn_6_8__C_1)
    eqn_6_8__C_2 = staticmethod(eqn_6_8__C_2)
    eqn_6_8__T_1 = staticmethod(eqn_6_8__T_1)
    eqn_6_8__T_2 = staticmethod(eqn_6_8__T_2)
    eqn_6_8__c_p = staticmethod(eqn_6_8__c_p)
    eqn_6_8__delta_h_c = staticmethod(eqn_6_8__delta_h_c)
    eqn_6_8__delta_h_v = staticmethod(eqn_6_8__delta_h_v)
    eqn_6_8__delta_t = staticmethod(eqn_6_8__delta_t)
    eqn_6_8__m_b = staticmethod(eqn_6_8__m_b)
    eqn_6_8__w_v = staticmethod(eqn_6_8__w_v)

    @kwasak_static
    def eqn_6_8(C_1=None, C_2=None, T_1=None, T_2=None, c_p=None, delta_h_c=None, delta_h_v=None, delta_t=None, m_b=None, w_v=None, **kwargs):
        return
