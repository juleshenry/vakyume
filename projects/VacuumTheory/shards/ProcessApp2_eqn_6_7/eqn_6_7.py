from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_6_7__C_1_cap import eqn_6_7__C_1
from .eqn_6_7__C_2_cap import eqn_6_7__C_2
from .eqn_6_7__T_1_cap import eqn_6_7__T_1
from .eqn_6_7__T_2_cap import eqn_6_7__T_2
from .eqn_6_7__c_p import eqn_6_7__c_p
from .eqn_6_7__delta_h_c import eqn_6_7__delta_h_c
from .eqn_6_7__delta_h_v import eqn_6_7__delta_h_v
from .eqn_6_7__m_b import eqn_6_7__m_b
from .eqn_6_7__m_v import eqn_6_7__m_v

class ProcessApp2:
    eqn_6_7__C_1 = eqn_6_7__C_1
    eqn_6_7__C_2 = eqn_6_7__C_2
    eqn_6_7__T_1 = eqn_6_7__T_1
    eqn_6_7__T_2 = eqn_6_7__T_2
    eqn_6_7__c_p = eqn_6_7__c_p
    eqn_6_7__delta_h_c = eqn_6_7__delta_h_c
    eqn_6_7__delta_h_v = eqn_6_7__delta_h_v
    eqn_6_7__m_b = eqn_6_7__m_b
    eqn_6_7__m_v = eqn_6_7__m_v

    @kwasak_static
    def eqn_6_7(self, C_1=None, C_2=None, T_1=None, T_2=None, c_p=None, delta_h_c=None, delta_h_v=None, m_b=None, m_v=None):
        return
