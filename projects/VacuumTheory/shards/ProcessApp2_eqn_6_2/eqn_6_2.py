from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_6_2__Q_v_cap import eqn_6_2__Q_v
from .eqn_6_2__T_1_cap import eqn_6_2__T_1
from .eqn_6_2__T_2_cap import eqn_6_2__T_2
from .eqn_6_2__T_R_cap import eqn_6_2__T_R
from .eqn_6_2__c_p import eqn_6_2__c_p
from .eqn_6_2__w_1 import eqn_6_2__w_1
from .eqn_6_2__w_2 import eqn_6_2__w_2

class ProcessApp2:
    eqn_6_2__Q_v = eqn_6_2__Q_v
    eqn_6_2__T_1 = eqn_6_2__T_1
    eqn_6_2__T_2 = eqn_6_2__T_2
    eqn_6_2__T_R = eqn_6_2__T_R
    eqn_6_2__c_p = eqn_6_2__c_p
    eqn_6_2__w_1 = eqn_6_2__w_1
    eqn_6_2__w_2 = eqn_6_2__w_2

    @kwasak_static
    def eqn_6_2(self, Q_v=None, T_1=None, T_2=None, T_R=None, c_p=None, w_1=None, w_2=None):
        return
