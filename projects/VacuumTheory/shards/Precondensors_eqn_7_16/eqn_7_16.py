from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_7_16__D_0_cap import eqn_7_16__D_0
from .eqn_7_16__D_LM_cap import eqn_7_16__D_LM
from .eqn_7_16__D_i_cap import eqn_7_16__D_i
from .eqn_7_16__R_f_0_cap import eqn_7_16__R_f_0
from .eqn_7_16__R_fi_cap import eqn_7_16__R_fi
from .eqn_7_16__U_0_cap import eqn_7_16__U_0
from .eqn_7_16__h_0 import eqn_7_16__h_0
from .eqn_7_16__h_i import eqn_7_16__h_i
from .eqn_7_16__k_w import eqn_7_16__k_w
from .eqn_7_16__x_w import eqn_7_16__x_w

class Precondensors:
    eqn_7_16__D_0 = eqn_7_16__D_0
    eqn_7_16__D_LM = eqn_7_16__D_LM
    eqn_7_16__D_i = eqn_7_16__D_i
    eqn_7_16__R_f_0 = eqn_7_16__R_f_0
    eqn_7_16__R_fi = eqn_7_16__R_fi
    eqn_7_16__U_0 = eqn_7_16__U_0
    eqn_7_16__h_0 = eqn_7_16__h_0
    eqn_7_16__h_i = eqn_7_16__h_i
    eqn_7_16__k_w = eqn_7_16__k_w
    eqn_7_16__x_w = eqn_7_16__x_w

    @kwasak_static
    def eqn_7_16(self, D_0=None, D_LM=None, D_i=None, R_f_0=None, R_fi=None, U_0=None, h_0=None, h_i=None, k_w=None, x_w=None):
        return
