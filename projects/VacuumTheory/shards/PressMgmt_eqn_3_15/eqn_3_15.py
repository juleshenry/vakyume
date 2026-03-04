from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_3_15__V_cap_P_capM_capI_capN_cap import eqn_3_15__V_PMIN


class PressMgmt:
    eqn_3_15__V_PMIN = eqn_3_15__V_PMIN

    @kwasak
    def eqn_3_15(self, V_PMIN=None):
        """
        V_PMIN := `PRACTICAL MIN, 1982`
        """
        return
