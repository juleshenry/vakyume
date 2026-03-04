from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_20_4__I_cap import eqn_20_4__I
from .eqn_20_4__R_capeff import eqn_20_4__Reff
from .eqn_20_4__V_capvoltmeter import eqn_20_4__Vvoltmeter


class ElectricCircuits:
    eqn_20_4__I = eqn_20_4__I
    eqn_20_4__Reff = eqn_20_4__Reff
    eqn_20_4__Vvoltmeter = eqn_20_4__Vvoltmeter

    @kwasak
    def eqn_20_4(self, I=None, Reff=None, Vvoltmeter=None):
        """
        Vvoltmeter := voltage measured by voltmeter
        Reff := effective resistance
        I := current
        RV := voltmeter resistance
        """
        return
