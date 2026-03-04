from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_8_5__E_capff import eqn_8_5__Eff
from .eqn_8_5__actual_brake_horsepower import eqn_8_5__actual_brake_horsepower
from .eqn_8_5__theoretical_adiabatic_horsepower import (
    eqn_8_5__theoretical_adiabatic_horsepower,
)


class SelectingPump:
    eqn_8_5__Eff = eqn_8_5__Eff
    eqn_8_5__actual_brake_horsepower = eqn_8_5__actual_brake_horsepower
    eqn_8_5__theoretical_adiabatic_horsepower = (
        eqn_8_5__theoretical_adiabatic_horsepower
    )

    @kwasak
    def eqn_8_5(
        self,
        Eff=None,
        actual_brake_horsepower=None,
        theoretical_adiabatic_horsepower=None,
    ):
        """
        Eff:= thermal efficiency
        """
        return
