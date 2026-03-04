from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_8_1__SCON(self, NC: float, NS: float, installation_cost: float, **kwargs):
    # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
    # Solve for SCON:
    # Step 1: (SCON / 1000) ** 0.35 = installation_cost / (16000 * (NS + 2 * NC))
    # Step 2: SCON / 1000 = (installation_cost / (16000 * (NS + 2 * NC))) ** (1.0 / 0.35)
    # Step 3: SCON = 1000 * (installation_cost / (16000 * (NS + 2 * NC))) ** (1.0 / 0.35)
    SCON = 1000 * (installation_cost / (16000 * (NS + 2 * NC))) ** (1.0 / 0.35)
    return [SCON]
