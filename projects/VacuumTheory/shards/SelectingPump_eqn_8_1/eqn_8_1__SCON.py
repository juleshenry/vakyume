from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_8_1__SCON(self, NC, NS, installation_cost, **kwargs):
    # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
    SCON = (installation_cost / (NC * 2) / 16000 - NS) / log(
        (NS + self.__class__.eneral_constants.log(NC + 2 * NC)) / 10 ** (7 / 35)
    )
    result = []
    SCON_value = eval(SCON)
    result.append(SCON_value)
    return [result]
