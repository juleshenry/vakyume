from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_8_3__hp(self, installed_costs, **kwargs):
    # [.pyeqn] installed_costs = 38000 * (hp / 10) ** 0.45
    result = []
    hp = 6.64818534458494e-10 * installed_costs ** (20 / 9)
    result.append(hp)
    hp = (
        -0.326678815618226 * installed_costs**0.111111111111111
        - 0.118901365050371 * I * installed_costs**0.111111111111111
    ) ** 20
    result.append(hp)
    hp = (
        -0.326678815618226 * installed_costs**0.111111111111111
        + 0.118901365050371 * I * installed_costs**0.111111111111111
    ) ** 20
    result.append(hp)
    hp = (
        -0.173822167159837 * installed_costs**0.111111111111111
        - 0.301068825002567 * I * installed_costs**0.111111111111111
    ) ** 20
    result.append(hp)
    hp = (
        -0.173822167159837 * installed_costs**0.111111111111111
        + 0.301068825002567 * I * installed_costs**0.111111111111111
    ) ** 20
    result.append(hp)
    hp = (
        0.0603678051308443 * installed_costs**0.111111111111111
        - 0.342362835728782 * I * installed_costs**0.111111111111111
    ) ** 20
    result.append(hp)
    hp = (
        0.0603678051308443 * installed_costs**0.111111111111111
        + 0.342362835728782 * I * installed_costs**0.111111111111111
    ) ** 20
    result.append(hp)
    hp = (
        0.266311010487382 * installed_costs**0.111111111111111
        - 0.223461470678411 * I * installed_costs**0.111111111111111
    ) ** 20
    result.append(hp)
    hp = (
        0.266311010487382 * installed_costs**0.111111111111111
        + 0.223461470678411 * I * installed_costs**0.111111111111111
    ) ** 20
    result.append(hp)
    return result
