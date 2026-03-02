from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_2_26__D_cap import eqn_2_26__D
from .eqn_2_26__L_cap import eqn_2_26__L
from .eqn_2_26__P_downstream_cap import eqn_2_26__P_downstream
from .eqn_2_26__P_p_cap import eqn_2_26__P_p
from .eqn_2_26__P_upstream_cap import eqn_2_26__P_upstream
from .eqn_2_26__mu import eqn_2_26__mu
from .eqn_2_26__q import eqn_2_26__q

class FluidFlowVacuumLines:
    eqn_2_26__D = eqn_2_26__D
    eqn_2_26__L = eqn_2_26__L
    eqn_2_26__P_downstream = eqn_2_26__P_downstream
    eqn_2_26__P_p = eqn_2_26__P_p
    eqn_2_26__P_upstream = eqn_2_26__P_upstream
    eqn_2_26__mu = eqn_2_26__mu
    eqn_2_26__q = eqn_2_26__q

    @kwasak_static
    def eqn_2_26(self, D=None, L=None, P_downstream=None, P_p=None, P_upstream=None, mu=None, q=None):
        return
