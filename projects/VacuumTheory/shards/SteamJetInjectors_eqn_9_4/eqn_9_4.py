from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_9_4__A_capE_capL_cap import eqn_9_4__AEL
from .eqn_9_4__S_capC_cap import eqn_9_4__SC
from .eqn_9_4__r import eqn_9_4__r
from .eqn_9_4__w_s import eqn_9_4__w_s

class SteamJetInjectors:
    eqn_9_4__AEL = eqn_9_4__AEL
    eqn_9_4__SC = eqn_9_4__SC
    eqn_9_4__r = eqn_9_4__r
    eqn_9_4__w_s = eqn_9_4__w_s

    @kwasak
    def eqn_9_4(self, AEL=None, SC=None, r=None, w_s=None):
        """
        w_s:= motive steam requirement
        r := pounds of steam required to compress 1 lb air from ejector suction pressure P_s to discharge pressure P_d
        SC := size correction factor
        """
        return
