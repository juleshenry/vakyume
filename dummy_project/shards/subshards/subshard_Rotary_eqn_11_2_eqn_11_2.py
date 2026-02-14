import math
from math import *
import numpy as np
import sympy
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
class UnsolvedException(Exception): pass
from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
import numpy as np



def eqn_11_2(Q=None, Q_0=None, Q_external_gas_throughput=None, SP_1=None, SP_2=None, S_vol_pump_speed=None, V=None, t=None, **kwargs):
    return