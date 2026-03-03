from cmath import *
from math import e, pi
import numpy as np

def check_harmony(_beta, mu, vel_grad, **kwargs):
    return (_beta) - (mu * vel_grad)
