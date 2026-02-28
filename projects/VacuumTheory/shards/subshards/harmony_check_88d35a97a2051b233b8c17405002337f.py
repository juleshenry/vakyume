from math import *
import numpy as np

def check_harmony(hp, installed_costs, **kwargs):
    return (installed_costs) - (38000 * (hp / 10) ** 0.45)
