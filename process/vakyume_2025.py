from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton, fsolve
from kwasak import kwasak_static
import pandas as pd
import numpy as np
from config import *
from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton, fsolve
from kwasak import kwasak_static
import pandas as pd
import numpy as np
from config import *

import tru
y = {}
for u,o in enumerate(filter(lambda o:str(o)[0].isalpha() and str(o)[0].capitalize()==str(o)[0] and str(o) not in map(lambda a:a.strip(),'I, Piecewise, LambertW, Eq, symbols'.split(',')),dir())):
    print(f'@@@{u+1}.',o, type(o))
    truth = False
    for tempt in range(budget:=5):
        try:
            truth = truth or tru.Verify(vars()[o]).verify() 
        except ValueError as ve:
            if (m:="math domain error") in str(ve):pass
            print("[ERROR]"+":"*99,m)
    print("+"*8*8,*((truth,) if (b:=isinstance(truth,bool)) else (truth.items())),sep=('\n\t'if not b else ''))
    y[o] = truth
print(*[yo for yo in y.items()],sep=('\n'))
def export_unfinished():
    return y
