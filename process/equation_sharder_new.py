from vakyume_2025_222 import export_unfinished as eu
from make_vakyume import Solver
import os
import shutil

o = eu()
S = Solver()


class Oowoo:
    def __init__(s):
        s.file_line_list = []

    def next_end(se, x):
        """proxima clase"""
        i = 0
        for i, l in enumerate(se.file_line_list[x + 1 :]):
            if "class" in l:
                return i + x + 2
            if (
                "import" in l
            ):  # yeah this is due to metaprogramming, pushing the script onto the end of the full-class dump
                return i + x + 2
        return len(se.file_line_list) - 1

    def fugly(se, v, v_eqn):
        with open("vakyume_2025_222.py") as f:
            se.file_line_list = list(f.readlines())
            aaa = se.get_def(v, v_eqn)
            es_clase = 0
            for il in enumerate(se.file_line_list):
                i, l = il
                if "class" in l:
                    es_clase = i
                if aaa in l:
                    return (es_clase, se.next_end(es_clase))

    def get_def(se, v, v_eqn):
        return f'def {v_eqn.replace("-","_")}__{v}'


ooo = Oowoo()


def process_equations(equations):
    """
    Processes a dictionary of equations grouped by class name, extracting equation names and their locations.

    Args:
        equations (dict): A dictionary where each key is a class name and each value is either:
            - a dictionary mapping variable names to equation strings (e.g., "eq1 x + y = z"), or
            - a boolean value (which will be ignored).

    Returns:
        tuple:
            - equation_names (dict): Maps each class name to a set of unique equation names extracted from the equations.
            - equation_locations (dict): Maps each class name to a location tuple as returned by `ooo.fugly(var_name, equation_name)`.

    Notes:
        - The function assumes that the equation name is the first word in the equation string.
        - Only classes with a dictionary of equations are processed; boolean values are skipped.
        - For each class, only the last computed location tuple is stored in `equation_locations`.
    """
    equation_names = {}  # Maps class name to set of equation names
    equation_locations = {}  # Maps class name to location tuples

    for class_name, equations_dict in equations.items():
        if not isinstance(equations_dict, bool):
            for var_name, equation in sorted(equations_dict.items()):
                equation_name = equation.split(" ")[0]
                location = ooo.fugly(var_name, equation_name)

                # Add equation name to set for this class
                if class_name not in equation_names:
                    equation_names[class_name] = set()
                equation_names[class_name].add(equation_name)

                # Store location tuple for this class
                equation_locations[class_name] = location

    return equation_names, equation_locations


wpp, wppp = process_equations(o)

print(wpp, wppp, "*" * 88, sep="\n")

raw_clipz = {k: "".join(ooo.file_line_list[v[0] : v[1] - 1]) for k, v in wppp.items()}


def murda(cli, base):
    nucli = ""
    mal = True
    clisp = cli.split("\n")
    gudz = wpp[base]
    print(gudz)
    prefix = "_".join(list(gudz)[0].split("_")[:2])
    q = lambda n: "".join(filter(lambda d: d.isnumeric(), n))
    rep = lambda n: q(n.split("__")[0].split("_")[-1])
    rep2 = lambda n: q(n.split("_")[2].split("(")[0])
    # below: bug; but now ignoring
    gudz_max = int(rep(max(gudz, key=lambda a: int(rep(a)))))
    badz = [
        prefix + "_" + ("0" + str(s) if int(s) < 9 else str(s))
        for s in range(1, gudz_max)
        if prefix + "_" + ("0" + str(s) if int(s) < 9 else str(s)) not in gudz
    ]
    for i, l in enumerate(clisp):
        # gonna assume to lookahead hear at cap at max line
        if i == len(clisp) - 1:
            break

        if any([j in clisp[i + 1] for j in gudz]):
            mal = False
            nucli += l + "\n"
        elif any([j in clisp[i + 1] for j in badz]) or (
            prefix in l and int(rep2(l)) > gudz_max
        ):
            mal = True

        elif not mal:
            nucli += l + "\n"
    return "\n".join(nucli.split("\n")[:-2])


ded_clipz = {}
for o in raw_clipz:
    ded_clipz[o] = murda(raw_clipz[o], o)


# Create shard directory if it doesn't exist
if not os.path.exists("shard"):
    os.makedirs("shard")
    shutil.copy("kwasak.py", "shard/kwasak.py")

# reshard a shard
n = "\n"


def f(s, k):

    inx = "😍"
    # print(k)
    wooo = []
    for i, l in enumerate(mama := k.split(n)):
        # print(l in  wpp[s].union({inx}))
        if any([j in l for j in wpp[s].union({inx})]):
            if inx not in l:
                inx = list(filter(lambda a: a in l, wpp[s]))[0]
                wooo += [i]
                # print(*["$$$$$"*5]*5,sep='\n')
            # else:
            #     inx = wpp[s].pop(s)
            # print(l)
    # print(wooo)
    # for a,b in zip(wooo, wooo[1:]+[len(mama)]):
    #     # print(a,b+1)
    #     print(*["@"*54]*5,sep='\n')
    #     print(w:=(n.join(  ['@kwasak_static'if i-1==0 or 0==len(mama)-1-i else mama[b]]+mama[a:b-1] )))
    return {
        s: [
            n.join(
                [
                    (
                        "    @kwasak_static"
                        if i - 1 == 0 or 0 == len(mama) - 1 - i
                        else mama[b]
                    )
                ]
                + mama[a : b - 1]
            )
            for a, b in zip(wooo, wooo[1:] + [len(mama)])
        ]
    }


# cut up clips into shards -> class header + eqn_X_Y__* .py
for s, k in ded_clipz.items():
    print(s)
    p = f(s, k)

    # s,k=list(ded_clipz.items())[0]
    # p = (f(s,k))
    for a_aa in p.items():
        a, aa = a_aa
        "class " + a + ":\n\n"
        for ii, aaa in enumerate(aa):
            z = sorted(wpp[s])[ii]
            # print(aaa[:93])
            print(z)
            with open(f"shard/{s}{z.replace("eqn","")}.py", "w") as ff:
                ff.write("from math import log, sqrt, exp, pow, e\n")
                ff.write("from sympy import I, Piecewise, LambertW, Eq, symbols, solve\n")
                ff.write("from scipy.optimize import newton\n")
                ff.write("from kwasak import kwasak_static\n")
                ff.write("import pandas as pd\n")
                ff.write("import numpy as np\n")
                ff.write("from kwasak import kwasak_static\n")
                ff.write("from suck_consts import *\n")
                ff.write("class " + a + ":\n\n")
                ff.write(aaa)


# save shard to shard fold
# run verify on every shard for budget until complete, allowing for transcendental,unsolvable ...



# wpp, wppp = {}, {}
# ooo = Oowoo()
# for a in o.items():
#     if type(a[1]) != bool:
#         # print(a[0]+'\n\t', '\n\t'.join(map(str,a[1].items())))
#         # print(a[0])
#         for aa in sorted(a[1].items()):
#             clip = ooo.fugly(aa[0], aa[1].split(" ")[0])
#             # wpp.update({aa[0]:})
#             wpp.update({a[0]: wpp.get(a[0], set()).union({aa[1].split(" ")[0]})})
#             wppp.update({a[0]: clip})