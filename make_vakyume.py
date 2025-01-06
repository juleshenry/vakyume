import os
import re
import timeout_decorator
from sympy import Symbol, solve, log

TAB = "    "
TYPE = ": float"
STD = 1
OUTFILE = "vakyume_2025.py"
MAX_COMP_TIME_SECONDS = 5 * 60


def stdout(s):
    if STD:
        with open(OUTFILE, "a+") as o:
            o.write(s + "\n")
    else:
        print(s)


class Solver:

    @timeout_decorator.timeout(MAX_COMP_TIME_SECONDS, timeout_exception=StopIteration)
    def get_solutions(s, nf, symb):
        try:
            return solve(nf, symb)
        except NotImplementedError:
            return []  # Unable to solve at present moment

    def valid_toke(s, t):
        # print(t)
        valid = t.isidentifier() or t.split("**")[0].strip().isidentifier()
        # print(t.split('**')[0].strip().isidentifier(),'~',valid)
        return valid

    def clean_t(s, t):
        d = t.strip().replace("(", "").replace(")", "")
        o = d.split("**")
        if len(o) > 1 and o[1]:
            d = f"{t.split('**')[0]} ** {''.join(o[1:])}"
            # print(d)
        return d

    def get_tokes(s, eqn):
        """
        Assumes equations are symbol-separated by spaces
        Anything like math.sin(a**2) must be treated via conversion to Sympy syntax
        B**2 is treated a token by space separation.
        `clean toke` cleans v**2 to v ** 2 for processing by Symbol
        """
        tokes = set()
        for t in eqn.split(" "):
            clean = s.clean_t(t)
            # print(clean,clean.isidentifier())
            if s.valid_toke(clean) and t not in {"ln", "log"}:
                tokes.add(clean)
            # else:
            #     print("invalid toke",clean)
        tokes = list(tokes)
        # print('tokes',tokes)
        return tokes

    def unexponentiate(s, t):
        return t.split("**")[0].strip()

    def permute(s, eqn, eqn_n):
        """1. gets tokes
        2. yields normal form
        """
        tokes = s.get_tokes(eqn)
        normal_form = (
            eqn.split("=")[1].strip().split("#")[0]
            + f" - ({eqn.split("=")[0].strip()})"
        )
        print("!!!!", normal_form)
        tokes = list(map(s.unexponentiate, tokes))
        for token in tokes:
            # print("investigating",f"{token}")
            # print(normal_form)
            args = sorted(filter(lambda x: x != token, tokes))
            typed_args = str(f"{TYPE}, ").join(args)
            if typed_args:
                typed_args += TYPE
            stdout(f"\n{TAB}@staticmethod")
            stdout(f'{TAB}def eqn_{eqn_n.replace("-","_")}__{token}({typed_args}):')
            stdout(
                f"{TAB*2}# [.pyeqn] {eqn.strip().replace('#','')}"
            )  # original text contains #-comment for units
            try:
                solns = s.get_solutions(normal_form, Symbol(token))
            except RuntimeError as unable_to_compute:
                print(unable_to_compute)
                solns = []
            # print('solns',solns)
            if not len(solns):
                stdout(f"{TAB*2}pass # unable to solve")
                continue
            stdout(TAB * 2 + "result = []")
            for soln in solns:
                stdout(f"{TAB*2}{token} = {soln}")
                stdout(f"{TAB*2}result.append({token})")
            stdout(TAB * 2 + f"return {token}")
            # except:
            #     stdout(f"{TAB*2}pass #NotImplementedError")

    def analyze(s, i):
        """1. opens a file in the chapters folder
        2. reads lines until equation is found.
        3. permutes the equation
        """
        root_dir = os.getcwd() + "/chapters/"
        get = list(filter(lambda x: i in x, os.listdir(root_dir)))[0]
        with open(root_dir + get) as file:
            eqn_number = ""
            for l in file.readlines():
                if x := re.compile(r"\d{1,2}-\d{1,2}\w{,2}").findall(l):
                    eqn_number = x[0]
                if " = " in l:
                    print("[DEBUG]", eqn_number)
                    print(l)
                    s.permute(l, eqn_number)


class SetupMethods:
    # These were used to convert my notes to code and check formatting
    # Universally, they are the 1st step in the process
    @staticmethod
    def reveal_blank_eqn_names():
        ix = 1
        for o in os.listdir(os.getcwd() + "/chapters"):
            if o[0] == "_":
                continue
            with open(os.getcwd() + "/chapters/" + o) as s:
                for l in s.readlines():
                    if x := re.compile(r"\d{1,2}-\d{1,2}\w").findall(l):
                        eqn_number = x[0]
                        if len(l) < 10:
                            ix += 1
                            # print(ix, l.strip(), "needs name!")

    def see_which_notes_are_valid_Python(s):
        for o in os.listdir(os.getcwd() + "/chapters"):
            if o[0] == "_":
                continue
            with open(os.getcwd() + "/chapters/" + o) as s:
                eqn_number = ""
                for l in s.readlines():
                    if x := re.compile(r"\d{1,2}-\d{1,2}").findall(l):
                        eqn_number = x[0]
                    if "=" in l and ":" not in l.split("=")[0]:
                        s.parse_eqn(l.strip().split("#")[0])

    def parse_eqn(s, l):
        try:
            eval(l)
        except SyntaxError as se:
            to_tokens = se.text.split("=")[0]
            y = re.compile(r"\w[A-Za-z0-9\-\_]+").findall(to_tokens)
            # for o in y:
            #     print(o, "=", 1.1)
            # print(l)
        except NameError as ne:
            pass


if __name__ == "__main__":
    X = Solver()
    stdout("from math import log, sqrt, exp")
    for modules in sorted(os.listdir(os.getcwd() + "/chapters")):
        if modules[2].isalpha():
            continue  # __.* files
        chap, mods = modules.split("_")[0], modules.split("_")[1:]
        cls_name = "".join(x[0].upper() + x[1:] for x in mods)[:-3]
        print(
            chap,
            mods,
            cls_name,
        )
        stdout(f"\n\nclass {cls_name}:")
        X.analyze(chap)
