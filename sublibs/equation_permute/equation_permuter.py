import os
import re
import time
from sympy import Symbol, solve, log
import subprocess

_DEBUG = False


def debug(s: str):
    if _DEBUG:
        print("[DEBUG]", s)


class Solver:
    def __init__(self, tab, typehint, outfile):
        self.tab, self.type, self.outfile = tab, typehint, outfile

    def stdout(self, s: str):
        if _DEBUG:
            print(s)
        else:
            with open(self.outfile, "a+") as o:
                o.write(s + "\n")

    def valid_toke(self, t):
        # print(t)
        valid = t.isidentifier() | t.split("**")[0].strip().isidentifier()
        # print(t.split('**')[0].strip().isidentifier(),'~',valid)
        return valid

    def clean_t(self, t):
        d = t.strip().replace("(", "").replace(")", "")
        o = d.split("**")
        if len(o) > 1 and o[1]:
            d = f"{t.split('**')[0]} ** {''.join(o[1:])}"
            # print(d)
        return d

    def get_tokes(self, eqn):
        """
        Assumes equations are symbol-separated by spaces
        Anything like math.sin(a**2) must be treated via conversion to Sympy syntax
        B**2 is treated a token by space separation.
        `clean toke` cleans v**2 to v ** 2 for processing by Symbol
        """
        tokes = set()
        for t in eqn.split(" "):
            clean = self.clean_t(t)
            # print(clean,clean.isidentifier())
            if self.valid_toke(clean) and t not in {"ln", "log"}:
                tokes.add(clean)
            # else:
            #     print("invalid toke",clean)
        tokes = list(tokes)
        # print('tokes',tokes)
        return tokes

    def permute(self, eqn: str, eqn_n: str):
        tokens = self.get_tokes(eqn)
        normal_form = eqn.split("=")[1].strip() + " - " + eqn.split("=")[0].strip()
        for token in tokens:
            # print("investigating",t)
            args = sorted(filter(lambda x: x != token, tokens))
            args = list(map(lambda x: x.split("**")[0].strip(), args))
            typed_args = str(f"{self.type}, ").join(args)
            if typed_args:
                typed_args += self.type
            # print(typed_args)
            token = token.split("**")[0]
            self.stdout(f"{self.tab}@staticmethod")
            self.stdout(
                f'{self.tab}def eqn_{eqn_n.replace("-","_")}__{token}({typed_args}):'
            )
            self.stdout(f"{self.tab}# {eqn.strip().replace('#','')}")
            solns = solve(normal_form, (token))
            if not len(solns):
                print("FAILED on:", rf"{normal_form}", rf"{token}", sep="\n")
                self.stdout(f"{self.tab*2}pass # unable  to solve")
                continue
            self.stdout(self.tab * 2 + "result = []")
            for soln in solns:
                self.stdout(f"{self.tab*2}{token} = {soln}")
                self.stdout(f"{self.tab*2}result.append({token})")
            self.stdout(self.tab * 2 + f"return {token}")

    def analyze(self, infile: str):
        full_file_path = os.getcwd() + "/" + infile
        with open(full_file_path) as file:
            eqn_number = ""
            for l in file.readlines():
                if x := re.compile("\d{1,2}-\d{1,2}\w{,2}").findall(l):
                    eqn_number = x[0]
                if " = " in l:
                    # print("[DEBUG]",eqn_number)
                    self.permute(l, eqn_number)


def get_class_name(pyeqn_file_path: str):
    # ends in .pyeqn, length 6
    # print(pyeqn_file_path)
    pyeqn_file = pyeqn_file_path
    return "".join(x[0].upper() + x[1:] for x in pyeqn_file)[:-6].split('/')[-1]

def empty_file(file_path):
    with open(file_path, 'w'):
        pass

def main(tab, typehint, infile, outfile):
    # Check if the file is empty
    try:
        if os.path.getsize(os.getcwd() + '/'+outfile) != 0:
            print(f"The file {outfile} already exists. Overwrite? y/n")
            x = input()
            if x.lower() == 'y':
                empty_file(outfile)
            else:
                raise ValueError("File exists already.")
    except FileNotFoundError:
        print('Creating '+outfile)
        pass
    
    s = Solver(tab, typehint, outfile)
    # standard imports
    s.stdout("from scipy import sqrt")
    s.stdout("from sympy import Piecewise,Eqn")
    if not infile.endswith(".pyeqn"):
        raise ValueError(f"Bad Name : {infile}; Must be .pyeqn file")
    pyeqn_file = infile.split("_")
    cls_name = get_class_name(pyeqn_file)
    s.stdout(f"class {cls_name}:")
    s.analyze(infile)
    # Run black on the file
    subprocess.run(["black", outfile])


if __name__ == "__main__":
    TAB = "    "
    TYPE = ": float"
    INFILE = "fluid_flow_vacuum_lines.pyeqn"
    OUTFILE = "out.py"
    main(TAB, TYPE, INFILE, OUTFILE)
