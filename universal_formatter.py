import os
import re
import time


# def see_which_notes_are_valid_Python():
# 	for o in os.listdir(os.getcwd()+'/chapters'):
# 		if o[0]=='_':continue
# 		with open(os.getcwd()+'/chapters/'+o) as s:
# 			eqn_number = ""
# 			for l in s.readlines():
# 				if (x:=re.compile('\d{1,2}-\d{1,2}').findall(l)):
# 					eqn_number = x[0]
# 				if '=' in l and ':' not in l.split('=')[0]:
# 					parse_eqn(l.strip().split('#')[0])
# def parse_eqn(l):
# 	try:
# 		eval(l)
# 	except SyntaxError as se:
# 		to_tokens = se.text.split('=')[0]
# 		y = re.compile('\w[A-Za-z0-9\-\_]+').findall(to_tokens)
# 		for o in y:
# 			print(o,'=',1.1)
# 		print(l)
# 	except NameError as ne:
# 		pass
from sympy import Symbol, solve, log

TAB = "    "
TYPE = ": float"
STD = "WRITE"
OUTFILE = "vakyume_lib.py"


def stdout(s):
    if STD == "WRITE":
        with open(OUTFILE, "a+") as o:
            o.write(s + "\n")
    else:
        print(s)


class Solver:
    def get_tokes(s, eqn):
        tokes = set()
        for t in eqn.split(" "):
            clean_t = t.strip().replace("(", "").replace(")", "")
            # print(clean_t,clean_t.isidentifier())
            if clean_t.isidentifier() and t not in {"ln", "log"}:
                tokes.add(clean_t)
        tokes = list(tokes)
        return tokes

    def permute(s, eqn, eqn_n):
        tokes = s.get_tokes(eqn)
        normal_form = eqn.split("=")[1].strip() + " - " + eqn.split("=")[0].strip()
        stdout(f"# {eqn.strip().replace('#','')}")
        for t in tokes:
            args = sorted(filter(lambda x: x != t, tokes))
            typed_args = str(f"{TYPE}, ").join(args)
            if typed_args:
                typed_args += TYPE
            stdout(f'{TAB}def eqn_{eqn_n.replace("-","_")}__{t}({typed_args}):')
            try:
                solns = solve(normal_form, Symbol(t))
                if not len(solns):
                    stdout(f"{TAB*2}pass #failed to solve")
                    continue
                for soln in solns:
                    stdout(f"{TAB*2}{t} = {soln}\n{TAB*2}return {t}")
            except:
                stdout(f"{TAB*2}pass #NotImplementedError")

    def analyze(s, i):
        root_dir = os.getcwd() + "/chapters/"
        get = list(filter(lambda x: i in x, os.listdir(root_dir)))[0]
        with open(root_dir + get) as file:
            eqn_number = ""
            for l in file.readlines():
                if x := re.compile("\d{1,2}-\d{1,2}\w{,2}").findall(l):
                    eqn_number = x[0]
                if " = " in l:
                    print(eqn_number)
                    s.permute(l, eqn_number)


def reveal_blank_eqn_names():
    ix = 1
    for o in os.listdir(os.getcwd() + "/chapters"):
        if o[0] == "_":
            continue
        with open(os.getcwd() + "/chapters/" + o) as s:
            for l in s.readlines():
                if x := re.compile("\d{1,2}-\d{1,2}\w").findall(l):
                    eqn_number = x[0]
                    if len(l) < 10:
                        ix += 1
                        print(ix, l.strip(), "needs name!")


if __name__ == "__main__":
    X = Solver()
    # o = X.get_tokes("S_pump_speed = (S_p * C) / (S_p + C)")
    # print(o)

    for modules in sorted(os.listdir(os.getcwd() + "/chapters")):
        if modules[2].isalpha():continue #__.* files
        chap, mods = modules.split('_')[0], modules.split('_')[1:]
        cls_name = ''.join(x[0].upper() + x[1:] for x in mods)[:-3]
        stdout(f"\n\nclass {cls_name}:")
        X.analyze(chap)









