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


class Solver:
    def permute(s, eqn, eqn_n):
        tokes = []
        for t in eqn.split(" "):
            t = t.strip()
            if t.isidentifier() and t not in {"ln", "log"}:
                tokes.append(t)
        normal_form = eqn.split("=")[1].strip() + " - " + eqn.split("=")[0].strip()
        print("ATTEMPT SOLVE", "`0 = ", normal_form, "`")

        for t in tokes:
            args = str(f"{TYPE}, ").join(sorted(filter(lambda x: x != t, tokes)))
            print(f'def eqn_{eqn_n.replace("-","_")}__{t}(' + args + f"{TYPE}):")
            solns = solve(normal_form, Symbol(t))
            if not len(solns):
                print("failed to solve")
                continue
            for soln in solns:
                print(f"{TAB}{t} = {soln}\n{TAB}return {t}")

    def analyze(s, i):
        root_dir = os.getcwd() + "/chapters/"
        get = list(filter(lambda x: i in x, os.listdir(root_dir)))[0]
        with open(root_dir + get) as file:
            eqn_number = ""
            for l in file.readlines():
                if x := re.compile("\d{1,2}-\d{1,2}").findall(l):
                    eqn_number = x[0]
                if " = " in l:
                    if "bhp_0" in l:
                        print("skipping", l)
                        continue
                    s.permute(l, eqn_number)


def reveal_blank_eqn_names():
    ix = 1
    for o in os.listdir(os.getcwd() + "/chapters"):
        if o[0] == "_":
            continue
        with open(os.getcwd() + "/chapters/" + o) as s:
            for l in s.readlines():
                if x := re.compile("\d{1,2}-\d{1,2}").findall(l):
                    eqn_number = x[0]
                    if len(l) < 10:
                        ix += 1
                        print(ix, l.strip(), "needs name!")


if __name__ == "__main__":
    X = Solver()
    for i in range(1,12):
        chap = ("0" if not len(str(i)) - 1 else "") + str(i)
        print(chap)
        try:
            X.analyze(chap)
        except:
            print('NotImplementedError')
    # reveal_blank_eqn_names()
    # except:
    # print(l)
    # Solver().permute("SS = S_Th * (P - p_s) / P ")
    # tokes = []
    # for t in tokes:
    # 	solve(,)

    # see_which_notes_are_valid_Python()
