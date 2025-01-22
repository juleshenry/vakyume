import os
import re
import timeout_decorator
from sympy import Symbol, solve, Eq
from llm import *
TAB = "    "
TYPE = ": float"
STD = 1
OUTFILE = "vakyume_2025.py"
MAX_COMP_TIME_SECONDS = 1 * 10
FUN = "*()/-+"


def stdout(s):
    if STD:
        with open(OUTFILE, "a+") as o:
            o.write(s + "\n")
    else:
        print(s)


class Solver:

    @timeout_decorator.timeout(MAX_COMP_TIME_SECONDS, timeout_exception=StopIteration)
    def get_solns_vanilla_nf(s, nf: str, symb: Symbol):
        try:
            return solve(nf, symb)
        except NotImplementedError:
            return []  # Unable to solve at present moment
        
    @timeout_decorator.timeout(MAX_COMP_TIME_SECONDS, timeout_exception=StopIteration)
    def get_solns_in_situ_nf(s, nf: str, symb: Symbol):
        try:
            #TODO: define in-situ the symbols and 
            #TODO: must be normal form as defined by python . perhaps requires another metaprogram to solve? 
            equation = Eq(0, nf)
            closed_form_solution = solve(equation, rho)
            return closed_form_solution[0] if closed_form_solution else None
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

    def tokenize(s, eqn, malos):
        # if any log\d+(), replace as log().... remove spaces trailing operator
        for m in malos:
            eqn = f"{m}".join([o.strip() for o in eqn.split(m)])
        # dilate functors
        eqn = eqn.split("#")[0]
        eqn = "".join(
            [
                (
                    f" {f} "
                    if 0 < i < len(eqn) - 1
                    and eqn[i + 1] != "*"
                    and eqn[i - 1] != "*"
                    and f in FUN
                    else f
                )
                for i, f in enumerate(eqn)
            ]
        )
        print("[dilated tokens]", eqn)
        return eqn.split(" ")

    def get_tokes(s, eqn):
        """
        Assumes equations are symbol-separated by spaces
        Anything like math.sin(a**2) must be treated via conversion to Sympy syntax
        B**2 is treated a token by space separation.
        `clean toke` cleans v**2 to v ** 2 for processing by Symbol
        """
        tokes = set()
        malos = {"ln", "log"}
        for t in s.tokenize(eqn, malos):
            clean = s.clean_t(t)
            # print(clean,clean.isidentifier())
            if s.valid_toke(clean) and not any(op in clean for op in malos):
                tokes.add(clean)
            elif any(op in clean for op in malos):  # TODO: approximate
                for m in malos:
                    clean = clean.replace(m, "")
                tokes.add(clean)

        def purge(t):
            for c in FUN:
                t = t.replace(c, "")
            return t

        tokes = set(filter(lambda a: a, map(purge, map(s.unexponentiate, set(tokes)))))
        print("tokens: ", tokes)
        return tokes

    def unexponentiate(s, t):
        return t.split("**")[0].strip()

    def sympy_backup_2(s, eqn_header, normal_form, token):
        stdout(TAB*2 + "# [Sympy Failover]")
        print(
            ans1 := escribir_codigo(
                eqn="0 = " + normal_form,
                single_variable=token,
                p1_i=4,
                p2_i=0,
            )
        )
        print(
            ans1 := escribir_codigo(
                eqn="",
                single_variable=token,
                header = eqn_header,
                pipin=ans1,
                p1_i=-1,
                p2_i=3,
            )
        )
        
        print(and1:=extract_code(ans1))
        and2 = s.fine_tune_extracted(and1)
        print("??$$"*6)
        print(and2)
        print("??$$"*6)
        for l in and2.split('\n'):
            stdout(l)

    @staticmethod
    def fine_tune_extracted(and1):
        """fte(extract_code) -> printable"""
        ans = ""
        tableau = lambda a:a.startswith(' '*4)or a.startswith('\t')
        for iii, line in enumerate(spl:=and1.split('\n')):
            # llm hacky
            line = line.replace("math.", "").replace('^','**') #LLm hacky
            if tableau(line):
                # attempt catch the naked, list-less return at the end of the solving code
                if iii + 1 == len(spl) and 'return' in line and not ('[' in line and ']' in line):
                    line = f'{TAB}return [{line.split('return')[1]} ]'
                ans += (TAB + line) + '\n'
        return ans
    def sympy_backup(s, eqn_header, normal_form, token):
        stdout(
            f"{TAB*2}# [FAILED TO PARSE VIA SYMPY]"
        )
        print(eqn_header)
        print(
            ans1 := escribir_codigo(
                eqn="0 = " + normal_form,
                lang="Python",
                single_variable=token,
                header = eqn_header,
                p1_i=1+1,
                p2_i=1,
            )
        )

        extract_0 = extract_code(ans1)
        print('!'*8**2)
        print(extract_0)
        print('!'*8**2)
        # extract_1 = extract_method_solving_for_(extract_0, token)
        # print(extract_1)
        # print('*'*8**2)

        for line in extract_0.split('\n'):
            # llm hacky
            line = line.replace("math.", "").replace('^','**') #LLm hacky
            if 0 and not line.strip() or any(m in line for m in ('import' , 'def' , '#', )):
                continue
            else:
                if 'return' in line and not ('[' in line and ']' in line):
                    #attempt to wrap formula
                    line = f'return [{line.split('return')[1]} ]'
                stdout(f"{TAB*2}{line.strip()}")
        # stdout(
        #     f"{TAB*2}pass # {'double parens issue. see reverse polish' if not ('** 0.' in normal_form or '**0.' in normal_form) else 'will not solve float exponential'}"
        # )

    def permute_and_print(s, eqn, eqn_n, comment=None):
        """1. gets tokes
        2. yields normal form
        3. injects comment
        4. forms method out of them

        """ 
        tokes = s.get_tokes(eqn)
        normal_form = (
            eqn.split("=")[1].strip().split("#")[0]
            + f" - ({eqn.split("=")[0].strip()})"
        )
        print("normal_form", normal_form)
        for token in tokes:
            # print("investigating",f"{token}")
            # print(normal_form)
            args = sorted(filter(lambda x: x != token, tokes))
            typed_args = str(f"{TYPE}, ").join(args)
            if typed_args:
                typed_args += TYPE
            stdout(f"\n{TAB}@staticmethod")
            stdout((eqn_header:=f'{TAB}def eqn_{eqn_n.replace("-","_")}__{token}({typed_args}):'))
            stdout(
                f"{TAB*2}# [.pyeqn] {eqn.strip().replace('#','')}"
            )  # original text contains #-comment for units
            '''TODO:
            if comment:
                stdout(f'{TAB*2}"""')
                for l in comment.split('\n')[1:-1]:
                    stdout(f"{TAB*2}{l}")
                stdout(f'{TAB*2}"""')
            '''
            try:
                solns = s.get_solns_vanilla_nf(normal_form, Symbol(token))
            except:
                solns = []
            if not len(solns):
                s.sympy_backup_2(eqn_header, normal_form, token)
                continue
            stdout(TAB * 2 + "result = []")
            for soln in solns:
                stdout(f"{TAB*2}{token} = {soln}")
                stdout(f"{TAB*2}result.append({token})")
            stdout(TAB * 2 + f"return result")

    def analyze(s, i):
        """1. opens a file in the chapters folder
        2. reads lines until equation is found.
        3. permutes the equation
        """
        root_dir = os.getcwd() + "/chapters/"
        get = list(filter(lambda x: i in x, os.listdir(root_dir)))[0]
        with open(root_dir + get) as file:
            eqn_number = comment = ""
            is_in_comment = 0
            # fmt: off
            for l in file.readlines():
                if '"""' in l:
                    is_in_comment =~ is_in_comment
                if is_in_comment:
                    comment = None #+= l
                else:
                    if x := re.compile(r"\d{1,2}-\d{1,2}\w{,2}").findall(l):
                        eqn_number = x[0]
                        comment = None #""
                    if " = " in l:
                        print("[DEBUG]", eqn_number,'\n',l)
                        # kludge:
                        # k = ''.join(comment.split('"""')[-2:]).split('\n')[1:-1]
                        comment = None # TODO:
                        s.permute_and_print(l, eqn_number, comment=comment)
                        # comment = ""
            # fmt: on


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
    stdout(
        "from math import log, sqrt, exp, pow, e\nfrom sympy import I, Piecewise, LambertW, Eq, symbols, solve"
    )
    for modules in sorted(os.listdir(os.getcwd() + "/chapters")):
        if modules[2].isalpha():
            continue  # __.* files
        chap, mods = modules.split("_")[0], modules.split("_")[1:]
        if int(chap) < 8:continue
        cls_name = "".join(x[0].upper() + x[1:] for x in mods)[:-3]
        print(
            chap,
            mods,
            cls_name,
        )
        stdout(f"\n\nclass {cls_name}:")
        X.analyze(chap)
