import os
import re
import timeout_decorator
import httpx
from sympy import Symbol, solve, Eq
from llm import *
from suck_consts import *


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
    

    def valid_toke(s, t):
        valid = t.isidentifier() or t.split("**")[0].strip().isidentifier()
        return valid

    def clean_t(s, t):
        d = t.strip().replace("(", "").replace(")", "")
        o = d.split("**")
        if len(o) > 1 and o[1]:
            d = f"{t.split('**')[0]} ** {''.join(o[1:])}"
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
                    and f in FUNKTORZ
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


        # print(clean,clean.isidentifier())
        """
        tokes = set()
        malos = {"ln", "log"}
        for t in s.tokenize(eqn, malos):
            clean = s.clean_t(t)
            if s.valid_toke(clean) and not any(op in clean for op in malos):
                tokes.add(clean)
            elif any(op in clean for op in malos):  # TODO: approximate
                for m in malos:
                    clean = clean.replace(m, "")
                tokes.add(clean)

        def purge(t):
            for c in FUNKTORZ:
                t = t.replace(c, "")
            return t

        tokes = set(filter(lambda a: a, map(purge, map(s.unexponentiate, set(tokes)))))
        print("tokens: ", tokes)
        return tokes

    def unexponentiate(s, t):
        return t.split("**")[0].strip()

    def sympy_backup_2(s, eqn_header, normal_form, token):
        stdout(TAB*2 + "# [Sympy Failover]")
        try:
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
        except httpx.ConnectError as s:
            print(s)or stdout(TAB * 2 + 'pass # Ollama offline') 
            return
        if any(a in ans1.lower() for a in ('transcendental','numerical','not be possible',)):
            # llm hacky
            stdout(TAB * 2 + 'pass # no closed form solution')
            return 
        # print("BEFORE:::")
        # ans1 = make_sure_python_annotated(ans1)
        # print("AFTER:::")
        xtracted_code = extract_code(ans1)
        if not xtracted_code.strip(): # hacky, can come back un-annotated and llm is quantum noisy these days
            xtracted_code = ans1
        print("xtractd()()()")
        print(xtracted_code)
        print("xtractd()()()")
        and2 = s.fine_tune_extracted(xtracted_code, eqn_header)
        print("??$$"*6)
        print(eqn_header)
        print(and2)
        print("??$$"*6)
        for l in and2.split('\n'):
            stdout(l)

    # def doit():
    #     class VacuumTheory:pass
    #     print(sak_funx:=filter(lambda a:a.startswith('eqn') and '__'not in a,dir(VacuumTheory)))
    #     print(funx:=filter(lambda a:a.startswith('eqn') and '__' in a,dir(VacuumTheory)))

    #     # iterate all methods and fill with dummy values.
    #     import inspect
    #     a,b = map(lambda a:inspect.signature(getattr(VacuumTheory,a)),list(funx)[:2])
    #     woa = lambda d:str(d).replace(')','').replace('(','').replace(', ','').split(TYPE)
    #     tokes = set([u for u in woa(a)+woa(b)if u])

    #     for o in sak_funx:
    #         print(o)
    #         print(inspect.signature(getattr(VacuumTheory,o)))
    
    @staticmethod
    def fine_tune_extracted(and1, header):
        """fte(extract_code) -> printable"""
        ans = ""
        tableau = lambda a:a.startswith(TAB) or a.startswith('\t')
        print("&"*88)
        print(and1)
        print("&"*88)
        
        def find_last_return(code_block) -> int:
            for i,o in enumerate(reversed(w:=code_block.split('\n'))):
                if o.lstrip().startswith('return'):
                    return len(w) - i - 1
        lix = find_last_return(and1)
        for iii, line in enumerate(spl:=and1.split('\n')):
            line = line.replace("math.", "").replace('^','**') #LLm hacky
            if tableau(line):
                # attempt catch the naked, list-less return at the end of the solving code; technically should be index of last return
                if iii == lix and not ('[' in line and ']' in line):
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
        # herin, g
        stdout(f"{n+TAB}@kwasak_static")
        stdout(f"{TAB}def eqn_{eqn_n.replace("-","_")}({', '.join(f'{a}{TYPE} = None'for a in tokes)},**kwargs):")
        stdout(f"{TAB*2}return{n}")
        for token in tokes:
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
                        subnum = x[0].split('-')[1]
                        w=int(''.join([s for s in subnum if not s.isalpha()]))
                        # there is a bug in Sun Jan 26 20:42:05 CST 2025 keeping inspect signature from processing with kwasak. kwasak is flawed because mro() is unintuitive
                        eqn_number = x[0].split('-')[0] + '-' + ('0' + str(w) if w<=9 else str(w)) +  ''.join([s for s in subnum if s.isalpha()])
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
        # detect pyeqn files 
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

def make():
    X = Solver()
    stdout(
        "from math import log, sqrt, exp, pow, e"
    )
    stdout("from sympy import I, Piecewise, LambertW, Eq, symbols, solve")
    stdout("from scipy.optimize import newton")
    stdout("from kwasak import kwasak_static")
    for modules in sorted(os.listdir(os.getcwd() + "/chapters")):
        if modules[2].isalpha():
            continue  # __.* files
        chap, mods = modules.split("_")[0], modules.split("_")[1:]
        # if int(chap) != 8:continue
        cls_name = "".join(x[0].upper() + x[1:] for x in mods)[:-3]
        print(
            chap,
            mods,
            cls_name,
        )
        stdout(f"\n\nclass {cls_name}:")
        X.analyze(chap)

truify = r"""""
import tru
y = {}
for u,o in enumerate(filter(lambda o:str(o)[0].isalpha() and str(o)[0].capitalize()==str(o)[0] and str(o) not in map(lambda a:a.strip(),'I, Piecewise, LambertW, Eq, symbols'.split(',')),dir())):
    print(f'@@@{u+1}.',o, type(o))
    # try:
    truth = False
    for tempt in range(budget:=5):
        try:
            truth = truth or tru.Verify(vars()[o]).verify() 
        except ValueError as ve:
            if (m:="math domain error") in str(ve):pass
            # elif(m:=)
            print("[ERROR]"+":"*99,m)
            # print(str(ve));1/0
    print("+"*8*8,*((truth,) if (b:=isinstance(truth,bool)) else (truth.items())),sep=('\n\t'if not b else ''))
    y[o] = truth
print(*[yo for yo in y.items()],sep=('\n\t'))

    """
import subprocess
if __name__ == "__main__":
    make()
    with open(OUTFILE, "a") as f:
        f.write(truify)
    def run_outfile():
        result = subprocess.run(['python3', OUTFILE], capture_output=True, text=True)
        return result.stdout

    output = run_outfile()
    print(output)