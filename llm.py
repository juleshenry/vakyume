import ollama
import secrets

n = "\n"
null_str = ""


def escribir_codigo(eqn, lang, p1_i=None, p2_i=None):
    int(sum([p1_i, p2_i,]))
    # print(f"solving { eqn }, buddy.")
    fluff = secrets.choice(["C++", "Erlang", "Haskell", "Rust", "Lisp"])
    response = ollama.chat(
        model="llama3:latest",
        messages=[
            {
                "role": "user",
                "content": (
                    p1 := (
                        f"""
                        You are a {lang} / {fluff} wizard giga chad developer who wakes up peer-reviewing arxiv.org
                          for News.YCombinator.com threads and sleeps coding {lang} at 3am.
                        Solve this equation: {eqn} in {lang}
                        """,
                        f"""
                        Hi Doctor Math, I need help solving this {eqn}  Give back the closed form expression for T
                        
                        """,
                        f"""identify the line that is a one liner for a closed-form expression for T:\n{eqn}""",
                    )[p1_i]
                )
                + (
                    p2 := (
                        null_str,
                        """print the python code syntax for an equation solving the formula for every value ablated ( taken one at a time out)
                    example: a**2 = b**2 + c**2


                    winning correct output:  a = [math.sqrt(b**2+c**2),math.sqrt(b**2+c**2),]; b = [math.sqrt(a**2-c**2),math.sqrt(a**2-c**2),];c = [math.sqrt(a**2-b**2),math.sqrt(a**2-b**2),];""",
                        """Give your answer as one line.""",
                    )[p2_i]
                ),
            },
        ],
    )
    dummy = lambda aa: print((n + "*" * 88 + n).join(["", aa, ""]))
    dummy(p1)
    dummy(p2)
    return n.join(filter(lambda a: a, response["message"]["content"].split(n)))

def eval_wrap(ans):
    max_bounds =[ 0,len(o:=ans.split(n))]
    for i,j in (s:=enumerate(len(o))):
        for ii,jj in s:
            try:
                eval('\n'.join(o[i:ii]))
                print('CODE',o[i:ii])
                if ii - i > max_bounds[1]-max_bounds[0]:
                    max_bounds = [i,ii]
            except:
                pass
            #O(n^2)
    print(n.join(o[max_bounds[0]:max_bounds[1]]))
print(
    ans1 := escribir_codigo(
        eqn="q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)",
        lang="Python",
        p1_i=1,
        p2_i=1,
    )
)
print(ans2 := escribir_codigo(eqn=ans1, lang=None, p1_i=2, p2_i=2))  # or 2
eval_wrap(ans)