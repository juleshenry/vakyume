import ollama
import secrets

n = "\n"
null_str = ""
dummy = lambda aa: print((n + "*" * 88 + n).join(["", *aa, ""]))


def escribir_codigo(eqn, lang, single_variable=None, p1_i=None, p2_i=None):
    int(
        sum(
            [
                p1_i,
                p2_i,
            ]
        )
    )
    print(f"solving { eqn }, queerly enough.")
    fluff = secrets.choice(["C++", "Erlang", "Haskell", "Rust", "Lisp"])
    mnms = [
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
                        Hi Doctor Math, I need help solving this {eqn}  Give back the closed form expression for {single_variable or 'T'}
                        
                        """,
                    f"""identify the line that is a one liner for a closed-form expression for {single_variable or 'T'}:\n{eqn}""",
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
            )
            + (
                p3 := (
                    "give back only the methods",
                    "no additional commentary",
                )[0]
                + p3[1]
            ),
        },
    ]
    print(mnms)
    1 / 0
    response = ollama.chat(
        model="llama3:latest",
        messages=mnms,
    )

    dummy(p1)
    dummy(p2)
    return n.join(filter(lambda a: a, response["message"]["content"].split(n)))




def llm_is_true(o):
    # print(o)
    response = ollama.chat(
        model="llama3:latest",
        messages=[
            {
                "role": "user",
                "content": [
                    f"Return only back the Python format code, signified by ```python. no comments or formatting {o}",
                    "Give the python code to extract text from string only from '```python' to nearest '```' ",
                ][1],
            }
        ],
    )
    return response["message"]["content"]



def extract_code(text):
    # dummy([text])
    maxxx = 0
    maxxx_str = ""
    max_str = ""
    for i, ii in enumerate(t := text.split("\n")):
        # print(i,ii)
        if "```python" in ii:
            max_str += n
        elif "```" in ii:
            if len(max_str.split(n)) > maxxx:
                maxxx, maxxx_str, max_str = len(max_str.split(n)), max_str, ""
                #  = 
                # max_str = ""
        if (a := len(max_str.split(n))) > 1:
            max_str += (ii if a - 2 else '') + n
    return maxxx_str


def eval_wrap(ans):
    max_bounds = [0, 0]
    print("evaluatin", len(ans), max_bounds)
    for i, j in (s := enumerate(o := ans.split(n))):
        print("$", j)
        for ii, jj in enumerate(o):
            try:
                # print('CODE'*8+'\n','\n'.join(o[i:ii]))
                f = "\n".join(map(lambda l: l.lstrip(), o[i:ii]))
                # dummy([f])
                # TODO: lmao; you have to eval with the variables already pre-initialized.
                # Your concept of "pyeqn formula applies"; or just reuse normal form code written? lmao
                # dummy(extract_code(f))
                if ii - i > max_bounds[1] - max_bounds[0]:
                    max_bounds = [i, ii]
                    print(f)
            except:
                pass
            # O(n^2)
    print("$#" * 88)
    return n.join(o[max_bounds[0] : max_bounds[1]])

# print(
#     ans1 := escribir_codigo(
#         eqn="q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)",
#         lang="Python",
#         p1_i=1,
#         p2_i=1,
#     )
# )
# print(eval_wrap(test_ans1))
# print(ans2 := escribir_codigo(eqn=ans1, lang=None, p1_i=2, p2_i=2))  # or 2
# eval_wrap(ans2)
