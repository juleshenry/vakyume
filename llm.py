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
        ],
    )

    dummy(p1)
    dummy(p2)
    return n.join(filter(lambda a: a, response["message"]["content"].split(n)))


# globe = 0
# def l_strip(s):


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


import re


def extract_code(text):
    dummy([text])
    # pattern = re.compile(r'```\s*python\s*(.*?)\s*```', re.DOTALL)
    # pattern = re.compile(r'```python\n(.*?)```', re.DOTALL)
    # print(pattern)
    # match = pattern.match(text)
    # if match:
    #     return match.group(1).strip()
    # else:
    #     return None
    maxxx = 0
    maxxx_str = ""
    max_str = ""
    for i, ii in enumerate(t := text.split("\n")):
        # print(i,ii)
        if "```python" in ii:
            max_str += ii + n
        elif "```" in ii:
            max_str += ii + n
            if len(max_str.split(n)) > maxxx:
                maxxx = len(max_str.split(n))
                maxxx_str = max_str
                max_str = ""
        if (a := len(max_str.split(n))) > 1:
            # print('->>>',a,maxxx)
            max_str += ii + n
    # print(maxxx_str)
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
test_ans1 = """
P = 30  # Replace with your value
q = 40  # Replace with your value
T = (492 * q) / (W * 5.99 / M * P)
print(T)
```
You can replace the values of W, M, P, and q with any values you want to solve for T.
If you want to automate this process for all possible combinations of values, you'll need to use a loop or a library like NumPy or SciPy. Here's an example using a simple loop:
```python
import math
# Define the variables and their ranges
Ws = [1, 2, 3]  # Replace with your range
Ms = [10, 20, 30]  # Replace with your range
Ps = [5, 15, 25]  # Replace with your range
qs = [50, 75, 100]  # Replace with your range
for W in Ws:
    for M in Ms:
        for P in Ps:
            for q in qs:
                T = (492 * q) / (W * 5.99 / M * P)
                print(f"W={W}, M={M}, P={P}, q={q}, T={T}")
```
This code will iterate over all possible combinations of values for W, M, P, and q, and print the corresponding value of T.
Note that you'll need to adjust the ranges and values of the variables to match your specific proble
"""
print(extract_code(test_ans1))
# print(eval_wrap(test_ans1))
# print(ans2 := escribir_codigo(eqn=ans1, lang=None, p1_i=2, p2_i=2))  # or 2
# eval_wrap(ans2)
