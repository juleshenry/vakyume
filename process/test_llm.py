from llm import *


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


def do_extract_code_1():
    # EVAL WRAP DOES NOT WORK.
    print(
        ans1 := escribir_codigo(
            eqn="q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)",
            lang="Python",
            p1_i=1 + 1,
            p2_i=1,
        )
    )
    print("###" * 19)
    print(extract_code(ans1))
    u = extract_code(test_ans1)
    # print(ans2 := escribir_codigo(eqn=ans1, lang=None, p1_i=2, p2_i=2))  # or 2
    # print(eval_wrap(ans2))
    print(u)


def do_extract_code_2():
    # Tue Jan 21 17:03:09 CST 2025
    print(
        ans1 := escribir_codigo(
            eqn="q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)",
            single_variable="T",
            p1_i=4,
            p2_i=0,
        )
    )

    print(
        ans1 := escribir_codigo(
            eqn="",
            single_variable="T",
            header="eqn_1_11__T(M: float, P: float, W: float, q: float):",
            pipin=ans1,
            p1_i=-1,
            p2_i=3,
        )
    )
    print(and1 := extract_code(ans1))
    # .below may be necessary ymmv
    # print(ans1:=(make_sure_python_annotated(ans1)))
    # print(and1:=extract_code(ans1))

do_extract_code_1()