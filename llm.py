import ollama
import secrets
n='\n';null_str=""

def escribir_codigo(eqn, lang,p1=1,p2=0):
    print(f"solving { eqn }, buddy.")
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
                        f"""identify the line that is a one liner for a closed-form expression for T:\n{eqn}"""
                    )[1]
                )
                + (
                    p2 := (
                        """print the python code syntax for an equation solving the formula for every value ablated ( taken one at a time out)
                    example: a**2 = b**2 + c**2


                    winning correct output:  a = [math.sqrt(b**2+c**2),math.sqrt(b**2+c**2),]; b = [math.sqrt(a**2-c**2),math.sqrt(a**2-c**2),];c = [math.sqrt(a**2-b**2),math.sqrt(a**2-b**2),];""",
                        null_str,
                        """Give your answer as one line."""
                    )[1]
                ),
            },
        ],
    )
    dummy = lambda aa: print((n + "*" * 88 + n).join(["", aa, ""]))
    dummy(p1)
    dummy(p2)
    return n.join(filter(lambda a: a, response["message"]["content"].split(n)))


print(
    escribir_codigo(
        eqn="q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)", lang="Python"
    )
)
# def identifica_importante(headlines):
#     response = ollama.chat(
#         model="llama3:latest",
#         messages=[
#             {
#                 "role": "user",
#                 "content": f"Order these headlines by most important {headlines}. Strictly use one list, ordinal 1,2,3",
#             }
#         ],
#     )
#     # Parse response
#     i, importantes = 1, []
#     for line in response["message"]["content"].split('\n'):
#         if line.startswith(f"{i}. "):
#             importantes.append(line.replace(f"{i}. ", "").replace('"', '').replace("'", ""))
#             i += 1
#     return importantes
