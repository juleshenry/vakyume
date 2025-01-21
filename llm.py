import ollama
import secrets

n = "\n"
null_str = ""
dummy = lambda aa: print((n + "*" * 88 + n).join(["", *aa, ""]))


def escribir_codigo(eqn: str, lang:str="Python", single_variable=None, header = None, pipin="", p1_i=None, p2_i=None, p3_i=0):
    """Walrus-routed method for writing code with prompts via ollama"""
    int(
        sum(
            [
                p1_i,
                p2_i,
            ]
        )
    )
    print(f"solving { eqn }, queerly enough.")
    print(f"ROUTING: {(p1_i,p2_i,p3_i,)}")
    fluff = secrets.choice(["C++", "Erlang", "Haskell", "Rust", "Lisp"])
    mnms = [
        {
            "role": "user",
            "content": (
                p1 := (
                    null_str,
                    f"""
                        You are a {lang} / {fluff} wizard giga chad developer who wakes up peer-reviewing arxiv.org
                          for News.YCombinator.com threads and sleeps coding {lang} at 3am.
                        Solve this equation: {eqn} in {lang}
                        """,
                    f"""
                        Hi Doctor Math, I need help solving this {eqn}  Give back the closed form expression for {single_variable or 'T'}
                        
                        """,
                    f"""identify the line that is a one liner for a closed-form expression for {single_variable or 'T'}:{n+eqn}""",
                    f"""{eqn};solve for {single_variable}""",
                    pipin + n
                )[p1_i]
            ) 
            + (
                p2 := (
                    null_str,
                    """print the python code syntax for an equation solving the formula the given variable
                    example: a**2 = b**2 + c**2
                    
                    closed-form for "a"
                    (mark and start your code with '```python')

                    winning correct output:
                    ```python
                    def eqn_9_9__a(b, c):
                        return = [math.sqrt(b**2+c**2),math.sqrt(-(b**2+c**2)),]
                    ```
                    {n}""",
                    """Give your answer as one line.""",
                    f"""Hello, You will write a specific Python function to return {single_variable}; Make sure variables are in terms of {"" if header else header.split('(')[1].replace(')','')}\n The method header will be:{n}{header}. """
                )[p2_i]
            )
            + (
                p3 := (
                    null_str,
                    "give back only the methods without additional commentary",
                    f"the equation header will be {header}",
                )
               
            )[p3_i],
        },
    ]
    dummy([mnms[0]['content']])
    response = ollama.chat(
        model="phi3:latest" or "llama3:latest",
        messages=mnms,
    )
    return n.join(filter(lambda a: a, response["message"]["content"].split(n)))




def make_sure_python_annotated(code_block):
    print("makin' sure iz annotated" + '...'*8)
    response = ollama.chat(
        model="llama3:latest",
        messages=[
            {
                "role": "user",
                "content": [
                    f"(Return the following code in Python, converting where necessary. \n \
                       (Make sure to annotate your answer with ```python and ```\n",
                    f"{code_block}",
                ][1],
            }
        ],
    )
    return response["message"]["content"]



def extract_code(text):
    """ gets actual code snippet (attempts)"""
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
                maxxx, maxxx_str, max_str = (len(max_str.split(n)), max_str, "",)
        if (a := len(max_str.split(n))) > 1:
            max_str += (ii if a - 2 else '') + n
    print("$#" * 32 + "extract_code" + "#$" * 32)
    return maxxx_str


def eval_wrap(ans):
    # EVAL WRAP DOES NOT WORK.  b/c you would have to know variables before
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
                    # print(f)
            except:
                pass
            # O(n^2)
    return n.join(o[max_bounds[0] : max_bounds[1]])

