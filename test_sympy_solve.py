from sympy import Symbol, solve

def test():
    nf = "(x + y) - (z)"
    symb = Symbol('z')
    try:
        solns = solve(nf, symb)
        print(f"Solutions: {solns}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    test()
