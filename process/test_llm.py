from llm import *


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
    function_text = """
    def eqn_1_11__T_test(M: float, P: float, W: float, q: float) -> float:
        T = (q * P * 492) / (W * (359/M) * 760)
        return T
    """
    eval(function_text)
    eval(and1)
    assert eqn_1_11__T(1,2,3,4) == eqn_1_11__T_test(1,2,3,4)
