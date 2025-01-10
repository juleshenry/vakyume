import vakyume_2025 as vakyume_lib
from inspect import signature
from typing import get_type_hints


def make_test_case(o):
    stro = str(o)
    if stro[0].isalpha():
        # print('sublib',stro)
        sublib = getattr(vakyume_lib, stro)
        for a in dir(sublib):
            astro = str(a)
            if stro[0].isalpha():
                x = vars(vakyume_lib)[stro]
                method = getattr(x, astro)
                if "<class 'function'>" in str(type(method)) and "eqn" in str(method):
                    test_args = [
                        i * 0.01 + 0.1 for i in range(len(signature(method).parameters))
                    ]
                    try:

                        method(*test_args)
                    except Exception as e:
                        print(e)
                        print(
                            f"{''.join(str(e).split(' ')[1:-2])}\n{method}({test_args})"
                        )


# TODO: jan 5 2025, solve issues of broken tests
if __name__ == "__main__":
    for o in dir(vakyume_lib):
        try:
            make_test_case(o)
        except Exception as e:
            print(e)

# 안
