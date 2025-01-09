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
                if "<class 'function'>" in str(type(method)):
                    test_args = [
                        i + 1 for i in range(len(signature(method).parameters))
                    ]
                    try:
                        # print(method, )
                        method(*test_args)
                    except Exception as e:
                        raise ValueError(f"{e}\nmethod failed {method}")


# TODO: jan 5 2025, solve issues of broken tests
if __name__ == "__main__":
    for o in dir(vakyume_lib):
        try:
            make_test_case(o)
        except Exception as e:
            print(e)
