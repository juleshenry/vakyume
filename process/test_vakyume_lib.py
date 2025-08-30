import shard.LiquidRing_10_10 as vakyume_lib
from inspect import signature


def make_test_case(o):
    cnt = 0
    stro = str(o)
    if stro[0].isalpha():
        # print('sublib',stro)
        sublib = getattr(vakyume_lib, stro)
        for a in dir(sublib):
            astro = str(a)
            if stro[0].isalpha():
                method = getattr(vars(vakyume_lib)[stro], astro)
                if "<class 'function'>" in str(type(method)) and "eqn" in str(method):
                    cnt += 1
                    # print(cnt)
                    test_args = [
                        1 + i * 0.5 for i in range(len(signature(method).parameters))
                    ]
                    try:
                        print("->" * 22)
                        print(f"{method}({test_args})")
                        j= method(*test_args)
                        print(j)
                    except Exception as e:
                        print("->" * 88)
                        print(e)


# TODO: jan 5 2025, solve issues of broken tests
if __name__ == "__main__":
    for s, o in enumerate(dir(vakyume_lib)):
        try:
            make_test_case(o)
        except Exception as e:
            print(e)