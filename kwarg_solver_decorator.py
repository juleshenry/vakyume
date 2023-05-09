# -*- coding: utf-8 -*-
"""
Created on Mon May  8 10:56:16 CDT 2023

@author: Julian Henry

Allows for arbitrary missing variable resolution
"""
import inspect


def kwarg_solver(func):
    """
    Intention:
    base method is `vanilla`... that means no `__syntax`

    and can calculate necessary formula to solve for a single variable

    Assumption: method is not static. why? we abuse inspect to inspect the
    [1:-1] kwargs in order to resolve the missing parameter
    """

    def wrapper(self, *args, **kwargs):
        method_name = func.__name__
        all_parameters = list(inspect.signature(func).parameters)[1:-1]
        missing_arg = [kw for kw in all_parameters if kw not in kwargs]
        if not len(missing_arg) == 1:
            raise ValueError(
                "Must have exactly one missing variable for which to solve."
            )
        correct_method = method_name + "__" + missing_arg[0]
        correct_args = [x[1] for x in sorted(kwargs.items(), key=lambda kv: kv[0])]
        return getattr(self, correct_method)(*correct_args)

    return wrapper


class Einstein:
    @kwarg_solver
    def einstein(s, e: float = None, m: float = None, **kwargs):
        return  # decorator skips return

    def einstein__m(s, e: float):
        return e / 8.98755179e16

    def einstein__e(s, m: float) -> float:
        return m * 8.98755179e16


e = Einstein()
ans = e.einstein(e=1000)  # returns m, (1000 / 8.98755179 e16), ~1.11265 e -14
print(ans)
ans = e.einstein(m=1000)  # returns e, 1000 * 8.98755179 e16, ~8.98755179 e19
print(ans)


class Pythagoras:
    @kwarg_solver
    def pythagorean(s, a: float = None, b: float = None, c: float = None, **kwargs):
        return

    def pythagorean__a(s, b: float, c: float):
        return (c ** 2 - b ** 2) ** 0.5

    def pythagorean__b(s, a: float, c: float):
        return (c ** 2 - a ** 2) ** 0.5

    def pythagorean__c(s, a: float, b: float):
        return (a ** 2 + b ** 2) ** 0.5


p = Pythagoras()

ans = p.pythagorean(a=12, c=13)  # 5
print(ans)
ans = p.pythagorean(a=12, b=5)  # 12
print(ans)
ans = p.pythagorean(b=12, c=13)  # 13
print(ans)
