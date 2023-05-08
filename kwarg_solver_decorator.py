# -*- coding: utf-8 -*-
"""
Created on Mon May  8 10:56:16 CDT 2023

@author: Julian Henry

Allows for arbitrary missing variable resolution
"""
def kwarg_solver(func):
    """
    Intention: 
    base method is `vanilla`
    To calculate necessary 
    
    formula to solve for a single variable 
    """
    def wrapper(self, *args, **kwargs):
        method_name = func.__name__
        missing_arg = list(filter(lambda x:x, kwargs))
        if not len(missing_arg) == 1:
            raise IllegalArgumentException
        correct_method = method_name + '__' + missing_arg[0] 
        correct_args = [x[1] for x in sorted(kwargs.items(), key=lambda kv:kv[0])]
        return getattr(self ,correct_method)(*correct_args)
    return wrapper


class Einstein:
    @kwarg_solver
    def einstein_1(s, **kwargs):
        return # decorator skips return

    def einstein_1__m(s, e:float):
        return e / 8.98755179e16


    def einstein_1__e(s, m: float) -> float:
        return m * 8.98755179e16

e = Einstein()
ans = e.einstein_1(e = 1000) # returns m, (1000 / 8.98755179 e16), ~1.11265 e -14
print(ans)
ans = e.einstein_1(m = 1000) # returns e, 1000 * 8.98755179 e16, ~8.98755179 e19
print(ans)


class Pythagoras:
    @kwarg_solver
    def pythagorean(s, **kwargs):
        return

    def pythagorean__a(s, b: float, c:float):
        return (c**2 - b**2) **.5
    def pythagorean__b(s, a: float, c:float):
        return (c**2 - a**2) **.5
    def pythagorean__c(s, a: float, b:float):
        return (a**2 + b**2) **.5


