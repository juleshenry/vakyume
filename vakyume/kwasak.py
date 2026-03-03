# -*- coding: utf-8 -*-
"""
kwasak -- arbitrary missing-variable resolution decorator.

Created on Mon May  8 10:56:16 CDT 2023
@author: Julian Henry

Given a class with solver methods named ``<eqn>__<variable>``, the
``@kwasak`` decorator routes a call with one missing keyword argument
to the appropriate solver automatically.
"""

import inspect


def kwasak(func):
    def wrapper(self, **kw):
        params = list(inspect.signature(func).parameters)
        # Skip 'self' so instance methods work correctly, and skip 'kwargs' if it exists so that it can be used to pass through extra variables to the underlying function
        params = params[params[0] == "self" : len(params) - (params[-1] == "kwargs")]
        if len(m := [w for w in params if w not in kw]) - 1:
            # Apparently, there is a bug in inspect.siganature mro() causing `X_1` to be preferred to match `X_10` when calling func
            raise ValueError(
                "Must have exactly one missing variable for which to solve."
            )
        return getattr(self, func.__name__ + "__" + m[0])(**kw)

    return wrapper
