"""Vakyume kwasak decorator: resolves missing-variable errors by algebraic rearrangement.

Given a class with solver methods named ``<eqn>__<variable>``, the
``@kwasak`` decorator routes a call with one missing keyword argument
to the appropriate solver automatically.

Created on Mon May  8 10:56:16 CDT 2023
@author: Julian Henry
"""

import functools
import inspect


def kwasak(func):
    """Dispatch ``func(**kw)`` to ``func__<missing>(**kw)``.

    Same rules as upstream juleshenry/kwasak: arguments are passed by name,
    ``x=None`` counts as missing, and unknown names are an error unless the
    stub declares ``**kwargs`` (then they pass through). A leading ``self``
    parameter is skipped, so static-style stubs work too.
    """
    params = list(inspect.signature(func).parameters.values())
    if params and params[0].name == "self":
        params = params[1:]
    variables = [
        p.name
        for p in params
        if p.kind not in (inspect.Parameter.VAR_POSITIONAL, inspect.Parameter.VAR_KEYWORD)
    ]
    extra_ok = any(p.kind is inspect.Parameter.VAR_KEYWORD for p in params)

    @functools.wraps(func)
    def wrapper(self, **kw):
        unknown = [k for k in kw if k not in variables]
        if unknown and not extra_ok:
            raise TypeError(f"{func.__name__}() got unknown variable(s): {', '.join(unknown)}")
        given = {k: v for k, v in kw.items() if v is not None}
        missing = [v for v in variables if v not in given]
        if len(missing) != 1:
            raise ValueError("Must have exactly one missing variable for which to solve.")
        return getattr(self, func.__name__ + "__" + missing[0])(**given)

    return wrapper
