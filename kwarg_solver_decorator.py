def kwarg_solver(func):
    def wrapper(self, *args, **kwargs):
        print(kwargs)
        method_name = func.__name__
        missing_arg = list(filter(lambda x:x, kwargs))
        if not len(missing_arg) == 1:
            raise IllegalArgumentException
        correct_method = method_name + '__' + missing_arg[0] 
        print(correct_method)
        correct_args = [x[1] for x in sorted(kwargs.items(), key=lambda kv:kv[0])]
        return getattr(self ,correct_method)(*correct_args)
    return wrapper

class Einstein:

    @kwarg_solver
    def einstein_1(s, **kwargs):
        return 0

    def einstein_1__m(s, e:float):
        return e / 8.98755179e16


    def einstein_1__e(s, m: float) -> float:
        return m * 8.98755179e16

e = Einstein()
ans = e.einstein_1(e = 1000) # returns m, (1000 / 8.98755179 e16), ~1.11265 e -14
print(ans)
ans = e.einstein_1(m = 1000) # returns e, 1000 * 8.98755179 e16, ~8.98755179 e19
print(ans)

# class A:
#     def o(s):
#         for k in dir(s):
#             print(k)
# A().o()