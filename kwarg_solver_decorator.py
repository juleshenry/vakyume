import inspect

def kwarg_solver(func):
    def wrapper(*args, **kwargs):
        print("Something is happening before the function is called.")
        func(*args, **kwargs)
        print("Something is happening after the function is called.")
    return wrapper

class Einstein:

    def einstein_1(s,**kwargs):
        methname = inspect.stack()[0][3]
        missing_arg = list(filter(lambda x:x,kwargs))
        if not len(missing_arg) == 1:
            raise IllegalArgumentException
        correct_method = methname + '__'+missing_arg[0] 
        for self_attr in dir(s):
            if self_attr == correct_method:
                args = [x[1] for x in sorted(kwargs.items(),key=lambda kv:kv[0])]
                return getattr(s,correct_method)(*args)
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