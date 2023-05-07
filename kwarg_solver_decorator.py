def kwarg_solver(func):
    def solver_decorator():
        print("Something is happening before the function is called.")
        func()
        print("Something is happening after the function is called.")
    return wrapper

class Einstein:

    @kwarg_solver
    def einstein(e=None, m=None):
        
        return 

    def einstein_

    def einstein_e(m: float) -> float:
        return m * 8.98755179 e16

einstein(e = 1000) # returns m, (1000 / 8.98755179 e16), ~1.11265 e -14
einstein(m = 1000) # returns e, 1000 * 8.98755179 e16, ~8.98755179 e19