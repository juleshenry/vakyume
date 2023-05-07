

solver_decorator

def my_decorator(func):
    def solver_decorator():
        print("Something is happening before the function is called.")
        func()
        print("Something is happening after the function is called.")
    return wrapper