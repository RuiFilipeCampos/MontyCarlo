from functools import wraps

def timer(func):
    @wraps(func) #copies metadata, id, docstring etc
    def wrapper(*args, **kwargs):
        import time
        t0 = time.perf_counter()
        result = func(*args, **kwargs)
        tf = time.perf_counter()
        
        print(func.__name__, tf-t0)
        
        return result
    return wrapper
