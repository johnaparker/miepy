from timeit import default_timer as timer


def time_function(func, runtime=0.1):
    """Time a function by running it repeatedly for at least 'runtime' seconds."""
    timer()
    t = 0
    count = 0

    while t < runtime:
        t0 = timer()
        func()
        tf = timer()
        t += tf - t0

        count += 1

    return t / count
