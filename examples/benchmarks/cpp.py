from functools import partial

import numpy as np
from termcolor import colored
from timer import time_function

import miepy


def display(name, Tcpp, Tpy):
    """Display name and times of C++/python runtimes."""
    print(colored(name, color="white", attrs=["underline"]))
    print(f"    cpp: {Tcpp * 1e6:>7.2f} μs")
    print(f"    py:  {Tpy * 1e6:>7.2f} μs")

    if Tcpp < Tpy:
        print(colored(f"C++   {Tpy / Tcpp:.2f}x", color="red", attrs=["bold"]))
    else:
        print(colored(f"Python   {Tcpp / Tpy:.2f}x", color="green", attrs=["bold"]))
    print()


def associated_legendre():
    n = 2
    m = 1
    x = np.linspace(-1, 1, 1000)

    f = partial(miepy.cpp.special.associated_legendre, n, m, x)
    Tcpp = time_function(f)

    f = partial(miepy.vsh.old_special.associated_legendre(n, m), x)
    Tpy = time_function(f)

    display(f"Associated Legendre (n = {n}, m = {m})", Tcpp, Tpy)


def hankel():
    n = 4
    x = np.linspace(0.3, 0.5, 1000)

    f = partial(miepy.cpp.special.spherical_hn, n, x)
    Tcpp = time_function(f)

    f = partial(miepy.vsh.old_special.spherical_hn, n, x)
    Tpy = time_function(f)

    display(f"Spherical Hankel, small x (n = {n})", Tcpp, Tpy)

    x = np.linspace(10.0, 12.0, 1000)

    f = partial(miepy.cpp.special.spherical_hn, n, x)
    Tcpp = time_function(f)

    f = partial(miepy.vsh.old_special.spherical_hn, n, x)
    Tpy = time_function(f)

    display(f"Spherical Hankel, large x (n = {n})", Tcpp, Tpy)


associated_legendre()
hankel()
