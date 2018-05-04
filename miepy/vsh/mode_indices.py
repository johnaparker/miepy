"""
Functions related to the (m,n) vsh mode indices
"""

def rmax_to_Lmax(rmax):
    """obtain Lmax from rmax"""
    Lmax = int(-1 + (1+rmax)**0.5)
    return Lmax

def Lmax_to_rmax(Lmax):
    """obtain rmax from Lmax"""
    rmax = Lmax*(Lmax + 2)
    return rmax

class mode_indices:
    """iterable to produce (r, n, m) mode indices in the order expected by MiePy"""

    def __init__(self, *args, **kwargs):
        """
        mode_indices(Lmax)
        mode_indices(Lmin, Lmax, m_start=m_start, m_stop=m_stop)

        Returns an object that produces a sequence of (r, n, m) indices with the following ranges:

            n: [Lmin:Lmax+1]     (default: Lmin=1)
            m: [m_start:m_stop]  (default: [-n:n+1])

        r is a counting index useful for indexing expansion coefficients of shape [2,rmax]
        """

        if len(args) == 1:
            n_start, n_stop = 1, args[0]
        elif len(args) == 2:
            n_start, n_stop = args
        else:
            raise TypeError(f"mode_indices expected at most 2 positional arguments, got {len(args)}")

        try:
            n_start, n_stop = int(n_start), int(n_stop)
        except ValueError:
            raise TypeError(f'integer arguments are required')

        if n_start < 1:
            raise ValueError('n_start must be greater than 1')

        self.n_start = n_start
        self.n_stop = n_stop
        self.m_start = kwargs.get('m_start', None)
        self.m_stop = kwargs.get('m_stop', None)
        self.idx_start = Lmax_to_rmax(n_start-1)

    def __repr__(self):
        if self.n_start == 1 and self.m_start is None and self.m_stop is None:
            return f'mode_indices(Lmax={self.n_stop})'
        elif self.m_start is None and self.m_stop is None:
            return f'mode_indices(Lmin={self.n_start}, Lmax={self.n_stop})'
        else:
            return f'mode_indices(Lmin={self.n_start}, Lmax={self.n_stop}, Mmin={self.m_start}, Mmax={self.m_stop})'

    def __iter__(self):
        idx = self.idx_start

        if self.m_start is None and self.m_stop is None:
            for n in range(self.n_start, self.n_stop+1):
                for m in range(-n,n+1):
                    yield idx, n, m
                    idx += 1
        else:
            for n in range(self.n_start, self.n_stop+1):
                m_start = -n if (self.m_start is None or self.m_start < -n) else self.m_start
                m_stop = n if (self.m_stop is None or self.m_stop > n) else self.m_stop
                idx = Lmax_to_rmax(n-1) + m_start + n
                for m in range(m_start, m_stop+1):
                    yield idx, n, m
                    idx += 1

    def __eq__(self, other):
        return isinstance(other, mode_indices) and \
               self.n_start == other.n_start   and \
               self.n_stop  == other.n_stop    and \
               self.m_start == other.m_start   and \
               self.m_stop  == other.m_stop

    def __len__(self):
        return self.n_stop - self.n_start + 1
