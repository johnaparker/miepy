"""
Functions related to the (m,n) vsh mode indices
"""

def rmax_to_lmax(rmax):
    """obtain lmax from rmax"""
    lmax = int(-1 + (1+rmax)**0.5)
    return lmax

def lmax_to_rmax(lmax):
    """obtain rmax from lmax"""
    rmax = lmax*(lmax + 2)
    return rmax

def reduced_index(n, m):
    r = n*(n+2) - n - m - 1;
    return r

class mode_indices:
    """iterable to produce (r, n, m) mode indices in the order expected by MiePy"""

    def __init__(self, *args, **kwargs):
        """
        mode_indices(lmax)
        mode_indices(Lmin, lmax, m_start=m_start, m_stop=m_stop)

        Returns an object that produces a sequence of (r, n, m) indices with the following ranges:

            n: [Lmin:lmax+1]     (default: Lmin=1)
            m: [m_start:m_stop]  (default: [-n:n+1])

        r is a counting index useful for indexing expansion coefficients of shape [2,rmax]
        """

        if len(args) == 1:
            n_start, n_stop = 1, args[0]
        elif len(args) == 2:
            n_start, n_stop = args
        else:
            raise TypeError("mode_indices expected at most 2 positional arguments, got {}".format(len(args)))

        try:
            n_start, n_stop = int(n_start), int(n_stop)
        except ValueError:
            raise TypeError('integer arguments are required')

        if n_start < 1:
            raise ValueError('n_start must be greater than 1')

        self.n_start = n_start
        self.n_stop = n_stop
        self.m_start = kwargs.get('m_start', None)
        self.m_stop = kwargs.get('m_stop', None)
        self.idx_start = lmax_to_rmax(n_start-1)

    def __repr__(self):
        if self.n_start == 1 and self.m_start is None and self.m_stop is None:
            return 'mode_indices(lmax={})'.format(self.n_stop)
        elif self.m_start is None and self.m_stop is None:
            return 'mode_indices(Lmin={}, lmax={})'.format(self.n_start, self.n_stop)
        else:
            return 'mode_indices(Lmin={}, lmax={}, Mmin={}, Mmax={})'.format(
                       self.n_start, self.n_stop, self.m_start, self.m_stop)

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
                idx = lmax_to_rmax(n-1) + m_start + n
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
