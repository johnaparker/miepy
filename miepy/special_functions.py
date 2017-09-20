"""
special_functions defines any additional special functions needed by MiePy
"""

import numpy as np
from scipy import special

def spherical_hn(n, z, derivative=False):
    return special.spherical_jn(n,z,derivative) + 1j*special.spherical_yn(n,z,derivative)

def riccati_1(nmax,x):
    """Riccati bessel function of the 1st kind

       returns (r1, r1'), n=0,1,...,nmax"""

    jn,jnp = special.sph_jn(nmax,x)

    r0 = x*jn
    r1 = jn + x*jnp
    return np.array([r0,r1])

def riccati_2(nmax,x):
    """Riccati bessel function of the 2nd kind

       returns (r2, r2'), n=0,1,...,nmax"""

    jn,jnp,yn,ynp = special.sph_jnyn(nmax,x)
    hn = jn + 1j*yn
    hnp = jnp + 1j*ynp

    r0 = x*hn
    r1 = hn + x*hnp
    return np.array([r0,r1])

def riccati_3(nmax,x):
    """Riccati bessel function of the 3rd kind

       returns (r3, r3'), n=0,1,...,nmax"""

    yn,ynp = special.sph_yn(nmax,x)

    r0 = x*yn
    r1 = yn + x*ynp
    return np.array([r0,r1])

def pi_tau_func(n):
    # if np.sin(theta) == 0: return 0
    lpn = special.legendre(n)
    lpn_p = lpn.deriv()
    lpn_p2 = lpn_p.deriv()

    def pi_func(theta):
        return -1*lpn_p(np.cos(theta))
        # with np.errstate(divide='ignore', invalid='ignore'):
            # val = lpn(np.cos(theta))/np.sin(theta)
            # val[val == np.inf] = 0
            # val = np.nan_to_num(val)
            # return val

    def tau_func(theta):
        # val = -1*np.sin(theta)*lpn_p(np.cos(theta))
        val = -1*np.cos(theta)*lpn_p(np.cos(theta)) + np.sin(theta)**2*lpn_p2(np.cos(theta))
        return val

    return pi_func, tau_func 

class vector_spherical_harmonics:
    def __init__(self, n, superscript=3):
        self.pi_func, self.tau_func = pi_tau_func(n)
        self.n = n

        if superscript == 1:
            self.z_func = lambda x: special.spherical_jn(n,x)
            self.zp_func = lambda x: special.spherical_jn(n,x, derivative=True)
        elif superscript == 3:
            self.z_func = lambda x: spherical_hn(n,x)
            self.zp_func = lambda x: spherical_hn(n,x, derivative=True)

    def M_o1n(self, k):
        def f(r, theta, phi):
            theta_comp = np.cos(phi)*self.pi_func(theta)*self.z_func(k*r)
            phi_comp = -1*np.sin(phi)*self.tau_func(theta)*self.z_func(k*r)
            r_comp = np.zeros(shape = theta.shape, dtype=np.complex)
            return np.array([r_comp, theta_comp, phi_comp])
        return f

    def M_e1n(self, k):
        def f(r, theta, phi):
            theta_comp = -1*np.sin(phi)*self.pi_func(theta)*self.z_func(k*r)
            phi_comp = -1*np.cos(phi)*self.tau_func(theta)*self.z_func(k*r)
            r_comp = np.zeros(shape = theta.shape, dtype=np.complex)
            return np.array([r_comp, theta_comp, phi_comp])
        return f

    def N_o1n(self, k):
        def f(r, theta, phi):
            p = k*r
            theta_comp = np.sin(phi)*self.tau_func(theta)*(self.z_func(p) + p*self.zp_func(p))/p
            phi_comp = np.cos(phi)*self.pi_func(theta)*(self.z_func(p) + p*self.zp_func(p))/p
            r_comp = np.sin(phi)*self.n*(self.n+1)*np.sin(theta)*self.pi_func(theta)*self.z_func(p)/p
            return np.array([r_comp, theta_comp, phi_comp])
        return f

    def N_e1n(self, k):
        def f(r, theta, phi):
            p = k*r
            theta_comp = np.cos(phi)*self.tau_func(theta)*(self.z_func(p) + p*self.zp_func(p))/p
            phi_comp = -1*np.sin(phi)*self.pi_func(theta)*(self.z_func(p) + p*self.zp_func(p))/p
            r_comp = np.cos(phi)*self.n*(self.n+1)*np.sin(theta)*self.pi_func(theta)*self.z_func(p)/p
            return np.array([r_comp, theta_comp, phi_comp])
        return f

def riccati_1_single(n,x):
    """Riccati_1, but only a single n value"""
    pre = (np.pi*x/2)**.5
    jn = pre*special.jv(n+0.5,x)
    jnp = jn/(2*x) + pre*special.jvp(n+0.5,x)

    return np.array([jn,jnp])

def riccati_2_single(n,x):
    """Riccati_2, but only a single n value"""
    pre = (np.pi*x/2)**.5
    hn = pre*special.hankel1(n+0.5,x)
    hnp = hn/(2*x) + pre*special.h1vp(n+0.5,x)

    return np.array([hn,hnp])

def riccati_3_single(n,x):
    """Riccati_3, but only a single n value"""
    # pre = (np.pi*x/2)**.5
    # yn = pre*special.yv(n+0.5,x)
    # ynp = yn/(2*x) + pre*special.yvp(n+0.5,x)

    # return np.array([yn,ynp])
    return riccati_2_single(n,x) - riccati_1_single(n,x)


if __name__ == "__main__":
    VCS = vector_spherical_harmonics(3, superscript=3)
    dipole = VCS.N_e1n(1)

    eps = 0.01
    theta = np.linspace(eps, np.pi-eps,60)
    phi = np.linspace(eps, 2*np.pi-eps,90)
    r = np.array([500])

    R, THETA, PHI = np.meshgrid(r, theta, phi, indexing='xy')
    print(R.shape)
    print(THETA.shape)
    print(PHI.shape)
    E_far = dipole(R, THETA, PHI)
    E_far = np.squeeze(E_far)
    I = np.sum(np.abs(E_far)**2, axis=0)

    import matplotlib.pyplot as plt
    plt.figure(1)
    plt.pcolormesh(np.squeeze(THETA),np.squeeze(PHI),I)
    plt.xlabel('theta')
    plt.ylabel('phi')
    plt.colorbar()

    plt.figure(2)
    plt.plot(theta, VCS.pi_func(theta))

    plt.show()
