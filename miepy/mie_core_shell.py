"""
mie_core_shell calculates the scattering coefficients of a core-shell structure using Mie theory
"""

import numpy as np
import miepy
from miepy.special_functions import riccati_1_single,riccati_2_single,riccati_3_single
from miepy.material_functions import constant_material
R1 = riccati_1_single
R2 = riccati_2_single
R3 = riccati_3_single

def M_matrix(m1,m2,x,y,mu,mu1,mu2,n):
    """M matrix in core_shell solver"""
    M = np.zeros([8,8,len(m1)], dtype=np.complex)
    z = np.zeros(len(m1))
    M[0] = np.array([z,z, -m2*R1(n,m1*x)[0],z, m1*R1(n,m2*x)[0],z, -m1*R3(n,m2*x)[0],z])
    M[1] = np.array([z,z,z, m2*R1(n,m1*x)[1],z, -m1*R1(n,m2*x)[1],z, m1*R3(n,m2*x)[1]])
    M[2] = np.array([z,z, mu2*R1(n,m1*x)[1],z, -mu1*R1(n,m2*x)[1],z, mu1*R3(n,m2*x)[1],z])
    M[3] = np.array([z,z,z, -mu2*R1(n,m1*x)[0],z, mu1*R1(n,m2*x)[0],z, -mu1*R3(n,m2*x)[0]])
    M[4] = np.array([-m2*R2(n,y)[1],z,z,z,z, -R1(n,m2*y)[1],z, R3(n,m2*y)[1]])
    M[5] = np.array([z, m2*R2(n,y)[0],z,z, R1(n,m2*y)[0],z, -R3(n,m2*y)[0],z])
    M[6] = np.array([-mu2*R2(n,y)[0],z,z,z,z, -mu*R1(n,m2*y)[0],z, mu*R3(n,m2*y)[0]])
    M[7] = np.array([z, mu2*R2(n,y)[1],z,z, mu*R1(n,m2*y)[1],z, -mu*R3(n,m2*y)[1],z])

    return np.transpose(M, (2,0,1))

def c_values(m2,y,mu2,n):
    """c array in core_shell solver"""
    z = np.zeros(len(m2))
    c = np.zeros([8,len(m2)], dtype=np.complex)
    c = np.array([z,z,z,z, -m2*R1(n,y)[1], m2*R1(n,y)[0], -mu2*R1(n,y)[0], mu2*R1(n,y)[1]])
    return np.transpose(c)

class single_mie_core_shell:
    def __init__(self, radius_in, radius_out, material_in, material_out, wavelength, Lmax, medium=None):
        """Solve traditional Mie theory: a single cores-shell in x-polarized plane wave illumination
               radius_in        core radius
               radius_out       core+shell radius
               material_in      core material
               material_out     shell material
               wavelength[N]    wavelength(s) to solve the system at
               Lmax             maximum number of orders to use in angular momentum expansion
               medium           material medium (must be non-absorbing; defaults to vacuum)
        """

        self.radius_in = radius_in
        self.radius_out = radius_out
        self.material_in = material_in
        self.material_out = material_out

        self.wavelength = np.asarray(np.atleast_1d(wavelength), dtype=float)
        self.Lmax = Lmax
        if medium is None:
            self.medium = constant_material(1.0, 1.0)
        else:
            self.medium = medium
            if (self.medium.eps(self.wavelength).imag != 0).any()  \
                         or (self.medium.mu(self.wavelength).imag != 0).any():
                raise ValueError('medium must be non-absorbing')

        self.Nfreq = len(self.wavelength)

        self.material_data = {}
        self.material_data['wavelength'] = self.wavelength
        self.material_data['eps_in']     = self.material_in.eps(self.wavelength)
        self.material_data['mu_in']      = self.material_in.mu(self.wavelength)
        self.material_data['n_in']       = np.sqrt(self.material_data['eps_in']*self.material_data['mu_in'])
        self.material_data['eps_out']    = self.material_out.eps(self.wavelength)
        self.material_data['mu_out']     = self.material_out.mu(self.wavelength)
        self.material_data['n_out']      = np.sqrt(self.material_data['eps_out']*self.material_data['mu_out'])
        self.material_data['eps_b']      = self.medium.eps(self.wavelength)
        self.material_data['mu_b']       = self.medium.mu(self.wavelength)
        self.material_data['n_b']        = np.sqrt(self.material_data['eps_b']*self.material_data['mu_b'])
        self.material_data['k']          = 2*np.pi*self.material_data['n_b']/self.wavelength
               
        self.an = np.zeros((self.Nfreq, self.Lmax), dtype=np.complex)
        self.bn = np.zeros((self.Nfreq, self.Lmax), dtype=np.complex)
        self.cn = np.zeros((self.Nfreq, self.Lmax), dtype=np.complex)
        self.dn = np.zeros((self.Nfreq, self.Lmax), dtype=np.complex)

        self.scattering_properties = (self.an, self.bn, self.material_data['k'])

        self.computed = False
    
    def solve(self):
        """solve the system"""
        mat = self.material_data

        m1 = mat['n_in']/mat['n_b']
        m2 = mat['n_out']/mat['n_b']
        xvals = mat['k']*self.radius_in
        yvals = mat['k']*self.radius_out

        for n in range(self.Lmax):
            M = M_matrix(m1, m2, xvals, yvals, mat['mu_b'], mat['mu_in'], mat['mu_out'], n+1)
            c = c_values(m2, yvals, mat['mu_out'], n+1)
            sol = np.linalg.solve(M,c)
            self.an[:,n] = sol[:,0]
            self.bn[:,n] = sol[:,1]
            self.cn[:,n] = sol[:,2]
            self.dn[:,n] = sol[:,3]
        
        self.an = np.nan_to_num(self.an)
        self.bn = np.nan_to_num(self.bn)

        self.computed = True
        return self.an, self.bn

    def cross_sections(self):
        """Return the 3 cross-sections: (Scattering, Absorbption, Extinction)"""
        if not self.computed: self.solve()
        return miepy.cross_sections(*self.scattering_properties)

    def cross_sections_per_multipole(self):
        """Return the 3 cross-sections per multipole: (Scattering, Absorbption, Extinction)"""
        if not self.computed: self.solve()
        return (miepy.scattering_per_multipole(*self.scattering_properties),
                miepy.absorbption_per_multipole(*self.scattering_properties),
                miepy.extinction_per_multipole(*self.scattering_properties))