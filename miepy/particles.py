"""
Particle class for reprententing more general excitations
"""

import numpy as np
from miepy.mie_sphere import single_mie_sphere
from mpl_toolkits.mplot3d import Axes3D
from my_pytools.my_numpy.integrate import simps_2d
from my_pytools.my_numpy.indices import levi_civita

levi = levi_civita()
Ntheta = 51
Nphi = 31

BUFFER_DEFAULT = 1

class particle:
    """Sphere with a position"""

    def __init__(self, sp, center, freq_index, source):
        self.center = np.asarray(center)
        self.radius = sp.r
        self.solution = sp
        self.k = sp.mat.k[freq_index]*(sp.eps_b*sp.mu_b)**0.5
        self.E_func = sp.E_field(freq_index)
        self.H_func = sp.H_field(freq_index)
        self.source = source
        self.amp = source.E(center, self.k)

    def E(self, X, Y, Z, inc=True):
        Xs = X - self.center[0]
        Ys = Y - self.center[1]
        Zs = Z - self.center[2]

        R = np.sqrt(Xs**2 + Ys**2 + Zs**2)
        THETA = np.arccos(Zs/R)
        PHI = np.arctan2(Ys,Xs)
        # PHI[PHI<0] += 2*np.pi

        # incident field
        xhat = np.array([np.sin(THETA)*np.cos(PHI), np.cos(THETA)*np.cos(PHI), -np.sin(PHI)])
        yhat = np.array([np.sin(THETA)*np.sin(PHI), np.cos(THETA)*np.sin(PHI), np.cos(PHI)])

        rhat = np.array([np.sin(THETA)*np.cos(PHI), np.sin(THETA)*np.sin(PHI), np.cos(THETA)])
        that = np.array([np.cos(THETA)*np.cos(PHI), np.cos(THETA)*np.sin(PHI), -np.sin(THETA)])
        phat = np.array([-np.sin(PHI), np.cos(PHI), np.zeros_like(THETA)])

        if inc:
            Einc = self.source.E(np.array([Xs,Ys,Zs]), self.k)
        else:
            Einc = 0
        Ax, Ay = self.amp[:2]

        Escat = Ax*self.E_func(R,THETA,PHI) + Ay*self.E_func(R,THETA,PHI-np.pi/2)
        
        # convert to cartesian
        Etot = Escat[0]*rhat + Escat[1]*that + Escat[2]*phat + Einc
        return Etot

    def H(self, X, Y, Z, inc=True):
        Xs = X - self.center[0]
        Ys = Y - self.center[1]
        Zs = Z - self.center[2]

        R = np.sqrt(Xs**2 + Ys**2 + Zs**2)
        THETA = np.arccos(Zs/R)
        PHI = np.arctan2(Ys,Xs)
        # PHI[PHI<0] += 2*np.pi

        # incident field
        xhat = np.array([np.sin(THETA)*np.cos(PHI), np.cos(THETA)*np.cos(PHI), -np.sin(PHI)])
        yhat = np.array([np.sin(THETA)*np.sin(PHI), np.cos(THETA)*np.sin(PHI), np.cos(PHI)])
        rhat = np.array([np.sin(THETA)*np.cos(PHI), np.sin(THETA)*np.sin(PHI), np.cos(THETA)])
        that = np.array([np.cos(THETA)*np.cos(PHI), np.cos(THETA)*np.sin(PHI), -np.sin(THETA)])
        phat = np.array([-np.sin(PHI), np.cos(PHI), np.zeros_like(THETA)])

        if inc:
            Hinc = self.source.H(np.array([Xs,Ys,Zs]), self.k)*(self.solution.eps_b/self.solution.mu_b)**0.5
        else:
            Hinc = 0
        Ax, Ay = self.amp[:2]

        Hscat = Ax*self.H_func(R,THETA,PHI) + Ay*self.H_func(R,THETA,PHI-np.pi/2)
        
        # convert to cartesian
        Htot = Hscat[0]*rhat + Hscat[1]*that + Hscat[2]*phat + Hinc
        return Htot

    def force(self, other_particles = [], buffer = BUFFER_DEFAULT, inc=True):
        r = np.array([self.radius + buffer])
        tau = np.linspace(-1,1, Ntheta) 
        theta = np.pi - np.arccos(tau)
        phi = np.linspace(0, 2*np.pi, Nphi)
        R, THETA, PHI = np.meshgrid(r,theta,phi, indexing='ij')

        X = self.center[0] + R*np.sin(THETA)*np.cos(PHI)
        Y = self.center[1] + R*np.sin(THETA)*np.sin(PHI) 
        Z = self.center[2] + R*np.cos(THETA)

        # E and H fields
        E = self.E(X,Y,Z, inc=inc)
        H = self.H(X,Y,Z, inc=inc)

        for p in other_particles:
            E += p.E(X,Y,Z, inc=False)
            H += p.H(X,Y,Z, inc=False)

        E = np.squeeze(E)
        H = np.squeeze(H)
        THETA = np.squeeze(THETA)
        PHI = np.squeeze(PHI)

        # cartesian unit vectors
        xhat = np.array([np.sin(THETA)*np.cos(PHI), np.cos(THETA)*np.cos(PHI), -np.sin(PHI)])
        yhat = np.array([np.sin(THETA)*np.sin(PHI), np.cos(THETA)*np.sin(PHI), np.cos(PHI)])
        zhat = np.array([np.cos(THETA), -np.sin(THETA), np.zeros_like(THETA)])

        # spherical unit vectors
        rhat = np.array([np.sin(THETA)*np.cos(PHI), np.sin(THETA)*np.sin(PHI), np.cos(THETA)])
        that = np.array([np.cos(THETA)*np.cos(PHI), np.cos(THETA)*np.sin(PHI), -np.sin(THETA)])
        phat = np.array([-np.sin(PHI), np.cos(PHI), np.zeros_like(THETA)])

        eps_b = self.solution.eps_b
        mu_b = self.solution.mu_b
        sigma = eps_b*np.einsum('ixy,jxy->ijxy', E, np.conj(E)) \
                + mu_b*np.einsum('ixy,jxy->ijxy', H, np.conj(H)) \
                - 0.5*np.einsum('ij,xy->ijxy', np.identity(3), eps_b*np.sum(np.abs(E)**2, axis=0)) \
                - 0.5*np.einsum('ij,xy->ijxy', np.identity(3), mu_b*np.sum(np.abs(H)**2, axis=0))

        # compute F
        dA = r[0]**2
        integrand = np.einsum('ijxy,jxy->ixy', sigma, rhat)*dA
        F = np.array([simps_2d(tau, phi, integrand[i].real) for i in range(3)])

        # compute T
        integrand = np.einsum('imn,mxy,njxy,jxy->ixy', levi, r[0]*rhat, sigma, rhat)*dA
        T = np.array([simps_2d(tau, phi, integrand[i].real) for i in range(3)])
        return F,T
    
    def flux(self, other_particles = [], buffer = BUFFER_DEFAULT, inc=False):
        r = np.array([self.radius + buffer])
        tau = np.linspace(-1,1, Ntheta) 
        theta = np.pi - np.arccos(tau)
        phi = np.linspace(0, 2*np.pi, Nphi)
        R, THETA, PHI = np.meshgrid(r,theta,phi, indexing='ij')

        X = self.center[0] + R*np.sin(THETA)*np.cos(PHI)
        Y = self.center[1] + R*np.sin(THETA)*np.sin(PHI) 
        Z = self.center[2] + R*np.cos(THETA)

        # E and H fields
        E = self.E(X,Y,Z, inc=inc)
        H = self.H(X,Y,Z, inc=inc)

        for p in other_particles:
            E += p.E(X,Y,Z, inc=False)
            H += p.H(X,Y,Z, inc=False)

        E = np.squeeze(E)
        H = np.squeeze(H)
        THETA = np.squeeze(THETA)
        PHI = np.squeeze(PHI)

        # cartesian unit vectors
        xhat = np.array([np.sin(THETA)*np.cos(PHI), np.cos(THETA)*np.cos(PHI), -np.sin(PHI)])
        yhat = np.array([np.sin(THETA)*np.sin(PHI), np.cos(THETA)*np.sin(PHI), np.cos(PHI)])
        zhat = np.array([np.cos(THETA), -np.sin(THETA), np.zeros_like(THETA)])

        # spherical unit vectors
        rhat = np.array([np.sin(THETA)*np.cos(PHI), np.sin(THETA)*np.sin(PHI), np.cos(THETA)])
        that = np.array([np.cos(THETA)*np.cos(PHI), np.cos(THETA)*np.sin(PHI), -np.sin(THETA)])
        phat = np.array([-np.sin(PHI), np.cos(PHI), np.zeros_like(THETA)])

        eps_b = self.solution.eps_b
        mu_b = self.solution.mu_b
        S = np.real(np.einsum('ijk,jxy,kxy->ixy', levi, E, np.conj(H)))

        # compute Flux
        dA = r[0]**2
        integrand = np.einsum('ixy,ixy->xy', S, rhat)*dA
        Flux = simps_2d(tau, phi, integrand)
        return Flux


class particle_system:
    def __init__(self, bodies, source, eps_b=1, mu_b=1, interactions=True):
        self.particles = [particle(single_mie_sphere(body['Nmax'], body['material'], body['radius'],
                          eps_b, mu_b), body['position'], 0, source) for body in bodies]
        self.source = source
        self.eps_b = eps_b
        self.mu_b = mu_b
        self.k = self.particles[0].k
        self.Nparticles = len(self.particles)

        if (interactions):
            amps = self.solve_interactions()
            for i in range(self.Nparticles):
                self.particles[i].amp = amps[:,i]

    def E(self,X,Y,Z, inc = True):
        Efield = self.particles[0].E(X,Y,Z, inc=inc).squeeze()
        Efield += sum([p.E(X,Y,Z, inc = False).squeeze() for p in self.particles[1:]])
        return Efield 

    def H(self,X,Y,Z, inc = True):
        Hfield = self.particles[0].H(X,Y,Z, inc=inc).squeeze()
        Hfield += sum([p.H(X,Y,Z, inc = False).squeeze() for p in self.particles[1:]])
        return Hfield 

    def solve_interactions(self):
        pos = np.array([p.center for p in self.particles]).T
        Einc = self.source.E(pos,self.k)
        Einc = Einc[:2,:]
        
        identity = np.zeros(shape = (2, self.Nparticles, 2, self.Nparticles), dtype=np.complex)
        np.einsum('xixi->xi', identity)[...] = 1
        
        MieMatrix = np.zeros(shape = (2, self.Nparticles, 2, self.Nparticles), dtype=np.complex)
        
        for i in range(self.Nparticles):
            for j in range(self.Nparticles):
                if i == j: continue
                pi = self.particles[i].center
                pj = self.particles[j].center
                dji = pi -  pj
                r_ji = np.linalg.norm(dji)
                theta_ji = np.arccos(dji[2]/r_ji)
                phi_ji = np.arctan2(dji[1], dji[0])
                
                rhat = np.array([np.sin(theta_ji)*np.cos(phi_ji), np.sin(theta_ji)*np.sin(phi_ji), np.cos(theta_ji)])
                that = np.array([np.cos(theta_ji)*np.cos(phi_ji), np.cos(theta_ji)*np.sin(phi_ji), -np.sin(theta_ji)])
                phat = np.array([-np.sin(phi_ji), np.cos(phi_ji), np.zeros_like(theta_ji)])
                
                xsol = self.particles[j].E_func(np.array([r_ji]), np.array([theta_ji]), np.array([phi_ji])).squeeze()
                ysol = self.particles[j].E_func(np.array([r_ji]), np.array([theta_ji]), np.array([phi_ji - np.pi/2])).squeeze()
                xsol = xsol[0]*rhat + xsol[1]*that + xsol[2]*phat
                ysol = ysol[0]*rhat + ysol[1]*that + ysol[2]*phat
                
                MieMatrix[:,i,0,j] = xsol[:2]
                MieMatrix[:,i,1,j] = ysol[:2]

        A = identity - MieMatrix
        sol = np.linalg.solve(A.reshape(2*self.Nparticles, 2*self.Nparticles), Einc.reshape(2*self.Nparticles)).reshape(2, self.Nparticles)
        
        return sol

    def particle_flux(self, i, buffer = BUFFER_DEFAULT, inc = False):
        other_particles = (self.particles[j] for j in range(self.Nparticles) if j != i)
        return self.particles[i].flux(other_particles, buffer=buffer, inc=inc)

    def particle_force(self, i, buffer = BUFFER_DEFAULT, inc = True):
        other_particles = (self.particles[j] for j in range(self.Nparticles) if j != i)
        return self.particles[i].force(other_particles, buffer=buffer, inc=inc)

    def forces(self, buffer = BUFFER_DEFAULT, inc = True):
        all_forces  = np.zeros([self.Nparticles, 3])
        all_torques = np.zeros([self.Nparticles, 3])
        for i in range(self.Nparticles):
            F,T = self.particle_force(i, buffer = buffer, inc = inc)
            all_forces[i]  = F
            all_torques[i] = T

        return all_forces, all_torques

    def fluxes(self, buffer = BUFFER_DEFAULT, inc = False):
        all_fluxes  = np.zeros(self.Nparticles)
        for i in range(self.Nparticles):
            F = self.particle_flux(i, buffer = buffer, inc = inc)
            all_fluxes[i]  = F

        return all_fluxes

    def net_force(self,center = None, radius = None):
        if center is None:
            center = self.center_of_mass()
        center = np.asarray(center)

        if radius is None:
            radius = 0
            for p in self.particles:
                radius = max(radius, np.linalg.norm(p.center - center) + p.radius)
            radius += 1

        r = np.array([radius])
        tau = np.linspace(-1,1, Ntheta) 
        theta = np.pi - np.arccos(tau)
        phi = np.linspace(0, 2*np.pi, Nphi)
        R, THETA, PHI = np.meshgrid(r,theta,phi, indexing='ij')

        X = center[0] + R*np.sin(THETA)*np.cos(PHI)
        Y = center[1] + R*np.sin(THETA)*np.sin(PHI) 
        Z = center[2] + R*np.cos(THETA)

        # E and H fields
        E = self.E(X,Y,Z)
        H = self.H(X,Y,Z)
        E = np.squeeze(E)
        H = np.squeeze(H)
        THETA = np.squeeze(THETA)
        PHI = np.squeeze(PHI)

        # cartesian unit vectors
        xhat = np.array([np.sin(THETA)*np.cos(PHI), np.cos(THETA)*np.cos(PHI), -np.sin(PHI)])
        yhat = np.array([np.sin(THETA)*np.sin(PHI), np.cos(THETA)*np.sin(PHI), np.cos(PHI)])
        zhat = np.array([np.cos(THETA), -np.sin(THETA), np.zeros_like(THETA)])

        # spherical unit vectors
        rhat = np.array([np.sin(THETA)*np.cos(PHI), np.sin(THETA)*np.sin(PHI), np.cos(THETA)])
        that = np.array([np.cos(THETA)*np.cos(PHI), np.cos(THETA)*np.sin(PHI), -np.sin(THETA)])
        phat = np.array([-np.sin(PHI), np.cos(PHI), np.zeros_like(THETA)])

        eps_b = self.eps_b
        mu_b = self.mu_b
        sigma = eps_b*np.einsum('ixy,jxy->ijxy', E, np.conj(E)) \
                + mu_b*np.einsum('ixy,jxy->ijxy', H, np.conj(H)) \
                - 0.5*np.einsum('ij,xy->ijxy', np.identity(3), eps_b*np.sum(np.abs(E)**2, axis=0)) \
                - 0.5*np.einsum('ij,xy->ijxy', np.identity(3), mu_b*np.sum(np.abs(H)**2, axis=0))

        # compute F
        dA = (4*np.pi*r[0]**2/(len(theta)*len(phi)))
        integrand = np.einsum('ijxy,jxy->ixy', sigma, rhat)*dA
        F = np.array([simps_2d(tau, phi, integrand[i].real) for i in range(3)])

        # compute T
        integrand = np.einsum('imn,mxy,njxy,jxy->ixy', levi, r[0]*rhat, sigma, rhat)*dA
        T = np.array([simps_2d(tau, phi, integrand[i].real) for i in range(3)])
        return F,T
