import numpy as np
import miepy
from miepy.sources import beam

class phase_only_slm(beam):
    """Phase-modulate a beam with a given phase function"""

    def __init__(self, beam, slm):
        """
        Arguments:
            beam      miepy.sources.beam object
            slm       slm(theta,phi)->(phase) function for modulated phase
        """
        super().__init__(power=beam.power, theta_max=beam.theta_max,
                phase=beam.phase, center=beam.center, theta=beam.theta, phi=beam.phi, standing=beam.standing)
        self.beam = beam
        self.slm = slm

    def __repr__(self):
        return 'phase-only SLM modified ' + str(self.beam)

    def E0(self, k):
        return self.beam.E0(k)

    def angular_spectrum(self, theta, phi, k):
        weights = np.exp(1j*self.slm(theta, phi))
        return self.beam.angular_spectrum(theta, phi, k)*weights

    def theta_cutoff(self, k, cutoff=1e-6, tol=None):
        return self.beam.theta_cutoff(k, cutoff, tol)
