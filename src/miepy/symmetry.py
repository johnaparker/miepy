import numpy as np

class periodic_lattice:
    def __init__(self, generate_func):
        self.generate_func = generate_func

    def generate(self, N):
        return self.generate_func(N)

def generator_1d(N, k):
    N_half = N//2
    
    pos = np.array([-k*(i+1) for i in range(N_half)], dtype=float)
    pos = np.concatenate([pos, np.array([k*(i+1) for i in range(N_half)], dtype=float)]) 

    return pos

def x_translation_1d(k):
    def generate(N):
        xpos = generator_1d(N, k)
        ypos = np.zeros_like(xpos)
        zpos = np.zeros_like(xpos)

        return xpos, ypos, zpos

    return periodic_lattice(generate)

def y_translation_1d(k):
    def generate(N):
        ypos = generator_1d(N, k)
        xpos = np.zeros_like(xpos)
        zpos = np.zeros_like(xpos)

        return xpos, ypos, zpos

    return periodic_lattice(generate)

def z_translation_1d(k):
    def generate(N):
        zpos = generator_1d(N, k)
        xpos = np.zeros_like(xpos)
        ypos = np.zeros_like(xpos)

        return xpos, ypos, zpos

    return periodic_lattice(generate)

def square_lattice_2d(k):
    def generate(N):
        Nx = int(np.sqrt(N))//2*2 + 1
        Ny = Nx
        
        x = np.arange(Nx)*k
        y = np.arange(Ny)*k
        x -= np.average(x)
        y -= np.average(y)

        X, Y = np.meshgrid(x, y)
        xpos = X.flatten()
        ypos = Y.flatten()

        rad = xpos**2 + ypos**2
        idx = np.argsort(rad)

        xpos = xpos[idx][1:]
        ypos = ypos[idx][1:]

        zpos = np.zeros_like(xpos)

        return xpos, ypos, zpos

    return periodic_lattice(generate)
