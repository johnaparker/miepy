import os
import subprocess
import tempfile 

import numpy as np
import miepy
from functools import namedtuple
from required_files import main_input_file, sct_input_file
from axisymmetric_file import axisymmetric_file


tmatrix_input = namedtuple('TmatrixInput', 'number, input_function')
class tmatrix_solvers:
    axisymmetric = tmatrix_input(number=1, input_function=axisymmetric_file)


def nfmds_solver(lmax, input_kwargs, solver=tmatrix_solvers.axisymmetric, extended_precision=False):
    """Return the T-matrix using the Null-Field Method with discrete sources (NFM-DS)
       
    Arguments:
        lmax           maximum number of multipoles
        input_kwargs   keyword arguments forwarded to solver.input_function
        solver         type of solver to use (default: axisymmetric)
        extended_precision (bool)    whether to use extended precision (default: False)
    """
    rmax = miepy.vsh.lmax_to_rmax(lmax)
    input_kwargs.update(dict(Nrank=lmax))

    ### create temporary directory tree
    with tempfile.TemporaryDirectory() as direc:
        ### create 4 sub-directories
        input_files_dir = f'{direc}/INPUTFILES'
        os.makedirs(input_files_dir)

        out_dir = f'{direc}/OUTPUTFILES'
        os.makedirs(out_dir)

        tmatrix_output_dir = f'{direc}/TMATFILES'
        os.makedirs(tmatrix_output_dir)

        sources_dir = f'{direc}/TMATSOURCES'
        os.makedirs(sources_dir)

        ### write 3 input files
        with open(f'{input_files_dir}/Input.dat', 'w') as f:
            f.write(main_input_file())

        with open(f'{input_files_dir}/InputSCT.dat', 'w') as f:
            f.write(sct_input_file())

        with open(f'{input_files_dir}/InputAXSYM.dat', 'w') as f:
            f.write((solver.input_function(**input_kwargs)))

        ### execute program and communicate
        if extended_precision:
            command = f'{miepy.__path__[0]}/bin/tmatrix_extended'
        else:
            command = f'{miepy.__path__[0]}/bin/tmatrix'

        proc = subprocess.Popen([command], cwd=sources_dir, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
        proc.communicate(f'{solver.number}'.encode())
        proc.wait()
        
        ### read T-matrix dimensions
        with open(f'{tmatrix_output_dir}/Infotmatrix.dat', 'r') as f:
            info_lines = f.readlines()

        for line in info_lines:
            if line.split()[0:4] == ['-', 'maximum', 'expansion', 'order,']:
                n_rank = int(line.split()[-1][0:-1])

            if line.split()[0:5] == ['-', 'number', 'of', 'azimuthal', 'modes,']:
                m_rank = int(line.split()[-1][0:-1])

        ### read T-matrix output
        with open(f'{tmatrix_output_dir}/tmatrix.dat') as f:
            tmat_lines = f.readlines()

        tmatrix_u = [[]]
        column_index = 0
        for line in tmat_lines[3:]:
            split_line = line.split()
            for i_entry in range(int(len(split_line) / 2)):
                if column_index == 2 * n_rank:
                    tmatrix_u.append([])
                    column_index = 0
                tmatrix_u[-1].append(complex(split_line[2 * i_entry]) + 1j * complex(split_line[2 * i_entry + 1]))
                column_index += 1

        ### restructure T-matrix
        T = np.zeros((2*rmax, 2*rmax), dtype=complex)

        m_max = lmax
        for m in range(-m_max, m_max + 1):
            n_max_nfmds = n_rank - max(1, abs(m)) + 1
            for tau1 in range(2):
                for l1 in range(max(1, abs(m)), lmax + 1):
                    n1 = l1*(l1+1) + m -1 + tau1*rmax
                    l1_nfmds = l1 - max(1, abs(m))
                    n1_nfmds = 2 * n_rank * abs(m) + tau1 * n_max_nfmds + l1_nfmds
                    for tau2 in range(2):
                        for l2 in range(max(1, abs(m)), lmax + 1):
                            n2 = l2*(l2+1) + m -1 + tau2*rmax
                            l2_nfmds = l2 - max(1, abs(m))
                            n2_nfmds = tau2 * n_max_nfmds + l2_nfmds
                            if abs(m) <= m_rank:
                                if m >= 0:
                                    T[n1, n2] = tmatrix_u[n1_nfmds][n2_nfmds]
                                else:
                                    T[n1, n2] = tmatrix_u[n1_nfmds][n2_nfmds] * (-1) ** (tau1 + tau2)

        ### restructure T-matrix again
        T = np.reshape(T, [2, rmax, 2, rmax])
        tmatrix = np.empty_like(T)
        tmatrix[0,:,0,:] = -T[1,:,1,:]
        tmatrix[1,:,1,:] = -T[0,:,0,:]
        tmatrix[0,:,1,:] = -T[1,:,0,:]
        tmatrix[1,:,0,:] = -T[0,:,1,:]

        for r1,n1,m1 in miepy.mode_indices(lmax):
            for r2,n2,m2 in miepy.mode_indices(lmax):
                tmatrix[:,r1,:,r2] *= 1j**(n2-n1)

        return tmatrix
