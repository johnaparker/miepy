import os
import subprocess
import tempfile 

import numpy as np
import miepy
import pandas
from functools import namedtuple
from .required_files import main_input_file, sct_input_file
from .axisymmetric_file import axisymmetric_file


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

    ### create temporary directory tree
    with tempfile.TemporaryDirectory() as direc:
        ### create 4 sub-directories
        input_files_dir = '{direc}/INPUTFILES'.format(direc=direc)
        os.makedirs(input_files_dir)

        out_dir = '{direc}/OUTPUTFILES'.format(direc=direc)
        os.makedirs(out_dir)

        tmatrix_output_dir = '{direc}/TMATFILES'.format(direc=direc)
        os.makedirs(tmatrix_output_dir)

        sources_dir = '{direc}/TMATSOURCES'.format(direc=direc)
        os.makedirs(sources_dir)

        ### write 3 input files
        with open('{input_files_dir}/Input.dat'.format(input_files_dir=input_files_dir), 'w') as f:
            f.write(main_input_file())

        with open('{input_files_dir}/InputSCT.dat'.format(input_files_dir=input_files_dir), 'w') as f:
            f.write(sct_input_file())

        with open('{input_files_dir}/InputAXSYM.dat'.format(input_files_dir=input_files_dir), 'w') as f:
            f.write((solver.input_function(Nrank=lmax, **input_kwargs)))

        ### execute program and communicate
        install_path = miepy.__path__[0]
        if extended_precision:
            command = '{install_path}/bin/tmatrix_extended'.format(install_path=install_path)
        else:
            command = '{install_path}/bin/tmatrix'.format(install_path=install_path)

        proc = subprocess.Popen([command], cwd=sources_dir, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
        proc.communicate(str(solver.number).encode())
        proc.wait()

        ### read T-matrix dimensions
        with open('{tmatrix_output_dir}/Infotmatrix.dat'.format(tmatrix_output_dir=tmatrix_output_dir), 'r') as f:
            lines = f.readlines()

            ### last entries in the final two lines give Nrank, Mrank, respectively
            m_rank_str = lines[-1].split()[-1]
            n_rank_str = lines[-2].split()[-1]

            ### Remove the pesky comma and period from the string
            m_rank = int(m_rank_str[:-1])
            n_rank = int(n_rank_str[:-1])

        ### read T-matrix output
        tmatrix_file = '{tmatrix_output_dir}/tmatrix.dat'.format(tmatrix_output_dir=tmatrix_output_dir)
        data = pandas.read_csv(tmatrix_file, skiprows=3,
                  delim_whitespace=True, header=None).values.flatten()  # read as flat array
        data = data[~np.isnan(data)]  # throw out the NaNs
        data_real = data[::2]  # every other element is the real part
        data_imag = data[1::2]
        T_nfmds = np.reshape(data_real, (-1, 2*n_rank)) \
                  + 1j*np.reshape(data_imag, (-1, 2*n_rank))  # reshape the final result to have 2*n_rank columns

        ### restructure T-matrix
        T = np.zeros((2, rmax, 2, rmax), dtype=complex)
        for r1,n1,m1 in miepy.mode_indices(lmax, m_start=-m_rank, m_stop=m_rank):
            for r2,n2,m2 in miepy.mode_indices(lmax, m_start=-m_rank, m_stop=m_rank):
                if m1 != m2:
                    continue

                n_max = n_rank - max(1, abs(m1)) + 1
                l1 = n1 - max(1, abs(m1)) 
                l2 = n2 - max(1, abs(m2)) 

                x = 2*n_rank*abs(m1) + l1
                y = l2

                factor = -1j**(n2-n1)
                T[1,r1,1,r2] = T_nfmds[x, y]*factor
                T[0,r1,0,r2] = T_nfmds[x+n_max, y+n_max]*factor
                T[0,r1,1,r2] = T_nfmds[x+n_max, y]*factor*np.sign(m1)
                T[1,r1,0,r2] = T_nfmds[x, y+n_max]*factor*np.sign(m2)

        return T
