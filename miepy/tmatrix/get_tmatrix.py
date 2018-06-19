import os
import subprocess
import tempfile 

import numpy as np
import miepy
from functools import namedtuple

tmatrix_input = namedtuple('TmatrixInput', 'number, input_function')
class tmatrix_solvers:
    axisymmetric = tmatrix_input(number=1, input_function=None)


def nfmds_solver(lmax, input_kwargs, solver=tmatrix_solvers.axisymmetric, extended_precision=False):
    rmax = miepy.vsh.lmax_to_rmax(lmax)
    input_kwargs.update(dict(Nrank=lmax+2))

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
            f.write(input_file_template)

        with open(f'{input_files_dir}/InputSCT.dat', 'w') as f:
            f.write(input_sct_template)

        with open(f'{input_files_dir}/InputAXSYM.dat', 'w') as f:
            f.write(file_str_template.format(**input_kwargs))

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


file_str_template = """OptProp
6.0
1.0
(2.0,0.0)
Variables:
 - wavelength - wavelength of the incident light in vacuo.
 - ind_refMed - refractive index of the ambient medium.
 - ind_refRel - relative refractive index of the particle.  

MatProp
.false.
.false.
0.1
 Variables:
 - perfectcond - if perfectcond = t, the particle is perfectly conducting.
 - chiral      - if chiral = t, the particle is optical active (chiral).
 - kb          - parameter of chirality.

GeomProp
.false.
'../GEOMFILES/prolate.fem'
1
2
1.5
1.0
1
1.0
1.0
.false.
 Variables:
 - FileGeom - if FileGeom = t, the particle geometry is supplied by the 
              input file FileFEM. 
 - FileFEM  - name of the file containing the particle geometry. 
 - TypeGeom - parameter specifying the type of the particle geometry.
 - Nsurf	   - number of surface parameters. 
 - surf(1)  - surface parameter.
 - ...  
 - surf(Nsurf  
 - Nparam   - number of smooth curves forming the generatrix curve.    
 - anorm    - characteristic length of the particle which is used to 
              normalize the differential scattering cross sections.	 
 - Rcirc    - characteristic length of the particle for computing Nrank. 
 - miror    - if miror = t, the particle is mirror symmetric.	            
 NOTE: FOR CHIRAL PARTICLES AND DISTRIBUTED SOURCES SET miror = f.

ConvTest
.false.
.false.
 Variables:
 - DoConvTest   - if DoConvTest = t, the interactive convergence tests 
                  over Nint and Nrank are performed.   
 - MishConvTest - if MishConvTest = t, estimates of Nint and Nrank are  
                  computed with the convergence criterion proposed by 
                  Mishchenko.        
 NOTE: IF THE PARTICLE IS OPTICAL ACTIVE (chiral = t) OR THE PARTICLE
 GEOMETRY IS SUPPLIED BY THE FILE FileFEM (FileGeom = t), THE CODE SETS
 MishConvTest = f. IN FACT, MISHCHENKOS CONVERGENCE TEST WILL BE 
 PERFORMED IF (DS = f AND DoConvTest = t AND chiral = f AND FileGeom = f), 
 OR (DS = t AND autGenDS = t AND DoConvTest = t AND chiral = f AND 
 FileGeom = f).   

Sources
.true.
.true.
 Variables:
 - DS       - if DS = t, distributed sources are used for T-matrix 
              calculation. 	
 - autGenDS - if autGenDS = t, the coordinates of the distributed sources
              are generated by the code.
 NOTE: IF THE PARTICLE GEOMETRY IS READ FROM FILE (FileGeom = t),
 THE CODE SETS autgenDS = f.                                 

SourcePosAut
.true.
0.95
 Variables: 
 - ComplexPlane - if ComplexPlane = t, the distributed sources are placed
                  in the complex plane.
 - EpsZReIm     - parameter controlling the distribution of the discrete 
                  sources.
 NOTE: THESE VARIABLES MUST BE PROVIDED IF (DS = t AND autgenDS = t).

NintNrank
{Nint}
{Nrank}
 Variables: 
 - Nint  - number of integration points in computing integrals over the 
           generatrix curve.
 - Nrank - maximum expansion order.  
 NOTE: THESE VARIABLES MUST BE PROVIDED IF ((DoConvTest = f) OR 
 (DS = t AND autgenDS = f)).                  

Errors
5.e-2
5.e-2
1.e-2
4
50
 Variables:
 - epsNint    - error tolerance for the integration test.    
 - epsNrank   - error tolerance for the expansion order test.  
 - epsMrank   - error tolerance for the azimuthal order test.  
 - dNint	     - number of division points for the integration test 
                and Mishchenkos convergence test.   
 - dNintMrank - number of division points for azimuthal mode 
                calculation.

Tmat
'../TMATFILES/tmatrix.dat'
 Variable:
 - FileTmat - name of the file to which the T matrix is written.  

PrintProgress
.false.
 Variable:
 - PrnProgress - if PrnProgress = t, the progress of calculation 
                 is printed. 
 Comment
 This file was generated by the routine smuthi.nfmds_wrappers.taxsym_write_input_spheroid 
"""

input_file_template = """MatrixSolver 
'LU2'      
 100       
' '       
 6         
 1.e-6       			             
Variables:
- TypeMatrSolv - specifies the method for solving linear algebraic 
                 equations. The permissive values are: 'LU1', 'LU2', 
                 'LU3' and 'BICG'. 
- itmax        - maximum number of iterations for the 'LU3' method or 
                'BICG' method.        
- TypePrecond  - specifies the type of preconditioning for the BICG 
                 method. The permissive values are: 'NEUMANN', 'SILU' 
                 and ' '. 
- NPrecOrder   - truncation order of the Neumann series.  
- epsilon      - tolerance for the BICG method.     


Interpolation
'LINEAR' 
Variables:
- TypeInterp - specifies the type of interpolation for scattering matrix 
               calculation. The permissive values are: 'LINEAR', 'SPLINE' 
               and 'HERMITE'.


Integration
'MET1'             
 1.e-15     
 1.e-10    
Variables:
- TypeIntegr  - specifies the source of the numerical integration routines.
                The permissive values are: 'MET1' (routines from Numerical 
                Recipes) and 'MET2' (modified routines from Slatec library).                     
- epsGauss    - tolerance for computing the roots of the Legendre polynoms.  
- epsLaguerre - tolerance for computing the roots of the Laguerre polynoms. 


RandomNumbers
'SLAT'
Variable:
- TypeRND - specifies the type of the random number generator.
            The permissive values are: 
            'SLAT' (modified routine from SLATEC library), 
            'LPCK' (modified routine from LAPACK library) and 
            'ZIGG' (Ziggurat method of Marsaglia and Tsang)
"""

input_sct_template = """OptProp       
 0.628318530717959                           
Variable:
- wavelength - wavelength in the surrounding medium.
NOTE: THIS VARIABLE MUST BE PROVIDED IF THE ROUTINE SCT.90 IS CALLED.


Tmat          
'../TMATFILES/TProlA10B5.dat'  
 17                  
  7                
.true.           
.false.           
.false.             
Variables:
- FileTmat - name of the file containing the T matrix.   
- Nrank    - maximum expansion order.   
- Mrank    - number of azimuthal modes.  
- axsym    - if axsym = t, the scatterer is a rotationally symmetric 
             particle.  
- sphere   - if sphere = t, the scatterer is a spherical particle.  
- chiral   - if chiral = t, the scatterer is an optical activ particle. 
NOTE: THESE VARIABLES MUST BE PROVIDED IF THE ROUTINE SCT.90 IS CALLED.


TypePDF      
.false.  
Variable:
- RandomOrientation - if RandomOrientation = t, the PDF is complete. 
 
 
***************************    complete PDF     *************************** 
CompletePDF - TypeMedium
.true.
Variable:
- MirorSym - if MirorSym = t, the scattering medium is isotropic and
             mirror-symmetric. 
  
  
CompletePDF - AvrgMtrSS
.false.        
181        
Variables:
- DoNumAvrg - if DoNumAvrg = t, the <SS*> matrix is computed by using the 
              numerical averaging procedure.
- NthetaGS  - number of scattering angles at which the <SS*> matrix is
              computed. 
NOTE: NthetaGS MUST BE AN ODD NUMBER.     
  
  
CompletePDF - AnalytComputAvrgMtrSS  
.true.         
 1.e-5 
Variables:   
- ReducedOrder - if ReducedOrder = t, the analytical averaging procedure 
                 (DoNumAvrg = f) computes the average matrix <SS*> for 
                  reduced values of Mrank and Nrank. 
- deltaOrder   - error tolerance for the convergence test over the 
                 extinction and scattering cross sections. 
  
  
CompletePDF - NumComputAvrgMtrSS
.true.
31
61
31
Variables:   
- UseSimpson -  if UseSimpson = t, the integration over the orientation 
                angle beta is transformed into an integration over 
                cos(beta) and the Simpson rule is used for numerical 
                integration (uniformly sample in cos(beta)).
- Nalpha      - number of division points for averaging <SS*> over alpha.   
- Nbeta       - number of division points for averaging <SS*> over beta. 
- Ngamma      - number of division points for averaging <SS*> over gamma.
NOTE: Nalpha and Ngamma MUST BE ODD NUMBERS. IF UseSimpson = t, NGAMMA 
MUST BE ALSO AN ODD NUMBER.                        


CompletePDF - NormConst
1.0
Variable:  
- anorm - characteristic length of the particle which is used to 
          normalize the average differential scattering cross sections 
          and the optical cross sections.
  
  
CompletePDF - DSCS  
.true. 
(1.0,0.0)
(1.0,0.0)         
 0.0                               
.true.          
'../OUTPUTFILES/DscsRandom.dat'     			   
Variables:
- ComputeDSCS - if ComputeDSCS = t, the average differential scattering 
                cross sections are computed in the azimuthal plane phiGS.
- EI_betaGI   - complex amplitude in the beta-direction of the incident 
                plane wave (parallel component).
- EI_alphaGI  - complex amplitude in the alpha-direction of the incident
                plane wave (perpendicular component).   
- phiGS       - azimuthal angle specifying the scattering plane at which 
                the DSCS are computed.     
- normalized  - if normalized = t, the average differential scattering 
                cross sections are normalized by Pi * anorm**2.  
- FileDSCS    - name of the file to which the average differential 
                scattering cross sections are written.  
NOTE: THE DSCS ARE COMPUTED IN THE AZIMUTHAL PLANE phiGS AT NthetaGS 
SCATTERING ANGLES, I.E., AT THE SAME ANGLES AT WHICH THE AVERAGE MATRIX
<SS*> IS COMPUTED.  
  
  
CompletePDF - ScatPars  
.true.  
 91
  0.0
180.0  
'../OUTPUTFILES/ScatRandom.dat'
 10
 11
 22
 33
 44
 21
 32
 43
 31
 42
 41    
Variables:
- ComputeScatPar - if ComputeScatPar = t, the scattering characteristics 
                   are computed at specified scattering directions. 
                   Specifically, the scattering matrix will be computed 
                   at NthetaRND zenith angles uniformly spaced in the 
                   interval (thetaminRND, thetamaxRND).  
- NthetaRND      - number of zenith angle values at which the scattering 
                   matrix is computed.     
- thetaminRND    - minimal value of the zenith angle at which the scattering 
                   matrix is computed.
- thetamaxRND    - maximal value of the zenith angle at which the scattering 
                   matrix is computed.  
- FileScat       - name of the file to which the scattering characteristics 
                   are written.
- Nelem          - number of scattering matrix elements to be printed out 
                   in the output file FileSCAT.
- MatrixElem(1)  - matrix element to be printed.
- ...
- MatrixElem(Nelem)  
  
  			     
**************************    incomplete PDF     **************************   
IncompletePDF - TypeIntegration
.true.
Variable:  
- UseSimpson - if UseSimpson = t, the integration over the orientation 
               angle beta is transformed into an integration over cos(beta) 
               and Simpson's rule is used for numerical integration 
              (uniformly sample in cos(beta)). In addition, if 
               UseSimpson = t, the integration over theta for computing the 
               mean direction of propagation of the scattered field is also 
               performed with Simpson's rule.  
    
  
IncompletePDF - OrientationAngles      
 45.0              
 45.0              
  1                
 45.0              
 45.0              
  1                
  0.0              
  0.0              
  1                  
Variables:
- alphamin - minimal value of the Euler orientation angle alpha.
- alphamax - maximal value of the Euler orientation angle alpha. 
- Nalpha   - number of division points for averaging the scattering 
             characteristics over alpha.   
- betamin  - minimal value of the Euler orientation angle beta. 
- betamax  - maximal value of the Euler orientation angle beta. 
- Nbeta    - number of division points for averaging the scattering 
             characteristics over beta.    
- gammamin - minimal value of the Euler orientation angle gamma.
- gammamax - maximal value of the Euler orientation angle gamma.
- Ngamma   - number of division points for averaging the scattering 
             characteristics over gamma.
NOTE: Nalpha and Ngamma MUST BE ODD NUMBERS. IF UseSimpson = t, NGAMMA 
MUST BE ALSO AN ODD NUMBER.    


IncompletePDF - NormConst
1.0
Variable:  
- anorm - characteristic length of the particle which is used to normalize 
          the average differential scattering cross sections and the 
          optical cross sections.
                      

IncompletePDF - IncWave
'PLANE'           
  0.0              
  0.0                               
  0.0              
  0.0              
  0.0              
  2.5       
Variables:
- TypeExcit - specifies the type of the external excitation. 
              TypeExcit = 'PLANE', for a plane wave excitation and 
              TypeExcit = 'GAUSS' for a Gaussian beam excitation.
- betaGI    - zenith angle of the incident direction in the global 
              coordinate system.
- alphaGI   - azimuthal angle of the incident direction in the global 
              coordinate system.
- x0        - x-coordinate of the middle of the Gaussian beam waist.     
- y0        - y-coordinate of the middle of the Gaussian beam waist.    
- z0        - z-coordinate of the middle of the Gaussian beam waist.
- w0        - waist radius of the Gaussian beam.   


IncompletePDF - DSCS          
.true.     
(1.0,0.0)
(1.0,0.0)
 45.0 
  0.0           
181             
.false.         
.true.          
'../OUTPUTFILES/Dscs.dat'     			   
Variables:
- ComputeDSCS - if ComputeDSCS = t, the average differential scattering 
                cross sections are computed in the azimuthal plane phiGS.
- EI_betaGI   - complex amplitude in the beta-direction of the incident
                plane wave (parallel component).
- EI_alphaGI  - complex amplitude in the alpha-direction of the incident 
                plane wave (perpendicular component).    
- alphapGauss - polarization angle of the Gaussian beam.   
- phiGS       - azimuthal angle specifying the scattering plane.     
- NthetaGS    - number of scattering angles at which the average 
                differential scattering cross sections are computed.   
- ExtThetaDom - if ExtThetaDom = t the DSCS are computed for scattering 
                angles ranging from 0° to 180° in the azimuthal plane phiGS, 
                and from 180° to 0° in the azimuthal plane phiGS + 180°.
- normalized  - if normalized = t, the average differential scattering 
                cross sections are normalized by Pi * anorm**2.  
- FileDSCS    - name of the file to which the average differential 
                scattering cross sections are written. 

  
IncompletePDF - ScatPars
.true.  
  2  
 45.0
  3  
 30.0 
150.0 
225.0  
  3 
 30.0
150.0
'../OUTPUTFILES/Scat.dat'
 10
 11
 22
 33
 44
 21
 32
 43
 31
 42
 41
Variables:
- ComputeScatPar - if ComputeScatPar = t, the scattering characteristics 
                   are computed at specified scattering directions. 
                   Specifically, the average phase matrix is computed at 
                   Nphi scattering planes (characterized by the azimuthal 
                   angles phi(1),...,phi(Nphi))and at Ntheta(1),...,Ntheta(Nphi) 
                   zenith angles (in the intervals (thetamin(1),thetamax(1)),
                   ...,(thetamin(Nphi),thetamax(Nphi)).		        
- Nphi           - number of scattering planes at which the phase matrix is 
                   computed.
- phi(1)         - azimuthal angles specifying the scattering planes at 
                   which the phase matrix is computed.
- Ntheta(1)      - number of zenith angle values in each scattering plane 
                   at which the phase matrix is computed.
- thetamin(1)    - specifies the minimal values of the zenith angle in each 
                   scattering plane at which the phase matrix is computed.
- thetamax(1)    - specifies the maximal values of the zenith angle in each 
                   scattering plane at which the phase matrix is computed.         
- ...
- phi(Nphi)
- Ntheta(Nphi)         
- thetamin(Nphi)  
- thetamax(Nphi)                          
- FileScat       - name of the file to which the scattering characteristics 
                   are written.
- Nelem          - number of phase matrix elements to be printed out in the 
                   output file FileSCAT.
- MatrixElem(1)  - matrix element to be printed.
- ...
- MatrixElem(Nelem)  
  
  
IncompletePDF - MeanDirPropagatScatWave 
.true.
 61              
 31              
Variables:
- ComputeAsymPar - if ComputeAsymPar = t, the mean direction of 
                   propagation of the scattered wave is computed.
- NthetaAsym     - number of integration points over theta for computing 
                   the mean direction of propagation of the scattered wave.
- NphiAsym       - number of integration points over phi for computing 
                   the mean direction of propagation of the scattered wave.
***************************************************************************


PrintInfo     
.true.      
.true.       
Variables:
- PrnProgress    - if PrnProgress = t, the progress of calculation is 
                   printed.  
- WriteInputInfo - if WriteInputInfo = t, the input parameters of the 
                   scattering problem are written to the output files 
                   FileDSCS and FileSCAT.  

  
Comment  
As provided, the input file is setup to calculate the scattering 
characteristics for an incomplete PDF. The incident wave is a 
plane wave propagating along the Z-axis, and a single orientation of 
the particle (characterized by alpha = beta = 45 °, gamma = 0° ) is 
considered. The DSCS is computed for a linearly polarized vector plane 
wave in the azimuthal plane phiGS = 0° and at NthetaGS = 181 scattering 
angles. The phase matrix is calculated at Nphi = 2 scattering planes
characterized by the azimuthal angles phi(1) = 45° and phi(2) = 225°.
In the first and second scattering plane, 3 zenith angles uniformly 
spaced in the intervals [30°,150°] and [30°,150°], respectively, 
are considered.         

If the routine SCT.90 is called, the input file is setup to calculate 
scattering characteristics of a prolate spheroid. The T matrix of the 
spheroid is stored in the file '../TMATFILES/TProlA10B5.dat' and 
Nrank = 20 and Mrank = 7.    
"""

nfmds_solver(3, dict(Nint=200), extended_precision=False)
