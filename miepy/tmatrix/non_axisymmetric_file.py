def non_axisymmetric_file(geometry_type, geometry_parameters, Nrank, Mrank, 
        wavelength, index, index_m, R_symmetry=0, kb=None, conducting=False,
        Nparam=1, Nint_1=50, Nint_2=50):
    """Create input file for non-axisymmetric particles

    Arguments:
        geometry_type (int)          choose from 1 (ellipsoid), 2 (square-prism), 3 (N-prism)
        geometry_parameters (list)   geometric parameters
        Nrank (int)                  maximum number of multipoles
        Mrank (int)                  maximum number of multipoles (azimuthal)
        wavelength (float)           wavelength of incident light
        index (complex)              index of refraction of the particle
        index_m (float)              index of refraction of the medium
        R_symmetry (int)             discrete rotational symmetry (default: 0)
        kb (float)                   parameter of chirality (default: None [no chirality])
        conducting (bool)            if True, particle is conducting (default: False)
        Nparam (int)                 number of smooth curves used in approximate surface (default: 1)
        Nint_1 (int)                 number of points used in first integration (default: 200)
        Nint_2 (int)                 number of points used in second integration (default: 200)
    """
    scaled_geometry_parameters = [p/wavelength for p in geometry_parameters]
    geo_str = ''.join([str(p) + '\n' for p in scaled_geometry_parameters])[:-1]
    Nsurf = len(geometry_parameters)
    wavelength = 1

    if kb is None:
        chiral = False
        kb = 1
    else:
        chiral = True

    file_str_template = """OptProp 
{wavelength}
{index_m.real}
({index.real}, {index.imag})
Variables:
 - wavelength - wavelength of the incident light in vacuo.
 - ind_refMed - refractive index of the ambient medium.
 - ind_refRel - relative refractive index of the particle.  

MatProp
.{conducting}.
.false.
.{chiral}.
{kb}
Variables:
- perfectcond - if perfectcond = t, the particle is perfectly conducting.
- anisotropic - if anisotropic = t, the particle is a uniaxial anisotropic 
                crystal. 
- chiral      - if chiral = t, the particle is optical active.   
- kb          - parameter of chirality.

         
GeomProp
.false.
'../GEOMFILES/prolate.fem'
{geometry_type}
{Nsurf}
{geometry_parameters}
{Nparam}
1.0
1.0
.true.
{R_symmetry}
Variables:
- FileGeom   - if FileGeom = t, the particle geometry is supplied by the 
               input file FileFEM.
- FileFEM    - name of the file containing the particle geometry.			        
- TypeGeom   - parameter specifying the type of the particle geometry.
- Nsurf      - number of surface parameters. 
- surf(1)    - surface parameter.
- ...
- surf(Nsurf)		       
- Nparam     - number of smooth surfaces forming the particle surface.
- anorm      - characteristic length of the particle which is used to 
               normalize the differential scattering cross sections.
- Rcirc      - characteristic length of the particle which is used to 
               compute an estimate of Nrank.
- miror      - if miror = t, the particle is mirror symmetric.
- Nazimutsym - number of azimuthal symmetric sections.
NOTE: FOR CHIRAL AND ANISOTROPIC PARTICLES SET miror = f. FURTHERMORE,
FOR ANISOTROPIC PARTICLES, A REGULAR N-HEDRAL PRISM (TypeGeom = 3) 
IS NOT A VALID GEOMETRY, AND THEREFORE SET Nazimutsym = 0. 


ConvTest      
.false.    
.false.        
Variables: 
- DoConvTest  - if DoConvTest = t, the interactive convergence tests over 
                Nint, Nrank and Mrank are performed.
- ExtThetaDom - if ExtThetaDom = t the DSCS are computed for scattering 
                angles ranging from 0° to 180° in the azimuthal plane 
                phiGS = 0°, and from 180° to 0° in the azimuthal plane 
                phiGS = 180°.


AnSVWF   
(1.7,0.0)       
0.0
0.0
181               
Variables:
- ind_refRelZ - the second relative refractive index of the uniaxial 
                anisotropic particle.
- alphaPR     - azimuthal angle specifying the orientation of the principal 
                coordinate system with respect to the particle coordinate 
                system.
- betaPR      - zenith angle specifying the orientation of the principal 
                coordinate system with respect to the particle coordinate 
                system.
- Nbeta       - number of integration points for computing the vector 
                quasispherical wave functions. 
NOTE: THESE VARIABLES MUST BE PROVIDED IF anisotropic = t.  


NrankMrank    
{Nrank}
{Mrank}
Variables:
- Nrank - maximum expansion order.
- Mrank - maximum azimuthal order.
NOTE: THESE VARIABLES MUST BE PROVIDED IF DoConvTest = f.


Nint          
{Nint_1}
{Nint_2}
Variables: 
- Nint1 - first number of integration points in computing integrals over 
          the particle surface.
- Nint2 - second number of integration points in computing integrals over 
          the particle surface.
NOTE: THESE VARIABLES MUST BE PROVIDED IF 
(DoConvTest = f AND FileGeom = f).  


Errors        
 1.e-2            
 1.e-2             
 1.e-2             
 4               
 4                       
Variables: 
- epsNint  - error tolerance for the integration test.
- epsNrank - error tolerance for the expansion order test.
- epsMrank - error tolerance for the azimuthal order test.
- dNint1   - first number of division points for the integration test.
- dNint2   - second number of division points for the integration test.


Tmat          
'../TMATFILES/tmatrix.dat'    
Variable:
- FileTmat - name of the file to which the T matrix is written.


PrintProgress 
.true.        
Variable: 
- PrnProgress - if PrnProgress = t, the progress of calculation is 
                printed. 
"""
    return file_str_template.format(geometry_type=geometry_type, geometry_parameters=geo_str,
             Nsurf=Nsurf, Nrank=Nrank, Mrank=Mrank, wavelength=wavelength, index=index/index_m, 
             index_m=index_m, chiral=str(chiral).lower(), R_symmetry=R_symmetry, kb=kb,
             conducting=str(conducting).lower(), Nparam=Nparam, Nint_1=Nint_1, Nint_2=Nint_2)
