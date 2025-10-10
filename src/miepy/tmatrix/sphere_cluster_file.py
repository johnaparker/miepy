import numpy as np

def single_paricle_string(pos, radius, Nrank, index):
    return f"""
TmatPart           
{radius}
{Nrank}
{Nrank}
({index.real}, {index.imag})

GeomPartProp       
{pos[0]}
{pos[1]}
{pos[2]}
""".strip()

def sphere_cluster_file(pos, radii, Nrank_particles, Nrank, wavelength, index, index_m):
    """Create input file for axisymmetric particles

    Arguments:
        Nrank (int)                  maximum number of multipoles
        wavelength (float)           wavelength of incident light
        index (complex)              index of refraction of the particle
        index_m (float)              index of refraction of the medium
    """
    pos = pos/wavelength
    radii = radii/wavelength
    index = index/index_m
    wavelength = 1

    Nparticles = len(pos)
    radius_cluster = np.max(np.linalg.norm(pos, axis=-1))
    identical = str(np.all(radii == radii[0])).lower()

    particle_string = ''
    for i in range(Nparticles):
        particle_string += single_paricle_string(pos[i], radii[i], Nrank_particles[i], index[i]) + '\n\n'

    return f"""OptProp
{wavelength}
{index_m.real}
Variables:
 - wavelength - wavelength of the incident light in vacuo.
 - ind_refMed - refractive index of the ambient medium.

GenProp            
{Nparticles}
{radius_cluster}
{radius_cluster}
.{identical}.       
Variables:
- Npart     - number of particles.
- anorm     - characteristic length of the cluster which is used 
              to normalize the differential scattering cross sections.
- Rcirc     - characteristic length of the cluster which is used to 
              compute an estimate of Nrank.
- identical - if identical = true, the spherical particles have the 
              same radius.

{particle_string}

ConvTest           
.false.   
.false.     
Variables:
- DoConvTest  - if DoConvTest = t, the interactive convergence tests 
                over Nrank and Mrank are performed.  
- ExtThetaDom - if ExtThetaDom = t the DSCS is computed for scattering 
                angles ranging from 0° to 180° in the azimuthal plane 
                phiGS = 0° and from 180° to 0° in the azimuthal plane 
                phiGS = 180°.


NrankMrankCluster  
{Nrank}
{Nrank}
Variables:
- Nrank - maximum expansion order for the cluster.
- Mrank	- maximum azimuthal order for the cluster.
NOTE: THESE VARIABLES MUST BR PROVIDED IF DoConvTest = f.

						
Errors             
 5.e-2        
 5.e-2        
Variables:
- epsNrank - error tolerance for the expansion order test. 
- epsMrank - error tolerance for the azimuthal order test.       

		    
Tmat               
'../TMATFILES/tmatrix.dat'
Variables:
- FileTmat - name of the file to which the T matrix of the cluster 
             is written.   

 
PrintProgress      
.false.      
Variables:
- PrnProgress - if PrnProgress = t, the progress of calculation 
                is printed.
"""
