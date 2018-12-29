subroutine TNONAXSYM
!------------------------------------------------------------------------------------
! 1. General Considerations                                                         !
! --------------------------                                                        !
! TNONAXSYM is a routine for computing the T matrix and the scattering              !
! characteristics of homogeneous, dielectric (isotropic, uniaxial anisotropic,      !
! chiral) and perfectly conducting, nonaxisymmetric particles.                      !
!                                                                                   !
! For electrically anisotropic media, the permittivity is a tensor and in general,  !
! the permittivity tensor of a crystal is symmetric. Since there always exists a    !
! coordinate transformation that transforms a symmetric matrix into a diagonal      !
! matrix, we can take this coordinate system as reference frame. This coordinate    !
! system is called the principal coordinate system, and in the principal coordinate !
! system the permittivity tensor is a diagonal matrix with elements eps_x, eps_y    !
! and eps_z. The three coordinate axes are the principal axes of the crystal, and   !
! if eps = eps_x = eps_y and eps_z /= eps, the  medium is uniaxial anisotropic.     ! 
!                                                                                   !
! For anisotropic particles, the internal fields are approximated by the system of  !
! vector quasispherical wave functions. These vector functions can be regarded as a !
! generalization of the regular vector spherical wave functions for uniaxial        !
! anisotropic media, and can be expressed as integrals of plane waves over the      !
! spherical angles alpha and beta. The integration over the azimuth angle alpha can !
! be analytically performed and essentially, the vector quasispherical wave         !
! functions are expressed as simple integrals over the zenith angle beta.           !
!                                                                                   !
! The integrals over the particle surface are evaluated by means of the             !
! Gauss-Legendre quadrature method. For particles with a plane of symmetry          !
! perpendicular to the axis of rotation, i.e., mirror symmetric particles with the  !
! surface parametrization                                                           !
!                                                                                   !
!                      r(theta,phi) = r(Pi - theta,phi),                            !
!                                                                                   !
! the T matrix can be computed by integrating theta over the interval [0,Pi/2].     !
! For particles with azimuthal symmetry, i.e., particles with the surface           !       
! parametrization                                                                   !
!                                                                                   !   
!                      r(theta,phi) = r(theta,phi + 2*Pi / Nazimutsym),             !
!                                                                                   !                  
! where Nazimutsym >= 2, the T matrix can be computed by integrating phi over the   !
! interval [0,2*Pi/Nazimutsym]. The numerical stability and accuracy of T-matrix    !
! computations for an unsmooth surface can be enhanced by dividing the surface      !
! into piecewise smooth surfaces and applying separate quadrature formulas to each  !
! surface. Note that no mirror and azimuthal symmetries are considered for          !
! anisotropic particles.                                                            !
!                                                                                   !
! Particle geometries currently supported include: ellipsoids, quadratic prisms     !
! and regular polyhedral prisms (regular N-hedral finite prisms with azimuthal      !
! and mirror symmetries). A regular N-hedral finite prism refers to a shape with    !
! N + 2 facets, that is, to a prism with top and bottom facets with the shape of a  !
! regular polygon with N corners and with N equal-sized rectangular side facets.    !
! Certainly, quadratic prisms also possess azimuthal and mirror symmetries, but in  !
! this specific case, the T matrix is computed by disregarding the surface          !
! symmetries. Quadratic prisms have been included for anisotropic particles and     !
! testing purposes. As mentioned before, for anisotropic particles, no mirror and   !
! azimuthal symmetries are taken into account. Therefore, the regular N-hedral      !
! finite prisms with azimuthal and mirror symmetries is NOT a valid geometry for    !
! anisotropic particles. The following routines (from the file "GeomLib.f90")       !
! provide the required geometry parameters:                                         !
! - the routine "interpolation_list3D" computes the integration nodes and weights   !
!   for a nonaxisymmmetric surface,                                                 !
! - the routine "elem_geom3D" provides the geometry parameters (position vector,    !
!   polar and azimuthal angles and normal unit vector in spherical coordinates)     !
!   at each integration node, and                                                   !
! - the routine "SurfaceElemNint" computes some parameters of the discretized       !
!   surface: the mean, maximum and minimum areas of the surface elements and the    !
!   number of elements per unit surface (number density).                           !
! The user can easily modify these routines to generate particles with other        !
! geometries. Note that the list of parameters must be maintained and only          !
! geometries with an analytical description of the surface can be implemented.      !
! Alternatively, the option FileGeom  = .true., causes the code to read the         !
! particle geometry information from the file FileFEM. The FileFEM is read by the   !
! routine "read_FileFEM" from the file "InputOutput.f90". The user can customize    !
! the routine "read_FileFEM" as needed to conform to the manner in which the        !
! particle geometry is stored in the file FileFEM. However, as supplied, the        !
! routine "read_FileFEM" expects the file FileFEM to have the following structure:  !
! - one line providing the number of smooth surfaces Nfaces; the format for         !
!   reading Nfaces is i7, i.e.,                                                     !
!                 read (iFEM, "(i7)", iostat = ios)  Nfaces;                        ! 
! - a set of Nfaces-sequences of data specifying the geometry parameters on each    !
!   smooth surface; each sequence includes                                          !
!   - one line containing the number of surface elements (vertices) NVvr; the       !
!     format for reading NVvr is i7, i.e.,                                          !
!                 read (iFEM, "(i7)", iostat = ios) NVvr;                           !                                         
!   - NVvr lines containing the current index ielem, the Cartesian coordinates      !   
!     of the surface element center: x = r(1), y = r(2) and z = r(3), the Cartesian ! 
!     components of the unit normal vector at the surface element center:           !
!     nx = n(1), ny = n(2) and nz = n(3), and the area of the surface element;      !
!     the format for reading the data is:                                           !
!                 read (iFEM, "(i7,2x,7(e15.7,2x))", iostat = ios)                  !
!                       ielem, r(1), r(2), r(3), n(1), n(2), n(3), area.            ! 
! As provided, the file FileFEM is contained in the directory "GEOMFILES".          !
!                                                                                   !
! 2. Convergence Test                                                               !
! -------------------                                                               !
! Convergence tests precede the scattering characteristics calculation. The         !
! following parameters control the T-matrix computation:                            !
! - the numbers of integration points for computing integrals over the particle     !
!   surface Nint1 and Nint2,                                                        !
! - the maximum expansion order Nrank and                                           !
! - the maximum azimuthal order Mrank.                                              ! 
! Nint1 and Nint2 are the global numbers of integration points on the particle      !
! surface. They depend on the surface parametrization and have the following        ! 
! significance.                                                                     !
! a. For an ellipsoid with semi-axes a, b and c, the surface parametrization is     !
! (explicit form)                                                                   !
!                                                                                   !
!                 x = a * sin(u) * cos(v),                                          !
!                 y = b * sin(u) * sin(v),                                          !
!                 z = c * cos(u),                                                   !
!                                                                                   !
! where u belongs to [0,Pi] and v belongs to [0,2Pi]. Nint1 and Nint2 are the       !
! number of division points for the integration over u and v, respectively. If      ! 
! miror = .true., u is integrated over [0,Pi/2] and the number of integration       !
! points in the interval [0,Pi/2] is Nint1 / 2.                                     !
! b. The surface parameters describing regular polyhedral prism are the prism length!
! and the side of the regular polygon basis. In this specific case it is appropriate! 
! to describe each face with respect to a local (two-dimensional) Cartesian         !
! coordinate system. For regular polyhedral prisms, Nint1 is the number of          !
! integration points along the prism length (for the integration over the lateral   !
! rectangular faces), and  Nint2 is the number of integration points along the side !
! of the regular polygon basis (for the integration over the lateral rectangular    !
! faces, and over the top and bottom faces). Roughly speaking, the number of        !
! integration points on a lateral rectangular face is Nint1*Nint2, while the number !
! of integration points on the regular polygon basis is Nint2*Nint2.                !
! If miror = .true., the integration over a lateral rectangular face is performed   !
! over the half-face, and the number of integration points along the prism length   !
! is Nint1 / 2.  For a regular N-hedral finite prims with azimuthal symmetry, the   !
! top and bottom integration faces are triangles; the triangle basis is the side of !
! the regular polygon. The integration is performed by passing to the polar         !
! coordinates ro and phi, and the numbers of devision points for the integration over
! ro and phi are both equal to Nint2.                                               !
!                                                                                   !
! For convergence tests, the incident wave is assumed to be a plane wave traveling  !
! along the Z-axis of the global coordinate system and the scattering               !
! characteristics are computed in the azimuthal plane phi = 0°. The particle is     !
! placed in a general orientation (characterized by the Euler rotation angles       !
! alpha = beta = 45° and gamma = 0°), and interactive convergence tests over Nint,  !
! Nrank and Mrank are performed.                                                    !
!                                                                                   !
! For the integration test, the scattering problem is solved for the pairs          !
! (Nint1, Nint2) and (Nint1 + dNint1, Nint2 + dNint2), while Nrank and Mrank are    !
! keep constant. Note that no integration test will be performed if the particle    !
! geometry is generated from input file. For the convergence test over the expansion! 
! order, the scattering problem is solved for the pairs (Nrank, Mrank) and          !
! (Nrank - 1, Mrank), while for the azimuthal order test the cases (Nrank, Mrank)   !
! and (Nrank, Mrank - 1) are considered. The normalized differential scattering     !
! cross section will be checked at 20° increments for convergence within epsX       !
! (epsNint, epsNrank or epsMrank) tolerance. If the calculated results are converged! 
! within this tolerance at 80% of the scattering angles, then convergence is        !
! achieved, and the T matrix is stored for later use by other programs. The values  !
! of Nrank and Mrank are printed to the screen and to the T-matrix information file.!  
! These values together with the T matrix serve as INPUT PARAMETERS for other       !
! programs.                                                                         ! 
!                                                                                   !
! At each calculation, the extinction and scattering cross sections Cext and Cscat  !
! are also computed. Although the convergence of Cext and Cscat does not guarantee  !
! that the DSCS converges, the divergence of Cext and Cscat implies the divergence  !
! of the DSCS. The results of the convergence analysis are written to the file      !
! "/OUTPUTFILES/Output.dat". The above convergence tests (based on the analysis     !
! of the differential scattering cross section) have been proposed by Barber and    !
! Hill [P. W. Barber and S. C. Hill, Light scattering by particles: Computational   !
! Methods, World Scientific, Singapore, 1990].                                      !
!                                                                                   !
! 3. Estimates of Nint1, Nint2, Nrank and Mrank                                     !
! ----------------------------------------------                                    !
! The above convergence tests require estimates of Nint1, Nint2, Nrank and Mrank.   !
! These estimates must be supplied by the user, and for this purpose, Wiscombe's    !
! truncation limit criterion [W. J. Wiscombe, Improved Mie scattering algorithms,   !
! Applied Optics, 19, 1505-1509, 1980] can be used.                                 !
!                                                                                   !
! The truncation limit criterion proposed by Wiscombe provides the estimate         !
!                                                                                   !
!                NrankW = int(x + 4.05 * x**0.33 + 2),                              !
!                                                                                   !
! where x is the size parameter, x = k * Rcirc, and Rcirc is the radius of the      !
! smallest circumscribing sphere. For almost spherical particles, Mrank can be      !
! chosen as Nrank - 2,..,Nrank, while for less spherical particles Mrank can be     !
! smaller than Nrank - 2.                                                           !
!                                                                                   !
! The values of Nint1 and Nint2 depend on the the size parameter, the particle      !
! shape and the relative refractive index of the particle. For each input pair      !
! (Nint1,Nint2), the code computes some parameters of the discretized surface: the  !
! mean, maximum and minimum areas of the surface elements and the number of         !
! elements per unit surface (number density). Large values of the mean area (or     !
! small values of the number density) mean that the mesh (the discretization)       !
! is too rough. Large deviations of the maximum and minimum areas from the mean     !
! value show that the mesh is not uniform. In this case, the surface                !
! parametrization must be changed, i.e., the routines "interpolation_list3D" and    !
! "elem_geom3D" from the file "GeomLib.f90" have to be modified. The code also      !
! computes specific elements of the matrix Q31 for increasing Nint1- and            !
! Nint2-values, i.e., (Nint1 + 1,Nint2 + 1),...,(Nint1 + 4,Nint2 + 4). The          !
! computed quantities are:                                                          !
! - MAX {                                                                           ! 
!   A11(m,n,m1,n1) for m = Mrank, m1 = Mrank, n = Mrank, Nrank, n1 = Mrank, Nrank,  !
!   A22(m,n,m1,n1) for m = Mrank, m1 = Mrank, n = Mrank, Nrank, n1 = Mrank, Nrank } !
! - MIN {                                                                           !
!   A11(m,n,m1,n1) for m = Mrank, m1 = Mrank, n = Mrank, Nrank, n1 = Mrank, Nrank,  !
!   A22(m,n,m1,n1) for m = Mrank, m1 = Mrank, n = Mrank, Nrank, n1 = Mrank, Nrank } !  
! - MAX {                                                                           !
!   A12(m,n,m1,n1) for m = Mrank, m1 = Mrank, n = Mrank, Nrank, n1 = Mrank, Nrank,  ! 
!   A21(m,n,m1,n1) for m = Mrank, m1 = Mrank, n = Mrank, Nrank, n1 = Mrank, Nrank } !
! - MIN {                                                                           !
!   A12(m,n,m1,n1) for m = Mrank, m1 = Mrank, n = Mrank, Nrank, n1 = Mrank, Nrank,  !
!   A21(m,n,m1,n1) for m = Mrank, m1 = Mrank, n = Mrank, Nrank, n1 = Mrank, Nrank } !
! and                                                                               !
! - MAX {                                                                           !
!   A11(m,n,m1,n1) for m = 1, m1 = Mrank, n = 1, Nrank, n1 = Mrank, Nrank,          ! 
!   A22(m,n,m1,n1) for m = 1, m1 = Mrank, n = 1, Nrank, n1 = Mrank, Nrank }         !
! - MIN {                                                                           !
!   A11(m,n,m1,n1) for m = 1, m1 = Mrank, n = 1, Nrank, n1 = Mrank, Nrank,          !
!   A22(m,n,m1,n1) for m = 1, m1 = Mrank, n = 1, Nrank, n1 = Mrank, Nrank }         !
! - MAX {                                                                           ! 
!   A12(m,n,m1,n1) for m = 1, m1 = Mrank, n = 1, Nrank, n1 = Mrank, Nrank,          !
!   A21(m,n,m1,n1) for m = 1, m1 = Mrank, n = 1, Nrank, n1 = Mrank, Nrank }         !
! - MIN {                                                                           !
!   A12(m,n,m1,n1) for m = 1, m1 = Mrank, n = 1, Nrank, n1 = Mrank, Nrank,          ! 
!   A21(m,n,m1,n1) for m = 1, m1 = Mrank, n = 1, Nrank, n1 = Mrank, Nrank }         !
! where A11, A12, A21 and A22 are the four block matrices of Q31. The submatrices   !
! corresponding to m = Mrank and m1 = Mrank are diagonal block matrices,            !
! while the submatrices corresponding to m = 1 and m1 = Mrank are off-diagonal      !
! block matrices (of Aij, i = 1,2, j = 1,2). The variation of the computed          !
! quantities give an impression about the required numbers of integration points.   !
! Note that no convergence can be expected for small values of the matrix elements  !
! and the convergence over the number of integration points does not imply the      !
! convergence of the T-matrix calculation.                                          !      
!                                                                                   !
! 4. Strategy for Performing Convergence Tests                                      !
! ----------------------------------------------                                    !
! For nonaxisymmetric particles, a semi-convergent behaviour is usually attended:   !
! the errors decrease with increasing the maximum expansion order, attain a         !
! relative constant level and afterwards increase. The region of stability can be   !
! localized by using the following technique:                                       !
! 1. choose a value of Nrank close to the value predicted by Wiscombe's truncation  !
!    limit criterion (or smaller) and set Mrank = Nrank - 2 (or smaller);           !
! 2. find reliable Nint1- and Nit2-values such that the relative errors in          !
!    computing the matrix elements are sufficiently small;                          !
! 3. perform the convergence tests over Nrank and Mrank;                            ! 
! 4. if convergence is not achieved, set Nrank = Nrank + 1, maintain the relation   !
!    Mrank = Nrank - 2 and go to Step 2;                                            !
! 5. if convergence is achieved, perform several integration tests with the         !
!    initial Nint1- and Nint2-values found at step 2, until the DSCS converges.     !
!                                                                                   !
! For particles which are too extreme in size and/or aspect ratio, the errors of    !
! the extinction and scattering cross sections (instead of the differential         !
! scattering cross section) can be analyzed.                                        ! 
!                                                                                   !
! 5. Additional Comments                                                            !
! -----------------------                                                           !
! If the geometry of the particle is provided by the file FileFEM, only the         !
! integration test over Nrank and Mrank can be performed.                           !
!                                                                                   !
! The convergence tests can be switched off by setting the logical variable         !
! DoConvTest to false. In this case, the values of Nint1, Nint2 (if FileGeom = f),  !
! Nrank and Mrank must be specified in the input file.                              !
!                                                                                   !
! 6. Input Parameters                                                               !
! --------------------                                                              !
! The parameters specified in the input file "/INPUTFILES/InputNONAXSYM.dat" are    !
! listed below.                                                                     !
!                                                                                   !
! - wavelength (real) - wavelength of the incident light in vacuo.                  !
!                                                                                   !
! - ind_refMed (real) - refractive index of the ambient medium (nonabsorbing        !
!   dielectric medium).                                                             !
!                                                                                   !
! - ind_refRel (complex) - relative refractive index of the particle with respect   !
!   to the ambient medium. The imaginary part of the relative refractive index      !
!   must be zero for nonabsorbing particles and positive for absorbing particles.   !
!                                                                                   !  
! - perfectcond (logical) - if perfectcond = t, the particle is perfectly           !
!   conducting.                                                                     !
!                                                                                   !
! - anisotropic (logical) - if anisotropic = t, the particle is a uniaxial          !
!   anisotropic crystal.                                                            !  
!                                                                                   !
! - chiral (logical) - if chiral = t, the particle is optical active (chiral).      !
!                                                                                   !
!   The permissive values of the logical variables perfectcond and chiral are       !
!   summarized below.                                                               !
!                                                                                   !
!       particle              perfectcond            chiral       anisotropic       !
!      dielectric                  f                   f               f            !
!   perfectly conducting           t                   f               f            !
!      anisotropic                 f                   f               t            !
!        chiral                    f                   t               f            !
!                                                                                   !
! - kb (real) - parameter of chirality. This parameter is used for chiral           !
!   particles. The wave numbers for the left- and right-handed circularly           !
!   polarized waves are given by:                                                   !
!                kil = ki / (1.0 - kb)  and   kir = ki / (1.0 + kb)                 !
!   where ki is the wave number in the chiral medium.                               !
!                                                                                   !
! - FileGeom (logical) - if FileGeom = t, the particle geometry is supplied by the  !
!   input file FileFEM.                                                             ! 
!                                                                                   !
! - TypeGeom (integer) - parameter specifying the type of the particle geometry.    !
!                                                                                   !
! - FileFEM (character(80)) - name of the file containing the particle geometry.    !
!                                                                                   !
! - Nsurf (integer) - number of surface parameters.                                 !
!                                                                                   !
! - surf (real array: surf(1), surf(2),...,surf(Nsurf)) - surface parameters        !
!   specifying the shape of the particle. The dimension of the array surf is        !
!   NsurfPD. The integer parameter NsurfPD is specified in the routine              !
!   "Parameters.f90" and has the value NsurfPD = 10. If Nsurf > NsurfPD, the        !
!   execution is automatically terminated.                                          !
!                                                                                   !
! - Nparam (integer) - number of smooth surfaces forming the particle surface.      !
!                                                                                   !
!   The permissive values of the surface parameters (for the supported geometries)  !
!   are given below.                                                                !
!                                                                                   !
!   Particle      TypeGeom   Nsurf   Nparam                surf                     !
!   ellipsoid        1         3       1         surf(1) - length of the semi-      !
!                                                          axis along the x-axis    !         
!                                                surf(2) - length of the semi-      !
!                                                          axis along the y-axis    !
!                                                surf(3) - length of the semi-      !
!                                                          axis along the z-axis    !
!                                                                                   !  
!   quadratic        2         2       6         surf(1) - half-length of           !
!     prism                                                the prism                ! 
!                                                surf(2) - half-length of the       !
!                                                          square side              !
!                                                                                   !
!  regular N-hedral  3         2       2         surf(1) - half-length of           !
!      prism                       if miror = t;           the prism                !
!  (with azimuthal                     3         surf(2) - half-length of the       ! 
!   symmetry and                   if miror = f            side of the regular      !
!   OPTIONALLY with                                        basis                    !
!   mirror symmetry                                                                 !
!                                                                                   !
! NOTE: FOR ANISOTROPIC PARTICLES, A REGULAR N-HEDRAL PRISM IS NOT A VALID          !
! GEOMETRY.                                                                         !
!                                                                                   !
! - anorm (real) - characteristic length of the particle which is used to           !
!   normalize the differential scattering cross sections.                           !
!                                                                                   !
! - Rcirc (real) - characteristic length of the particle (usually the radius of     !
!   the smallest circumscribing sphere) which is used to compute an estimate of the !
!   maximum expansion order by using Wiscombe's truncation limit criterion (the     !
!   size parameter is x = k * Rcirc, where k is the wave number in the ambient      !
!   medium). Alternatively, Rcirc can be chosen as the equal-volume sphere radius   !
!   or the surface-equivalent-sphere radius. This parameter is used if the          !
!   interactive convergence tests are performed (DoConvTest = t).                   ! 
!                                                                                   !
! - miror (logical) - if miror = t, the particle is mirror symmetric (the plane of  !
!   symmetry or the plane of reflection is perpendicular to the axis of rotation).  !
!                                                                                   !
! NOTE: FOR CHIRAL AND ANISOTROPIC PARTICLES SET miror = f.                         !
!                                                                                   !
! - Nazimutsym (integer) - number of azimuthal symmetric sections. The T-matrix     !
!   calculation accounts for the azimuthal symmetry if Nazimutsym >= 2. If          !
!   Nazimutsym = 0, no azimuthal symmetry is considered.                            !
!                                                                                   !
! NOTE: FOR ANISOTROPIC PARTICLES SET Nazimutsym = 0.                               ! 
!                                                                                   !
!   For dielectric and perfectly conducting particles, the permissive values of     !
!   the variables miror and Nazimutsym are summarized below.                        !
!             TypeGeom           miror           Nazimutsym                         !  
!                1               t or f               0                             !
!                2                 f                  0                             !
!                3               t or f              >=3                            !
!   For chiral particles, the permissive values of the variables miror and          !
!   Nazimutsym are given below.                                                     !
!             TypeGeom           miror           Nazimutsym                         !  
!                1                 f                  0                             !
!                2                 f                  0                             !
!                3                 f                 >=3                            ! 
!                                                                                   !
! - DoConvTest (logical) - if DoConvTest = t, the interactive convergence tests     !
!   over Nint, Nrank and Mrank are invoked. An estimates of Nrank is given by       !
!   Wiscombe's truncation limit criterion. If DoConvTest = f, the values of Nint1,  !
!   Nint2, Nrank and Mrank must be supplied in the input file.                      !
!                                                                                   !
! - ExtThetaDom (logical) - if ExtThetaDom = t, the DSCS is computed for            !
!   scattering angles ranging from 0° to 180° in the azimuthal plane phiGS = 0°,    !
!   and from 180° to 0° in the azimuthal plane phiGS = 180°. The total number of    !
!   scattering angles is 10. If ExtThetaDom = f, the DSCS is computed for           !
!   scattering angles ranging from 0° to 180° in the azimuthal plane phiGS = 0°.    !
!                                                                                   !
! - ind_refRelZ (complex) - the second relative refractive index of the uniaxial    !
!   anisotropic particle. The imaginary part must be positive.                      !
!   If ind_refRelZ > ind_refRel, the crystal is positive anisotropic and if         !
!   ind_refRelZ < ind_refRel, the crystal is negative anisotropic.                  !    
!                                                                                   !
! - alphaPR, betaPR (real variables) - Euler angles specifying the orientation of   !
!   the principal coordinate system with respect to the particle coordinate system. !
!                                                                                   !
! - Nbeta (integer) - number of integration points for computing the vector         !
!   quasispherical wave functions.                                                  !
!                                                                                   !
! - Nrank (integer) - maximum expansion order. This parameter is used if the        !
!   convergence tests are not performed (DoConvTest = f).                           !
!                                                                                   !
! - Mrank (integer) - maximum azimuthal order. This parameter is used if the        !
!   convergence tests are not performed (DoConvTest = f).                           !
!                                                                                   !
! - Nint1, Nint2 (integer variables) - numbers of integration points in computing   !
!   integrals over the particle surface. These parameters are used if the           !
!   convergence tests are not performed (DoConvTest = f) and the particle geometry  !
!   is not supplied by the file FileFEM (FileGeom = f).                             !
!                                                                                   !
!   The following table explain the significance of the variables Nint1, Nint2,     !
!   Nrank and Mrank.                                                                !
!     DoConvTest    FileGeom    Nrank and Mrank   Nint1 and Nint2   TypeConvTest    !
!         t            f         console input     console input       1 or 2       !
!         t            t         console input      not required         2          !
!         f            f          file input         file input          0          ! 
!         f            t          file input        not required         0          !
!   The significance of the variable TypeConvTest is given below.                   !
!           TypeConvTest                    Significance                            !
!               0                        no convergence test                        !  
!               1                 convergence test over Nint1 and Nint2             ! 
!               2                 convergence test over Nrank and Mrank             !
!                                                                                   !
! - epsNint (real) - error tolerance for the integration test.                      !
!                                                                                   !
! - epsNrank (real) - error tolerance for the expansion order test.                 !
!                                                                                   !
! - epsMrank (real) - error tolerance for the azimuthal order test.                 !
!                                                                                   !
! - dNint1, dNint2 (integer variables) - numbers of division points for the         !
!   integration test. Note that the scattering problem is solved for the pairs      !
!   (Nint1,Nint2) and (Nint1 + dNint1,Nint2 + dNint2).                              ! 
!                                                                                   !
! - FileTmat (character(80)) - name of the file to which the T matrix is written.   !
!                                                                                   !
! - PrnProgress (logical) - if PrnProgress = t, the progress of calculation         !
!   is printed.                                                                     !
!                                                                                   !
! Note that all input parameters specifying lengths must be provided in the same    !
! units of lengths.                                                                 !
!                                                                                   !
! 7. Logical Scheme                                                                 !
! -------------------                                                               !  
! The organization of the code is as follows:                                       !
! 1. < read input data >                                                            !
! 2. < T-matrix calculation >                                                       !
!                                                                                   !
!    if ( do convergence test ) then                                                !
!                                                                                   !
!       if ( .not. read geometry from file) then                                    !        
!          < the code computes an estimate of Nrank by using                        !
!            Wiscombe's truncation criterion >                                      !
!          < the code prompts for the estimated values of Nrank and Mrank >         !
!          < the code prompts for the estimated values of Nint1 and Nint2,          !
!            computes the surface parameters for these input values and performs    ! 
!            a convergence test over some elements of the matrix Q31 >              !
!          < the code prompts for the type of convergence test:                     !
!            1 - over Nint, 2 - over Nrank and Mrank >                              !
!       else if ( read geometry from file) then                                     !
!          < the code computes an estimate of Nrank by using                        !
!            Wiscombe's truncation criterion >                                      !
!          < the code prompts for the estimated values of Nrank and Mrank >         !
!           type of convergence test = 2 (convergence test over Nrank and Mrank)    !
!       end if                                                                      !
!                                                                                   !
!    else if ( .not. do convergence test ) then                                     !
!                                                                                   !
!       if ( .not. read geometry from file) then                                    !        
!          < the input file provides the values of Nrank and Mrank >                !
!          < the input file provides the values of Nint1 and Nint2 >                !
!          < the code computes the surface parameters for Nint1 and Nint2           !
!            and performs a convergence test over some elements of the matrix Q31 > !
!          type of convergence test = 0 (no convergence test)                       !
!       else if ( read geometry from file) then                                     !
!          < the input file provides the values of Nrank and Mrank >                !
!          type of convergence test = 0 (no convergence test)                       !
!    end if                                                                         !
!    if ( do integration test, i.e., type of convergence test = 1 ) then            !
!       < the code computes the DSCS for (Nint1,Nint2) and                          !
!        (Nint1 + dNint1,Nint2 + dNint2) and write the results                      !
!         to the file "Output.dat" >                                                !
!       < the execution is terminated >                                             !
!    end if                                                                         !
!    if ( do expansion and azimuthal order tests,                                   !
!         i.e., type of convergence test = 2 ) then                                 !
!       < the code computes the DSCS for (Nrank,Mrank), (Nrank - 1,Mrank) and       ! 
!         (Nrank, Mrank - 1) and write the results to the file "Output.dat" >       !
!       < the T matrix is stored in the file FileTmat and the T-matrix information  !
!         file is created >                                                         !
!       < Mrank and Nrank are written to the screen and to the T-matrix             !
!         information file >                                                        !
!       < go to Step 3 >                                                            !
!    end if                                                                         !
!    if ( .not. do convergence test, i.e., type of convergence test = 0 ) then      !
!       < the T matrix is computed and stored in the file FileTmat >                !
!       < the T-matrix information file is created >                                ! 
!       < Mrank and Nrank are written to the screen and to the T-matrix             !
!         information file >                                                        !
!       < go to Step 3 >                                                            !
!    end if                                                                         !
! 3. < scattering characteristics calculation >                                     !    
!------------------------------------------------------------------------------------
  use parameters  
  implicit none 
  integer       :: TypeGeom, Nsurf, Nface, Nparam, Nazimutsym, TypeConvTest, Mrank, &
                   Nrank, Nbeta, Nint1, Nint2, dNint1, dNint2                     
  real(O)       :: k, ind_refMed, wavelength, anorm, surf(NsurfPD), snorm,          &
                   kb, epsNint, epsNrank, epsMrank, Rcirc, rp(3,NfacePD),           &
                   np(3,NfacePD), area(NfacePD), alphaPR, betaPR, gammaPR
  complex(O)    :: ind_refRel, ind_RefRelZ
  logical       :: FileGeom, miror, perfectcond, chiral, anisotropic, DoConvTest,   &
                   ExtThetaDom, PrnProgress
  character(80) :: FileTmat, FileFEM  
! -----------------------------------------------------------------------------------
!                              Read the input file                                  ! 
! -----------------------------------------------------------------------------------  
  call readinputNONAXSYM (wavelength, ind_refMed, ind_refRel, ind_RefRelZ,          &
       alphaPR, betaPR, perfectcond, anisotropic, chiral, kb, FileGeom,             &
       TypeGeom, FileFEM, Nsurf, surf, Nparam, anorm, Rcirc, miror, Nazimutsym,     &   
       DoConvTest, ExtThetaDom, Nbeta, Nint1, Nint2, Nrank, Mrank, epsNint,         &
       epsNrank, epsMrank, dNint1, dNint2, FileTmat, PrnProgress, k, gammaPR,       &
       snorm, rp, np, area, Nface, TypeConvTest )    
! -----------------------------------------------------------------------------------
!                                      Main                                         !
! -----------------------------------------------------------------------------------                            
  if (DoConvTest) then
    open (unit = iOutput, file = FileOutput, status = "replace")
    call printinputNONAXSYM (FileGeom, TypeGeom, FileFEM, Nsurf, Nparam, Nazimutsym,&
         dNint1, dNint2, ind_refMed, wavelength, anorm, Rcirc, surf, kb, epsNint,   &
         epsNrank, epsMrank,  ind_refRel, ind_refRelZ, alphaPR, betaPR, anisotropic,&
         miror, perfectcond, chiral) 
    if (TypeConvTest == 1) then  
      if (.not. anisotropic) then             
        call convergence_NintNONAXSYM (FileGeom, TypeGeom, k, ind_refRel, snorm,    &
             Nsurf, surf, rp, np, area, Nface, Nparam, Mrank, Nrank, Nint1, Nint2,  &
             dNint1, dNint2, miror, Nazimutsym, perfectcond, chiral, kb, epsNint,   &
             ExtThetaDom, PrnProgress)      
      else
        call convergence_NintAnis (FileGeom, TypeGeom, k, ind_refRel, ind_refRelZ,  &
             alphaPR, betaPR, gammaPR, snorm, Nsurf, surf, rp, np, area, Nface,     &
             Nparam, Mrank, Nrank, Nbeta, Nint1, Nint2, dNint1, dNint2, epsNint,    &
             ExtThetaDom, PrnProgress) 
      end if
    else if (TypeConvTest == 2) then      
      if (.not. anisotropic) then
        call convergence_Nrank_MrankNONAXSYM (FileGeom, TypeGeom, k, ind_refRel,    &
             snorm, Nsurf, surf, rp, np, area, Nface, Nparam, Mrank, Nrank, Nint1,  &
             Nint2, miror, Nazimutsym, perfectcond, chiral, kb, epsNrank, epsMrank, &
             ExtThetaDom, FileTmat, PrnProgress)             
      else
        call convergence_Nrank_MrankAnis (FileGeom, TypeGeom, k, ind_refRel,        &
             ind_refRelZ, alphaPR, betaPR, gammaPR, snorm, Nsurf, surf, rp, np,     &
             area, Nface, Nparam, Mrank, Nrank, Nbeta, Nint1, Nint2, epsNrank,      &
             epsMrank, ExtThetaDom, FileTmat, PrnProgress)  
        end if
    end if
    close (unit = iOutput)           
  else       
    if (.not. anisotropic) then          
      call TMatrix_Nrank_MrankNONAXSYM (FileGeom, TypeGeom, k, ind_refRel, Nsurf,   &
           surf, rp, np, area, Nface, Nparam, Mrank, Nrank, Nint1, Nint2, miror,    &
           Nazimutsym, perfectcond, chiral, kb, FileTmat, PrnProgress)                                                     
    else
      call TMatrix_Nrank_MrankAnis (FileGeom, TypeGeom, k, ind_refRel, ind_refRelZ, &
           alphaPR, betaPR, gammaPR, Nsurf, surf, rp, np, area, Nface, Nparam,      &
           Mrank, Nrank, Nbeta, Nint1, Nint2, FileTmat, PrnProgress) 
    end if 
  end if    
end subroutine TNONAXSYM
!***********************************************************************************
subroutine readinputNONAXSYM (wavelength, ind_refMed, ind_refRel, ind_RefRelZ,      &
           alphaPR, betaPR, perfectcond, anisotropic, chiral, kb, FileGeom,         &
           TypeGeom, FileFEM, Nsurf, surf, Nparam, anorm, Rcirc, miror, Nazimutsym, &   
           DoConvTest, ExtThetaDom, Nbeta, Nint1, Nint2, Nrank, Mrank, epsNint,     &
           epsNrank, epsMrank, dNint1, dNint2, FileTmat, PrnProgress, k, gammaPR,   &
           snorm, rp, np, area, Nface, TypeConvTest )
  use parameters
  use derived_parameters
  implicit none 
  integer       :: TypeGeom, Nsurf, Nface, Nparam, Nazimutsym, TypeConvTest, Mrank, &
                   Nrank, NrankW, Nbeta, Nint1, Nint2, dNint1, dNint2, i, j, ios,   &
                   icall                     
  real(O)       :: k, ind_refMed, wavelength, anorm, surf(NsurfPD), xpart, snorm,   &
                   kb, epsNint, epsNrank, epsMrank, Rcirc, x, rp(3,NfacePD),        &
                   np(3,NfacePD), area(NfacePD), dp, alphaPR, betaPR, gammaPR, grd
  complex(O)    :: ind_refRel, ind_RefRelZ
  logical       :: FileGeom, miror, perfectcond, chiral, anisotropic, DoConvTest,   &
                   ExtThetaDom, PrnProgress, more, continuare, IntTest, XFindPar
  character(80) :: FileTmat, FileFEM, string    
! -----------------------------------------------------------------------------------
!                        Read the input file FileInputNONAXSYM                      ! 
! ----------------------------------------------------------------------------------- 
  call DrvParameters   
  open (unit = iInputNONAXSYM, file = FileInputNONAXSYM, status = "old",            &
        position = "rewind")   
  wavelength  = 0.1_O * 2._O * Pi
  ind_refMed  = 1._O
  ind_refRel  = (1.5_O,0._O)  
  string  = 'OptProp'
  if (XFindPar (iInputNONAXSYM, string)) then
    read (iInputNONAXSYM, *, iostat = ios) wavelength
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable wavelength;')"
      stop
    end if
    read (iInputNONAXSYM, *, iostat = ios) ind_refMed
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refMed;')"
      stop
    end if    
    read (iInputNONAXSYM, *, iostat = ios) ind_refRel
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refRel;')"
      stop
    end if        
  else
    print "(/,2x,'Group name OptProp not found;')"
    stop  
  end if                          
  k   = 2._O * Pi * ind_refMed / wavelength
  grd = Pi / 180._O 
!
  perfectcond = .false.
  anisotropic = .false.  
  chiral      = .false.
  kb     = 0._O
  string = 'MatProp'    
  if (XFindPar (iInputNONAXSYM, string)) then
    read (iInputNONAXSYM, *, iostat = ios) perfectcond
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable perfectcond;')"
      stop
    end if
    read (iInputNONAXSYM, *, iostat = ios) anisotropic
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable anisotropic;')"
      stop
    end if    
    read (iInputNONAXSYM, *, iostat = ios) chiral
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable chiral;')"
      stop
    end if    
    read (iInputNONAXSYM, *, iostat = ios) kb
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable kb;')"
      stop
    end if              
  else
    print "(/,2x,'Group name MatProp not found;')"
    stop  
  end if   
  call check_MatPropNONAXSYM (perfectcond, anisotropic, chiral, kb) 
  if (chiral) call check_chirality (kb) 
!
  FileGeom = .false.  
  FileFEM  = ' '
  TypeGeom = 1  
  Nsurf = 3
  do i = 1, NsurfPD
    surf(i) = 1._O
  end do
  Nparam = 1
  anorm  = 1._O
  Rcirc  = 1._O
  miror  = .true.
  Nazimutsym = 0
  string = 'GeomProp'    
  if (XFindPar (iInputNONAXSYM, string)) then
    read (iInputNONAXSYM, *, iostat = ios) FileGeom
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable FileGeom;')"
      stop
    end if      
    read (iInputNONAXSYM, *, iostat = ios) FileFEM
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable FileFEM;')"
      stop
    end if
    read (iInputNONAXSYM, *, iostat = ios) TypeGeom
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable TypeGeom;')"
      stop
    end if          
    read (iInputNONAXSYM, *, iostat = ios) Nsurf
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nsurf;')"
      stop
    end if
    if (Nsurf > NsurfPD) then
      print "(/,2x,'Input error: Nsurf exceeds NsurfPD;')"                                    
      stop
    end if
    do i = 1, Nsurf
      read (iInputNONAXSYM, *, iostat = ios) surf(i)
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable surf;')"
        stop
      end if
    end do 
    read (iInputNONAXSYM, *, iostat = ios) Nparam
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nparam;')"
      stop
    end if
    read (iInputNONAXSYM, *, iostat = ios) anorm
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable anorm;')"
      stop
    end if
    read (iInputNONAXSYM, *, iostat = ios) Rcirc
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Rcirc;')"
      stop
    end if
    read (iInputNONAXSYM, *, iostat = ios) miror
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable miror;')"
      stop
    end if
    read (iInputNONAXSYM, *, iostat = ios) Nazimutsym
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nazimutsym;')"
      stop
    end if
  else
    print "(/,2x,'Group name GeomProp not found;')"
    stop  
  end if
  call check_inputNONAXSYM (FileGeom, FileFEM, miror, chiral, anisotropic,          &
       Nazimutsym)
  call check_geom3D (TypeGeom, Nsurf, Nparam, anisotropic, miror, Nazimutsym)    
  call check_anorm (anorm)
  xpart = k * anorm
  snorm = Pi * xpart * xpart 
  Nface = 1
  do i = 1, NfacePD
    do j = 1, 3
      rp(j,i) = 0._O
      np(j,i) = 0._O
    end do
    area(i) = 0._O
  end do 
  if (FileGeom) then    
    call read_FileFEM (FileFEM, Nface, rp, np, area) 
    Rcirc = 0._O
    do i = 1, Nface
      dp = sqrt(rp(1,i)**2 + rp(2,i)**2 + rp(3,i)**2)
      if (dp > Rcirc) Rcirc = dp
    end do            
  end if                    
!
  DoConvTest   = .true.
  ExtThetaDom  = .true.
  string = 'ConvTest'    
  if (XFindPar (iInputNONAXSYM, string)) then
    read (iInputNONAXSYM, *, iostat = ios) DoConvTest
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable DoConvTest;')"
      stop
    end if
    read (iInputNONAXSYM, *, iostat = ios) ExtThetaDom
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ExtThetaDom;')"
      stop
    end if                     
  else
    print "(/,2x,'Group name ConvTest not found;')"
    stop  
  end if     
!
  ind_refRelZ = (1.5_O,0._O)
  alphaPR = 0._O
  betaPR  = 0._O
  Nbeta   = 60
  if (anisotropic) then
    string = 'AnSVWF'    
    if (XFindPar (iInputNONAXSYM, string)) then
      read (iInputNONAXSYM, *, iostat = ios) ind_refRelZ
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable ind_refRelZ;')"
        stop
      end if    
      read (iInputNONAXSYM, *, iostat = ios) alphaPR
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable alphaPR;')"
        stop
      end if    
      read (iInputNONAXSYM, *, iostat = ios) betaPR
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable betaPR;')"
        stop
      end if            
      read (iInputNONAXSYM, *, iostat = ios) Nbeta
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Nbeta;')"
        stop
      end if
    else
      print "(/,2x,'Group name AnSVWF not found;')"
      stop  
    end if      
  end if
  if (anisotropic) call check_ind_ref3 (ind_refRel, ind_refRelZ)  
  alphaPR = alphaPR * grd
  betaPR  = betaPR  * grd
  gammaPR = 0._O  
!  
  if (DoConvTest) then
    if (.not. FileGeom) then
      print "(/,2x,'Convergence Test for a Nonaxisymmetric Particle')"
      print "(  2x,'-----------------------------------------------')"      
    else
      print "(/,2x, a)",                                                            &
     'Convergence Test for a Nonaxisymmetric Particle over Nrank and Mrank'
      print "(  2x, a)",                                                            &
     '--------------------------------------------------------------------'             
    end if      
  else
    print "(/,2x,'T-Matrix Computation for a Nonaxisymmetric Particle')"
    print "(  2x,'---------------------------------------------------')"
  end if
!  
  x = k * Rcirc
  NrankW = int(x + 4.05_O * x**0.33_O + 2._O)
  Nrank  = 16
  Mrank  = 14
  if (.not. DoConvTest) then
    string = 'NrankMrank'    
    if (XFindPar (iInputNONAXSYM, string)) then
      read (iInputNONAXSYM, *, iostat = ios) Nrank
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Nrank;')"
        stop
      end if
      read (iInputNONAXSYM, *, iostat = ios) Mrank
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Mrank;')"
        stop
      end if
    else
      print "(/,2x,'Group name NrankMrank not found;')"
      stop  
    end if
    print "(/,2x,'Nrank and Mrank input values:')"
    print "(  2x,'the input values of Nrank and Mrank are ', i3,' and ', i3, ',')",  &
              Nrank, Mrank
    print "(  2x, a, i3, a)",                                                        &
   'while the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';' 
  else  
    print "(/,2x,'Nrank estimate:')"  
    print "(  2x, a, i3, a)" ,                                                       &
   'the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';' 
    print "(/,2x,'- enter the estimated values of Nrank and Mrank,')"
    print "(  2x, a)",                                                               &
   '  where for almost spherical particles Mrank = Nrank - 2,...,Nrank,'
    print "(  2x, a)",                                                               &
   '  while for less spherical particles Mrank can be smaller than Nrank - 2;'
    call read_integer2 (Nrank, Mrank)                                                               
  end if 
  call check_MrankNrank (Mrank, Nrank)  
! 
  Nint1 = 12
  Nint2 = 12  
  if (.not. FileGeom) then
    if (DoConvTest) then
      more  = .true.
      icall = 0
      do while (more)
        icall = icall + 1     
        print "(/,2x,'- enter the estimated values of Nint1 and Nint2;')"
        call read_integer2 (Nint1, Nint2)
        if (Nint1 < 5) Nint1 = 5
        if (Nint2 < 5) Nint2 = 5
        if (.not. anisotropic) then          
          call SurfaceElemNint (k, TypeGeom, Nparam, Nsurf, surf, miror,            &
               Nazimutsym, Nint1, Nint2)
          if (icall == 1) then
            print "(/,2x, a)",                                                      &
           '- enter true to perform an integration test over specific matrix '	                             
            print "(  2x,'elements and false otherwise;')"
            call read_logical (IntTest)
          end if           
          if (IntTest) then            
            call ConvergenceMatrixElem (3, 1, TypeGeom, Nsurf, Nparam, Mrank,       &
                 Mrank, Nrank, 5, Nint1, Nint2, Nazimutsym, k, surf, kb,            &
                 ind_refRel, miror, perfectcond, chiral)
            call ConvergenceMatrixElem (3, 1, TypeGeom, Nsurf, Nparam, 1, Mrank,    &
                 Nrank, 5, Nint1, Nint2, Nazimutsym, k, surf, kb, ind_refRel,       &
                 miror, perfectcond, chiral)
          end if
        else
          call SurfaceElemNint (k, TypeGeom, Nparam, Nsurf, surf, .false., 0,       &
          Nint1, Nint2)
          if (icall == 1) then
            print "(/,2x, a)",                                                      &
           '- enter true to perform an integration test over specific matrix '	                      
            print "(  2x,'elements and false otherwise;')"
            call read_logical (IntTest)
          end if           
          if (IntTest) then                  
            call ConvergenceMatrixElemAnis (TypeGeom, 3, 1, k, ind_refRel,          &
                 ind_refRelZ, alphaPR, betaPR, gammaPR, Nsurf, surf, Mrank, Mrank,  &
                 Nrank, 5, Nint1, Nint2, Nbeta, Nparam)
            call ConvergenceMatrixElemAnis (TypeGeom, 3, 1, k, ind_refRel,          &
                 ind_refRelZ, alphaPR, betaPR, gammaPR, Nsurf, surf, 1, Mrank,      &
                 Nrank, 5, Nint1, Nint2, Nbeta, Nparam)     
          end if
        end if 
        print "(2x,'- enter true for a new input of number of integration points')"
        print "(2x,'or false to continue;')"        
        call read_logical (continuare)
        if (.not. continuare) more = .false.
      end do
    else
      string = 'Nint'    
      if (XFindPar (iInputNONAXSYM, string)) then
        read (iInputNONAXSYM, *, iostat = ios) Nint1
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable Nint1;')"
          stop
        end if
        read (iInputNONAXSYM, *, iostat = ios) Nint2
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable Nint2;')"
          stop
        end if
      else
        print "(/,2x,'Group name Nint not found;')"
        stop  
      end if
      if (Nint1 < 5) Nint1 = 5
      if (Nint2 < 5) Nint2 = 5
      print "(/,2x,'Nint1 and Nint2 input values:')"
      print "(  2x,'the input values of Nint1 and Nint2 are ', i4,' and',i4,';')",  &
                Nint1, Nint2
      if (.not. anisotropic) then
        call SurfaceElemNint (k, TypeGeom, Nparam, Nsurf, surf, miror, Nazimutsym,  &
             Nint1, Nint2)        
        print "(/,2x,'- enter true to perform an integration test over')"                  
        print "(  2x,'specific matrix elements and false otherwise;')"
        call read_logical (IntTest)                        
        if (IntTest) then       
          call ConvergenceMatrixElem (3, 1, TypeGeom, Nsurf, Nparam, Mrank, Mrank,  &
               Nrank, 5, Nint1, Nint2, Nazimutsym, k, surf, kb, ind_refRel, miror,  &
               perfectcond, chiral )
          call ConvergenceMatrixElem (3, 1, TypeGeom, Nsurf, Nparam, 1, Mrank,      &
               Nrank, 5, Nint1, Nint2, Nazimutsym, k, surf, kb, ind_refRel, miror,  &
               perfectcond, chiral )                                                         
        end if
      else
        call SurfaceElemNint (k, TypeGeom, Nparam, Nsurf, surf, .false., 0,         &
             Nint1, Nint2)
        print "(/,2x,'- enter true to perform an integration test over')"                  
        print "(  2x,'specific matrix elements and false otherwise;')"
        call read_logical (IntTest)                        
        if (IntTest) then                          
          call ConvergenceMatrixElemAnis (TypeGeom, 3, 1, k, ind_refRel,            &
               ind_refRelZ, alphaPR, betaPR, gammaPR, Nsurf, surf, Mrank, Mrank,    &
               Nrank, 5, Nint1, Nint2, Nbeta, Nparam)
          call ConvergenceMatrixElemAnis (TypeGeom, 3, 1, k, ind_refRel,            &
               ind_refRelZ, alphaPR, betaPR, gammaPR, Nsurf, surf, 1, Mrank,        &
               Nrank, 5, Nint1, Nint2, Nbeta, Nparam)       
        end if
      end if
    end if
  end if
!  
  if (DoConvTest) then   
    if (.not. FileGeom) then   
      print "(/,2x, a)",                                                            &
     '- enter the type of convergence test: 1 - Nint, 2 - Nrank and Mrank;'
      call read_integerbound (TypeConvTest, 1, 2)              
    else    
      TypeConvTest = 2 
    end if  
  else
    TypeConvTest = 0
  end if
!         
  epsNint  = 5.e-2_O
  epsNrank = 5.e-2_O
  epsMrank = 5.e-2_O
  dNint1 = 4
  dNint2 = 4
  string   = 'Errors'
  if (XFindPar (iInputNONAXSYM, string)) then
    read (iInputNONAXSYM, *, iostat = ios) epsNint
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsNint;')"
      stop
    end if
    read (iInputNONAXSYM, *, iostat = ios) epsNrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsNrank;')"
      stop
    end if
    read (iInputNONAXSYM, *, iostat = ios) epsMrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsMrank;')"
      stop
    end if 
    read (iInputNONAXSYM, *, iostat = ios) dNint1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable dNint1;')"
      stop
    end if
    read (iInputNONAXSYM, *, iostat = ios) dNint2
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable dNint2;')"
      stop
    end if         
  else
    print "(/,2x,'Group name Errors not found;')"
    stop  
  end if  
!
  FileTmat = '../TMATFILES/T.dat'
  string   = 'Tmat' 
  if (XFindPar (iInputNONAXSYM, string)) then
    read (iInputNONAXSYM, *, iostat = ios) FileTmat
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable FileTmat;')"
      stop
    end if             
  else
    print "(/,2x,'Group name Tmat not found;')"
    stop  
  end if   
!
  PrnProgress = .true.
  string   = 'PrintProgress' 
  if (XFindPar (iInputNONAXSYM, string)) then
    read (iInputNONAXSYM, *, iostat = ios) PrnProgress
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable PrnProgress;')"
      stop
    end if             
  else
    print "(/,2x,'Group name PrintProgress not found;')"
    stop  
  end if     
  close (unit = iInputNONAXSYM)  
end subroutine readinputNONAXSYM                            
!***********************************************************************************
subroutine printinputNONAXSYM (FileGeom, TypeGeom, FileFEM, Nsurf, Nparam,          &
           Nazimutsym, dNint1, dNint2, ind_refMed, wavelength, anorm, Rcirc, surf,  &
           kb, epsNint, epsNrank, epsMrank, ind_refRel, ind_refRelZ, alphaPR,       &
           betaPR, anisotropic, miror, perfectcond, chiral)
  use parameters
  implicit none
  integer        :: TypeGeom, Nsurf, Nparam, Nazimutsym, dNint1, dNint2, i, LenString                     
  real(O)        :: ind_refMed, wavelength, anorm, surf(Nsurf), kb, epsNint,        &
                    epsNrank, epsMrank, Rcirc, alphaPR, betaPR
  complex(O)     :: ind_refRel, ind_refRelZ
  character(80)  :: FileFEM, FileFEMWrite
  logical        :: FileGeom, miror, perfectcond, chiral, anisotropic
!
  write (iOutput,"(/,2x,'Input Parameters of the Scattering Problem:',/)")
  write (iOutput,"(2x,'wavelength of the free space, wavelength = ',1pe13.4,';')")  &
         wavelength
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'refractive index of the ambient medium,ind_refMed = ', ind_refMed, ';'
  if (.not. anisotropic) then
    write (iOutput,"(2x, a, 1pe10.3, ',', 1pe10.3, a)")                             &
   'relative refractive index of the particle, ind_refRel = (', ind_refRel, ');'
  else
    write (iOutput,"(2x,'relative refractive indices of the uniaxial particle:')")
    write (iOutput,"(2x, a, 1pe10.3, ',', 1pe10.3, a, 1pe10.3, ',', 1pe10.3, a)")   &
   'ind_refRel = (', ind_refRel, '), ind_refRelZ = (', ind_refRelZ, ');' 
  end if
  write (iOutput,*)
  if (FileGeom) then
    FileFEMWrite = FileFEM(14:LenString(FileFEM))
    write (iOutput,"(2x, a, a30)")                                                  &
   'name of the file containing the particle geometry, FileFEM = ',     FileFEMWrite
  else
    write (iOutput,"(2x,'index of the particle geometry, TypeGeom = ',i2,';')")     &
           TypeGeom
    if (TypeGeom == 1) then
      write (iOutput,"(2x,'ellipsoid;')")       
    else if (TypeGeom == 2) then
      write (iOutput,"(2x,'quadratic prism;')") 
    else if (TypeGeom == 3) then
      write (iOutput,"(2x,'regular N-hedral prism;')")
    end if
    write (iOutput,"(2x,'number of surface parameters, Nsurf = ',i2,';')") Nsurf
    write (iOutput,"(2x,'surface parameters:')")
    do i = 1, Nsurf
      write (iOutput,"(2x,'surf(',i2,') = ',1pe10.3,',')") i, surf(i)
    end do       
    if (TypeGeom == 1) then
      write (iOutput,"(2x,'where surf(1) is the semi-axis along the x-axis,')")
      write (iOutput,"(2x,'surf(2) is the semi-axis along the y-axis and')")
      write (iOutput,"(2x,'surf(3) is the semi-axis along the z-axis;')")
    else if (TypeGeom == 2 .or. TypeGeom == 3) then
      write (iOutput,"(2x,'where surf(1) is the half-length of the prism,')")
      write (iOutput,"(2x,'and   surf(2) is the half-length of the basis side;')")
    end if 
    write (iOutput,"(2x,'number of integration surfaces, Nparam = ',i2,';')")       &
           Nparam
  end if
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'characteristic length of the particle, anorm = ', anorm, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'radius of the circumscribing sphere, Rcirc = ', Rcirc, ';'
  if (anisotropic) then
    write (iOutput,"(2x, a)")                                                       &
   'Euler angles specifying the orientation of the principal coordinate system:'
    write (iOutput,"(2x, a, 1pe10.3, a, 1pe10.3, a, 1pe10.3, a)")                   &
   'alpha = ', alphaPR * 180._O / Pi, ', beta = ', betaPR * 180._O / Pi,            &
   ', gamma = ', 0._O, ';'        
  end if
  if (miror) write (iOutput,"(2x,'mirror symmetric particle;')")
  write (iOutput,"(2x, a, i2, a)")                                                  &
 'number of azimuthal symmetric sections, Nazimutsym = ', Nazimutsym, ';'
  if (Nazimutsym >=2) then
    write (iOutput,"(2x,'azimuthal symmetry is taken into account;')")
  else
    write (iOutput,"(2x,'azimuthal symmetry is disregarded;')")
  end if
  write (iOutput,*)
  if (perfectcond) then
    write (iOutput,"(2x,'perfectly conducting particle;')")
  else if (chiral) then
    write (iOutput,"(2x,'chiral particle;')")
    write (iOutput,"(2x,'characteristic of chirality, kb = ',1pe10.3,';')") kb  
  else if (anisotropic) then
    write (iOutput,"(2x,'uniaxial anisotropic particle;')")
  else 
    write (iOutput,"(2x,'dielectric particle;')")
  end if 
  write (iOutput,*)    
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum expansion order tolerance, epsNrank = ', epsNrank, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum azimuthal order tolerance, epsMrank = ', epsMrank, ';'
  if (.not. FileGeom) then               
    write (iOutput,"(2x,'integration tolerance, epsNint = ',1pe10.3,';')") epsNint                         
    write (iOutput,"(2x,'integration steps, dNint1 = ',i4,', dNint2 = ',i4,'.')")   &
           dNint1, dNint2 
  end if
  write (iOutput,"(/)")                       
end subroutine printinputNONAXSYM
!***********************************************************************************
subroutine convergence_NintNONAXSYM (FileGeom, TypeGeom, k, ind_ref, snorm, Nsurf,  &
           surf, rp, np, area, Nface, Nparam, Mrank, Nrank, Nint1, Nint2, dNint1,   &
           dNint2, miror, Nazimutsym, perfectcond, chiral, kb, epsNint, ExtThetaDom,&
           PrnProgress)
  use parameters
  implicit none
  integer    :: TypeGeom, Nsurf, Nface, Nparam, Mrank, Nrank, Nint1, Nint2, dNint1, &
                dNint2, Nazimutsym                  
  real(O)    :: k, snorm, surf(Nsurf), rp(3,NfacePD), np(3,NfacePD), area(NfacePD), &
                kb, epsNint
  complex(O) :: ind_ref
  logical    :: FileGeom, miror, perfectcond, chiral, ExtThetaDom, PrnProgress
!      
  integer    :: Nmax, Nteta, i, NthetaConv, iNint, NintAL   
  real(O)    :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,        &
                Cext, Qext
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG1(:,:), paramG2(:,:), weightsG(:,:), h(:), v(:),  &
                            oldh(:), oldv(:)
  complex(O),allocatable :: a(:,:), c(:), b(:,:), c1(:)   
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = Pi / 4._O
  beta   = Pi / 4._O
  gama   = 0._O
  alfap  = Pi / 4._O
  Nmax   = Nrank + Mrank * (2 * Nrank - Mrank + 1) 
  call write_TypeConvHead (1)
  allocate (a(2*Nmax,2*Nmax), c(2*Nmax), b(2*Nmax,2*Nmax), c1(2*Nmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do
  if (PrnProgress) call write_progress (.true., 1, 7)
  do iNint = 1, 2
    NintAL = max(Nint1,Nint2)                            
    allocate (paramG1(Nparam,NintAL*NintAL), paramG2(Nparam,NintAL*NintAL),         &
              weightsG(Nparam,NintAL*NintAL))
    allocate (Nintparam(Nparam))
    call interpolation_list3D (TypeGeom, Nsurf, surf, Nint1, Nint2, NintAL, Nparam, &
         Nintparam, paramG1, paramG2, weightsG, miror, Nazimutsym) 
    call matrix_Q (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np, area, &
         Nface, Mrank, Nrank, Nmax, NintAL, Nparam, Nintparam, paramG1, paramG2,    &
         weightsG, miror, Nazimutsym, perfectcond, chiral, kb, a, Nmax, Nmax)  
    if (PrnProgress) call write_progress (.false., 2+3*(iNint-1), 7)
    call matrix_Q (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np, area, &
         Nface, Mrank, Nrank, Nmax, NintAL, Nparam, Nintparam, paramG1, paramG2,    &
         weightsG, miror, Nazimutsym, perfectcond, chiral, kb, b, Nmax, Nmax)
    if (PrnProgress) call write_progress (.false., 3+3*(iNint-1), 7)
    call LU_SYSTEM (a, 2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, 2*Nmax)
    if (PrnProgress) call write_progress (.false., 4+3*(iNint-1), 7)
    call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrank, Nrank,   &
         Nmax, c)    
    call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, c, c1)
    call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfa, beta, gama, k, snorm,    &
         ExtThetaDom,.true., h, v)
    call CQscat (c1, Mrank, Nrank, Nmax, k, snorm, Cscat, Qscat)
    call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfa, beta, gama,            &
         alfap, k, snorm, Cext, Qext)     
    call delta_DSCS (Nteta, h, v, oldh, oldv, epsNint, NthetaConv)  
    call write_4ConvParam (Nint1, Nint2, Nrank, Mrank)
    call write_DSCS (Nteta, ExtThetaDom, h, v)
    call write_Effic (Qscat, Qext)
    Nint1 = Nint1 + dNint1
    Nint2 = Nint2 + dNint2
    deallocate (paramG1, paramG2, weightsG, Nintparam)
  end do
  call write_NintConvRes (NthetaConv, Nteta, epsNint)
  deallocate (a, b, c, c1, h, v, oldh, oldv)    
end subroutine convergence_NintNONAXSYM 
!***********************************************************************************
subroutine convergence_Nrank_MrankNONAXSYM (FileGeom, TypeGeom, k, ind_ref, snorm,  &
           Nsurf, surf, rp, np, area, Nface, Nparam, Mrank, Nrank, Nint1, Nint2,    &
           miror, Nazimutsym, perfectcond, chiral, kb, epsNrank, epsMrank,          &
           ExtThetaDom, FileTmat, PrnProgress)
  use parameters
  implicit none
  integer       :: TypeGeom, Nsurf, Nface, Nparam, Mrank, Nrank, Nint1, Nint2,      &
                   Nazimutsym                   
  real(O)       :: k, snorm, surf(Nsurf), rp(3,NfacePD), np(3,NfacePD),             &
                   area(NfacePD), kb, epsNrank, epsMrank
  complex(O)    :: ind_ref
  logical       :: FileGeom, miror, perfectcond, chiral, ExtThetaDom, PrnProgress 
  character(80) :: FileTmat
!      
  integer       :: Nmax, Nteta, i, j, NthetaConvN, NthetaConvM, NintAL
  real(O)       :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,     &
                   Cext, Qext
  integer,allocatable    :: Nintparam(:) 
  real(O),allocatable    :: paramG1(:,:), paramG2(:,:), weightsG(:,:), h(:), v(:),  &
                            oldh(:), oldv(:), oldh0(:), oldv0(:)
  complex(O),allocatable :: a(:,:), c(:), b(:,:), c1(:), a1(:,:), b1(:,:)   
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = Pi / 4._O
  beta   = Pi / 4._O
  gama   = 0._O
  alfap  = Pi / 4._O
  Nmax   = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  NintAL = max(Nint1,Nint2) 
  call write_TypeConvHead (4)  
  open (unit = iTmat, file = FileTmat, status = 'replace')
  call write_HeadFileTmat (Nmax, Nmax)     
  allocate (a1(2*Nmax,2*Nmax), b1(2*Nmax,2*Nmax))  
  allocate (a(2*Nmax,2*Nmax), c(2*Nmax), b(2*Nmax,2*Nmax), c1(2*Nmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta), oldh0(Nteta),             &
            oldv0(Nteta))                        
  allocate (paramG1(Nparam,NintAL*NintAL), paramG2(Nparam,NintAL*NintAL),           &
            weightsG(Nparam,NintAL*NintAL))
  allocate (Nintparam(Nparam))
  if (.not. FileGeom) then    
    call interpolation_list3D (TypeGeom, Nsurf, surf, Nint1, Nint2, NintAL, Nparam, &
         Nintparam, paramG1, paramG2, weightsG, miror, Nazimutsym)                            
  else
    do i = 1, Nparam
      do j = 1, NintAL*NintAL
        paramG1(i,j)  = 0._O
        paramG2(i,j)  = 0._O
        weightsG(i,j) = 0._O
      end do
      Nintparam(i) = 1
    end do   
  end if
  if (PrnProgress) call write_progress (.true., 1, 4)
  call matrix_Q (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np, area,   &
       Nface, Mrank, Nrank, Nmax, NintAL, Nparam, Nintparam, paramG1, paramG2,      &
       weightsG, miror, Nazimutsym, perfectcond, chiral, kb, a, Nmax, Nmax)
  if (PrnProgress) call write_progress (.false., 2, 4)    
  call copy_matrix (2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax, a1, 2*Nmax, 2*Nmax)  
  call matrix_Q (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np, area,   &
       Nface, Mrank, Nrank, Nmax, NintAL, Nparam, Nintparam, paramG1, paramG2,      &
       weightsG, miror, Nazimutsym, perfectcond, chiral, kb, b, Nmax, Nmax)
  if (PrnProgress) call write_progress (.false., 3, 4)    
  call copy_matrix (2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, b1, 2*Nmax, 2*Nmax) 
  call LU_SYSTEM (a, 2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, 2*Nmax)
  if (PrnProgress) call write_progress (.false., 4, 4)
  call write_FileTmat (Nmax, Nmax, b)
  call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrank, Nrank,     &
       Nmax, c)    
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, c, c1)
  call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfa, beta, gama, k, snorm,      &
       ExtThetaDom,.true., h, v)
  call CQscat (c1, Mrank, Nrank, Nmax, k, snorm, Cscat, Qscat)
  call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfa, beta, gama,              &
       alfap, k, snorm, Cext, Qext)  
  if (.not. FileGeom) then             
    call write_4ConvParam (Nint1, Nint2, Nrank, Mrank)                  
  else
    call write_2ConvParam (Nrank, Mrank)
  end if
  call write_DSCS (Nteta, ExtThetaDom, h, v)
  call write_Effic (Qscat, Qext)
  deallocate (paramG1, paramG2, weightsG, Nintparam)
  do i = 1, Nteta
    oldh(i)  = h(i)
    oldv(i)  = v(i)
    oldh0(i) = h(i)
    oldv0(i) = v(i)
  end do
  close (unit = iTmat)  
! --- (Nrank - 1) configuration ---
  if (PrnProgress) call write_progress_low  
  call copy_matrix (2*Nmax, 2*Nmax, a1, 2*Nmax, 2*Nmax, a, 2*Nmax,2*Nmax)    
  call matrix_Nrank_1_left (Mrank, Nrank, Nmax, a, Nmax, Nmax)  
  call copy_matrix (2*Nmax, 2*Nmax, b1, 2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax) 
  call matrix_Nrank_1_right (Mrank, Nrank, Nmax, b, Nmax, Nmax)
  call LU_SYSTEM (a,2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, 2*Nmax)          
  call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrank, Nrank,     &
       Nmax, c)    
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, c, c1)
  call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfa, beta, gama, k, snorm,      &
       ExtThetaDom,.true., h, v)
  call CQscat (c1, Mrank, Nrank, Nmax, k, snorm, Cscat, Qscat)
  call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfa, beta, gama,              &
       alfap, k, snorm, Cext, Qext)    
  call delta_DSCS (Nteta, h, v, oldh, oldv, epsNrank, NthetaConvN)
  if (.not. FileGeom) then    
    call write_4ConvParam (Nint1, Nint2, Nrank - 1, Mrank)
  else    
    call write_2ConvParam (Nrank - 1, Mrank)
  end if
  call write_DSCS (Nteta, ExtThetaDom, h, v)
  call write_Effic (Qscat, Qext)
  call write_NrankConvRes (NthetaConvN, Nteta, epsNrank)
! --- (Mrank - 1) configuration ---  
  call copy_matrix (2*Nmax, 2*Nmax, a1, 2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax)   
  call matrix_Mrank_1_left (Mrank, Nrank, Nmax, a, Nmax, Nmax)  
  call copy_matrix (2*Nmax, 2*Nmax, b1, 2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax) 
  call matrix_Mrank_1_right (Mrank, Nrank, Nmax, b, Nmax, Nmax)
  call LU_SYSTEM (a, 2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, 2*Nmax)         
  call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrank, Nrank,     &
       Nmax, c)    
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, c, c1)
  call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfa, beta, gama, k, snorm,      &
       ExtThetaDom,.true., h, v)
  call CQscat (c1, Mrank, Nrank, Nmax, k, snorm, Cscat, Qscat)
  call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfa, beta, gama,              &
       alfap, k, snorm, Cext, Qext)    
  call delta_DSCS (Nteta, h, v, oldh0, oldv0, epsMrank, NthetaConvM)
  if (.not. FileGeom) then
    call write_4ConvParam (Nint1, Nint2, Nrank, Mrank - 1)
  else
    call write_2ConvParam (Nrank, Mrank - 1)
  end if
  call write_DSCS (Nteta, ExtThetaDom, h, v) 
  call write_Effic (Qscat, Qext)
  call write_MrankConvRes (NthetaConvM, epsMrank)
  if (NthetaConvN >= int(0.8*Nteta) .and. NthetaConvM >= int(0.8*Nteta)) then
    print "(/,2x,'Convergence criteria for Nrank and Mrank are satisfied;')"                                
  else
    print "(/,2x,'Convergence criteria for Nrank and Mrank are not satisfied;')"
  end if
  call write_InfoFileTmat (FileTmat, Mrank, Nrank, .false., .false., chiral)
  call ScatCharact (k, FileTmat, Mrank, Nrank, .false., .false., chiral)
  print "(/,2x,'T matrix is stored in ',a50)", FileTmat
  print "(  2x,'The dimensions of the T matrix are given by:')"
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", Nrank
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank    
  deallocate (a1, b1)  
  deallocate (a, b, c, c1, h, v, oldh, oldv, oldh0, oldv0)      
end subroutine convergence_Nrank_MrankNONAXSYM 
!***********************************************************************************
subroutine TMatrix_Nrank_MrankNONAXSYM (FileGeom, TypeGeom, k, ind_ref, Nsurf, surf,&
           rp, np, area, Nface, Nparam, Mrank, Nrank, Nint1, Nint2, miror,          &
           Nazimutsym, perfectcond, chiral, kb, FileTmat, PrnProgress)
  use parameters
  implicit none
  integer       :: TypeGeom, Nsurf, Nface, Nparam, Mrank, Nrank, Nint1, Nint2,      &
                   Nazimutsym                   
  real(O)       :: k, surf(Nsurf), rp(3,NfacePD), np(3,NfacePD), area(NfacePD), kb
  complex(O)    :: ind_ref
  logical       :: FileGeom, miror, perfectcond, chiral, PrnProgress 
  character(80) :: FileTmat
!      
  integer       :: Nmax, i, j, NintAL
  integer,allocatable    :: Nintparam(:) 
  real(O),allocatable    :: paramG1(:,:), paramG2(:,:), weightsG(:,:)
  complex(O),allocatable :: a(:,:), b(:,:)
!  
  Nmax   = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  NintAL = max(Nint1,Nint2)  
  open (unit = iTmat, file = FileTmat, status = 'replace')
  call write_HeadFileTmat (Nmax, Nmax)   
  allocate (a(2*Nmax,2*Nmax), b(2*Nmax,2*Nmax))
  allocate (paramG1(Nparam,NintAL*NintAL), paramG2(Nparam,NintAL*NintAL),           &
            weightsG(Nparam,NintAL*NintAL))
  allocate (Nintparam(Nparam))
  if (.not. FileGeom) then    
    call interpolation_list3D (TypeGeom, Nsurf, surf, Nint1, Nint2, NintAL, Nparam, &
         Nintparam, paramG1, paramG2, weightsG, miror, Nazimutsym)                            
  else
    do i = 1, Nparam
      do j = 1, NintAL*NintAL
        paramG1(i,j)  = 0._O
        paramG2(i,j)  = 0._O
        weightsG(i,j) = 0._O
      end do
      Nintparam(i) = 1
    end do   
  end if
  if (PrnProgress) call write_progress (.true., 1, 4)
  call matrix_Q (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np, area,   &
       Nface, Mrank, Nrank, Nmax, NintAL, Nparam, Nintparam, paramG1, paramG2,      &
       weightsG, miror, Nazimutsym, perfectcond, chiral, kb, a, Nmax, Nmax)
  if (PrnProgress) call write_progress (.false., 2, 4)  
  call matrix_Q (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np, area,   &
       Nface, Mrank, Nrank, Nmax, NintAL, Nparam, Nintparam, paramG1, paramG2,      &
       weightsG, miror, Nazimutsym, perfectcond, chiral, kb, b, Nmax, Nmax)
  if (PrnProgress) call write_progress (.false., 3, 4)  
  call LU_SYSTEM (a, 2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, 2*Nmax)
  if (PrnProgress) call write_progress (.false., 4, 4)
  call write_FileTmat (Nmax, Nmax, b) 
  close (unit = iTmat)  
  call write_InfoFileTmat (FileTmat, Mrank, Nrank, .false., .false., chiral)
  call ScatCharact (k, FileTmat, Mrank, Nrank, .false., .false., chiral)
  print "(/,2x,'T matrix is stored in ',a50)", FileTmat
  print "(  2x,'The dimensions of the T matrix are given by:')"
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", Nrank
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank  
  deallocate (a, b, paramG1, paramG2, weightsG, Nintparam)      
end subroutine TMatrix_Nrank_MrankNONAXSYM 
!***********************************************************************************
!                   ROUTINES FOR UNIAXIAL ANISOTROPIC PARTICLES                    *
!***********************************************************************************
subroutine convergence_NintAnis (FileGeom, TypeGeom, k, ind_ref, ind_refZ, alfaPR,  &
           betaPR, gamaPR, snorm, Nsurf, surf, rp, np, area, Nface, Nparam, Mrank,  &
           Nrank, Nbeta, Nint1, Nint2, dNint1, dNint2, epsNint, ExtThetaDom,        &
           PrnProgress)
  use parameters
  implicit none
  integer    :: TypeGeom, Nsurf, Nface, Nparam, Mrank, Nrank, Nbeta, Nint1, Nint2,  &
                dNint1, dNint2
  real(O)    :: k, snorm, surf(Nsurf), rp(3,NfacePD), np(3,NfacePD), area(NfacePD), &
                epsNint, alfaPR, betaPR, gamaPR
  complex(O) :: ind_ref, ind_refZ  
  logical    :: FileGeom, ExtThetaDom, PrnProgress
!      
  integer    :: Nmax, Nteta, i, NthetaConv, iNint, NintAL  
  real(O)    :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,        &
                Cext, Qext
  integer,allocatable    :: Nintparam(:) 
  real(O),allocatable    :: paramG1(:,:), paramG2(:,:), weightsG(:,:), h(:), v(:),  &
                            oldh(:), oldv(:)
  complex(O),allocatable :: a(:,:), c(:), b(:,:), c1(:)   
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = Pi / 4._O
  beta   = Pi / 4._O
  gama   = 0._O
  alfap  = Pi / 4._O
  Nmax   = Nrank + Mrank * (2 * Nrank - Mrank + 1)         
  call write_TypeConvHead (1)
  allocate (a(2*Nmax,2*Nmax), c(2*Nmax), b(2*Nmax,2*Nmax), c1(2*Nmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do
  if (PrnProgress) call write_progress (.true., 1, 7)
  do iNint = 1, 2                       
    NintAL = max(Nint1,Nint2) 
    allocate (paramG1(Nparam,NintAL*NintAL), paramG2(Nparam,NintAL*NintAL),         &
              weightsG(Nparam,NintAL*NintAL))
    allocate (Nintparam(Nparam))
    call interpolation_list3D (TypeGeom, Nsurf, surf, Nint1, Nint2, NintAL, Nparam, &
         Nintparam, paramG1, paramG2, weightsG,.false., 1)   
    call matrix_Q_anis (FileGeom, TypeGeom, 3, 1, k, ind_ref, ind_refZ, alfaPR,     &
         betaPR, gamaPR, Nsurf, surf, rp, np, area, Nface, Mrank, Nrank, Nmax,      &
         Nbeta, NintAL, Nparam, Nintparam, paramG1, paramG2, weightsG, a, Nmax, Nmax)
    if (PrnProgress) call write_progress (.false., 2+3*(iNint-1), 7)                                                    
    call matrix_Q_anis (FileGeom, TypeGeom, 1, 1, k, ind_ref, ind_refZ, alfaPR,     &
         betaPR, gamaPR, Nsurf, surf, rp, np, area, Nface, Mrank, Nrank, Nmax,      &
         Nbeta, NintAL, Nparam, Nintparam, paramG1, paramG2, weightsG, b, Nmax, Nmax)
    if (PrnProgress) call write_progress (.false., 3+3*(iNint-1), 7)
    call LU_SYSTEM (a, 2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, 2*Nmax)
    if (PrnProgress) call write_progress (.false., 4+3*(iNint-1), 7)
    call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrank, Nrank,   &
         Nmax, c)    
    call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, c, c1)
    call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfa, beta, gama, k, snorm,    &
         ExtThetaDom,.true., h, v)
    call CQscat (c1, Mrank, Nrank, Nmax, k, snorm, Cscat, Qscat)
    call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfa, beta, gama,            &
         alfap, k, snorm, Cext, Qext)     
    call delta_DSCS (Nteta, h, v, oldh, oldv, epsNint, NthetaConv) 
    call write_4ConvParam (Nint1, Nint2, Nrank, Mrank)                  
    call write_DSCS (Nteta, ExtThetaDom, h, v)
    call write_Effic (Qscat, Qext)
    Nint1 = Nint1 + dNint1
    Nint2 = Nint2 + dNint2
    deallocate (paramG1, paramG2, weightsG, Nintparam)
  end do
  call write_NintConvRes (NthetaConv, Nteta, epsNint)                           
  deallocate (a, b, c, c1, h, v, oldh, oldv)    
end subroutine convergence_NintAnis
!***********************************************************************************
subroutine convergence_Nrank_MrankAnis (FileGeom, TypeGeom, k, ind_ref, ind_refZ,   &
           alfaPR, betaPR, gamaPR, snorm, Nsurf, surf, rp, np, area, Nface, Nparam, &
           Mrank, Nrank, Nbeta, Nint1, Nint2, epsNrank, epsMrank, ExtThetaDom,      &
           FileTmat, PrnProgress)
  use parameters
  implicit none
  integer       :: TypeGeom, Nsurf, Nface, Nparam, Mrank, Nrank, Nint1, Nint2, Nbeta
  real(O)       :: k, snorm, surf(Nsurf), rp(3,NfacePD), np(3,NfacePD),             &
                   area(NfacePD), epsNrank, epsMrank, alfaPR, betaPR, gamaPR
  complex(O)    :: ind_ref, ind_refZ
  character(80) :: FileTmat
  logical       :: FileGeom, ExtThetaDom, PrnProgress
!        
  integer       :: Nmax, Nteta, i, NthetaConvN, NthetaConvM, NintAL
  real(O)       :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,     &
                   Cext, Qext
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG1(:,:), paramG2(:,:), weightsG(:,:), h(:), v(:),  &
                            oldh(:), oldv(:), oldh0(:), oldv0(:)
  complex(O),allocatable :: a(:,:), c(:), b(:,:), c1(:), a1(:,:), b1(:,:)   
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = Pi / 4._O
  beta   = Pi / 4._O
  gama   = 0._O
  alfap  = Pi / 4._O
  Nmax   = Nrank + Mrank * (2 * Nrank - Mrank + 1)   
  NintAL = max(Nint1,Nint2)
  call write_TypeConvHead (4)  
  open (unit = iTmat, file = FileTmat, status = 'replace')
  call  write_HeadFileTmat (Nmax, Nmax)    
  allocate (a1(2*Nmax,2*Nmax), b1(2*Nmax,2*Nmax))  
  allocate (a(2*Nmax,2*Nmax), c(2*Nmax), b(2*Nmax,2*Nmax), c1(2*Nmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta), oldh0(Nteta),             &
            oldv0(Nteta))                        
  allocate (paramG1(Nparam,NintAL*NintAL), paramG2(Nparam,NintAL*NintAL),           &
            weightsG(Nparam,NintAL*NintAL))
  allocate (Nintparam(Nparam))
  call interpolation_list3D (TypeGeom, Nsurf, surf, Nint1, Nint2, NintAL, Nparam,   &
       Nintparam, paramG1, paramG2, weightsG,.false., 1)   
  if (PrnProgress) call write_progress (.true., 1, 4)
  call matrix_Q_anis (FileGeom, TypeGeom, 3, 1, k, ind_ref, ind_refZ, alfaPR,       &
       betaPR, gamaPR, Nsurf, surf, rp, np, area, Nface, Mrank, Nrank, Nmax, Nbeta, &
       NintAL, Nparam, Nintparam, paramG1, paramG2, weightsG, a, Nmax, Nmax)
  if (PrnProgress) call write_progress (.false., 2, 4)   
  call copy_matrix (2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax, a1, 2*Nmax, 2*Nmax)
  call matrix_Q_anis (FileGeom, TypeGeom, 1, 1, k, ind_ref, ind_refZ, alfaPR,       &
       betaPR, gamaPR, Nsurf, surf, rp, np, area, Nface, Mrank, Nrank, Nmax, Nbeta, &
       NintAL, Nparam, Nintparam, paramG1, paramG2, weightsG, b, Nmax, Nmax)
  if (PrnProgress) call write_progress (.false., 3, 4)  
  call copy_matrix (2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, b1, 2*Nmax, 2*Nmax)
  call LU_SYSTEM (a, 2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, 2*Nmax)
  if (PrnProgress) call write_progress (.false., 4, 4)
  call write_FileTmat (Nmax, Nmax, b)
  call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrank, Nrank,     &
       Nmax, c)    
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, c, c1)
  call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfa, beta, gama, k, snorm,      &
       ExtThetaDom,.true., h, v)
  call CQscat (c1, Mrank, Nrank, Nmax, k, snorm, Cscat, Qscat)
  call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfa, beta, gama,              &
       alfap, k, snorm, Cext, Qext)             
  call write_4ConvParam (Nint1, Nint2, Nrank, Mrank)
  call write_DSCS (Nteta, ExtThetaDom, h, v)
  call write_Effic (Qscat, Qext)
  deallocate (paramG1, paramG2, weightsG, Nintparam)
  do i = 1, Nteta
    oldh(i)  = h(i)
    oldv(i)  = v(i)
    oldh0(i) = h(i)
    oldv0(i) = v(i)
  end do
  close (unit = iTmat)  
! --- (Nrank - 1) configuration ---
  if (PrnProgress) call write_progress_low  
  call copy_matrix (2*Nmax, 2*Nmax, a1, 2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax) 
  call matrix_Nrank_1_left (Mrank, Nrank, Nmax, a, Nmax, Nmax)  
  call copy_matrix (2*Nmax, 2*Nmax, b1,  2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax)  
  call matrix_Nrank_1_right (Mrank, Nrank, Nmax, b, Nmax, Nmax)
  call LU_SYSTEM (a, 2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, 2*Nmax)         
  call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrank, Nrank,     &
       Nmax, c)    
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, c, c1)
  call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfa, beta, gama, k, snorm,      &
       ExtThetaDom,.true., h, v)
  call CQscat (c1, Mrank, Nrank, Nmax, k, snorm, Cscat, Qscat)
  call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfa, beta, gama,              &
       alfap, k, snorm, Cext, Qext)    
  call delta_DSCS (Nteta, h, v, oldh, oldv, epsNrank, NthetaConvN)
  call write_4ConvParam (Nint1, Nint2, Nrank - 1, Mrank)
  call write_DSCS (Nteta, ExtThetaDom, h, v)
  call write_Effic (Qscat, Qext)
  call write_NrankConvRes (NthetaConvN, Nteta, epsNrank)                                          
! --- (Mrank - 1) configuration ---  
  call copy_matrix (2*Nmax, 2*Nmax, a1, 2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax) 
  call matrix_Mrank_1_left (Mrank, Nrank, Nmax, a, Nmax, Nmax)  
  call copy_matrix (2*Nmax, 2*Nmax, b1, 2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax)  
  call matrix_Mrank_1_right (Mrank, Nrank, Nmax, b, Nmax, Nmax)
  call LU_SYSTEM (a, 2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, 2*Nmax)         
  call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrank, Nrank,     &
       Nmax, c)    
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, c, c1)
  call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfa, beta, gama, k, snorm,      &
       ExtThetaDom,.true., h, v)
  call CQscat (c1, Mrank, Nrank, Nmax, k, snorm, Cscat, Qscat)
  call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfa, beta, gama,              &
       alfap, k, snorm, Cext, Qext)    
  call delta_DSCS (Nteta, h, v, oldh0, oldv0, epsMrank, NthetaConvM)
  call write_4ConvParam (Nint1, Nint2, Nrank, Mrank - 1)
  call write_DSCS (Nteta, ExtThetaDom, h, v) 
  call write_Effic (Qscat, Qext)
  call write_MrankConvRes (NthetaConvM, epsMrank)
  if (NthetaConvN >= int(0.8*Nteta) .and. NthetaConvM >= int(0.8*Nteta)) then
    print "(/,2x,'Convergence criteria for Nrank and Mrank are satisfied;')"
  else
    print "(/,2x,'Convergence criteria for Nrank and Mrank are not satisfied;')"
  end if          
  deallocate (a1, b1)  
  call write_InfoFileTmat (FileTmat, Mrank, Nrank, .false., .false., .false.)
  call ScatCharact (k, FileTmat, Mrank, Nrank, .false., .false., .false.)
  print "(/,2x,'The T matrix is stored in ',a50)", FileTmat
  print "(  2x,'The dimensions of the T matrix are given by:')"
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,';')", Nrank
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank  
  deallocate (a, b, c, c1, h, v, oldh, oldv, oldh0, oldv0)      
end subroutine convergence_Nrank_MrankAnis
!***********************************************************************************
subroutine TMatrix_Nrank_MrankAnis (FileGeom, TypeGeom, k, ind_ref, ind_refZ,       &
           alfaPR, betaPR, gamaPR, Nsurf, surf, rp, np, area, Nface, Nparam, Mrank, &
           Nrank, Nbeta, Nint1, Nint2, FileTmat, PrnProgress)
  use parameters
  implicit none
  integer       :: TypeGeom, Nsurf, Nface, Nparam, Mrank, Nrank, Nint1, Nint2, Nbeta
  real(O)       :: k, surf(Nsurf), rp(3,NfacePD), np(3,NfacePD), area(NfacePD),     &
                   alfaPR, betaPR, gamaPR
  complex(O)    :: ind_ref, ind_refZ
  character(80) :: FileTmat
  logical       :: FileGeom, PrnProgress
!        
  integer       :: Nmax, NintAL  
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG1(:,:), paramG2(:,:), weightsG(:,:)
  complex(O),allocatable :: a(:,:), b(:,:) 
!  
  Nmax   = Nrank + Mrank * (2 * Nrank - Mrank + 1)   
  NintAL = max(Nint1,Nint2)  
  open (unit = iTmat, file = FileTmat, status = 'replace')
  call  write_HeadFileTmat (Nmax, Nmax)    
  allocate (a(2*Nmax,2*Nmax), b(2*Nmax,2*Nmax))                          
  allocate (paramG1(Nparam,NintAL*NintAL), paramG2(Nparam,NintAL*NintAL),           &
            weightsG(Nparam,NintAL*NintAL))
  allocate (Nintparam(Nparam))
  call interpolation_list3D (TypeGeom, Nsurf, surf, Nint1, Nint2, NintAL, Nparam,   &
       Nintparam, paramG1, paramG2, weightsG,.false., 1)   
  if (PrnProgress) call write_progress (.true., 1, 4)
  call matrix_Q_anis (FileGeom, TypeGeom, 3, 1, k, ind_ref, ind_refZ, alfaPR,       &
       betaPR, gamaPR, Nsurf, surf, rp, np, area, Nface, Mrank, Nrank, Nmax, Nbeta, &
       NintAL, Nparam, Nintparam, paramG1, paramG2, weightsG, a, Nmax, Nmax)
  if (PrnProgress) call write_progress (.false., 2, 4)   
  call matrix_Q_anis (FileGeom, TypeGeom, 1, 1, k, ind_ref, ind_refZ, alfaPR,       &
       betaPR, gamaPR, Nsurf, surf, rp, np, area, Nface, Mrank, Nrank, Nmax, Nbeta, &
       NintAL, Nparam, Nintparam, paramG1, paramG2, weightsG, b, Nmax, Nmax)
  if (PrnProgress) call write_progress (.false., 3, 4)  
  call LU_SYSTEM (a, 2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, 2*Nmax)
  if (PrnProgress) call write_progress (.false., 4, 4)
  call write_FileTmat (Nmax, Nmax, b)
  close (unit = iTmat)
  call write_InfoFileTmat (FileTmat, Mrank, Nrank, .false., .false., .false.)
  call ScatCharact (k, FileTmat, Mrank, Nrank, .false., .false., .false.)
  print "(/,2x,'The T matrix is stored in ',a50)", FileTmat
  print "(  2x,'The dimensions of the T matrix are given by:')"
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,';')", Nrank
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank  
  deallocate (a, b, paramG1, paramG2, weightsG, Nintparam)      
end subroutine TMatrix_Nrank_MrankAnis
