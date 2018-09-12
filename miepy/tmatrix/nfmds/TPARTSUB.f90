subroutine TPARTSUB
!------------------------------------------------------------------------------------
! 1. General Considerations                                                         !
! --------------------------                                                        !
! TPARTSUB is a routine for computing the scattering characteristics of a           !
! homogeneous, dielectric or perfectly conducting, axisymmetric particle on or      !
! near a plane surface.                                                             !
!                                                                                   !
! The axisymmetric particle is situated in the neighbourhood of a plane surface,    !
! so that its axis of symmetry coincides with the normal to the plane surface.      !
! The particle is illuminated by a vector pane wave traveling in the ambient medium,! 
! or by a vector plane wave traveling in the substrate (Figure 1). In the first     !
! case, the direction of the wave vector is (theta = beta, phi = Pi), while in the  !
! second case, the direction is (theta = Pi - beta, phi = 0), where beta is the     !
! incident angle.                                                                   !
!                                                                                   !
!                                                  beta ^ Z                         !
!                                                    !  | 		                      !                                                        			!
!                                                  * ...| 		                      ! 
!                                                   *   | 		                      !
!           Z ^ substrate              incident - -> *  |  substrate	              !
!             |    !                   direction      * |    !                      !
!             |    !                                   *|    !                      !
!     ================/=                         =======o=======/= --> x            ! 
!             |       !                                 |       !                   !
!          -------    ! z0                           -------    ! z0                ! 
!         |   |   |   !                             |   |   |   !                   !
!         |   |   |   !    X                        |   |   |   !    X              !
!      ---|---O---|---/---->                     ---|---O---|---/---->              ! 
!         |   |*  |                                 |   |   |                       !
!         |   | * |<- - particle                    |   |   |< - - particle         !
!          ------*                                   -------                        !
!             |   * <- - incident                       |                           !
!             |... *     direction                      |                           !
!             | !   *                                                               !
!              beta                                                                 !
!                                                                                   !
!            (a)                                        (b)                         !
!                                                                                   !
!  Figure 1. The particle is illuminated by: (a) a vector plane wave traveling in   !
!  the ambient medium and (b) a vector plane plane wave traveling in the substrate. !
!                                                                                   !
! The scattering problem is a multiple particles problem and the solution method    !
! is the separation of variables technique. To model the scattering problem in the  !
! framework of the separation of variables technique we must address how the        !
! radiation interacts with the particle. The incident field strikes the particle    !
! either directly or after interacting with the surface, while the fields           !
! emanating from the particle may also reflect off the surface and interact with    !
! the particle again. The transition matrix relating the incident and scattered     !
! field coefficients is computed in the framework of the null-field method.         !
! For axisymmetric particles, the T-matrix computation decouples over the           !
! azimuthal modes and the surface integrals reduce to line integrals over the       !
! generatrix. The integrals are evaluated by means of Gauss-Legendre quadrature     !
! method, and note that for particles with a plane of symmetry perpendicular to the !
! axis of rotation (mirror symmetric particles), the integration can be performed   !
! along the half-generatrix curve. The numerical stability and accuracy of T-matrix !
! computations for an unsmooth generatrix can be enhanced by dividing the           !
! generatrix into piecewise smooth curves and applying separate quadrature          !
! formulas to each curve. The reflection matrix characterizing the reflection of    !
! the scattered field by the surface is computed by using the integral              !
! representation for the vector spherical wave functions over plane waves. As a     !
! result, the elements of the reflection matrix are expressed as integrals over a   !
! complex angle.                                                                    !
!                                                                                   !
! Particle geometries currently supported include: spheroids, cylinders and         !
! rounded oblate cylinders. As for axisymmetric particles, the routines             !
! "interpolation_listAXSYM" and "elem_geomAXSYM" provide the required geometry      !
! parameters. For more details, consult the TAXSYM routine.                         !
!                                                                                   !
! IMPORTANT NOTE: For cylinders and rounded oblate cylinders, THE GEOMETRIC         !
! CONSTRAINT of the T-matrix method IS VIOLATED and convergence problems occur.     !
! Therefore, we recommend to use the code only for PROLATE SPHEROIDS.               !
!                                                                                   !
! The code computes the differential scattering cross sections for parallel and     !
! perpendicular polarizations in the azimuthal plane phiGS and at NthetaGS          !
! scattering angles. Denoting by                                                    !
!                                                                                   !
!            ES = [exp(jkR) / kR] * ESinf = [exp(jkR) / R] * FS                     !
!                                                                                   ! 
! the scattered field in the direction (theta,phi), where FS is the far-field       !
! pattern and ESinf = k * FS, and assuming the decomposition                        !
!                                                                                   !  
!            ESinf = ESinf_theta * e_theta + ESinf_phi * e_phi,                     !
!                                                                                   !
! where e_theta and e_phi are the spherical unit vectors, we define the             !
! differential scattering cross sections for parallel and perpendicular             !
! polarizations by                                                                  !
!                                                                                   ! 
!            h = |ESinf_theta|**2 / k**2 = |FS_theta|**2                            !
!                                                                                   !
! and                                                                               !
!            v = |ESinf_phi|**2 / k**2   = |FS_phi|**2,                             !
!                                                                                   !
! respectively. The normalized differential scattering cross sections for parallel  !
! and perpendicular polarizations are given by                                      ! 
!                                                                                   ! 
!     hn = |ESinf_theta|**2 / [Pi * (k*anorm)**2] = |FS_theta|**2 / [Pi * anorm**2] !
!                                                                                   !
! and                                                                               !
!                                                                                   !
!     vn = |ESinf_phi|**2 / [Pi * (k*anorm)**2]   = |FS_phi|**2 / [Pi * anorm**2],  !
!                                                                                   !
! where anorm is a characteristic dimension of the particle. For a more complex     !
! analysis, the code also computes the electromagnetic fields ESinf_theta and       !
! ES_inf_phi at Nphi scattering planes.                                             !
!                                                                                   !
! 2. Convergence Test                                                               !
! --------------------                                                              !
! Convergence tests precede the scattering characteristics calculation. The         !
! following parameters control the T-matrix computation:                            !
! - the number of integration points for computing integrals over the generatrix,   ! 
!   Nint,                                                                           !
! - the number of integration points for computing the elements of the reflection   !
!   matrix, NXint,                                                                  !
! - the maximum expansion order Nrank, and                                          !
! - the maximum azimuthal order Mrank.                                              ! 
! Nint is the global number of integration points on the generatrix curve and has   !
! the same significance as in the routine "TAXSYM".                                 ! 
! The integrals appearing in the expression of the reflection matrix are of the     !
! form                                                                              !
!                                                                                   !
!           |Pi/2 - j * inf                                                         !
!       I = |     f(cos(beta)) * exp{ 2 * j * q * cos(beta) } * sin(beta) dbeta.    !
!           |0                                                                      !
!                                                                                   !
! Changing variables from beta to x = - 2 * j * q * [cos(beta) - 1], we have        !
!                                                                                   !
!                                           |inf                                    !
!       I = [exp{2 * j* q} / (2 * j * q)] * |   f(1 - x/(2 * j * q)) * exp{ -x } dx.!
!                                           |0                                      !
!                                                                                   ! 
! and integrals of this type can be computed efficiently by using Laguerre          !
! polynomials. In this context, an important parameter of the convergence procedure !
! is the number of integration points NXint for the Laguerre quadrature method.     !
!                                                                                   !
! For convergence tests, the incident wave is assumed to be a vector plane wave and !
! the scattering characteristics are computed in the azimuthal plane  phi = 0°. The !
! convergence tests over Nint, NXint and Nrank are interactive, while the           !
! convergence test over Mrank is automatically performed.                           ! 
!                                                                                   !
! For the integration and expansion order test, the incident wave propagates along  !
! the axis of symmetry of the particle. In this case, all azimuthal modes are zero  !
! except for m = - 1 and m = 1. For the convergence test over Nint and NXint, the   !
! scattering problem is solved for the pairs (Nint,NXint) and                       !
! (Nint + dNint, NXint + dNint), while for the convergence test over Nrank, the     !
! scattering problem is solved for Nrank and Nrank - 1. The normalized              !
! differential scattering cross section (DSCS) will be checked at 20° increments    !
! for convergence within epsX (epsNint or epsNrank) tolerance. If the calculated    !
! results converge within this tolerance at 80% of the scattering angles, then      !
! convergence is achieved.                                                          !
!                                                                                   ! 
! After Nrank, Nint and NXint have been determined we pass to the azimuthal         !
! order test. The program automatically sets the incident angle to a more general   !
! orientation, i.e., beta = 45°, and solves the scattering problem  for increasing  !
! m values until convergence of the angular scattering is achieved.                 !
!                                                                                   !
! The convergence tests over Nint, NXint and Nrank can be switched off by setting   !
! the logical variable DoConvTest to false. In this case, the values of Nint,       !
! NXint and Nrank must be specified in the input file.                              !
!                                                                                   !
! 3. Estimates of Nint and Nrank                                                    !
! --------------------------------                                                  !
! The above convergence tests require estimates of Nint, NXint and Nrank.           ! 
! These estimates must be supplied by the user, and for this purpose, Wiscombe's    !
! truncation limit criterion [W. J. Wiscombe, Improved Mie scattering algorithms,   ! 
! Applied Optics, 19, 1505-1509, 1980] can be used.                                 !
!                                                                                   !
! The truncation limit criterion proposed by Wiscombe provides the estimate         !
!                                                                                   !
!                NrankW = int(x + 4.05 * x**0.33 + 2),                              !
!                                                                                   !
! where x is the size parameter, x = k * Rcirc, k is the wave number and Rcirc is   !
! the radius of the smallest circumscribing sphere.                                 !
!                                                                                   ! 
! The value of Nint depends on the size parameter, the particle shape and the       !
! relative refractive index. A conservative estimate of Nint is                     !
!                                                                                   !
!                        Nint = Ndgs * Nrank,                                       !
!                                                                                   !
! where Ndgs = 10, 15, 20, ....                                                     !
!                                                                                   !
! The value of NXint also depends on the size parameter of the particle. For size   !
! parameters smaller than 5, NXint = 60, 70,...,90.                                 !
!                                                                                   !
! 4. Input Parameters                                                               !
! --------------------                                                              !
! The parameters specified in the input file "/INPUTFILES/InputPARTSUB.dat" are     !
! listed below.                                                                     !
!                                                                                   !
! - wavelength (real) - wavelength of the incident light in the ambient medium.     !
!                                                                                   !
! - ind_refPART (complex) - relative refractive index of the particle with respect  !
!   to the ambient medium. The imaginary part of the relative refractive index      !
!   must be zero for nonabsorbing particles and positive for absorbing particles.   !
!                                                                                   !
! - ind_refSUB (complex) - relative refractive index of the substrate with respect  !
!   to the ambient medium. The imaginary part of the relative refractive index      !
!   must be zero or positive.                                                       !
!                                                                                   !
! - TypeScat (integer) - if TypeScat = 1, the particle is illuminated by a vector   !
!   pane wave traveling in the ambient medium, while for TypeScat = 2 the particle  !
!   is illuminated by a vector plane wave traveling in the substrate. In the first  !
!   case the direction of the wave vector is (theta = beta, phi = Pi), while in the !
!   second case the direction is (theta = Pi - beta, phi = 0), where beta is the    !
!   incident angle.                                                                 !
!                                                                                   !
! - TypeGeom (integer) - parameter specifying the type of the particle geometry.    !
!                                                                                   !
! - Nsurf (integer) - number of surface parameters.                                 !
!                                                                                   !
! - surf (real array: surf(1), surf(2),...,surf(Nsurf)) - surface parameters        !
!   specifying the shape of the particle. The dimension of the array surf is        !
!   NsurfPD. The integer parameter NsurfPD is specified in the routine              !
!   "Parameters.f90" and has the value NsurfPD = 10. If Nsurf > NsurfPD, the        !
!   execution is automatically terminated.                                          !
!                                                                                   !
! - Nparam (integer) - number of smooth curves forming the generatrix curve.        !
!                                                                                   !
!   The permissive values of the surface parameters (for the supported geometries)  !
!   are given below.                                                                !
!                                                                                   !
!   Particle      TypeGeom   Nsurf   Nparam                surf                     !
!   spheroid         1         2       1         surf(1) - length of the semi-      !
!                                                          axis along the           !
!                                                          symmetry axis            !
!                                                surf(2) - length of the second     !
!                                                          semi-axis                !
!                                                                                   !
!   cylinder         2         2       3         surf(1) - half-length of           !
!                                                          the cylinder             !
!                                                surf(2) - cylinder radius          !
!                                                                                   !
!   rounded          3         2       3         surf(1) - half-length of           !
!    oblate                                                the cylinder             !
!   cylinder                                     surf(2) - cylinder radius          ! 
!                                                          including the rounded    !
!                                                          part                     !
!                                                                                   !
!   NOTE: Due to convergence problems occurring for cylindrical particles, it is    !
!   recommended to use the code only for PROLATE SPHEROIDS.                         !
!                                                                                   !
! - z0 (real) - axial position of the substrate in the particle coordinate system.  !
!   Note that z0 must be a POSITIVE quantity.                                       !         
!                                                                                   !
! - anorm (real) - characteristic length of the particle which is used to           !
!   normalize the differential scattering cross sections.                           !
!                                                                                   !
! - Rcirc (real) - characteristic length of the particle (usually the radius of     !
!   the smallest circumscribing sphere) which is used to compute an estimate of the !
!   maximum expansion order by using Wiscombe's truncation limit criterion (the size!
!   parameter is x = k * Rcirc, where k is the wave number in the ambient medium).  !
!   Alternatively, Rcirc can be chosen as the equal-volume sphere radius or the     !
!   surface-equivalent-sphere radius. This parameter must be specified if the       !
!   interactive convergence tests over Nint, NXint and Nrank are used               !
!   (DoConvTest = t).                                                               ! 
!                                                                                   !
! - miror (logical) - if miror = t, the particle is mirror symmetric (the plane of  !
!   symmetry or the plane of reflection is perpendicular to the axis of rotation).  !
!                                                                                   !
! - perfectcond (logical) - if perfectcond = t, the particle is perfectly           !
!   conducting.                                                                     !
!                                                                                   !
! - DoConvTest (logical) - if DoConvTest = t, the interactive convergence tests     !
!   over Nint, NXint and Nrank are invoked. An estimate of Nrank is given by        !
!   Wiscombe's truncation limit criterion. If DoConvTest = f, the values of Nint,   !
!   NXint and Nrank must be supplied in the input file.                             !
!                                                                                   !
! - Nint (integer) - number of integration points in computing integrals over       !
!   the generatrix curve. This parameter is used if the convergence tests are not   !
!   performed (DoConvTest = f).                                                     !
!                                                                                   !
! - NXint (integer) - number of integration points for computing the elements of    !
!   the reflection matrix. This parameter is used if the convergence tests are      !
!   not performed (DoConvTest = f).                                                 ! 
!                                                                                   ! 
! - Nrank (integer) - maximum expansion order. This parameter is used if the        !
!   convergence tests are not performed (DoConvTest = f).                           !
!                                                                                   !
! - epsNint (real) - error tolerance for the integration test.                      !
!                                                                                   !
! - epsNrank (real) - error tolerance for the expansion order test.                 !
!                                                                                   !
! - epsMrank (real) - error tolerance for the azimuthal order test.                 !
!                                                                                   !
! - dNint (integer) - number of division points for the integration test. Note      !
!   that the scattering problem is solved for Nint and Nint + dNint.                ! 
!                                                                                   !
! The following parameters control the scattering characteristics calculation.      !
! - beta (real) - incident angle.                                                   !
!                                                                                   !
! - alphap (real) - polarization angle of the incident wave. If                     ! 
!  (e_r,e_theta,e_phi) are the spherical unit vectors of the incident               !
!   direction (thetaGI,phiGI), then alphap is the angle between the electric        ! 
!   field vector and the unit vector e_theta.                                       !
!                                                                                   !
! - ComputeDSCS (logical) - if ComputeDSCS = t, the differential scattering         !
!   cross sections are computed in the azimuthal plane phiGS at NthetaGS            !
!   scattering angles.                                                              !
!                                                                                   !
! - ComputeFields (logical) - if ComputeFields = t, the electromagnetic fields      !
!   are computed at Nphi scattering planes.                                         !
!                                                                                   !
! - phiGS (real) - azimuthal angle of the scattering plane in which the             !
!   differential scattering cross sections are computed.                            !
!                                                                                   !
! - NthetaGS (integer) - number of scattering angle at which the differential       !
!   scattering cross sections are computed. The scattering angles are uniformly     !
!   spaced.                                                                         !
!                                                                                   !
! - Nphi (integer) - number of scattering planes at which the electromagnetic       !
!   fields are computed.                                                            !
!                                                                                   !
! - phi (real array: phi(1), phi(2),..., phi(Nphi)) - azimuthal angles of the       !
!   scattering planes at which the electromagnetic fields are computed.             !
!                                                                                   !
! - Ntheta (integer array: Ntheta(1), Ntheta(2),..., Ntheta(Nphi)) - number of      !
!   zenith angle values in each scattering plane at which the electromagnetic       !
!   fields are computed.                                                            !
!                                                                                   !
! - thetamin, thetamax (real arrays: thetamin(1), thetamin(2),...,                  !
!   thetamin(Nphi) and thetamax(1), thetamax(2),..., thetamax(Nphi)) - specify      !
!   the minimal and maximal values of the zenith angle in each scattering plane     !
!   at which the electromagnetic fields are computed.                               !
!                                                                                   !
! - normalized (logical) - if normalized = t, the average differential              ! 
!   scattering cross sections are normalized by Pi * anorm**2, where anorm is the   !
!   characteristic length of the particle.                                          !
!                                                                                   !
! - FileDSCS (character(80)) - name of the file to which the differential           !
!   scattering cross sections are written.                                          !
!                                                                                   !
! - FileEMF (character(80)) - name of the file to which the electromagnetic         !
!   fields are written.                                                             !
!                                                                                   ! 
! - WriteInputInfo (logical) - if WriteInputInfo = t, the input parameters of       !
!   the scattering problem are written to the output files FileDSCS and             !
!   FileEMF.                                                                        !
!                                                                                   !
! - PrnProgress (logical) - if PrnProgress = t, the progress of calculation         !
!   is printed.                                                                     !
!                                                                                   !
! Note that all input parameters specifying lengths must be provided in the same    !
! units of lengths.                                                                 !
!                                                                                   !
! 5. Logical Scheme                                                                 !
! ------------------                                                                !  
! The organization of the code is as follows:                                       !
! 1. < read input data >                                                            !
! 2. < T-matrix calculation >                                                       !
!                                                                                   !
!    if ( do convergence test) then                                                 !                
!      < the code computes an estimate of Nrank accordingly to                      !
!         Wiscombe's truncation criterion >                                         !
!      < the code prompts for the estimated values of Nint, NXint and Nrank >       !
!      < the code prompts for the type of convergence test:                         !
!        1 - over Nint and NXint, 2 - over Nrank or 3 - over Mrank >                !
!    else if ( .not. do convergence test) then                                      !
!      < the input file provides the values of Nint, NXint and Nrank >              !
!      type of convergence test = 3 (convergence test over Mrank)                   !
!    end if                                                                         !
!                                                                                   !
!    if ( do integration test, i.e., type of convergence test = 1 ) then            !
!       < the code computes the DSCS for (Nint,NXint) and (Nint + dNint,            !
!         NXint + dNint) and write the results to the file "Output.dat" >           !
!       < the execution is terminated >                                             !
!    end if                                                                         !
!    if ( do expansion order test, i.e., type of convergence test = 2 ) then        !
!       < the code computes the DSCS for Nrank and Nrank - 1 and                    !
!         write the results to the file "Output.dat" >                              !
!       < the execution is terminated >                                             !
!    end if                                                                         !
!    if ( do azimuthal order test, i.e., type of convergence test = 3 ) then        !
!       < the code computes the DSCS for increasing m values and                    !
!         write the results to the file "Output.dat" >                              !
!       < go to Step 3 >                                                            !
!    end if                                                                         !
! 3. < scattering characteristics calculation >                                     !
!------------------------------------------------------------------------------------
  use parameters
  implicit none 
  integer       :: TypeGeom, TypeScat, Nsurf, Nparam, TypeConvTest, dNint, Nrank,   &
                   Nint, NXint                                   
  real(O)       :: wavelength, z0, anorm, surf(NsurfPD), snorm, epsNint, epsNrank,  &
                   epsMrank, k, Rcirc                                       
  complex(O)    :: ind_refPART, ind_refSUB
  logical       :: miror, perfectcond, DoConvTest, PrnProgress       
! -----------------------------------------------------------------------------------
!                                 Read the input file                               ! 
! -----------------------------------------------------------------------------------
  call readinputPARTSUB ( wavelength, ind_refPART, ind_refSUB, TypeScat,            &
       TypeGeom, Nsurf, surf, Nparam, z0, anorm, Rcirc, miror, perfectcond,         &
       DoConvTest, Nint, NXint, Nrank, epsNint, epsNrank, epsMrank, dNint,          &
       PrnProgress, TypeConvTest, k, snorm)  
! -----------------------------------------------------------------------------------
!                                      Main                                         !
! -----------------------------------------------------------------------------------   
  open (unit = iOutput, file = FileOutput, status = "replace")   
  call printinputPARTSUB (TypeGeom, TypeScat, Nsurf, Nparam, dNint, wavelength,     &
       anorm, Rcirc, z0, surf, epsNint, epsNrank, epsMrank, ind_refPART, ind_refSUB,&
       miror, perfectcond) 
  if (DoConvTest) then
    if (TypeConvTest == 1) then
      call convergence_NintPARTSUB (TypeGeom, TypeScat, k, ind_refSUB, ind_refPART, &
           z0, snorm, Nsurf, surf, Nparam, Nrank, Nint, NXint, dNint, miror,        &
           perfectcond, epsNint, PrnProgress)      
    else if (TypeConvTest == 2) then
      call convergence_NrankPARTSUB (TypeGeom, TypeScat, k, ind_refSUB, ind_refPART,&
           z0, snorm, Nsurf, surf, Nparam, Nrank, Nint, NXint, miror, perfectcond,  &
           epsNrank, PrnProgress)        
    else 
      call convergence_MrankPARTSUB (TypeGeom, TypeScat, k, ind_refSUB, ind_refPART,&
           z0, snorm, Nsurf, surf, Nparam, Nrank, Nint, NXint, miror, perfectcond,  &
           epsMrank, PrnProgress)         
    end if
  else
    call convergence_MrankPARTSUB (TypeGeom, TypeScat, k, ind_refSUB, ind_refPART,  &
         z0, snorm, Nsurf, surf, Nparam, Nrank, Nint, NXint, miror, perfectcond,    &
         epsMrank, PrnProgress) 
  end if
  close (unit = iOutput)           
end subroutine TPARTSUB
!***********************************************************************************
subroutine readinputPARTSUB ( wavelength, ind_refPART, ind_refSUB, TypeScat,        &
           TypeGeom, Nsurf, surf, Nparam, z0, anorm, Rcirc, miror, perfectcond,     &
           DoConvTest, Nint, NXint, Nrank, epsNint, epsNrank, epsMrank, dNint,      &
           PrnProgress, TypeConvTest, k, snorm)     
  use parameters
  use derived_parameters
  implicit none 
  integer       :: TypeGeom, TypeScat, Nsurf, Nparam, TypeConvTest, dNint, Nrank,   &
                   Nint, NXint, i, ios, NrankW                                    
  real(O)       :: wavelength, z0, anorm, surf(NsurfPD), xpart, snorm, epsNint,     &
                   epsNrank, epsMrank, k, Rcirc, x                                       
  complex(O)    :: ind_refPART, ind_refSUB
  character(80) :: string
  logical       :: miror, perfectcond, DoConvTest, PrnProgress, XFindPar     
! -----------------------------------------------------------------------------------
!                         Read the input file FileInputPARTSUB                      ! 
! -----------------------------------------------------------------------------------    
  call DrvParameters 
  open (unit = iInputPARTSUB, file = FileInputPARTSUB, status = "old",              &
        position = "rewind")   
  wavelength  = 0.1_O * 2._O * Pi  
  ind_refPART = (1.5_O,0._O)                                                                    
  ind_refSUB  = (1.3_O,0._O)
  string      = 'OptProp'
  if (XFindPar (iInputPARTSUB, string)) then
    read (iInputPARTSUB, *, iostat = ios) wavelength
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable wavelength;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) ind_refPART
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refPART;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) ind_refSUB
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refSUB;')"
      stop
    end if      
  else
    print "(/,2x,'Group name OptProp not found;')"
    stop  
  end if 
  call check_ind_ref2 (ind_refPART, ind_refSUB)  
  k = 2._O * Pi / wavelength 
!
  TypeScat = 1
  TypeGeom = 1 
  Nsurf = 2
  do i = 1, NsurfPD
    surf(i) = 0.1_O
  end do
  Nparam = 1
  z0     = 0.1_O
  anorm  = 0.1_O
  Rcirc  = 0.1_O
  miror  = .true.
  perfectcond = .false.
  string      = 'GeomProp'
  if (XFindPar (iInputPARTSUB, string)) then
    read (iInputPARTSUB, *, iostat = ios) TypeScat
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable TypeScat;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) TypeGeom
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable TypeGeom;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) Nsurf
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nsurf;')"
      stop
    end if
    if (Nsurf > NsurfPD) then
      print "(/,2x,'Input error: Nsurf exceeds NsurfPD;')"                                    
      stop
    end if          
    do i = 1, Nsurf
      read (iInputPARTSUB, *, iostat = ios) surf(i)
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable surf;')"
        stop
      end if
    end do
    read (iInputPARTSUB, *, iostat = ios) Nparam
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nparam;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) z0
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable z0;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) anorm
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable anorm;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) Rcirc
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Rcirc;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) miror
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable miror;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) perfectcond
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable perfectcond;')"
      stop
    end if
  else
    print "(/,2x,'Group name GeomProp not found;')"
    stop  
  end if
  call check_geomAXSYM (TypeGeom, Nsurf, Nparam)
  call check_geomAXSYMOblate (TypeGeom, Nsurf, surf) 
  call check_anorm (anorm)  
  xpart = k * anorm
  snorm = Pi * xpart * xpart
!
  DoConvTest = .true.
  string     = 'ConvTest'
  if (XFindPar (iInputPARTSUB, string)) then
    read (iInputPARTSUB, *, iostat = ios) DoConvTest
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable DoConvTest;')"
      stop
    end if        
  else
    print "(/,2x,'Group name ConvTest not found;')"
    stop  
  end if 
!
  if (DoConvTest) then
    print "(/,2x,'Convergence Test for a Particle on or Near a Plane Surface')"
    print "(  2x,'----------------------------------------------------------')"
  else
    print "(/,2x, a)",                                                              &
   'Convergence Test for a Particle on or Near a Plane Surface over Mrank'
    print "(  2x, a)",                                                              &
   '---------------------------------------------------------------------'
  end if   
! 
  x = k * Rcirc
  NrankW = int(x + 4.05_O * x**0.33_O + 2._O)
  if (.not. DoConvTest) then    
    Nint   = 60
    NXint  = 60 
    Nrank  = 10
    string = 'NintNrank'
    if (XFindPar (iInputPARTSUB, string)) then
      read (iInputPARTSUB, *, iostat = ios) Nint
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Nint;')"
        stop
      end if
      read (iInputPARTSUB, *, iostat = ios) NXint
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable NXint;')"
        stop
      end if
      read (iInputPARTSUB, *, iostat = ios) Nrank
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Nrank;')"
        stop
      end if
    else
      print "(/,2x,'Group name NintNrank not found;')"
      stop  
    end if 
    print "(/,2x,'Input values:')"
    print "(  2x, a, i4, a, i4, a, i3, a)",                                         &
   'the input values of Nint, NXint and Nrank are ', Nint, ', ', NXint, ' and ',    &
    Nrank, ',' 
    print "(  2x, a, i3, a)",                                                       &
   'while the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';'                                             
  else          
    print "(/,2x,'Nrank estimate:')"                                                    
    print "(  2x, a, i3, a)",                                                       &
   'the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';'
    print "(/,2x,'- enter the estimated values of Nint, NXint and Nrank, where  ')"
    print "(  2x,'  Nint = Ndgs * Nrank, Ndgs = 10, 15,..., and NXint = 60, 70, ...;')"         
    call read_integer3 (Nint, NXint, Nrank) 
  end if
!
  if (DoConvTest) then
    print "(/,2x, a)",                                                              &
   '- enter the type of convergence test: 1 - Nint/NXint, 2 - Nrank, 3 - Mrank;'
    call read_integerbound (TypeConvTest, 1, 3)
  else
    TypeConvTest = 3 
  end if
!
  epsNint  = 5.e-2_O
  epsNrank = 5.e-2_O
  epsMrank = 5.e-2_O
  dNint    = 4
  string   = 'Errors'
  if (XFindPar (iInputPARTSUB, string)) then
    read (iInputPARTSUB, *, iostat = ios) epsNint
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsNint;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) epsNrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsNrank;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) epsMrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsMrank;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) dNint
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable dNint;')"
      stop
    end if
  else
    print "(/,2x,'Group name Errors not found;')"
    stop  
  end if  
!
  PrnProgress = .true.
  string      = 'PrintProgress'
  if (XFindPar (iInputPARTSUB, string)) then
    read (iInputPARTSUB, *, iostat = ios) PrnProgress
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable PrnProgress;')"
      stop
    end if        
  else
    print "(/,2x,'Group name PrintProgress not found;')"
    stop  
  end if     
  close (unit = iInputPARTSUB)  
end subroutine readinputPARTSUB
!***********************************************************************************
subroutine printinputPARTSUB (TypeGeom, TypeScat, Nsurf, Nparam, dNint, wavelength, &
           anorm, Rcirc, z0, surf, epsNint, epsNrank, epsMrank, ind_refPART,        &
           ind_refSUB, miror, perfectcond)             
  use parameters
  implicit none
  integer    :: TypeGeom, TypeScat, Nsurf, Nparam, dNint, i
  real(O)    :: wavelength, anorm, Rcirc, surf(Nsurf), epsNint, epsNrank,           &
                epsMrank, z0
  complex(O) :: ind_refPART, ind_refSUB
  logical    :: miror, perfectcond
!
  write (iOutput,"(/,2x,'Input Parameters of the Scattering Problem:',/)")
  write (iOutput,"(2x, a, 1pe13.4, a)")                                             &
 'wavelength of the ambient medium, wavelength = ', wavelength, ';'
  write (iOutput,"(2x, a, 1pe10.3, ',', 1pe10.3, a)")                               &
 'relative refractive index of the substrate, ind_refSUB = (', ind_refSUB, ');' 
  write (iOutput,"(2x, a, 1pe10.3, ',', 1pe10.3, a)")                               &
 'relative refractive index of the particle, ind_refPART = (', ind_refPART, ');'
  write (iOutput,*)
  write (iOutput,"(2x,'index of the particle geometry, TypeGeom = ',i2,';')")       &
         TypeGeom
  if (TypeGeom == 1) then
    write (iOutput,"(2x,'spheroid;')")  
  else if (TypeGeom == 2) then
    write (iOutput,"(2x,'cylinder;')")  
  end if
  write (iOutput,"(2x,'number of surface parameters, Nsurf = ',i2,';')") Nsurf
  write (iOutput,"(2x,'surface parameters:')")
  do i = 1, Nsurf
    write (iOutput,"(2x,'surf(',i2,') = ',1pe10.3,',')") i, surf(i)
  end do
  if (TypeGeom == 1) then
    write (iOutput,"(2x,'where surf(1) is the half-axis along the symmetry axis')")
    write (iOutput,"(2x,'and   surf(2) is the second half-axis;')")
  else if (TypeGeom == 2) then
    write (iOutput,"(2x,'where surf(1) is the half-length of the cylinder')")
    write (iOutput,"(2x,'and   surf(2) is the cylinder radius;')")
  end if         
  write (iOutput,"(2x,'number of integration surfaces, Nparam = ',i2,';')") Nparam
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'characteristic length of the particle, anorm = ', anorm, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'radius of the circumscribing sphere, Rcirc = ', Rcirc, ';'
  write (iOutput,"(2x,'axial position of the particle, z0 = ',1pe10.3,';')") z0
  if (miror) write (iOutput,"(2x,'mirror symmetric particle;')")  
  if (perfectcond) then
    write (iOutput,"(2x,'perfectly conducting particle;')") 
  else
    write (iOutput,"(2x,'dielectric particle;')")
  end if
  if (TypeScat == 1)  then
    write (iOutput,"(2x, a)")                                                       &
   'illumination by a plane wave traveling in the ambient medium;'          
  else if (TypeScat == 2)  then
    write (iOutput,"(2x,'illumination by a plane wave traveling in the substrate;')")   
  end if  
  write (iOutput,*)  
  write (iOutput,"(2x,'integration tolerance, epsNint = ',1pe10.3,';')") epsNint
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum expansion order tolerance, epsNrank = ', epsNrank, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum azimuthal order tolerance, epsMrank = ', epsMrank, ';'          
  write (iOutput,"(2x,'integration step, dNint = ',i4,'.')") dNint
  write (iOutput,"(/)")                       
end subroutine printinputPARTSUB
! **********************************************************************************
subroutine convergence_NintPARTSUB (TypeGeom, TypeScat, k, ind_ref_s, ind_ref_p,    &
           z0, snorm, Nsurf, surf, Nparam, Nrank, Nint, NXint, dNint, miror,        &
           perfectcond, epsNint, PrnProgress) 
  use parameters             
  implicit none 
  integer    :: Nsurf, Nparam, Nrank, Nint, NXint, dNint, TypeGeom, TypeScat   
  real(O)    :: surf(Nsurf), snorm, k, epsNint, z0
  complex(O) :: ind_ref_s, ind_ref_p
  logical    :: miror, perfectcond, PrnProgress 
!  
  integer    :: iNint, Mstart, Nteta, Mrank, Nmax, Nmaxmax, i, m, j, NthetaConv,    &
                Nface
  real(O)    :: beta, alfap, phiGS, kb 
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG(:,:), weightsG(:,:), Xint(:), pondereX(:),       &
                            h(:), v(:), oldh(:), oldv(:), zRe(:), zIm(:),           &
                            rp(:,:), np(:,:), area(:)
  complex(O),allocatable :: a(:,:), b(:,:), t(:,:), c(:), c1(:), cc(:), em(:,:)
!
  phiGS  = 0._O
  Nteta  = 11       
  beta   = 0._O
  alfap  = Pi / 4._O
  Mstart = 1
  Mrank  = 1
  Nmax   = Nrank
  m  = 1
  write (iOutput,"(3x,'--- Convergence Test over Nint and NXint ---',/)")  
  allocate (a(2*Nrank,2*Nrank), b(2*Nrank,2*Nrank), t(2*Nrank,2*Nrank),             &
            c(2*Nrank), c1(2*Nrank), zRe(Nrank), zIm(Nrank),                        &
            rp(2,NfacePD), np(2,NfacePD), area(NfacePD))
  Nmaxmax = Nrank + Mrank*(2*Nrank - Mrank + 1)
  allocate (cc(2*Nmaxmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta), em(3,Nteta))
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do 
  do i = 1, Nrank
    zRe(i) = 0._O
    zIm(i) = 0._O
  end do
  Nface = 1
  do i = 1, NfacePD
    do j = 1, 2
      rp(j,i) = 0._O
      np(j,i) = 0._O
    end do
    area(i) = 0._O
  end do 
  kb = 0._O 
  if (PrnProgress) call write_progress (.true., 1, 11)   
  do iNint = 1, 2     
    allocate (paramG(Nparam,Nint), weightsG(Nparam,Nint), Nintparam(Nparam))
    call interpolation_listAXSYM (TypeGeom, Nsurf, surf, Nint, Nparam, Nintparam,   &
         paramG, weightsG, miror)      
    allocate (Xint(NXint), pondereX(NXint))
    call Laguerre (NXint, Xint, pondereX)       
    call matrix_Q_m (.false., TypeGeom, 3, 1, k, ind_ref_p, Nsurf, surf,            &
         rp, np, area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam,    &
         paramG, weightsG, miror, perfectcond, .false., .false., kb, a, Nrank, Nrank)
    if (PrnProgress) call write_progress (.false., 2+5*(iNint-1), 11)
    call matrix_Q_m (.false., TypeGeom, 1, 1, k, ind_ref_p, Nsurf, surf,            &
         rp, np, area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam,    &
         paramG, weightsG, miror, perfectcond, .false., .false., kb, b, Nrank, Nrank)
    if (PrnProgress) call write_progress (.false., 3+5*(iNint-1), 11)
    call LU_SYSTEM (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nmax)        
    if (PrnProgress) call write_progress (.false., 4+5*(iNint-1), 11) 
    call copy_matrix (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, t, 2*Nrank, 2*Nrank)
    call matrixPARTSUB (k, ind_ref_s, z0, m, Nrank, Nmax, NXint, Xint, pondereX,    &
         a, Nrank, Nrank)   
    call product_matrixSUB (b, 2*Nrank, 2*Nrank, a, 2*Nrank, 2*Nrank, 2*Nmax)  
    if (PrnProgress) call write_progress (.false., 5+5*(iNint-1), 11)
    call LU_SYSTEM_DIRECT (b, 2*Nrank, 2*Nrank, t, 2*Nrank, 2*Nrank, 2*Nmax, 2*Nmax)
    if (PrnProgress) call write_progress (.false., 6+5*(iNint-1), 11)
    if (TypeScat == 1) then                       
      call PWcoefficientsPARTSUB (beta, alfap, m, Nrank, Nmax, c)
      call PWcoefficientsPARTSUBrefl (beta, alfap, z0, k, ind_ref_s, m, Nrank,      &
           Nmax, c1)
      call sum_vectors (c, c1, 2*Nmax)    
    else  
      call PWcoefficientsPARTSUBtrans (beta, alfap, z0, k, ind_ref_s, m, Nrank,     &
           Nmax, c)
    end if
    call product_matrix_vector (2*Nmax, 2*Nmax, t, 2*Nrank, 2*Nrank, c, c1)
    call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)         
    call matrix_m_negativ (Nmax, Nmax, t, Nrank, Nrank) 
    if (TypeScat == 1) then   
      call PWcoefficientsPARTSUB (beta, alfap, -m, Nrank, Nmax, c)
      call PWcoefficientsPARTSUBrefl (beta, alfap, z0, k, ind_ref_s, -m, Nrank,     &
           Nmax, c1)
      call sum_vectors (c, c1, 2*Nmax)    
    else
      call PWcoefficientsPARTSUBtrans (beta, alfap, z0, k, ind_ref_s, -m, Nrank,    &
           Nmax, c)
    end if
    call product_matrix_vector (2*Nmax, 2*Nmax, t, 2*Nrank, 2*Nrank, c, c1)
    call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)                                              
    do j = 1, 3
      do i = 1, Nteta        
        em(j,i) = zero 
      end do
    end do               
    call DSCSPARTSUB (1, cc, k, z0, ind_ref_s, Mrank, Nrank, Nmaxmax, Nteta, phiGS, &
         snorm, em,.true., h, v)
    call DSCSPARTSUB (2, cc, k, z0, ind_ref_s, Mrank, Nrank, Nmaxmax, Nteta, phiGS, &
         snorm, em,.true., h, v)                       
    call delta_DSCSPARTSUB (m, Nteta, h, v, oldh, oldv, epsNint, NthetaConv)
    write (iOutput, "(1x, a, i4, a, i4, a, i3, a, i3, /)")                          &
   'Nint = ', Nint, ', NXint = ', NXint, ', Nrank = ', Nrank, ', Mrank = ', m    
    call write_DSCSPARTSUB (Nteta, h, v)
    Nint  = Nint + dNint
    NXint = NXint + dNint
    deallocate (paramG, weightsG, Nintparam, Xint, pondereX)
  end do  
  write (iOutput,"(3x, a, i2, a, /, 3x, a, 1f5.2, a, /)")                           &
 '--- the solution converges in ', NthetaConv, ' points ---',                       &
 '--- with an relative error of', 100*epsNint, ' %    ---' 
  if (NthetaConv >= int(0.8*(Nteta - 2))) then
    print "(/,2x,'Convergence criterion for Nint and NXint is satisfied;')" 
  else
    print "(/,2x,'Convergence criterion for Nint and NXint is not satisfied;')"                 
  end if
  deallocate (a, b, t, c, c1, cc, h, v, oldh, oldv, em, zRe, zIm, rp, np, area)  
end subroutine convergence_NintPARTSUB 
! **********************************************************************************
subroutine convergence_NrankPARTSUB (TypeGeom, TypeScat, k, ind_ref_s, ind_ref_p,   &
           z0, snorm, Nsurf, surf, Nparam, Nrank, Nint, NXint, miror, perfectcond,  &
           epsNrank, PrnProgress) 
  use parameters                 
  implicit none 
  integer    :: Nsurf, Nparam, Nrank, Nint, NXint, TypeGeom, TypeScat        
  real(O)    :: surf(Nsurf), snorm, k, epsNrank, z0
  complex(O) :: ind_ref_s, ind_ref_p
  logical    :: miror, perfectcond, PrnProgress 
!  
  integer    :: Mstart, Nteta, Mrank, Nmax, Nmaxmax, i, m, j, NthetaConv, Nface
  real(O)    :: beta, alfap, phiGS, kb
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG(:,:), weightsG(:,:), Xint(:), pondereX(:),       &
                            h(:), v(:), oldh(:), oldv(:), zRe(:), zIm(:),           &
                            rp(:,:), np(:,:), area(:)
  complex(O),allocatable :: a(:,:), b(:,:), t(:,:), c(:), c1(:), cc(:), em(:,:),    &
                            a1(:,:), b1(:,:), a0(:,:)
!
  phiGS  = 0._O
  Nteta  = 11 
  beta   = 0._O
  alfap  = Pi / 4._O
  Mstart = 1
  Mrank  = 1
  Nmax   = Nrank
  m  = 1  
  call write_TypeConvHead (2)
  allocate (a(2*Nrank,2*Nrank), b(2*Nrank,2*Nrank), t(2*Nrank,2*Nrank),             &
            c(2*Nrank), c1(2*Nrank), zRe(Nrank), zIm(Nrank),                        &
            rp(2,NfacePD), np(2,NfacePD), area(NfacePD))    
  allocate (a1(2*Nrank,2*Nrank), b1(2*Nrank,2*Nrank), a0(2*Nrank,2*Nrank))
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  allocate (cc(2*Nmaxmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta), em(3,Nteta))        
  allocate (paramG(Nparam,Nint), weightsG(Nparam,Nint), Nintparam(Nparam))
  call interpolation_listAXSYM (TypeGeom, Nsurf, surf, Nint, Nparam, Nintparam,     &
       paramG, weightsG, miror)      
  allocate (Xint(NXint), pondereX(NXint))
  call Laguerre (NXint, Xint, pondereX)           
  do i = 1, Nrank
    zRe(i) = 0._O
    zIm(i) = 0._O
  end do
  Nface = 1
  do i = 1, NfacePD
    do j = 1, 2
      rp(j,i) = 0._O
      np(j,i) = 0._O
    end do
    area(i) = 0._O
  end do
  kb = 0._O  
  if (PrnProgress) call write_progress (.true., 1, 10)
! --- Nrank configuration ---            
  call matrix_Q_m (.false., TypeGeom, 3, 1, k, ind_ref_p, Nsurf, surf,              &
       rp, np, area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam,      &
       paramG, weightsG, miror, perfectcond, .false., .false., kb, a, Nrank, Nrank)
  if (PrnProgress) call write_progress (.false., 2, 10)
  call copy_matrix (2*Nmax, 2*Nmax, a, 2*Nrank, 2*Nrank, a1, 2*Nrank, 2*Nrank)  
  call matrix_Q_m (.false., TypeGeom, 1, 1, k, ind_ref_p, Nsurf, surf,              &
       rp, np, area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam,      &
       paramG, weightsG, miror, perfectcond, .false., .false., kb, b, Nrank, Nrank)
  if (PrnProgress) call write_progress (.false., 3, 10)
  call copy_matrix (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, b1, 2*Nrank, 2*Nrank)
  call LU_SYSTEM (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nmax)          
  if (PrnProgress) call write_progress (.false., 4, 10)
  call copy_matrix (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, t, 2*Nrank, 2*Nrank)
  call matrixPARTSUB (k, ind_ref_s, z0, m, Nrank, Nmax, NXint, Xint, pondereX, a,   &
       Nrank, Nrank)                                                                                  
  call copy_matrix (2*Nmax, 2*Nmax, a, 2*Nrank, 2*Nrank, a0, 2*Nrank, 2*Nrank)
  call product_matrixSUB (b, 2*Nrank, 2*Nrank, a, 2*Nrank, 2*Nrank, 2*Nmax)      
  if (PrnProgress) call write_progress (.false., 5, 10)
  call LU_SYSTEM_DIRECT (b, 2*Nrank, 2*Nrank, t, 2*Nrank, 2*Nrank, 2*Nmax, 2*Nmax)
  if (PrnProgress) call write_progress (.false., 6, 10)
  if (TypeScat == 1) then                         
    call PWcoefficientsPARTSUB (beta, alfap, m, Nrank, Nmax, c)
    call PWcoefficientsPARTSUBrefl (beta, alfap, z0, k, ind_ref_s, m, Nrank,        &
         Nmax, c1)
    call sum_vectors (c, c1, 2*Nmax)   
  else
    call PWcoefficientsPARTSUBtrans (beta, alfap, z0, k, ind_ref_s, m, Nrank,       &
         Nmax, c)   
  end if
  call product_matrix_vector (2*Nmax, 2*Nmax, t, 2*Nrank, 2*Nrank, c, c1)
  call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)         
  call matrix_m_negativ (Nmax, Nmax, t, Nrank, Nrank)
  if (TypeScat == 1) then
    call PWcoefficientsPARTSUB (beta, alfap, -m, Nrank, Nmax, c)
    call PWcoefficientsPARTSUBrefl (beta, alfap, z0, k, ind_ref_s, -m, Nrank,       &
         Nmax, c1)
    call sum_vectors (c, c1, 2*Nmax)  
  else
    call PWcoefficientsPARTSUBtrans (beta, alfap, z0, k, ind_ref_s, -m, Nrank,      &
         Nmax, c)
  end if
  call product_matrix_vector (2*Nmax, 2*Nmax, t, 2*Nrank, 2*Nrank, c, c1)
  call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)                                            
  do j = 1, 3
    do i = 1, Nteta
      em(j,i) = zero
    end do
  end do          
  call DSCSPARTSUB (1, cc, k, z0, ind_ref_s, Mrank, Nrank, Nmaxmax, Nteta, phiGS,   &
       snorm, em,.true., h, v)
  call DSCSPARTSUB (2, cc, k, z0, ind_ref_s, Mrank, Nrank, Nmaxmax, Nteta, phiGS,   &
       snorm, em,.true., h, v)                   
  do i = 1, Nteta
    oldh(i) = h(i)
    oldv(i) = v(i)
  end do  
  write (iOutput, "(1x, a, i4, a, i4, a, i3, a, i3, /)")                            &
 'Nint = ', Nint, ', NXint = ', NXint, ', Nrank = ', Nrank, ', Mrank = ', m 
  call write_DSCSPARTSUB (Nteta, h, v)  
! --- (Nrank - 1) configuration ---            
  call copy_matrix (2*Nmax, 2*Nmax, a1, 2*Nrank, 2*Nrank, a, 2*Nrank, 2*Nrank)
  call matrix_Nrank_m_left (Nmax, a, Nrank, Nrank)      
  call copy_matrix (2*Nmax, 2*Nmax, b1, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank)
  call matrix_Nrank_m (Nmax, b, Nrank, Nrank)                                      
  if (PrnProgress) call write_progress (.false., 7, 10)
  call LU_SYSTEM (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nmax)        
  if (PrnProgress) call write_progress (.false., 8, 10)
  call copy_matrix (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, t, 2*Nrank, 2*Nrank)  
  call copy_matrix (2*Nmax, 2*Nmax, a0, 2*Nrank, 2*Nrank, a, 2*Nrank, 2*Nrank)
  call matrix_Nrank_m (Nmax, a, Nrank, Nrank)
  call product_matrixSUB (b, 2*Nrank, 2*Nrank, a, 2*Nrank, 2*Nrank, 2*Nmax)        
  if (PrnProgress) call write_progress (.false., 9, 10)
  call LU_SYSTEM_DIRECT (b, 2*Nrank, 2*Nrank, t, 2*Nrank, 2*Nrank, 2*Nmax, 2*Nmax)
  if (PrnProgress) call write_progress (.false., 10, 10)
  if (TypeScat == 1) then                 
    call PWcoefficientsPARTSUB (beta, alfap, m, Nrank, Nmax, c)
    call PWcoefficientsPARTSUBrefl (beta, alfap, z0, k, ind_ref_s, m, Nrank,        &
         Nmax, c1)
    call sum_vectors (c, c1, 2*Nmax)  
  else
    call PWcoefficientsPARTSUBtrans (beta, alfap, z0, k, ind_ref_s, m, Nrank,       &
         Nmax, c)
  end if
  call product_matrix_vector (2*Nmax, 2*Nmax, t, 2*Nrank, 2*Nrank, c, c1)
  call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)         
  call matrix_m_negativ (Nmax, Nmax, t, Nrank, Nrank)
  if (TypeScat == 1) then
    call PWcoefficientsPARTSUB (beta, alfap, -m, Nrank, Nmax, c)
    call PWcoefficientsPARTSUBrefl (beta, alfap, z0, k, ind_ref_s, -m, Nrank,       &
         Nmax, c1)
    call sum_vectors (c, c1, 2*Nmax)  
  else
    call PWcoefficientsPARTSUBtrans (beta, alfap, z0, k, ind_ref_s, -m, Nrank,      &
         Nmax, c)
  end if
  call product_matrix_vector (2*Nmax, 2*Nmax, t, 2*Nrank, 2*Nrank, c, c1)
  call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)                                              
  do j = 1, 3
    do i = 1, Nteta
      em(j,i) = zero
    end do
  end do                
  call DSCSPARTSUB (1, cc, k, z0, ind_ref_s, Mrank, Nrank, Nmaxmax, Nteta, phiGS,   &
       snorm, em,.true., h, v)
  call DSCSPARTSUB (2, cc, k, z0, ind_ref_s, Mrank, Nrank, Nmaxmax, Nteta, phiGS,   &
       snorm, em,.true., h, v)                   
  call delta_DSCSPARTSUB (m, Nteta, h, v, oldh, oldv, epsNrank, NthetaConv)    
  write (iOutput, "(1x, a, i4, a, i4, a, i3, a, i3, /)")                            &
 'Nint = ', Nint, ', NXint = ', NXint, ', Nrank = ', Nrank - 1, ', Mrank = ', m 
  call write_DSCSPARTSUB (Nteta, h, v)
  call write_NrankConvRes (NthetaConv, Nteta - 2, epsNrank)
  deallocate (a1, b1, a0)
  deallocate (paramG, weightsG, Nintparam)
  deallocate (Xint, pondereX)
  deallocate (a, b, t, c, c1, cc, h, v, oldh, oldv, em, zRe, zIm, rp, np, area)      
end subroutine convergence_NrankPARTSUB      
! **********************************************************************************
subroutine convergence_MrankPARTSUB (TypeGeom, TypeScat, k, ind_ref_s, ind_ref_p,   &
           z0, snorm, Nsurf, surf, Nparam, Nrank, Nint, NXint, miror, perfectcond,  &
           epsMrank, PrnProgress) 
  use parameters                 
  implicit none 
  integer       :: Nsurf, Nparam, Nrank, Nint, NXint, TypeGeom, TypeScat     
  real(O)       :: surf(Nsurf), snorm, k, epsMrank, z0
  complex(O)    :: ind_ref_s, ind_ref_p
  logical       :: miror, perfectcond, PrnProgress      
!  
  integer       :: Mstart, Nteta, Mrank, Nmax, Nmaxmax, i, m, j, NthetaConv, Nface
  real(O)       :: beta, alfap, phiGS, kb
  character(80) :: FileTmat
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG(:,:), weightsG(:,:), Xint(:), pondereX(:),       &
                            h(:), v(:), oldh(:), oldv(:), zRe(:), zIm(:),           &
                            rp(:,:), np(:,:), area(:)
  complex(O),allocatable :: a(:,:), b(:,:), t(:,:), c(:), c1(:), cc(:), em(:,:)               
!
  phiGS  = 0._O
  Nteta  = 11        
  beta   = Pi / 4._O
  alfap  = Pi / 4._O
  Mstart = 0
  Mrank  = Nrank  
  FileTmat = '../TMATFILES/TPartSub.dat'
  open (unit = iTmat, file = FileTmat, status = 'replace')
  call  write_HeadFileTmat (Nrank, Nrank) 
  call write_TypeConvHead (3)
  write (iOutput, "(1x, a, i4, a, i4, a, i3, /)")                                   &
 'Nint = ', Nint, ', NXint = ', NXint, ', Nrank = ', Nrank
  allocate (a(2*Nrank,2*Nrank), b(2*Nrank,2*Nrank), t(2*Nrank,2*Nrank),             &
            c(2*Nrank), c1(2*Nrank), zRe(Nrank), zIm(Nrank),                        &
            rp(2,NfacePD), np(2,NfacePD), area(NfacePD))
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  allocate (cc(2*Nmaxmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta), em(3,Nteta))   
  do i = 1, Nteta   
    oldh(i) = 0._O
    oldv(i) = 0._O          
  end do
  do i = 1, Nrank
    zRe(i) = 0._O
    zIm(i) = 0._O
  end do
  Nface = 1
  do i = 1, NfacePD
    do j = 1, 2
      rp(j,i) = 0._O
      np(j,i) = 0._O
    end do
    area(i) = 0._O
  end do
  kb = 0._O        
  allocate (paramG(Nparam,Nint), weightsG(Nparam,Nint), Nintparam(Nparam))
  call interpolation_listAXSYM (TypeGeom, Nsurf, surf, Nint, Nparam, Nintparam,     &
       paramG, weightsG, miror)      
  allocate (Xint(NXint), pondereX(NXint))
  call Laguerre (NXint, Xint, pondereX) 
  Mrank = - 1
  do m = Mstart, Nrank   
    call write_1ConvParam (m)         
    Mrank = Mrank + 1
    if (m == 0) then
      Nmax = Nrank
    else
      Nmax = Nrank - m + 1
    end if
    if (PrnProgress) call write_progress_m (.true., m, 1, 6)            
    call matrix_Q_m (.false., TypeGeom, 3, 1, k, ind_ref_p, Nsurf, surf,            &
         rp, np, area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam,    &
         paramG, weightsG, miror, perfectcond, .false., .false., kb, a, Nrank, Nrank)
    if (PrnProgress) call write_progress_m (.false., m, 2, 6)                                   
    call matrix_Q_m (.false., TypeGeom, 1, 1, k, ind_ref_p, Nsurf, surf,            &
         rp, np, area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam,    &
         paramG, weightsG, miror, perfectcond, .false., .false., kb, b, Nrank, Nrank)
    if (PrnProgress) call write_progress_m (.false., m, 3, 6)
    call LU_SYSTEM (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nmax)           
    if (PrnProgress) call write_progress_m (.false., m, 4, 6)
    call copy_matrix (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, t, 2*Nrank, 2*Nrank)
    call matrixPARTSUB (k, ind_ref_s, z0, m, Nrank, Nmax, NXint, Xint, pondereX,    &
         a, Nrank, Nrank)
    call product_matrixSUB (b, 2*Nrank, 2*Nrank, a, 2*Nrank, 2*Nrank, 2*Nmax)     
    if (PrnProgress) call write_progress_m (.false., m, 5, 6)
    call LU_SYSTEM_DIRECT (b, 2*Nrank, 2*Nrank, t, 2*Nrank, 2*Nrank, 2*Nmax, 2*Nmax)  
    if (PrnProgress) call write_progress_m (.false., m, 6, 6)
    call write_FileTmat (Nrank, Nrank, t)
    if (TypeScat == 1) then                      
      call PWcoefficientsPARTSUB (beta, alfap, m, Nrank, Nmax, c)
      call PWcoefficientsPARTSUBrefl (beta, alfap, z0, k, ind_ref_s, m, Nrank,      &
           Nmax, c1)
      call sum_vectors (c, c1, 2*Nmax)     
    else
      call PWcoefficientsPARTSUBtrans (beta, alfap, z0, k, ind_ref_s, m, Nrank,     &
           Nmax, c)
    end if
    call product_matrix_vector (2*Nmax, 2*Nmax, t, 2*Nrank, 2*Nrank, c, c1)
    call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)         
    if (m /= 0) then
      call matrix_m_negativ (Nmax, Nmax, t, Nrank, Nrank)  
      if (TypeScat == 1) then    
        call PWcoefficientsPARTSUB (beta, alfap, -m, Nrank, Nmax, c)
        call PWcoefficientsPARTSUBrefl (beta, alfap, z0, k, ind_ref_s, -m, Nrank,   &
             Nmax, c1)
        call sum_vectors (c, c1, 2*Nmax)      
      else
        call PWcoefficientsPARTSUBtrans (beta, alfap, z0, k, ind_ref_s, -m, Nrank,  &
             Nmax, c)
      end if
      call product_matrix_vector (2*Nmax, 2*Nmax, t, 2*Nrank, 2*Nrank, c, c1)
      call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)                                          
    end if
    do j = 1, 3
      do i = 1, Nteta
        em(j,i) = zero
      end do
    end do               
    call DSCSPARTSUB (1, cc, k, z0, ind_ref_s, Mrank, Nrank, Nmaxmax, Nteta, phiGS, &
         snorm, em,.true., h, v)
    call DSCSPARTSUB (2, cc, k, z0, ind_ref_s, Mrank, Nrank, Nmaxmax, Nteta, phiGS, &
         snorm, em,.true., h, v)                       
    call delta_DSCSPARTSUB (m, Nteta, h, v, oldh, oldv, epsMrank, NthetaConv)           
    call write_DSCSPARTSUB (Nteta, h, v)
    if (NthetaConv >= int(0.8*(Nteta - 2))) exit 
  end do
  close (unit = iTmat)
  call write_MrankConvRes (NthetaConv, epsMrank)
  if (NthetaConv >= int(0.8*(Nteta - 2))) then
    print "(/,2x,'Convergence criterion for Mrank is satisfied;')"                                  
  else
    print "(/,2x,'Convergence criterion for Mrank is not satisfied;')"
  end if
  call DiffScatCrossSectPARTSUB (TypeScat, k, ind_ref_s, z0, snorm, Nrank, Mrank,   &
       FileTmat)
  deallocate (paramG, weightsG, Nintparam)
  deallocate (Xint, pondereX)
  deallocate (a, b, t, c, c1, cc, h, v, oldh, oldv, em, zRe, zIm, rp, np, area)       
end subroutine convergence_MrankPARTSUB      
