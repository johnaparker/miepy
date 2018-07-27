subroutine TPARTSUBFILM
!------------------------------------------------------------------------------------
! 1. General Considerations                                                         !
! --------------------------                                                        !
! TPARTSUBFILM is a routine for computing the scattering characteristics of a       !
! homogeneous, dielectric naxisymmetric particle on or near a plane surface coated  !
! with a film.                                                                      !
!                                                                                   !
! The particle is illuminated by a vector pane wave traveling in the ambient medium,! 
! or by a vector plane wave traveling in the substrate (Figure 1). In the first     !
! case, the direction of the wave vector is (theta = beta, phi = Pi), while in the  !
! second case, the direction is (theta = Pi - beta, phi = 0), where beta is the     !
! incident angle.                                                                   !
!                                                                                   !
!                                                  beta ^ Z                         !
!                                                    !  | 		                    !                                                        			!
!                                                  * ...| 		                    ! 
!                                                   *   | 		                    !
!           Z ^ substrate              incident - -> *  |  substrate	            !
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
! 2. Input Parameters                                                               !
! --------------------                                                              !
! The parameters specified in the input file "/INPUTFILES/InputPARTSUBFILM.dat" are !
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
! - ind_refFILM (complex) - relative refractive index of the film with respect      !
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
!   NOTE: Due to convergence problems occuring for cylindrical particles, it is     !
!   recommended to use the code only for PROLATE SPHEROIDS.                         !
!                                                                                   !
! - z0 (real) - axial position of the film in the particle coordinate system.       !
!   Note that z0 must be a POSITIVE quantity.                                       !         
!                                                                                   !
! - d (real) - thickness of the film.                                               !
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
!------------------------------------------------------------------------------------
  use parameters
  implicit none 
  integer       :: TypeGeom, TypeScat, Nsurf, Nparam, TypeConvTest, dNint, Nrank,   &
                   Nint, NXint                                   
  real(O)       :: wavelength, z0, d, anorm, surf(NsurfPD), snorm, epsNint,         &
                   epsNrank, epsMrank, k, Rcirc                                       
  complex(O)    :: ind_refPART, ind_refSUB, ind_refFILM
  logical       :: miror, perfectcond, DoConvTest, PrnProgress       
! -----------------------------------------------------------------------------------
!                                 Read the input file                               ! 
! -----------------------------------------------------------------------------------
  call readinputPARTSUBFILM ( wavelength, ind_refPART, ind_refSUB, ind_refFILM,     &
       TypeGeom, TypeScat, Nsurf, surf, Nparam, z0, d, anorm, Rcirc, miror,         &
	   perfectcond, DoConvTest, Nint, NXint, Nrank, epsNint, epsNrank, epsMrank,    &
	   dNint, PrnProgress, TypeConvTest, k, snorm)  
! -----------------------------------------------------------------------------------
!                                      Main                                         !
! -----------------------------------------------------------------------------------   
  open (unit = iOutput, file = FileOutput, status = "replace")
  call printinputPARTSUBFILM (TypeGeom, TypeScat, Nsurf, Nparam, dNint, wavelength, &
       anorm, Rcirc, z0, d, surf, epsNint, epsNrank, epsMrank, ind_refPART,         &
	   ind_refSUB, ind_refFILM, miror, perfectcond)    
  if (DoConvTest) then
    if (TypeConvTest == 1) then
      call convergence_NintPARTSUBFILM (TypeGeom, TypeScat, k, ind_refSUB,          &
	       ind_refPART, ind_refFILM, z0, d, snorm, Nsurf, surf, Nparam, Nrank,      &
		   Nint, NXint, dNint, miror, perfectcond, epsNint, PrnProgress)      
    else if (TypeConvTest == 2) then
      call convergence_NrankPARTSUBFILM (TypeGeom, TypeScat, k, ind_refSUB,         &
	       ind_refPART, ind_refFILM, z0, d, snorm, Nsurf, surf, Nparam, Nrank,      &
		   Nint, NXint, miror, perfectcond, epsNrank, PrnProgress)        
    else 
      call convergence_MrankPARTSUBFILM (TypeGeom, TypeScat, k, ind_refSUB,         &
	       ind_refPART, ind_refFILM, z0, d, snorm, Nsurf, surf, Nparam, Nrank,      &
		   Nint, NXint, miror, perfectcond, epsMrank, PrnProgress)         
    end if
  else
    call convergence_MrankPARTSUBFILM (TypeGeom, TypeScat, k, ind_refSUB,           &
	     ind_refPART, ind_refFILM, z0, d, snorm, Nsurf, surf, Nparam, Nrank,        &
		 Nint, NXint, miror, perfectcond, epsMrank, PrnProgress) 
  end if
  close (unit = iOutput)           
end subroutine TPARTSUBFILM
!***********************************************************************************
subroutine readinputPARTSUBFILM ( wavelength, ind_refPART, ind_refSUB, ind_refFILM, &
           TypeGeom, TypeScat, Nsurf, surf, Nparam, z0, d, anorm, Rcirc, miror,     &
		   perfectcond, DoConvTest, Nint, NXint, Nrank, epsNint, epsNrank, epsMrank,& 
		   dNint, PrnProgress, TypeConvTest, k, snorm)     
  use parameters
  use derived_parameters
  implicit none 
  integer       :: TypeGeom, TypeScat, Nsurf, Nparam, TypeConvTest, dNint, Nrank,   &
                   Nint, NXint, i, ios, NrankW                                    
  real(O)       :: wavelength, z0, d, anorm, surf(NsurfPD), xpart, snorm, epsNint,  &
                   epsNrank, epsMrank, k, Rcirc, x                                       
  complex(O)    :: ind_refPART, ind_refSUB, ind_refFILM 
  character(80) :: string
  logical       :: miror, perfectcond, DoConvTest, PrnProgress, XFindPar     
! -----------------------------------------------------------------------------------
!                         Read the input file FileInputPARTSUBFILM                  ! 
! -----------------------------------------------------------------------------------    
  call DrvParameters 
  open (unit = iInputPARTSUB, file = FileInputPARTSUBFILM, status = "old",          &
        position = "rewind")   
  wavelength  = 0.1_O * 2._O * Pi  
  ind_refPART = (1.5_O,0._O)                                                                    
  ind_refSUB  = (1.3_O,0._O)
  ind_refFILM = (1.3_O,0._O)
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
	read (iInputPARTSUB, *, iostat = ios) ind_refFILM
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refFILM;')"
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
  d      = 0.1_O
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
	read (iInputPARTSUB, *, iostat = ios) d
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable d;')"
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
end subroutine readinputPARTSUBFILM
!***********************************************************************************
subroutine printinputPARTSUBFILM (TypeGeom, TypeScat, Nsurf, Nparam, dNint,         &
           wavelength, anorm, Rcirc, z0, d, surf, epsNint, epsNrank, epsMrank,      &
		   ind_refPART, ind_refSUB, ind_refFILM, miror, perfectcond)             
  use parameters
  implicit none
  integer    :: TypeGeom, TypeScat, Nsurf, Nparam, dNint, i
  real(O)    :: wavelength, anorm, Rcirc, surf(Nsurf), epsNint, epsNrank,           &
                epsMrank, z0, d
  complex(O) :: ind_refPART, ind_refSUB, ind_refFILM
  logical    :: miror, perfectcond
!
  write (iOutput,"(/,2x,'Input Parameters of the Scattering Problem:',/)")
  write (iOutput,"(2x, a, 1pe13.4, a)")                                             &
 'wavelength of the ambient medium, wavelength = ', wavelength, ';'
  write (iOutput,"(2x, a, 1pe10.3, ',', 1pe10.3, a)")                               &
 'relative refractive index of the substrate, ind_refSUB = (', ind_refSUB, ');' 
  write (iOutput,"(2x, a, 1pe10.3, ',', 1pe10.3, a)")                               &
 'relative refractive index of the film,     ind_refFILM = (', ind_refFILM, ');' 
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
  write (iOutput,"(2x,'thickness of the film, d = ',1pe10.3,';')") d
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
end subroutine printinputPARTSUBFILM
! **********************************************************************************
subroutine convergence_NintPARTSUBFILM (TypeGeom, TypeScat, k, ind_ref_s, ind_ref_p,&
           ind_ref_f, z0, d, snorm, Nsurf, surf, Nparam, Nrank, Nint, NXint, dNint, &
		   miror, perfectcond, epsNint, PrnProgress) 
  use parameters             
  implicit none 
  integer    :: Nsurf, Nparam, Nrank, Nint, NXint, dNint, TypeGeom, TypeScat
  real(O)    :: surf(Nsurf), snorm, k, epsNint, z0, d
  complex(O) :: ind_ref_s, ind_ref_p, ind_ref_f
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
    call matrixPARTSUBFILM (k, ind_ref_f, ind_ref_s, d, z0, m, Nrank, Nmax, NXint,  &
	     Xint, pondereX, a, Nrank, Nrank)   
    call product_matrixSUB (b, 2*Nrank, 2*Nrank, a, 2*Nrank, 2*Nrank, 2*Nmax)  
    if (PrnProgress) call write_progress (.false., 5+5*(iNint-1), 11)
    call LU_SYSTEM_DIRECT (b, 2*Nrank, 2*Nrank, t, 2*Nrank, 2*Nrank, 2*Nmax, 2*Nmax)
    if (PrnProgress) call write_progress (.false., 6+5*(iNint-1), 11)
!  
    if (TypeScat == 1) then                     
      call PWcoefficientsPARTSUB (beta, alfap, m, Nrank, Nmax, c)
      call PWcoefficientsPARTSUBreflFILM (beta, alfap, d, z0, k, ind_ref_f,         &
	       ind_ref_s, m, Nrank, Nmax, c1)
      call sum_vectors (c, c1, 2*Nmax)    
    else
      call PWcoefficientsPARTSUBtransFILM (beta, alfap, d, z0, k, ind_ref_f,        &
	       ind_ref_s, m, Nrank, Nmax, c)
	end if
!    
    call product_matrix_vector (2*Nmax, 2*Nmax, t, 2*Nrank, 2*Nrank, c, c1)
    call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)         
    call matrix_m_negativ (Nmax, Nmax, t, Nrank, Nrank) 
!
    if (TypeScat == 1) then
      call PWcoefficientsPARTSUB (beta, alfap, -m, Nrank, Nmax, c)
      call PWcoefficientsPARTSUBreflFILM (beta, alfap, d, z0, k, ind_ref_f,         &
	       ind_ref_s, -m, Nrank, Nmax, c1)
      call sum_vectors (c, c1, 2*Nmax)    
    else
      call PWcoefficientsPARTSUBtransFILM (beta, alfap, d, z0, k, ind_ref_f,        &
	       ind_ref_s, -m, Nrank, Nmax, c)
	end if
!
    call product_matrix_vector (2*Nmax, 2*Nmax, t, 2*Nrank, 2*Nrank, c, c1)
    call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)                                              
    do j = 1, 3
      do i = 1, Nteta        
        em(j,i) = zero 
      end do
    end do		               
    call DSCSPARTSUBFILM (1, cc, k, z0, d, ind_ref_f, ind_ref_s, Mrank, Nrank,      &
	     Nmaxmax, Nteta, phiGS, snorm, em,.true., h, v)
    call DSCSPARTSUBFILM (2, cc, k, z0, d, ind_ref_f, ind_ref_s, Mrank, Nrank,      &
	     Nmaxmax, Nteta, phiGS, snorm, em,.true., h, v)                       
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
end subroutine convergence_NintPARTSUBFILM 
! **********************************************************************************
subroutine convergence_NrankPARTSUBFILM (TypeGeom, TypeScat, k, ind_ref_s,          &
           ind_ref_p, ind_ref_f, z0, d, snorm, Nsurf, surf, Nparam, Nrank,          &
		   Nint, NXint, miror, perfectcond, epsNrank, PrnProgress) 
  use parameters                 
  implicit none 
  integer    :: Nsurf, Nparam, Nrank, Nint, NXint, TypeGeom, TypeScat
  real(O)    :: surf(Nsurf), snorm, k, epsNrank, z0, d
  complex(O) :: ind_ref_s, ind_ref_p, ind_ref_f
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
  call matrixPARTSUBFILM (k, ind_ref_f, ind_ref_s, d, z0, m, Nrank, Nmax, NXint,    &
       Xint, pondereX, a, Nrank, Nrank)                                                                                  
  call copy_matrix (2*Nmax, 2*Nmax, a, 2*Nrank, 2*Nrank, a0, 2*Nrank, 2*Nrank)
  call product_matrixSUB (b, 2*Nrank, 2*Nrank, a, 2*Nrank, 2*Nrank, 2*Nmax)      
  if (PrnProgress) call write_progress (.false., 5, 10)
  call LU_SYSTEM_DIRECT (b, 2*Nrank, 2*Nrank, t, 2*Nrank, 2*Nrank, 2*Nmax, 2*Nmax)
  if (PrnProgress) call write_progress (.false., 6, 10)
!   
  if (TypeScat == 1) then                    
    call PWcoefficientsPARTSUB (beta, alfap, m, Nrank, Nmax, c)
    call PWcoefficientsPARTSUBreflFILM (beta, alfap, d, z0, k, ind_ref_f,           &
	     ind_ref_s, m, Nrank, Nmax, c1)
    call sum_vectors (c, c1, 2*Nmax)   
  else
    call PWcoefficientsPARTSUBtransFILM (beta, alfap, d, z0, k, ind_ref_f,          &
	     ind_ref_s, m, Nrank, Nmax, c)
  end if
!
  call product_matrix_vector (2*Nmax, 2*Nmax, t, 2*Nrank, 2*Nrank, c, c1)
  call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)         
  call matrix_m_negativ (Nmax, Nmax, t, Nrank, Nrank)
!
  if (TypeScat == 1) then
    call PWcoefficientsPARTSUB (beta, alfap, -m, Nrank, Nmax, c)
    call PWcoefficientsPARTSUBreflFILM (beta, alfap, d, z0, k, ind_ref_f,           &
	     ind_ref_s, -m, Nrank, Nmax, c1)
    call sum_vectors (c, c1, 2*Nmax)  
  else
    call PWcoefficientsPARTSUBtransFILM (beta, alfap, d, z0, k, ind_ref_f,          &
	     ind_ref_s, -m, Nrank, Nmax, c)
	end if
!
  call product_matrix_vector (2*Nmax, 2*Nmax, t, 2*Nrank, 2*Nrank, c, c1)
  call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)                                            
  do j = 1, 3
    do i = 1, Nteta
      em(j,i) = zero
    end do
  end do          
  call DSCSPARTSUBFILM (1, cc, k, z0, d, ind_ref_f, ind_ref_s, Mrank, Nrank,        &
       Nmaxmax, Nteta, phiGS, snorm, em,.true., h, v)
  call DSCSPARTSUBFILM (2, cc, k, z0, d, ind_ref_f, ind_ref_s, Mrank, Nrank,        &
       Nmaxmax, Nteta, phiGS, snorm, em,.true., h, v)                   
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
!
  if (TypeScat == 1) then
    call PWcoefficientsPARTSUB (beta, alfap, m, Nrank, Nmax, c)
    call PWcoefficientsPARTSUBreflFILM (beta, alfap, d, z0, k, ind_ref_f,           &
	     ind_ref_s, m, Nrank, Nmax, c1)
    call sum_vectors (c, c1, 2*Nmax)  
  else
    call PWcoefficientsPARTSUBtransFILM (beta, alfap, d, z0, k, ind_ref_f,          &
	     ind_ref_s, m, Nrank, Nmax, c)
  end if
!
  call product_matrix_vector (2*Nmax, 2*Nmax, t, 2*Nrank, 2*Nrank, c, c1)
  call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)         
  call matrix_m_negativ (Nmax, Nmax, t, Nrank, Nrank)
!
  if (TypeScat == 1) then
    call PWcoefficientsPARTSUB (beta, alfap, -m, Nrank, Nmax, c)
    call PWcoefficientsPARTSUBreflFILM (beta, alfap, d, z0, k, ind_ref_f,           &
	     ind_ref_s, -m, Nrank, Nmax, c1)
    call sum_vectors (c, c1, 2*Nmax)  
  else
    call PWcoefficientsPARTSUBtransFILM (beta, alfap, d, z0, k, ind_ref_f,          &
	     ind_ref_s, -m, Nrank, Nmax, c)
  end if
!
  call product_matrix_vector (2*Nmax, 2*Nmax, t, 2*Nrank, 2*Nrank, c, c1)
  call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)                                              
  do j = 1, 3
    do i = 1, Nteta
      em(j,i) = zero
    end do
  end do                
  call DSCSPARTSUBFILM (1, cc, k, z0, d, ind_ref_f, ind_ref_s, Mrank, Nrank,        &
       Nmaxmax, Nteta, phiGS, snorm, em,.true., h, v)
  call DSCSPARTSUBFILM (2, cc, k, z0, d, ind_ref_f, ind_ref_s, Mrank, Nrank,        &
       Nmaxmax, Nteta, phiGS, snorm, em,.true., h, v)                   
  call delta_DSCSPARTSUB (m, Nteta, h, v, oldh, oldv, epsNrank, NthetaConv)    
  write (iOutput, "(1x, a, i4, a, i4, a, i3, a, i3, /)")                            &
 'Nint = ', Nint, ', NXint = ', NXint, ', Nrank = ', Nrank - 1, ', Mrank = ', m 
  call write_DSCSPARTSUB (Nteta, h, v)
  call write_NrankConvRes (NthetaConv, Nteta - 2, epsNrank)
  deallocate (a1, b1, a0)
  deallocate (paramG, weightsG, Nintparam)
  deallocate (Xint, pondereX)
  deallocate (a, b, t, c, c1, cc, h, v, oldh, oldv, em, zRe, zIm, rp, np, area)      
end subroutine convergence_NrankPARTSUBFILM      
! **********************************************************************************
subroutine convergence_MrankPARTSUBFILM (TypeGeom, TypeScat, k, ind_ref_s,          &
           ind_ref_p, ind_ref_f, z0, d, snorm, Nsurf, surf, Nparam, Nrank,          &
		   Nint, NXint, miror, perfectcond, epsMrank, PrnProgress) 
  use parameters                 
  implicit none 
  integer       :: Nsurf, Nparam, Nrank, Nint, NXint, TypeGeom, TypeScat
  real(O)       :: surf(Nsurf), snorm, k, epsMrank, z0, d
  complex(O)    :: ind_ref_s, ind_ref_p, ind_ref_f
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
    call matrixPARTSUBFILM (k, ind_ref_f, ind_ref_s, d, z0, m, Nrank, Nmax, NXint,  &
	     Xint, pondereX, a, Nrank, Nrank)
    call product_matrixSUB (b, 2*Nrank, 2*Nrank, a, 2*Nrank, 2*Nrank, 2*Nmax)     
    if (PrnProgress) call write_progress_m (.false., m, 5, 6)
    call LU_SYSTEM_DIRECT (b, 2*Nrank, 2*Nrank, t, 2*Nrank, 2*Nrank, 2*Nmax, 2*Nmax)  
    if (PrnProgress) call write_progress_m (.false., m, 6, 6)
    call write_FileTmat (Nrank, Nrank, t)
!                      
    if (TypeScat == 1) then
      call PWcoefficientsPARTSUB (beta, alfap, m, Nrank, Nmax, c)
      call PWcoefficientsPARTSUBreflFILM (beta, alfap, d, z0, k, ind_ref_f,         &
	       ind_ref_s, m, Nrank, Nmax, c1)
      call sum_vectors (c, c1, 2*Nmax)     
    else
      call PWcoefficientsPARTSUBtransFILM (beta, alfap, d, z0, k, ind_ref_f,        &
	       ind_ref_s, m, Nrank, Nmax, c)
	end if
!
    call product_matrix_vector (2*Nmax, 2*Nmax, t, 2*Nrank, 2*Nrank, c, c1)
    call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)         
    if (m /= 0) then
      call matrix_m_negativ (Nmax, Nmax, t, Nrank, Nrank)  
!
      if (TypeScat == 1) then
        call PWcoefficientsPARTSUB (beta, alfap, -m, Nrank, Nmax, c)
        call PWcoefficientsPARTSUBreflFILM (beta, alfap, d, z0, k, ind_ref_f,       &
	         ind_ref_s, -m, Nrank, Nmax, c1)
        call sum_vectors (c, c1, 2*Nmax)      
      else
        call PWcoefficientsPARTSUBtransFILM (beta, alfap, d, z0, k, ind_ref_f,      &
	         ind_ref_s, -m, Nrank, Nmax, c)
	  end if
!
      call product_matrix_vector (2*Nmax, 2*Nmax, t, 2*Nrank, 2*Nrank, c, c1)
      call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)                                          
    end if
    do j = 1, 3
      do i = 1, Nteta
        em(j,i) = zero
      end do
    end do               
    call DSCSPARTSUBFILM (1, cc, k, z0, d, ind_ref_f, ind_ref_s, Mrank, Nrank,      &
	     Nmaxmax, Nteta, phiGS, snorm, em,.true., h, v)
    call DSCSPARTSUBFILM (2, cc, k, z0, d, ind_ref_f, ind_ref_s, Mrank, Nrank,      &
	     Nmaxmax, Nteta, phiGS, snorm, em,.true., h, v)                       
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
  call DiffScatCrossSectPARTSUBFILM (TypeScat, k, ind_ref_f, ind_ref_s, d, z0,      &
       snorm, Nrank, Mrank, FileTmat)
  deallocate (paramG, weightsG, Nintparam)
  deallocate (Xint, pondereX)
  deallocate (a, b, t, c, c1, cc, h, v, oldh, oldv, em, zRe, zIm, rp, np, area)       
end subroutine convergence_MrankPARTSUBFILM 
! **********************************************************************************
!                These files must be stored in IncCoef.f90                         *
! **********************************************************************************
subroutine FresnelCoefFilmReflection ( cosb, k, d, z0, m1, m2, Rpar, Rperp)
! ----------------------------------------------------------------------------------
!   Significance of input parameters:                                              ! 
!   - m1 - refractive index of the film,                                           !
!   - m2 - refractive index of the substrate                                       !
!   - cosb - cosine of the incident angle of the wave propagating in the ambient   !                                                  !
! ----------------------------------------------------------------------------------
  use parameters
  implicit none
  real(O)    :: k, d, z0
  complex(O) :: cosb, m1, m2, Rpar, Rperp 
!  
  complex(O) :: sinb, sinb1, cosb1, sinb2, cosb2, argR, argR1, Rb12, Ra12, Rb01, Ra01
! 
  sinb = sqrt(1._O - cosb**2)
!
  sinb1  = sinb / m1
  cosb1  = sqrt(1._O - sinb1**2)
  if (aimag(m1 * cosb1) < 0._O) cosb1 = - cosb1
!
  sinb2  = sinb / m2
  cosb2  = sqrt(1._O - sinb2**2)
  if (aimag(m2 * cosb2) < 0._O) cosb2 = - cosb2  
!
  argR  = 2._O * im * k * z0 * cosb  
  argR1 = 2._O * im * m1 * k * d * cosb1
!
!
  Rb12 = (m2 * cosb1 - m1 * cosb2) / (m2 * cosb1 + m1 * cosb2)    
  Ra12 = (m1 * cosb1 - m2 * cosb2) / (m1 * cosb1 + m2 * cosb2)
!
  Rb01 = (m1 * cosb - cosb1) / (cosb1 + m1 * cosb)  
  Ra01 = (cosb - m1 * cosb1) / (m1 * cosb1 + cosb)
!
  Rpar  = ( Rb01 + Rb12 * exp(argR1) ) / ( 1._O + Rb01 * Rb12 * exp(argR1) )
  Rperp = ( Ra01 + Ra12 * exp(argR1) ) / ( 1._O + Ra01 * Ra12 * exp(argR1) )
!
  Rpar  = Rpar  * exp(argR)
  Rperp = Rperp * exp(argR)
end subroutine FresnelCoefFilmReflection
! **********************************************************************************
subroutine FresnelCoefFilmReflectionPhase ( cosb, k, d, z0, m1, m2, Rpar, Rperp)
  use parameters
  implicit none
  real(O)    :: k, d, z0
  complex(O) :: cosb, m1, m2, Rpar, Rperp 
!  
  complex(O) :: sinb, sinb1, cosb1, sinb2, cosb2, argR1, Rb12, Ra12, Rb01, Ra01
! 
  sinb = sqrt(1._O - cosb**2)
!
  sinb1  = sinb / m1
  cosb1  = sqrt(1._O - sinb1**2)
  if (aimag(m1 * cosb1) < 0._O) cosb1 = - cosb1
!
  sinb2  = sinb / m2
  cosb2  = sqrt(1._O - sinb2**2)
  if (aimag(m2 * cosb2) < 0._O) cosb2 = - cosb2  
!  
  argR1 = 2._O * im * m1 * k * d * cosb1
!
  Rb12 = (m2 * cosb1 - m1 * cosb2) / (m2 * cosb1 + m1 * cosb2)    
  Ra12 = (m1 * cosb1 - m2 * cosb2) / (m1 * cosb1 + m2 * cosb2)
!
  Rb01 = (m1 * cosb - cosb1) / (cosb1 + m1 * cosb)  
  Ra01 = (cosb - m1 * cosb1) / (m1 * cosb1 + cosb)
!
  Rpar  = ( Rb01 + Rb12 * exp(argR1) ) / ( 1._O + Rb01 * Rb12 * exp(argR1) )
  Rperp = ( Ra01 + Ra12 * exp(argR1) ) / ( 1._O + Ra01 * Ra12 * exp(argR1) )
end subroutine FresnelCoefFilmReflectionPhase
! **********************************************************************************
subroutine PWcoefficientsPARTSUBreflFILM (beta0, alphap, d, z0, wavenumber,         &
           ind_refFILM, ind_refSUB, m, Nrank, Nmax, c)
!-----------------------------------------------------------------------------------
! The routine computes the expansion coefficients of the reflected plane wave      !
! (traveling in the Oxz plane) for the azimuthal mode m.                           !
!----------------------------------------------------------------------------------- 
  use parameters
  implicit none
  integer    :: m, Nrank, Nmax
  real(O)    :: beta0, alphap, d, z0, wavenumber      
  complex(O) :: c(2*Nmax), ind_refFILM, ind_refSUB
!
  integer    :: k, n      
  real(O)    :: nm, betaM0, mr
  complex(O) :: Rpar, Rperp, et, ep, cosb0, factc, factp, factt
  real(O),allocatable :: Pnm(:), dPnm(:), pinm(:), taunm(:)
!
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))
  cosb0 = cmplx(cos(beta0),0._O)
  call FresnelCoefFilmReflection ( cosb0, wavenumber, d, z0, ind_refFILM,           &
       ind_refSUB, Rpar, Rperp)
  et     = cos(alphap) * Rpar 
  ep     = sin(alphap) * Rperp
  betaM0 = Pi - beta0
  call leg_normalized (betaM0, abs(m), Nrank, Pnm, dPnm, pinm, taunm)
  mr = real(m,O)  
  do k = 1, Nmax
    if (m == 0) then
      n = k
    else
      n = abs(m) + k - 1
    end if
    nm    = real(2 * n * (n + 1),O)         
    nm    = 1._O / sqrt(nm)     
    factc = 4._O * im**n * nm * cos(mr * Pi)        
    factp = factc * im * mr * pinm(n)
    factt = factc * taunm(n)
    c(k)      = - factp * et - factt * ep
    c(Nmax+k) = - im * (factt * et - factp * ep)
  end do      
  deallocate (Pnm, dPnm, pinm, taunm)      
end subroutine PWcoefficientsPARTSUBreflFILM
! **********************************************************************************
subroutine FresnelCoefFilmTransmission ( cosb, k, d, z0, m1, m, Tpar, Tperp,        &
           sinb2, cosb2)
! ----------------------------------------------------------------------------------
!   Significance of input parameters:                                              ! 
!   - m1 - refractive index of the film,                                           !
!   - m  - refractive index of the substrate                                       !
!   - cosb - cosine of the incident angle of the wave propagating in the substrate !
!   - cosb2 - cosine of the transmitted wave propagating in the ambient            !
!   - sinb2 - sine of the transmitted wave propagating in the ambient              !
! ----------------------------------------------------------------------------------
  use parameters
  implicit none
  real(O)    :: k, d, z0
  complex(O) :: cosb, m1, m, Tpar, Tperp, sinb2, cosb2 
!  
  complex(O) :: sinb, sinb1, cosb1, argT, argR1, Rb12, Ra12, Rb01, Ra01
! 
  sinb = sqrt(1._O - cosb**2)
!
  sinb1  = m * sinb / m1
  cosb1  = sqrt(1._O - sinb1**2)
  if (aimag(m1 * cosb1) < 0._O) cosb1 = - cosb1
!
  sinb2  = m * sinb 
  cosb2  = sqrt(1._O - sinb2**2)
  if (aimag(cosb2) < 0._O) cosb2 = - cosb2  
!
  argT  = im * k * ( d * (m1 * cosb1 - m * cosb) + z0 * (cosb2 - m * cosb) )   
  argR1 = 2._O * im * m1 * k * d * cosb1  
!
  Rb12 = (cosb1 - m1 * cosb2) / (cosb1 + m1 * cosb2)    
  Ra12 = (m1 * cosb1 - cosb2) / (m1 * cosb1 + cosb2)
!
  Rb01 = (m1 * cosb - m * cosb1) / (m * cosb1 + m1 * cosb)  
  Ra01 = (m * cosb - m1 * cosb1) / (m1 * cosb1 + m * cosb)
!
  Tpar  = m * ( 1._O + Rb12) * (1._O + Rb01) / ( 1._O + Rb01 * Rb12 * exp(argR1) )
  Tperp =     ( 1._O + Ra12) * (1._O + Ra01) / ( 1._O + Ra01 * Ra12 * exp(argR1) )
!
  Tpar  = Tpar  * exp(argT)
  Tperp = Tperp * exp(argT)
end subroutine FresnelCoefFilmTransmission
! **********************************************************************************
subroutine PWcoefficientsPARTSUBtransFILM (beta, alphap, d, z0, wavenumber,         &
           ind_refFILM, ind_refSUB, m, Nrank, Nmax, c)
!----------------------------------------------------------------------------------
! The routine computes the expansion coefficients of the transmitted plane wave   !
! (traveling in the Oxz plane) for the azimuthal mode m.                          !
! --------------------------------------------------------------------------------- 
  use parameters
  implicit none
  integer    :: m, Nrank, Nmax
  real(O)    :: beta, alphap, z0, d, wavenumber
  complex(O) :: c(2*Nmax), ind_refFILM, ind_refSUB     
!
  integer    :: k, n
  real(O)    :: nm, mr
  complex(O) :: Tpar, Tperp, et, ep, cosb, sinb0, cosb0, factc, factp, factt
  complex(O),allocatable :: Pnm(:), dPnm(:), pinm(:), taunm(:)
!
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))
  cosb = cmplx(cos(beta),0._O)         
  call FresnelCoefFilmTransmission ( cosb, wavenumber, d, z0, ind_refFILM,          &
       ind_refSUB, Tpar, Tperp, sinb0, cosb0)         
  et    = cos(alphap) * Tpar 
  ep    = sin(alphap) * Tperp 
  cosb0 = - cosb0 
  call Leg_normalized_complex (sinb0, cosb0, abs(m), Nrank, Pnm, dPnm, pinm, taunm)
  mr = real(m,O) 
  do k = 1, Nmax
    if (m == 0) then
      n = k
    else
      n = abs(m) + k - 1
    end if
    nm    = real(2 * n * (n + 1),O)         
    nm    = 1._O / sqrt(nm)     
    factc = 4._O * im**n * nm
    factp = factc * im * mr * pinm(n)
    factt = factc * taunm(n)
    c(k)      = - factp * et - factt * ep
    c(Nmax+k) = - im * (factt * et - factp * ep)
  end do   
  deallocate (Pnm, dPnm, pinm, taunm)    
end subroutine PWcoefficientsPARTSUBtransFILM
!***********************************************************************************
!                       This file must be stored in Proces3.f90                    * 
!***********************************************************************************
subroutine matrixPARTSUBFILM (k, ind_refFILM, ind_refSUB, dist, z0, m, Nrank, Nmax, &
           Nx, x, pond, A, na, ma)
!-----------------------------------------------------------------------------------
! The routine computes the reflection matrix for the scattering of a particle      !
! near or on a plane surface coated with a film.                                   ! 
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: m, Nrank, Nmax, Nx, na, ma
  real(O)    :: k, dist, z0, x(Nx), pond(Nx)
  complex(O) :: ind_refFILM, ind_refSUB, A(2*na,2*ma)
!
  integer    :: i, j, pint, n, n1
  real(O)    :: D, nm, n1m, mr
  complex(O) :: fact, t, Rpar, Rperp, q, cosbM0, sinbM0, sinb0, cosb0, f, fpp, ftp, &
                fpt, ftt
  complex(O),allocatable :: Pnm(:), dPnm(:), pinm(:), taunm(:),                     &
                            PnmM(:), dPnmM(:), pinmM(:), taunmM(:)      
!
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))
  allocate (PnmM(0:Nrank), dPnmM(0:Nrank), pinmM(0:Nrank), taunmM(0:Nrank))
  do i = 1, 2*Nmax
    do j = 1, 2*Nmax
      A(i,j) = zero
    end do
  end do                           
  q = cmplx(k * z0,0.0,O)
  q = 2._O * im * q   
  if (abs(q) < MachEps) q = cmplx(MachEps,MachEps,O)
  mr = real(m,O)
  do pint = 1, Nx
    t     = one - x(pint) / q
    fact  = 4._O * exp(q) * pond(pint) / q    
    cosb0 = t
    sinb0 = sqrt(one - t * t)
    call Leg_normalized_complex (sinb0, cosb0, abs(m), Nrank, Pnm, dPnm, pinm, taunm) 
	call FresnelCoefFilmReflectionPhase ( cosb0, k, dist, z0, ind_refFILM,          &
	     ind_refSUB, Rpar, Rperp)                        							     
    cosbM0 = - t
    sinbM0 = sqrt(one - t * t)
    call Leg_normalized_complex (sinbM0, cosbM0, abs(m), Nrank, PnmM, dPnmM,        &
         pinmM, taunmM)                                                 
    do i = 1, Nmax
      if (m == 0) then
        n1 = i
      else
        n1 = abs(m) + i - 1
      end if
      n1m = real(2 * n1 * (n1 + 1),O)
      n1m = 1._O / sqrt(n1m)
      do j = 1, Nmax
        if (m == 0) then
          n = j
        else
          n = abs(m) + j - 1
        end if
        nm  = real(2 * n * (n + 1),O)
        nm  = 1._O / sqrt(nm)
        D   = nm * n1m
        f   = fact * D * im**(n1 - n)
        fpp = f * pinm(n)  * pinmM(n1)
        ftt = f * taunm(n) * taunmM(n1)
        fpt = f * pinm(n)  * taunmM(n1)
        ftp = f * taunm(n) * pinmM(n1)        
        A(i,j) = A(i,j) + (mr**2 * fpp * Rpar + ftt * Rperp)        
        A(i+Nmax,j) = A(i+Nmax,j) + mr * (fpt * Rpar + ftp * Rperp)
        A(i,j+Nmax) = A(i,j+Nmax) + mr * (ftp * Rpar + fpt * Rperp)
        A(i+Nmax,j+Nmax) = A(i+Nmax,j+Nmax) + (mr**2 * fpp * Rperp + ftt * Rpar)
      end do
    end do
  end do           
  deallocate (Pnm, dPnm, pinm, taunm, PnmM, dPnmM, pinmM, taunmM)        
end subroutine matrixPARTSUBFILM
! **********************************************************************************
!                These files must be stored in PostProces1.f90                     *
! **********************************************************************************
subroutine DiffScatCrossSectPARTSUBFILM (TypeScat, k, ind_ref_f, ind_ref_s, d, z0,  &
           snorm, Nrank, Mrank, FileTmat)                           
  use parameters
  implicit none 
  integer       :: Nrank, Mrank, TypeScat
  real(O)       :: snorm, k, z0, d
  complex(O)    :: ind_ref_s, ind_ref_f
  character(80) :: FileTmat  
!  
  integer       :: NthetaGS, Nphi, Ntheta(NphiMax), Mstart, Nmax, Nmaxmax,          &
                   m, i, itheta, iphi, ntl, mtl, NphiAL, NthetaAL
  real(O)       :: beta, alphap, phiGS, phi(NphiMax), thetamin(Nphimax),            &
                   thetamax(NphiMax)
  character(80) :: FileDSCS, FileEMF
  logical       :: ComputeDSCS, ComputeFields, normalized, WriteInputInfo
  real(O),allocatable    :: h(:), v(:)
  complex(O),allocatable :: t(:,:), c(:), c1(:), cc(:), em(:,:), emf(:,:,:)
! -----------------------------------------------------------------------------------
!                                 Read the input file                               !
! -----------------------------------------------------------------------------------               
  call readinputPARTSUBFILM1 ( beta, alphap, ComputeDSCS, ComputeFields, NthetaGS,  &
       phiGS, Nphi, phi, Ntheta, thetamin, thetamax, normalized, FileDSCS, FileEMF, &
       WriteInputInfo )  
  if (.not. ComputeDSCS .and. .not. ComputeFields) return
! -----------------------------------------------------------------------------------
!                                       Main                                        !
! ----------------------------------------------------------------------------------- 
  print "(/,2x,'Scattering Characteristics of a Particle on or Near a Plane Surface')"  
  print "(  2x,'-------------------------------------------------------------------')"  
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)  
  allocate (c(2*Nrank), c1(2*Nrank), cc(2*Nmaxmax))        
  open (unit = iTmat, file = FileTmat, status = "old",  position = "rewind")   
  call read_HeadFileTmat (ntl, mtl)
  call check_dimensionMat (ntl, mtl, Nrank)             
  allocate (t(2*ntl, 2*mtl))
  Mstart = 0
  do m = Mstart, Mrank    
    if (m == 0) then
      Nmax = Nrank
    else
      Nmax = Nrank - m + 1
    end if    
    call read_FileTmat (ntl, mtl, t)
! 
    if (TypeScat == 1) then                        
      call PWcoefficientsPARTSUB (beta, alphap, m, Nrank, Nmax, c)
      call PWcoefficientsPARTSUBreflFILM (beta, alphap, d, z0, k, ind_ref_f,        &
	       ind_ref_s, m, Nrank, Nmax, c1)
      call sum_vectors (c, c1, 2*Nmax)     
    else
      call PWcoefficientsPARTSUBtransFILM (beta, alphap, d, z0, k, ind_ref_f,       &
	       ind_ref_s, m, Nrank, Nmax, c)
	end if
!
    call product_matrix_vector (2*Nmax, 2*Nmax, t, 2*ntl, 2*mtl, c, c1)
    call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)         
    if (m /= 0) then
      call matrix_m_negativ (Nmax, Nmax, t, ntl, mtl)  
!    
      if (TypeScat == 1) then
        call PWcoefficientsPARTSUB (beta, alphap, - m, Nrank, Nmax, c)
        call PWcoefficientsPARTSUBreflFILM (beta, alphap, d, z0, k, ind_ref_f,      &
	         ind_ref_s, - m, Nrank, Nmax, c1)
        call sum_vectors (c, c1, 2*Nmax)      
      else
        call PWcoefficientsPARTSUBtransFILM (beta, alphap, d, z0, k, ind_ref_f,     &
	         ind_ref_s, -m, Nrank, Nmax, c)
	  end if
!
      call product_matrix_vector (2*Nmax, 2*Nmax, t, 2*ntl, 2*mtl, c, c1)
      call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)                                          
    end if
  end do
  close (unit = iTmat)
  deallocate (t, c, c1) 
  if (ComputeDSCS) then
    open (unit = iDSCS, file = FileDSCS, status = "replace")
    allocate (h(NthetaGS), v(NthetaGS), em(3,NthetaGS)) 
    do i = 1, 3
      do itheta = 1, NthetaGS
        em(i,itheta) = zero
      end do
    end do         
    call DSCSPARTSUBFILM (1, cc, k, z0, d, ind_ref_f, ind_ref_s, Mrank, Nrank,      &
	     Nmaxmax, NthetaGS, phiGS, snorm, em, normalized, h, v)
    call DSCSPARTSUBFILM (2, cc, k, z0, d, ind_ref_f, ind_ref_s, Mrank, Nrank,      &
	     Nmaxmax, NthetaGS, phiGS, snorm, em, normalized, h, v)                                
    if (WriteInputInfo) call inputDSCS_EMFFILM (.true., TypeScat, k, ind_ref_f,     &
	                         ind_ref_s, d, z0, Mrank, Nrank, phiGS, beta, alphap,   &
							 snorm, normalized, FileTmat)
    call write_DSCSPARTSUB1 (NthetaGS, normalized, h, v)                    
    close (unit = iDSCS) 
    deallocate (h, v, em)
  end if
  if (ComputeFields) then
    open (unit = iSCAT, file = FileEMF,  status = "replace")
    NphiAL   = Nphi
    NthetaAL = 0
    do iphi = 1, Nphi
      if (Ntheta(iphi) > NthetaAL) NthetaAL = Ntheta(iphi)
    end do
    allocate (emf(3,NphiAL,NthetaAL))
    do i = 1, 3
      do iphi = 1, NphiAL
        do itheta = 1, NthetaAL
          emf(i,iphi,itheta) = zero
        end do
      end do
    end do
    call EMFPARTSUBFILM (1, cc, k, z0, d, ind_ref_f, ind_ref_s, Mrank, Nrank,       &
	     Nmaxmax, Nphi, phi, Ntheta, thetamin, thetamax, emf, NphiAL, NthetaAL)
    call EMFPARTSUBFILM (2, cc, k, z0, d, ind_ref_f, ind_ref_s, Mrank, Nrank,       &
	     Nmaxmax, Nphi, phi, Ntheta, thetamin, thetamax, emf, NphiAL, NthetaAL)
    if (WriteInputInfo) call inputDSCS_EMFFILM (.false., TypeScat, k, ind_ref_f,    &
	                         ind_ref_s, d, z0, Mrank, Nrank, phiGS, beta, alphap,   &
							 snorm, normalized, FileTmat)
    call write_EMF (Nphi, phi, Ntheta, thetamin, thetamax, emf, NphiAL, NthetaAL)       
    close (unit = iSCAT) 
    deallocate (emf)
  end if      
  deallocate (cc)       
end subroutine DiffScatCrossSectPARTSUBFILM 
! **********************************************************************************
subroutine DSCSPARTSUBFILM (TypeField, c, wavenumber, z0, d, ind_refFILM,           &
           ind_refSUB, Mrank, Nrank, Nmax, NthetaGS, phiAZIMUT, snorm, em,          &
		   normalized, h, v)
!------------------------------------------------------------------------------------     
! The routine computes the differential scattering cross sections of an axisymmetric!
! particle on a plane substrate coated with a film, in the azimuth plane phiAZIMUTH.!
! -----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: TypeField, Mrank, Nrank, Nmax, NthetaGS
  real(O)    :: phiAZIMUT, snorm, h(NthetaGS), v(NthetaGS), wavenumber, z0, d      
  complex(O) :: ind_refFILM, ind_refSUB, em(3,NthetaGS), c(2*Nmax)
  logical    :: normalized
!
  integer    :: itheta, m, k, N0, l
  real(O)    :: thetaGS, phiGS, fact
  complex(O) :: sum(3)
  complex(O),allocatable :: Minf(:,:), Ninf(:,:)      
!
  if (normalized) then
    fact = snorm
  else 
    fact = wavenumber * wavenumber
  end if
  fact = 1._O / fact
  allocate (Minf(3,Nmax), Ninf(3,Nmax)) 
  do itheta = 1, NthetaGS  
    thetaGS = Pi / 2._O + real(itheta - 1,O) * Pi / real(NthetaGS - 1,O)
    phiGS   = phiAZIMUT
    if (thetaGS > Pi) then
      thetaGS = 2._O * Pi - thetaGS
      phiGS   = phiAZIMUT + Pi
    end if          
    if (TypeField == 1) then               
      call MN_infinit_complete (thetaGS, phiGS, Mrank, Nrank, Nmax, .true.,         &
           Minf, Ninf)
    else if (TypeField == 2) then   
	  call MN_infinit_reflection_completeFILM (wavenumber, z0, d, ind_refFILM,      &
           ind_refSUB, thetaGS, phiGS, Mrank, Nrank, Nmax, Minf, Ninf)	                      
    end if                         
    sum(1) = zero
    sum(2) = zero
    sum(3) = zero
    do m = 0, Mrank
      if (m == 0) then
        do k = 1, Nrank
          sum(1) = sum(1) + (Minf(1,k) * c(k) + Ninf(1,k) * c(Nmax+k))
          sum(2) = sum(2) + (Minf(2,k) * c(k) + Ninf(2,k) * c(Nmax+k))
          sum(3) = sum(3) + (Minf(3,k) * c(k) + Ninf(3,k) * c(Nmax+k))
        end do
      else          
        N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
        do l = 1, 2
          do k = 1, Nrank - m + 1
            sum(1) = sum(1) + (Minf(1,N0+k) * c(N0+k) + Ninf(1,N0+k) * c(Nmax+N0+k))
            sum(2) = sum(2) + (Minf(2,N0+k) * c(N0+k) + Ninf(2,N0+k) * c(Nmax+N0+k))
            sum(3) = sum(3) + (Minf(3,N0+k) * c(N0+k) + Ninf(3,N0+k) * c(Nmax+N0+k))
          end do           
          N0 = N0 + Nrank - m + 1
        end do
      end if
    end do                          
    em(1,itheta) = em(1,itheta) + sum(1)
    em(2,itheta) = em(2,itheta) + sum(2)
    em(3,itheta) = em(3,itheta) + sum(3)                        
    h(itheta) = abs(em(2,itheta))**2 * fact
    v(itheta) = abs(em(3,itheta))**2 * fact      
  end do
  deallocate (Minf, Ninf)      
end subroutine DSCSPARTSUBFILM
! **********************************************************************************
subroutine EMFPARTSUBFILM (TypeField, c, wavenumber, z0, d, ind_refFILM, ind_refSUB, &
           Mrank, Nrank, Nmax, Nphi, phi, Ntheta, thetamin, thetamax, emf, NphiAL,   &
		   NthetaAL)
!------------------------------------------------------------------------------------     
! The routine computes the electromagnetic scattered field of an axisymmetric       !
! particle on a plane substrate coated with a film.                                 !
! -----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: TypeField, Mrank, Nrank, Nmax, Nphi, Ntheta(NphiMax), NphiAL,       &
                NthetaAL
  real(O)    :: phi(NphiMax), thetamin(NphiMax), thetamax(NphiMax), wavenumber, z0, d      
  complex(O) :: ind_refFILM, ind_refSUB, emf(3,NphiAL,NthetaAL), c(2*Nmax)
!
  integer    :: iphi, itheta, m, k, N0, l
  real(O)    :: thetaGS, phiGS, dtheta
  complex(O) :: sum(3)
  complex(O),allocatable :: Minf(:,:), Ninf(:,:)      
!
  allocate (Minf(3,Nmax), Ninf(3,Nmax)) 
  do iphi = 1, Nphi
    phiGS = phi(iphi)
    if (Ntheta(iphi) /= 1) then
      dtheta = (thetamax(iphi) - thetamin(iphi)) / (Ntheta(iphi) - 1)
    else
      dtheta = 0._O
    end if          
    do itheta = 1, Ntheta(iphi)
      thetaGS = thetamin(iphi) + (itheta - 1) * dtheta              
      if (TypeField == 1) then               
        call MN_infinit_complete (thetaGS, phiGS, Mrank, Nrank, Nmax, .true.,       &
             Minf, Ninf)
      else if (TypeField == 2) then            
        call MN_infinit_reflection_completeFILM (wavenumber, z0, d, ind_refFILM,    & 
           ind_refSUB, thetaGS, phiGS, Mrank, Nrank, Nmax, Minf, Ninf)      
      end if                       
      sum(1) = zero
      sum(2) = zero
      sum(3) = zero
      do m = 0, Mrank
        if (m == 0) then
          do k = 1, Nrank
            sum(1) = sum(1) + (Minf(1,k) * c(k) + Ninf(1,k) * c(Nmax+k))
            sum(2) = sum(2) + (Minf(2,k) * c(k) + Ninf(2,k) * c(Nmax+k))
            sum(3) = sum(3) + (Minf(3,k) * c(k) + Ninf(3,k) * c(Nmax+k))
          end do
        else          
          N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
          do l = 1, 2
            do k = 1, Nrank - m + 1
              sum(1) = sum(1) + (Minf(1,N0+k) * c(N0+k) + Ninf(1,N0+k) *            &
                                   c(Nmax+N0+k))
              sum(2) = sum(2) + (Minf(2,N0+k) * c(N0+k) + Ninf(2,N0+k) *            &
                                   c(Nmax+N0+k))
              sum(3) = sum(3) + (Minf(3,N0+k) * c(N0+k) + Ninf(3,N0+k) *            &
                                   c(Nmax+N0+k))
            end do           
            N0 = N0 + Nrank - m + 1
          end do
        end if
      end do                        
      emf(1,iphi, itheta) = emf(1,iphi, itheta) + sum(1)
      emf(2,iphi, itheta) = emf(2,iphi, itheta) + sum(2)
      emf(3,iphi, itheta) = emf(3,iphi, itheta) + sum(3)                                 
    end do
  end do
  deallocate (Minf, Ninf)      
end subroutine EMFPARTSUBFILM
! **********************************************************************************
subroutine readinputPARTSUBFILM1 ( beta, alphap, ComputeDSCS, ComputeFields,        &
           NthetaGS, phiGS, Nphi, phi, Ntheta, thetamin, thetamax, normalized,      &
		   FileDSCS, FileEMF, WriteInputInfo ) 
  use parameters
  implicit none     
  integer       :: NthetaGS, Nphi, Ntheta(NphiMax), iphi, ios
  real(O)       :: beta, alphap, phiGS, phi(NphiMax), thetamin(NphiMax),            &
                   thetamax(NphiMax), deg
  character(80) :: FileDSCS, FileEMF, string
  logical       :: ComputeDSCS, ComputeFields, normalized, WriteInputInfo, XFindPar  
! -----------------------------------------------------------------------------------
!                        Read the input file FileInputPARTSUBFILM                   !
! -----------------------------------------------------------------------------------             
  open (unit = iInputPARTSUB, file = FileInputPARTSUBFILM, status = "old",          &
        position = "rewind")                
  beta   = 0._O
  alphap = 0._O
  string = 'IncWave'
  if (XFindPar (iInputPARTSUB, string)) then
    read (iInputPARTSUB, *, iostat = ios) beta
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable beta;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) alphap
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable alphap;')"
      stop
    end if         
  else
    print "(/,2x,'Group name IncWave not found;')"
    stop  
  end if 
  call check_incdir_partsub (beta)
  call check_polarization_angle (alphap) 
  deg    = Pi / 180._O
  beta   = beta   * deg
  alphap = alphap * deg
! 
  ComputeDSCS   = .true.
  ComputeFields = .true.
  NthetaGS = 181 
  phiGS    = 0._O  
  Nphi = 1
  do iphi = 1, NphiMax
    phi(iphi)      = 0._O
    Ntheta(iphi)   = 1
    thetamin(iphi) = 0._O
    thetamax(iphi) = 0._O  
  end do 
  normalized = .true.
  FileDSCS   = '../OUTPUTFILES/DSCSpartsub.dat'
  FileEMF    = '../OUTPUTFILES/EMFields.dat'
  WriteInputInfo = .true.
  string     = 'DSCSEM'
  if (XFindPar (iInputPARTSUB, string)) then
    read (iInputPARTSUB, *, iostat = ios) ComputeDSCS 
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ComputeDSCS;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) ComputeFields
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ComputeFields;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) NthetaGS
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable NthetaGS;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) phiGS
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable phiGS;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) Nphi
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nphi;')"
      stop
    end if
    do iphi = 1, Nphi
      read (iInputPARTSUB, *, iostat = ios) phi(iphi)
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable phi;')"
        stop
      end if
      read (iInputPARTSUB, *, iostat = ios) Ntheta(iphi)
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Ntheta;')"
        stop
      end if
      read (iInputPARTSUB, *, iostat = ios) thetamin(iphi)
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable thetamin;')"
        stop
      end if      
      read (iInputPARTSUB, *, iostat = ios) thetamax(iphi)
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable thetamax;')"
        stop
      end if      
    end do    
    read (iInputPARTSUB, *, iostat = ios) normalized
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable normalized;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) FileDSCS
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable FileDSCS;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) FileEMF
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable FileEMF;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) WriteInputInfo
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable WriteInputInfo;')"
      stop
    end if
  else
    print "(/,2x,'Group name DSCSEM not found;')"
    stop  
  end if  
  if (ComputeDSCS) then      
    call check_azimuthal_plane (phiGS)
    phiGS = phiGS * deg                 
  end if
  if (ComputeFields) then
    call check_Nphi (Nphi)
    call check_teta_phiminmax (.false., Nphi, phi, Ntheta, thetamin, thetamax)
    do iphi = 1, Nphi    
      phi(iphi)      = phi(iphi)      * deg
      thetamin(iphi) = thetamin(iphi) * deg
      thetamax(iphi) = thetamax(iphi) * deg      
    end do
  end if
  close (unit = iInputPARTSUB)
end subroutine readinputPARTSUBFILM1
! **********************************************************************************
!                       This file must be stored in InputOutput.f90                *
! **********************************************************************************
subroutine inputDSCS_EMFFILM (ComputeDSCS, TypeScat, k, ind_ref_f, ind_ref_s, d,    &
           z0, Mrank, Nrank, phiGS, beta, alphap, snorm, normalized, FileTmat)                             
  use parameters
  implicit none
  integer       :: TypeScat, Mrank, Nrank, LenString, iOut
  real(O)       :: k, phiGS, beta, alphap, snorm, d, z0, wavelength, anorm, grd
  complex(O)    :: ind_ref_f, ind_ref_s                   
  character(80) :: FileTmat, FileTmatWrite
  logical       :: ComputeDSCS, normalized
!
  wavelength = 2._O * Pi / k
  anorm = sqrt(snorm / Pi) / k
  grd   = 180._O / Pi
  if (ComputeDSCS) then
    iOut = iDSCS
  else
    iOut = iSCAT
  end if
  write (iOut,"(/,2x,'Input Parameters of the Scattering Problem:',/)")
  write (iOut,"(2x,'wavelength of the ambient medium, wavelength = ',1pe13.4,';')") & 
         wavelength
  write (iOut,"(2x, a, 1pe10.3, ',', 1pe10.3, a)")                                  &
 'relative refractive index of the substrate, ind_refSUB  = (', ind_ref_s, ');' 
  write (iOut,"(2x, a, 1pe10.3, ',', 1pe10.3, a)")                                  &
 'relative refractive index of the film,      ind_refFILM = (', ind_ref_f, ');' 
  FileTmatWrite = FileTmat(14:LenString(FileTmat))                 
  write (iOut,"(2x,'name of the file containing the T matrix, FileTmat = ',a40)")   &
         FileTmatWrite
  write (iOut,"(2x,'maximum expansion order, Nrank = ',i3,';')") Nrank
  write (iOut,"(2x,'maximum azimuthal order, Mrank = ',i3,';')") Mrank           
  write (iOut,"(2x,'axial position of the particle, z0 = ',1pe10.3,';')") z0 
  write (iOut,"(2x,'thickness of the film, d = ',1pe10.3,';')") d 
  if (TypeScat == 1)  then
    write (iOut,"(2x, a)")                                                          &
   'illumination by a plane wave traveling in the ambient medium;'
    write (iOut,"(2x,'incident angle, beta = ',f7.2,';')") beta * grd                          
  else if (TypeScat == 2)  then
    write (iOut,"(2x,'illumination by a plane wave traveling in the substrate;')")
    write (iOut,"(2x,'incident angle, beta = ',f7.2,';')") beta * grd       
    write (iOut,"(2x, a, f7.2, a)")                                                 &
   'propagation angle of the incident wave, Pi - beta = ', 180._O - beta * grd, ';'        
  end if                           
  write (iOut,"(2x,'polarization angle of the incident wave, alphap = ',f7.2,';')") &
         alphap * grd
  if (ComputeDSCS) write (iOut,"(2x,'scattering plane, phiGS = ', f7.2,';')")       &
                          phiGS * grd  
  write (iOut,"(2x, a, 1pe10.3, a)")                                                & 
 'characteristic length of the scatterer, anorm = ', anorm, ';'   
  if (normalized) write (iOut,"(2x, a, 1pe13.4, a)")                                &
 'normalization constant, pi * anorm**2 = ', Pi * anorm * anorm, ';'
  write (iOut,*)        
end subroutine inputDSCS_EMFFILM
!***********************************************************************************
!                      This file must be stored in SVWF.f90                        * 
!***********************************************************************************
subroutine MN_infinit_reflection_completeFILM (wavenumber, z0, d, ind_refFILM,      &
           ind_refSUB, theta, phi, Mrank, Nrank, Nmax, Minf, Ninf)
!-----------------------------------------------------------------------------------
! The routine computes the reflected localized vector spherical wave functions     !
! at infinity (excepting the factor exp(j*k*R)/(k*R)) and the azimuthal modes      !
! m = 0,+1,-1,...,+Mrank,-Mrank and for a substrate coated with a film.            !
!                                                                                  !
! Input parameters:                                                                !
! - wavenumber (real) - wave number in the ambient medium.                         !
! - z0 (real) - axial coordinate of the substrate.                                 !
! - d (real) - thickness of the film.                                              !
! - ind_refFILM (complex) - relative refractive index of the film.                 !
! - ind_refSUB (complex) - relative refractive index of the substrate.             !
! - theta, phi (real variables) - polar angles.                                    !
! - Mrank (integer) - number of azimuthal modes.                                   !
! - Nrank (integer) - maximum expansion order.                                     !
! - Nmax (integer) - specifies the dimension of the vectors.                       !
!                                                                                  !
! Output parameters:                                                               !
! - Minf, Ninf (complex arrays) - vector spherical wave functions.                 !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax, k, m, N0, ml, l, n
  real(O)    :: theta, theta0, phi, wavenumber, z0, d, arga, mlr, nm
  complex(O) :: ind_refFILM, ind_refSUB, Rpar, Rperp, cost0, fact, facta, factr,    &
                factc, factt, factp, arg, Minf(3,Nmax), Ninf(3,Nmax)
  real(O),allocatable :: Pnm(:), dPnm(:), pinm(:), taunm(:)       
!
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))
  theta0 = Pi - theta
  cost0  = cmplx(cos(theta0),0.0,O)
  call FresnelCoefFilmReflectionPhase ( cost0, wavenumber, d, z0, ind_refFILM,      &
	   ind_refSUB, Rpar, Rperp) 
  arg   = 2._O * cost0 * wavenumber * z0 
  factr = exp(im * arg)       
  do m = 0, Mrank                  
    call leg_normalized (theta0, m, Nrank, Pnm, dPnm, pinm, taunm)
    if (m == 0) then
      do k = 1, Nrank
        n     = k               
        nm    = real(2 * n * (n + 1),O)         
        nm    = 1._O / sqrt(nm)             
        factc = factr * (-im)**(n + 1) * nm
        factt = factc * taunm(n)                        
        Minf(1,k) =   zero
        Minf(2,k) =   zero
        Minf(3,k) = - factt * Rperp
        Ninf(1,k) =   zero
        Ninf(2,k) =   im * factt * Rpar
        Ninf(3,k) =   zero
      end do    
    else
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      ml = m
      do l = 1, 2
        mlr   = real(ml,O)
        arga  = mlr * phi
        facta = exp(im * arga)
        fact  = facta * factr
        do k = 1, Nrank - m + 1
          n     = m + k - 1
          nm    = real(2 * n * (n + 1),O)               
          nm    = 1._O / sqrt(nm)            
          factc = fact * (-im)**(n + 1) * nm    
          factp = factc * mlr * pinm(n)
          factt = factc * taunm(n)
          Minf(1,N0+k) =   zero
          Minf(2,N0+k) =   im * factp * Rpar
          Minf(3,N0+k) = - factt * Rperp
          Ninf(1,N0+k) =   zero
          Ninf(2,N0+k) =   im * factt * Rpar
          Ninf(3,N0+k) = - factp * Rperp         
        end do
        N0 =   N0 + Nrank - m + 1
        ml = - ml       
      end do
    end if                         
  end do
  deallocate (Pnm, dPnm, pinm, taunm)      
end subroutine MN_infinit_reflection_completeFILM




     
