subroutine TPARTSUB3DFILM
!------------------------------------------------------------------------------------
! 1. General Considerations                                                         !
! --------------------------                                                        !
! TPARTSUB3D is a routine for computing the scattering characteristics of a         !
! homogeneous, dielectric nonaxisymmetric particle on or near a plane surface.      !
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
! The parameters specified in the input file "/INPUTFILES/InputPARTSUB3D.dat" are   !
! listed below.                                                                     !
!                                                                                   !
! - wavelength (real) - wavelength of the incident light in the ambient medium.     !
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
! - FileTmat (character(80)) - name of the file containing the T matrix of          !
!   the particle.                                                                   !
!                                                                                   !
! - axsym (logical) - if axsym = t, the scatterer is a rotationally symmetric       !
!   particle (axisymmetric particle).                                               !
!                                                                                   !
! - Nrank, Mrank (integer variables) - the maximum expansion and azimuthal orders   !
!   of the particle.                                                                !
!                                                                                   !
! - z0 (real) - axial position of the substrate in the particle coordinate system.  !
!   Note that z0 must be a POSITIVE quantity.                                       !         
!                                                                                   !
! - anorm (real) - characteristic length of the particle which is used to           !
!   normalize the differential scattering cross sections.                           !
!                                                                                   !
! - DoConvTest (logical) - if DoConvTest = t, the interactive convergence tests     !
!   over Nint, NXint and Nrank are invoked. An estimate of Nrank is given by        !
!   Wiscombe's truncation limit criterion. If DoConvTest = f, the values of Nint,   !
!   NXint and Nrank must be supplied in the input file.                             !
!                                                                                   !
! - NXint (integer) - number of integration points for computing the elements of    !
!   the reflection matrix. This parameter is used if the convergence tests are      !
!   not performed (DoConvTest = f).                                                 ! 
!                                                                                   !
! - epsNrank (real) - error tolerance for the expansion order test.                 !
!                                                                                   !
! - epsMrank (real) - error tolerance for the azimuthal order test.                 !
!                                                                                   !
!                                                                                   !
! The following parameters control the scattering characteristics calculation.      !
! - beta (real) - incident polar angle.                                             !
!                                                                                   !
! - alpha (real) - incident azimuthal angle.                                        !
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
  integer       :: TypeScat, Nrank, Mrank, NXint                                   
  real(O)       :: wavelength, z0, d, anorm, snorm, epsNrank, epsMrank, k                                       
  complex(O)    :: ind_refFILM, ind_refSUB
  logical       :: axsym, DoConvTest, PrnProgress  
  character(80) :: FileTmat     
! -----------------------------------------------------------------------------------
!                                 Read the input file                               ! 
! -----------------------------------------------------------------------------------       
  call readinputPARTSUB3DFILM ( TypeScat, wavelength, ind_refSUB, ind_refFILM,      &
       Nrank, Mrank, FileTmat, axsym, z0, d, anorm, DoConvTest, NXint, epsNrank,    &
	   epsMrank, PrnProgress, k, snorm)          
! -----------------------------------------------------------------------------------
!                                      Main                                         !
! -----------------------------------------------------------------------------------   
  open (unit = iOutput, file = FileOutput, status = "replace")   
  call  printinputPARTSUB3DFILM (TypeScat, wavelength, anorm, z0, d, epsNrank,      &
        epsMrank, ind_refFILM, ind_refSUB, FileTmat, Mrank, Nrank, axsym)       
  if (DoConvTest) then
    call convergence_Nrank_MrankPARTSUBFILM (TypeScat, k, ind_refFILM, ind_refSUB,  &
           z0, d, snorm, FileTmat, Mrank, Nrank, NXint, axsym, epsNrank, epsMrank,  &
		   PrnProgress) 
  else    
    call TmatPARTSUBFILM (TypeScat, k, ind_refFILM, ind_refSUB, z0, d, snorm,       &
         FileTmat, Mrank, Nrank, NXint, axsym, PrnProgress) 
  end if
  close (unit = iOutput)           
end subroutine TPARTSUB3DFILM
!***********************************************************************************
subroutine readinputPARTSUB3DFILM ( TypeScat, wavelength, ind_refSUB, ind_refFILM,  &
           Nrank, Mrank, FileTmat, axsym, z0, d, anorm, DoConvTest, NXint, epsNrank,&
		   epsMrank, PrnProgress, k, snorm)     
  use parameters
  use derived_parameters
  implicit none 
  integer       :: TypeScat, Mrank, Nrank, NXint, i, ios                                    
  real(O)       :: wavelength, z0, d, anorm, snorm, epsNrank, epsMrank, k,          &
                   Rcirc, xpart                                       
  complex(O)    :: ind_refFILM, ind_refSUB
  character(80) :: string, FileTmat
  logical       :: axsym, DoConvTest, PrnProgress, XFindPar     
! -----------------------------------------------------------------------------------
!                         Read the input file FileInputPARTSUB3DFILM                ! 
! -----------------------------------------------------------------------------------    
  call DrvParameters 
  open (unit = iInputPARTSUB, file = FileInputPARTSUB3DFILM, status = "old",        &
        position = "rewind")   
  wavelength  = 0.1_O * 2._O * Pi  
  ind_refSUB  = (1.3_O,0._O)
  ind_refFILM = (1.3_O,0._O)
  string      = 'OptProp'
  if (XFindPar (iInputPARTSUB, string)) then
    read (iInputPARTSUB, *, iostat = ios) wavelength
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable wavelength;')"
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
  k = 2._O * Pi / wavelength 
!        
  TypeScat = 1  
  FileTmat = '../TMATFILES/T.dat' 
  axsym  = .true.
  Nrank  = 6
  Mrank  = 4
  z0     = 0.1_O
  d      = 0.1_O
  anorm  = 0.1_O   
  string = 'TmatPart'
  if (XFindPar (iInputPARTSUB, string)) then
    read (iInputPARTSUB, *, iostat = ios) TypeScat
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable TypeScat;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) FileTmat
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable FileTmat;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) axsym
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable axsym;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) Nrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nrank;')"
      stop                  
    end if    
    read (iInputPARTSUB, *, iostat = ios) Mrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Mrank;')"
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
  else
    print "(/,2x,'Group name GeomProp not found;')"
    stop  
  end if 
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
  end if   
!   
  NXint  = 60
  string = 'NintX'
  if (XFindPar (iInputPARTSUB, string)) then
    read (iInputPARTSUB, *, iostat = ios) NXint
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable NXint;')"
      stop
    end if        
  else
    print "(/,2x,'Group name NintX not found;')"
    stop  
  end if 
!
  epsNrank = 5.e-2_O
  epsMrank = 5.e-2_O 
  string   = 'Errors'
  if (XFindPar (iInputPARTSUB, string)) then   
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
end subroutine readinputPARTSUB3DFILM
!***********************************************************************************
subroutine printinputPARTSUB3DFILM (TypeScat, wavelength, anorm, z0, d, epsNrank,   &
           epsMrank, ind_refFILM, ind_refSUB, FileTmat, Mrank, Nrank, axsym)             
  use parameters
  implicit none
  integer    :: TypeScat, Mrank, Nrank, LenString
  real(O)    :: wavelength, anorm, epsNrank, epsMrank, z0, d
  complex(O) :: ind_refFILM, ind_refSUB
  logical    :: axsym
  character(80) :: string, FileTmat, FileTmatWrite
!
  write (iOutput,"(/,2x,'Input Parameters of the Scattering Problem:',/)")
  write (iOutput,"(2x, a, 1pe13.4, a)")                                             &
 'wavelength of the ambient medium, wavelength = ', wavelength, ';'
  write (iOutput,"(2x, a, 1pe10.3, ',', 1pe10.3, a)")                               &
 'relative refractive index of the substrate, ind_refSUB = (', ind_refSUB, ');'  
  write (iOutput,"(2x, a, 1pe10.3, ',', 1pe10.3, a)")                               &
 'relative refractive index of the film,     ind_refFILM = (', ind_refFILM, ');' 
  write (iOutput,*)
  FileTmatWrite = FileTmat(14:LenString(FileTmat))
  write (iOutput,"(2x,'name of the file containing the T matrix = ',a)")            &
         FileTmatWrite
  if (axsym) then
    write (iOutput,"(2x,'axisymmetric particle;')")              
  else
    write (iOutput,"(2x,'nonaxisymmetric particle;')")             
  end if
  write (iOutput,"(2x,'maximum expansion order, Nrank = ',i3,';')") Nrank
  write (iOutput,"(2x,'maximum azimuthal order, Mrank = ',i3,';')") Mrank
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'characteristic length of the particle, anorm = ', anorm, ';'  
  write (iOutput,"(2x,'axial position of the particle, z0 = ',1pe10.3,';')") z0
  write (iOutput,"(2x,'thickness of the film, d = ',1pe10.3,';')") d 
  if (TypeScat == 1)  then
    write (iOutput,"(2x, a)")                                                       &
   'illumination by a plane wave traveling in the ambient medium;'          
  else if (TypeScat == 2)  then
    write (iOutput,"(2x,'illumination by a plane wave traveling in the substrate;')")   
  end if  
  write (iOutput,*)    
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum expansion order tolerance, epsNrank = ', epsNrank, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum azimuthal order tolerance, epsMrank = ', epsMrank, ';'            
  write (iOutput,"(/)")                       
end subroutine printinputPARTSUB3DFILM
! **********************************************************************************
subroutine convergence_Nrank_MrankPARTSUBFILM (TypeScat, k, ind_ref_f, ind_ref_s,   &
           z0, d, snorm, FileTmat, Mrank, Nrank, NXint, axsym, epsNrank, epsMrank,  &
		   PrnProgress) 
  use parameters                 
  implicit none 
  integer       :: TypeScat, Mrank, Nrank, NXint
  real(O)       :: snorm, k, epsNrank, epsMrank, z0, d
  complex(O)    :: ind_ref_f, ind_ref_s
  character(80) :: FileTmat   
  logical       :: axsym, PrnProgress 
!
  integer       :: Nmax, ntg, mtg, i, j, Nteta, NthetaConvN, NthetaConvM 
  real(O)       :: beta0, alpha0, alfap, phiGS
  character(80) :: FileTmatPARTSUB
  real(O),allocatable    :: Xint(:), pondereX(:), h(:), v(:), oldh(:), oldv(:),     &
                            oldh0(:), oldv0(:)
  complex(O),allocatable :: a(:,:), b(:,:), t(:,:), r(:,:), c(:), c1(:), em(:,:)
!  
  phiGS  = 0._O
  Nteta  = 11
  beta0  = Pi / 4._O
  alpha0 = 0._O
  alfap  = Pi / 4._O  
  Nmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  FileTmatPARTSUB = '../TMATFILES/TPartSub.dat'  
!
  allocate (Xint(NXint), pondereX(NXint))
  call Laguerre (NXint, Xint, pondereX)    
!  
  if (PrnProgress) call write_progress (.true., 1, 6)
  open (unit = iTmat, file = FileTmat, status = "old", position = "rewind")
  if (.not. axsym) then
    call read_HeadFileTmat (ntg, mtg)
    call check_dimensionMat (ntg, mtg, Nmax)            
    allocate (t(2*ntg,2*mtg))
    call read_FileTmat (ntg, mtg, t) 
  else
    ntg = Nmax
    mtg = Nmax
    allocate (t(2*ntg,2*mtg))
    call read_Tmatrix (.false., Nrank, Mrank, Nmax, t, ntg, mtg)         
  end if
  close (unit = iTmat)
!
  allocate (a(2*Nmax,2*Nmax), b(2*Nmax,2*Nmax),  r(2*Nmax,2*Nmax))
  allocate (c(2*Nmax), c1(2*Nmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta), oldh0(Nteta), oldv0(Nteta),& 
            em(3,Nteta))  
!
  if (PrnProgress) call write_progress (.false., 2, 6)	    
  call copy_matrix (2*Nmax, 2*Nmax, t, 2*ntg, 2*mtg, b, 2*Nmax, 2*Nmax)    ! b = t
  call matrixPARTSUB3DFILM (k, ind_ref_f, ind_ref_s, d, z0, Mrank, Nrank,            &
       Nmax, NXint, Xint, pondereX, r, Nmax, Nmax)                         ! r = refl.
  call copy_matrix (2*Nmax, 2*Nmax, r, 2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax)  ! a = r
  call product_matrixSUB (b, 2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax, 2*Nmax)    ! b = I - b*a   
!
  if (PrnProgress) call write_progress (.false., 3, 6)
  call copy_matrix (2*Nmax, 2*Nmax, t, 2*ntg, 2*mtg, a, 2*Nmax, 2*Nmax)    ! a = t
  call LU_SYSTEM_DIRECT (b, 2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax,                       &
       2*Nmax, 2*Nmax)                                                     ! a = b^-1*a       
!
  open (unit = iTmat, file = FileTmatPARTSUB, status = 'replace')
  call write_HeadFileTmat (Nmax, Nmax) 
  call write_FileTmat (Nmax, Nmax, a) 
  close (unit = iTmat)         
!
  if (PrnProgress) call write_progress (.false., 4, 6)
    if (TypeScat == 1) then
    call PWcoefficientsPARTSUB3D (beta0, alpha0, alfap, Mrank, Nrank, Nmax, c)
    call PWcoefficientsPARTSUBrefl3DFILM (beta0, alpha0, alfap, d, z0,              &
         k, ind_ref_f, ind_ref_s, Mrank, Nrank, Nmax, c1)
    call sum_vectors (c, c1, 2*Nmax)                                       ! c  = c + c1
  else
    call PWcoefficientsPARTSUBtrans3DFILM (beta0, alpha0, alfap, d, z0,              &
         k, ind_ref_f, ind_ref_s, Mrank, Nrank, Nmax, c)    
  end if
  call product_matrix_vector (2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax, c, c1)    ! c1 = a*c  
!
  if (PrnProgress) call write_progress (.false., 5, 6)
  do j = 1, 3
    do i = 1, Nteta
      em(j,i) = zero
    end do
  end do  
  call DSCSPARTSUBFILM (1, c1, k, z0, d, ind_ref_f, ind_ref_s, Mrank, Nrank,         &
	   Nmax, Nteta, phiGS, snorm, em, .true., h, v)
  call DSCSPARTSUBFILM (2, c1, k, z0, d, ind_ref_f, ind_ref_s, Mrank, Nrank,         &
	   Nmax, Nteta, phiGS, snorm, em, .true., h, v)                            
  do i = 1, Nteta
    oldh(i)  = h(i)
    oldv(i)  = v(i)
    oldh0(i) = h(i)
    oldv0(i) = v(i)    
  end do  
  write (iOutput, "(1x, a, i4, a, i3, a, i3, /)")                                     &
 'NXint = ', NXint, ', Nrank = ', Nrank, ', Mrank = ', Mrank 
  call write_DSCSPARTSUB (Nteta, h, v)    
!  
! --- (Nrank - 1) configuration --- 
  if (PrnProgress) call write_progress (.false., 6, 6)
  call copy_matrix (2*Nmax, 2*Nmax, t, 2*ntg, 2*mtg, b, 2*Nmax, 2*Nmax)    ! b = t
  call matrix_Nrank_1_right (Mrank, Nrank, Nmax, b, Nmax, Nmax)            ! b reduce to 0 
  call copy_matrix (2*Nmax, 2*Nmax, r, 2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax)  ! a = r
  call matrix_Nrank_1_right (Mrank, Nrank, Nmax, a, Nmax, Nmax)            ! a reduce to 0
  call product_matrixSUB (b, 2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax, 2*Nmax)    ! b = I - b*a
  call matrix_Nrank_1_left (Mrank, Nrank, Nmax, b, Nmax, Nmax)             ! b reduce to 0 and 1
!    
  call copy_matrix (2*Nmax, 2*Nmax, t, 2*ntg, 2*mtg, a, 2*Nmax, 2*Nmax)    ! a = t
  call matrix_Nrank_1_right (Mrank, Nrank, Nmax, a, Nmax, Nmax)            ! a reduce to 0
  call LU_SYSTEM_DIRECT (b, 2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax,                        &
       2*Nmax, 2*Nmax)                                                     ! a = b^-1*a     
!
  call product_matrix_vector (2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax, c, c1)    ! c1 = a*c  
!
  do j = 1, 3
    do i = 1, Nteta
      em(j,i) = zero
    end do
  end do  
  call DSCSPARTSUBFILM (1, c1, k, z0, d, ind_ref_f, ind_ref_s, Mrank, Nrank,         &
	   Nmax, Nteta, phiGS, snorm, em, .true., h, v)
  call DSCSPARTSUBFILM (2, c1, k, z0, d, ind_ref_f, ind_ref_s, Mrank, Nrank,         &
	   Nmax, Nteta, phiGS, snorm, em, .true., h, v)  
  call delta_DSCSPARTSUB (Mrank, Nteta, h, v, oldh, oldv, epsNrank, NthetaConvN)    
  write (iOutput, "(1x, a, i4, a, i3, a, i3, /)")                                    &
 'NXint = ', NXint, ', Nrank = ', Nrank - 1, ', Mrank = ', Mrank 
  call write_DSCSPARTSUB (Nteta, h, v)
  call write_NrankConvRes (NthetaConvN, Nteta - 2, epsNrank)  
!  
! --- (Mrank - 1) configuration --- 
  call copy_matrix (2*Nmax, 2*Nmax, t, 2*ntg, 2*mtg, b, 2*Nmax, 2*Nmax)    ! b = t
  call matrix_Mrank_1_right (Mrank, Nrank, Nmax, b, Nmax, Nmax)            ! b reduce to 0 
  call copy_matrix (2*Nmax, 2*Nmax, r, 2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax)  ! a = r
  call matrix_Mrank_1_right (Mrank, Nrank, Nmax, a, Nmax, Nmax)            ! a reduce to 0
  call product_matrixSUB (b, 2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax, 2*Nmax)    ! b = I - b*a
  call matrix_Mrank_1_left (Mrank, Nrank, Nmax, b, Nmax, Nmax)             ! b reduce to 0 and 1
!    
  call copy_matrix (2*Nmax, 2*Nmax, t, 2*ntg, 2*mtg, a, 2*Nmax, 2*Nmax)    ! a = t
  call matrix_Mrank_1_right (Mrank, Nrank, Nmax, a, Nmax, Nmax)            ! a reduce to 0
  call LU_SYSTEM_DIRECT (b, 2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax,                        &
       2*Nmax, 2*Nmax)                                                     ! a = b^-1*a     
!
  call product_matrix_vector (2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax, c, c1)    ! c1 = a*c  
!
  do j = 1, 3
    do i = 1, Nteta
      em(j,i) = zero
    end do
  end do 
  call DSCSPARTSUBFILM (1, c1, k, z0, d, ind_ref_f, ind_ref_s, Mrank, Nrank,         &
	   Nmax, Nteta, phiGS, snorm, em, .true., h, v)
  call DSCSPARTSUBFILM (2, c1, k, z0, d, ind_ref_f, ind_ref_s, Mrank, Nrank,         &
	   Nmax, Nteta, phiGS, snorm, em, .true., h, v)  	                   
  call delta_DSCSPARTSUB (Mrank, Nteta, h, v, oldh0, oldv0, epsMrank, NthetaConvM)    
  write (iOutput, "(1x, a, i4, a, i3, a, i3, /)")                                    &
 'NXint = ', NXint, ', Nrank = ', Nrank, ', Mrank = ', Mrank - 1
  call write_DSCSPARTSUB (Nteta, h, v)
  call write_MrankConvRes (NthetaConvM, epsMrank) 
!
!
  if (NthetaConvN >= int(0.8*Nteta) .and. NthetaConvM >= int(0.8*Nteta)) then
    print "(/,2x,'Convergence criteria for Nrank and Mrank are satisfied;')"                                
  else
    print "(/,2x,'Convergence criteria for Nrank and Mrank are not satisfied;')"
  end if
  call DiffScatCrossSectPARTSUB3DFILM (TypeScat, k, ind_ref_f, ind_ref_s, d, z0,     &
           snorm, Nrank, Mrank, FileTmatPARTSUB)     	        
  deallocate (Xint, pondereX, t, a, b, r, c, c1, h, v, oldh, oldv, oldh0, oldv0, em) 
end subroutine convergence_Nrank_MrankPARTSUBFILM
! ***********************************************************************************  
subroutine TmatPARTSUBFILM (TypeScat, k, ind_ref_f, ind_ref_s, z0, d, snorm,         &
           FileTmat, Mrank, Nrank, NXint, axsym, PrnProgress) 
  use parameters                 
  implicit none 
  integer       :: TypeScat, Mrank, Nrank, NXint
  real(O)       :: snorm, k, z0, d
  complex(O)    :: ind_ref_f, ind_ref_s
  character(80) :: FileTmat   
  logical       :: axsym, PrnProgress 
!
  integer       :: Nmax, ntg, mtg   
  character(80) :: FileTmatPARTSUB
  real(O),allocatable    :: Xint(:), pondereX(:)
  complex(O),allocatable :: a(:,:), b(:,:), t(:,:), r(:,:)
!   
  Nmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  FileTmatPARTSUB = '../TMATFILES/TPartSub.dat'  
!
  allocate (Xint(NXint), pondereX(NXint))
  call Laguerre (NXint, Xint, pondereX)    
!  
  if (PrnProgress) call write_progress (.true., 1, 3)
  open (unit = iTmat, file = FileTmat, status = "old", position = "rewind")
  if (.not. axsym) then
    call read_HeadFileTmat (ntg, mtg)
    call check_dimensionMat (ntg, mtg, Nmax)            
    allocate (t(2*ntg,2*mtg))
    call read_FileTmat (ntg, mtg, t) 
  else
    ntg = Nmax
    mtg = Nmax
    allocate (t(2*ntg,2*mtg))
    call read_Tmatrix (.false., Nrank, Mrank, Nmax, t, ntg, mtg)         
  end if
  close (unit = iTmat)
!
  allocate (a(2*Nmax,2*Nmax), b(2*Nmax,2*Nmax),  r(2*Nmax,2*Nmax))  
!	    
  if (PrnProgress) call write_progress (.false., 2, 3)
  call copy_matrix (2*Nmax, 2*Nmax, t, 2*ntg, 2*mtg, b, 2*Nmax, 2*Nmax)    ! b = t
  call matrixPARTSUB3DFILM (k, ind_ref_f, ind_ref_s, d, z0, Mrank, Nrank,            &
       Nmax, NXint, Xint, pondereX, r, Nmax, Nmax)                         ! r = refl.
  call copy_matrix (2*Nmax, 2*Nmax, r, 2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax)  ! a = r
  call product_matrixSUB (b, 2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax, 2*Nmax)    ! b = I - b*a   
!
  if (PrnProgress) call write_progress (.false., 3, 3)
  call copy_matrix (2*Nmax, 2*Nmax, t, 2*ntg, 2*mtg, a, 2*Nmax, 2*Nmax)    ! a = t
  call LU_SYSTEM_DIRECT (b, 2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax,                       &
       2*Nmax, 2*Nmax)                                                     ! a = b^-1*a       
!
  open (unit = iTmat, file = FileTmatPARTSUB, status = 'replace')
  call write_HeadFileTmat (Nmax, Nmax) 
  call write_FileTmat (Nmax, Nmax, a) 
  close (unit = iTmat)   
  call DiffScatCrossSectPARTSUB3DFILM (TypeScat, k, ind_ref_f, ind_ref_s, d, z0,     &
           snorm, Nrank, Mrank, FileTmatPARTSUB) 		           
  deallocate (Xint, pondereX, t, a, b, r) 
end subroutine TmatPARTSUBFILM
!***********************************************************************************
!                     These files must be stored in IncCoef.f90                    *
!***********************************************************************************
subroutine PWcoefficientsPARTSUBrefl3DFILM (beta0, alpha0, alphap, d, z0,           &
           wavenumber, ind_refFILM, ind_refSUB, Mrank, Nrank, Nmax, c)
!-----------------------------------------------------------------------------------
! The routine computes the expansion coefficients of the reflected plane wave      !
! traveling in the direction (Pi - beta0,alpha0).                                  !
!----------------------------------------------------------------------------------- 
  use parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax
  real(O)    :: beta0, alpha0, alphap, d, z0, wavenumber      
  complex(O) :: c(2*Nmax), ind_refFILM, ind_refSUB
!
  integer    :: k, m, n, N0, ml, l      
  real(O)    :: nm, betaM0, mlr, arg
  complex(O) :: Rpar, Rperp, et, ep, cosb0, factc, factp, factt, fact
  real(O),allocatable :: Pnm(:), dPnm(:), pinm(:), taunm(:)
!
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))
  cosb0 = cmplx(cos(beta0),0._O)
  call FresnelCoefFilmReflection ( cosb0, wavenumber, d, z0, ind_refFILM,           &
       ind_refSUB, Rpar, Rperp)
  et     = cos(alphap) * Rpar 
  ep     = sin(alphap) * Rperp
  betaM0 = Pi - beta0 
  do m = 0, Mrank
    call leg_normalized (betaM0, m, Nrank, Pnm, dPnm, pinm, taunm)    
	if (m == 0) then 
      do k = 1, Nrank
        n  = k
        nm = real(2 * n * (n + 1),O)         
        nm    = 1._O / sqrt(nm)     
        factc = 4._O * im**n * nm                 
        factt = factc * taunm(n)
        c(k)      = - factt * ep
        c(Nmax+k) = - im * factt * et
      end do
    else  	                  
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      ml = m        
      do l = 1, 2
        mlr  = real(ml,O)
        arg  = mlr * alpha0
        fact = exp(- im * arg)                  
        do k = 1, Nrank - m + 1
          n     = m + k - 1
          nm    = real(2 * n * (n + 1),O)               
          nm    = 1._O / sqrt(nm)             
          factc = 4._O * im**n * fact * nm
          factp = factc * im * mlr * pinm(n)
          factt = factc * taunm(n)          
          c(N0+k)      = - factp * et - factt * ep
          c(Nmax+N0+k) = - im * (factt * et - factp * ep) 
        end do
        N0 =   N0 + Nrank - m + 1
        ml = - ml
      end do                  
    end if
  end do
  deallocate (Pnm, dPnm, pinm, taunm)      
end subroutine PWcoefficientsPARTSUBrefl3DFILM
! **********************************************************************************
subroutine PWcoefficientsPARTSUBtrans3DFILM (beta, alpha0, alphap, d, z0,           &
           wavenumber, ind_refFILM, ind_refSUB, Mrank, Nrank, Nmax, c)
!-----------------------------------------------------------------------------------
! The routine computes the expansion coefficients of the transmitted plane wave    !
! traveling in the direction (Pi - beta,alpha0).                                   !
!----------------------------------------------------------------------------------- 
  use parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax
  real(O)    :: beta, alpha0, alphap, d, z0, wavenumber      
  complex(O) :: c(2*Nmax), ind_refFILM, ind_refSUB
!
  integer    :: k, m, n, N0, ml, l      
  real(O)    :: nm, mlr, arg
  complex(O) :: Tpar, Tperp, et, ep, cosb, sinb0, cosb0, factc, factp, factt, fact
  complex(O),allocatable :: Pnm(:), dPnm(:), pinm(:), taunm(:)
!
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))
  cosb = cmplx(cos(beta),0._O)
  call FresnelCoefFilmTransmission ( cosb, wavenumber, d, z0, ind_refFILM,          &
       ind_refSUB, Tpar, Tperp, sinb0, cosb0)         
  et    = cos(alphap) * Tpar 
  ep    = sin(alphap) * Tperp 
  cosb0 = - cosb0  
  do m = 0, Mrank
	call Leg_normalized_complex (sinb0, cosb0, m, Nrank, Pnm, dPnm, pinm, taunm) 		   
	if (m == 0) then 
      do k = 1, Nrank
        n  = k
        nm = real(2 * n * (n + 1),O)         
        nm    = 1._O / sqrt(nm)     
        factc = 4._O * im**n * nm                 
        factt = factc * taunm(n)
        c(k)      = - factt * ep
        c(Nmax+k) = - im * factt * et
      end do
    else  	                  
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      ml = m        
      do l = 1, 2
        mlr  = real(ml,O)
        arg  = mlr * alpha0
        fact = exp(- im * arg)                  
        do k = 1, Nrank - m + 1
          n     = m + k - 1
          nm    = real(2 * n * (n + 1),O)               
          nm    = 1._O / sqrt(nm)             
          factc = 4._O * im**n * fact * nm
          factp = factc * im * mlr * pinm(n)
          factt = factc * taunm(n)          
          c(N0+k)      = - factp * et - factt * ep
          c(Nmax+N0+k) = - im * (factt * et - factp * ep) 
        end do
        N0 =   N0 + Nrank - m + 1
        ml = - ml
      end do                  
    end if
  end do  
  deallocate (Pnm, dPnm, pinm, taunm)      
end subroutine PWcoefficientsPARTSUBtrans3DFILM
!***********************************************************************************
!                      This file must be stored in Proces3.f90                     *
!***********************************************************************************
subroutine matrixPARTSUB3DFILM (k, ind_refFILM, ind_refSUB, dist, z0, Mrank, Nrank, &
           Nmax, Nx, x, pond, A, na, ma)
!-----------------------------------------------------------------------------------
! The routine computes the reflection matrix for the scattering of a particle      !
! near or on a plane surface coated with a film.                                   ! 
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax, Nx, na, ma
  real(O)    :: k, dist, z0, x(Nx), pond(Nx)
  complex(O) :: ind_refFILM, ind_refSUB, A(2*na,2*ma)
!
  integer    :: i, j, pint, m, n, n1, N0, NN1
  real(O)    :: D, nm, n1m, mr
  complex(O) :: fact, t, Rpar, Rperp, q, cosbM0, sinbM0, sinb0, cosb0, f, fpp,      &
                ftp, fpt, ftt
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
  do pint = 1, Nx
    t     = one - x(pint) / q
    fact  = 4._O * exp(q) * pond(pint) / q    
    cosb0 = t
    sinb0 = sqrt(one - t * t)     	
	call FresnelCoefFilmReflectionPhase ( cosb0, k, dist, z0, ind_refFILM,          &
	     ind_refSUB, Rpar, Rperp) 	                       	     
    cosbM0 = - t
    sinbM0 = sqrt(one - t * t)
    do m = 0, Mrank     
      mr = real(m,O)
	  call Leg_normalized_complex (sinb0, cosb0, m, Nrank, Pnm, dPnm, pinm, taunm) 
      call Leg_normalized_complex (sinbM0, cosbM0, m, Nrank, PnmM, dPnmM,           &
	                               pinmM, taunmM)                                                 
      if (m == 0) then
        do i = 1, Nrank
          n1  = i
          n1m = real(2 * n1 * (n1 + 1),O)
          n1m = 1._O / sqrt(n1m)
		  do j = 1, Nrank
            n  = j
            nm = real(2 * n * (n + 1),O)
            nm = 1._O / sqrt(nm)
!
            D   = nm * n1m
            f   = fact * D * im**(n1 - n)
            fpp = f * pinm(n)  * pinmM(n1)
            ftt = f * taunm(n) * taunmM(n1)
            fpt = f * pinm(n)  * taunmM(n1)
            ftp = f * taunm(n) * pinmM(n1)        
            A(i,j) = A(i,j) +  ftt * Rperp            
            A(i+Nmax,j+Nmax) = A(i+Nmax,j+Nmax) + ftt * Rpar
          end do
        end do
      else
        N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
		do i = 1, Nrank - m + 1
          n1  = m + i - 1
		  n1m = real(2 * n1 * (n1 + 1),O)
          n1m = 1._O / sqrt(n1m)
          do j = 1, Nrank - m + 1
            n  = m + j - 1
            nm = real(2 * n * (n + 1),O)
            nm = 1._O / sqrt(nm)	
!
            D   = nm * n1m
            f   = fact * D * im**(n1 - n)
            fpp = f * pinm(n)  * pinmM(n1)
            ftt = f * taunm(n) * taunmM(n1)
            fpt = f * pinm(n)  * taunmM(n1)
            ftp = f * taunm(n) * pinmM(n1)        
            A(i+N0,j+N0) = A(i+N0,j+N0) + (mr**2 * fpp * Rpar + ftt * Rperp)        
            A(i+Nmax+N0,j+N0) = A(i+Nmax+N0,j+N0) + mr * (fpt * Rpar + ftp * Rperp)
            A(i+N0,j+Nmax+N0) = A(i+N0,j+Nmax+N0) + mr * (ftp * Rpar + fpt * Rperp)
            A(i+Nmax+N0,j+Nmax+N0) = A(i+Nmax+N0,j+Nmax+N0) +                       &
			                        (mr**2 * fpp * Rperp + ftt * Rpar)							  		 	    		    		    
          end do
		end do    
!
        NN1 = N0 + Nrank - m + 1
		do i = 1, Nrank - m + 1
          do j = 1, Nrank - m + 1			 
            A(i+NN1,j+NN1)           =   A(i+N0,j+N0) 
            A(i+Nmax+NN1,j+NN1)      = - A(i+Nmax+N0,j+N0)
            A(i+NN1,j+Nmax+NN1)      = - A(i+N0,j+Nmax+N0)
            A(i+Nmax+NN1,j+Nmax+NN1) =   A(i+Nmax+N0,j+Nmax+N0)												  
          end do
        end do  			    
      end if  
    end do
  end do         
  deallocate (Pnm, dPnm, pinm, taunm, PnmM, dPnmM, pinmM, taunmM)        
end subroutine matrixPARTSUB3DFILM  
! ***********************************************************************************
!             These files must be stored in PostProces1.f90                         *
! ***********************************************************************************  
subroutine DiffScatCrossSectPARTSUB3DFILM (TypeScat, k, ind_ref_f, ind_ref_s, d, z0, &
           snorm, Nrank, Mrank, FileTmat)                           
  use parameters
  implicit none 
  integer       :: TypeScat, Nrank, Mrank
  real(O)       :: snorm, k, d, z0
  complex(O)    :: ind_ref_f, ind_ref_s
  character(80) :: FileTmat  
!  
  integer       :: NthetaGS, Nphi, Ntheta(NphiMax), Nmax, i, itheta, iphi, ntg, mtg, &
                   NphiAL, NthetaAL
  real(O)       :: beta, alpha, alphap, phiGS, phi(NphiMax), thetamin(Nphimax),      &
                   thetamax(NphiMax)
  character(80) :: FileDSCS, FileEMF
  logical       :: ComputeDSCS, ComputeFields, normalized, WriteInputInfo
  real(O),allocatable    :: h(:), v(:)
  complex(O),allocatable :: t(:,:), c(:), c1(:), em(:,:), emf(:,:,:)
! -----------------------------------------------------------------------------------
!                                 Read the input file                               !
! -----------------------------------------------------------------------------------               
  call readinputPARTSUB3DFILM1 ( beta, alpha, alphap, ComputeDSCS, ComputeFields,   &
       NthetaGS, phiGS, Nphi, phi, Ntheta, thetamin, thetamax, normalized,          &
       FileDSCS, FileEMF, WriteInputInfo )  
  if (.not. ComputeDSCS .and. .not. ComputeFields) return
! -----------------------------------------------------------------------------------
!                                       Main                                        !
! ----------------------------------------------------------------------------------- 
  print "(/,2x,'Scattering Characteristics of a Particle on or Near a Plane Surface')"  
  print "(  2x,'-------------------------------------------------------------------')"  
  Nmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)  
  allocate (c(2*Nmax), c1(2*Nmax))       
  open (unit = iTmat, file = FileTmat, status = "old",  position = "rewind")   
  call read_HeadFileTmat (ntg, mtg)
  call check_dimensionMat (ntg, mtg, Nmax)            
  allocate (t(2*ntg,2*mtg))
  call read_FileTmat (ntg, mtg, t) 
  close (unit = iTmat)
!
  if (TypeScat == 1) then
    call PWcoefficientsPARTSUB3D (beta, alpha, alphap, Mrank, Nrank, Nmax, c)
    call PWcoefficientsPARTSUBrefl3DFILM (beta, alpha, alphap, d, z0,               &
         k, ind_ref_f, ind_ref_s, Mrank, Nrank, Nmax, c1)
    call sum_vectors (c, c1, 2*Nmax)                                       ! c  = c + c1
  else
    call PWcoefficientsPARTSUBtrans3DFILM (beta, alpha, alphap, d, z0,              &
         k, ind_ref_f, ind_ref_s, Mrank, Nrank, Nmax, c)    
  end if
  call product_matrix_vector (2*Nmax, 2*Nmax, t, 2*ntg, 2*mtg, c, c1)      ! c1 = t*c
  deallocate(t, c)                                  
!  
  if (ComputeDSCS) then
    open (unit = iDSCS, file = FileDSCS, status = "replace")
    allocate (h(NthetaGS), v(NthetaGS), em(3,NthetaGS)) 
    do i = 1, 3
      do itheta = 1, NthetaGS
        em(i,itheta) = zero
      end do
    end do 
	call DSCSPARTSUBFILM (1, c1, k, z0, d, ind_ref_f, ind_ref_s, Mrank, Nrank,      &
	     Nmax, NthetaGS, phiGS, snorm, em, normalized, h, v)
    call DSCSPARTSUBFILM (2, c1, k, z0, d, ind_ref_f, ind_ref_s, Mrank, Nrank,      &
	     Nmax, NthetaGS, phiGS, snorm, em, normalized, h, v)                                
    if (WriteInputInfo) call inputDSCS_EMF3DFILM (.true., TypeScat, k, ind_ref_f,   &
	                         ind_ref_s, d, z0, Mrank, Nrank, phiGS, beta, alpha,    &
							 alphap, snorm, normalized, FileTmat)
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
    call EMFPARTSUBFILM (1, c1, k, z0, d, ind_ref_f, ind_ref_s, Mrank, Nrank,       &
	     Nmax, Nphi, phi, Ntheta, thetamin, thetamax, emf, NphiAL, NthetaAL)
    call EMFPARTSUBFILM (2, c1, k, z0, d, ind_ref_f, ind_ref_s, Mrank, Nrank,       &
	     Nmax, Nphi, phi, Ntheta, thetamin, thetamax, emf, NphiAL, NthetaAL)
    if (WriteInputInfo) call inputDSCS_EMF3DFILM (.false., TypeScat, k, ind_ref_f,  &
	                         ind_ref_s, d, z0, Mrank, Nrank, phiGS, beta, alpha,    &
							 alphap, snorm, normalized, FileTmat)					                         
    call write_EMF (Nphi, phi, Ntheta, thetamin, thetamax, emf, NphiAL, NthetaAL)       
    close (unit = iSCAT) 
    deallocate (emf)
  end if      
  deallocate (c1)       
end subroutine DiffScatCrossSectPARTSUB3DFILM 
! **********************************************************************************
subroutine readinputPARTSUB3DFILM1 ( beta, alpha, alphap, ComputeDSCS,              &
           ComputeFields, NthetaGS, phiGS, Nphi, phi, Ntheta, thetamin, thetamax,   &
		   normalized, FileDSCS, FileEMF, WriteInputInfo ) 
  use parameters
  implicit none     
  integer       :: NthetaGS, Nphi, Ntheta(NphiMax), iphi, ios
  real(O)       :: beta, alpha, alphap, phiGS, phi(NphiMax), thetamin(NphiMax),      &
                   thetamax(NphiMax), deg
  character(80) :: FileDSCS, FileEMF, string
  logical       :: ComputeDSCS, ComputeFields, normalized, WriteInputInfo, XFindPar  
! -----------------------------------------------------------------------------------
!                        Read the input file FileInputPARTSUB3DFILM                 !
! -----------------------------------------------------------------------------------             
  open (unit = iInputPARTSUB, file = FileInputPARTSUB3DFILM, status = "old",         &
        position = "rewind")                
  beta   = 0._O
  alpha  = 0._O  
  alphap = 0._O
  string = 'IncWave'
  if (XFindPar (iInputPARTSUB, string)) then
    read (iInputPARTSUB, *, iostat = ios) beta
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable beta;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) alpha
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable alpha;')"
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
  alpha  = alpha  * deg   
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
end subroutine readinputPARTSUB3DFILM1
! **********************************************************************************
!                       This file must be stored in InputOutput.f90                *
! **********************************************************************************
subroutine inputDSCS_EMF3DFILM (ComputeDSCS, TypeScat, k, ind_ref_f, ind_ref_s, d,  &
           z0, Mrank, Nrank, phiGS, beta, alpha, alphap, snorm, normalized, FileTmat)                             
  use parameters
  implicit none
  integer       :: TypeScat, Mrank, Nrank, LenString, iOut
  real(O)       :: k, phiGS, beta, alpha, alphap, snorm, d, z0, wavelength,         &
                   anorm, grd
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
    write (iOut,"(2x,'polar incident angle, beta = ',f7.2,';')") beta * grd                          
  else if (TypeScat == 2)  then
    write (iOut,"(2x,'illumination by a plane wave traveling in the substrate;')")
    write (iOut,"(2x,'polar incident angle, beta = ',f7.2,';')") beta * grd       
    write (iOut,"(2x, a, f7.2, a)")                                                 &
   'polar propagation angle of the incident wave, Pi - beta = ', 180._O - beta * grd, ';'        
  end if 
  write (iOut,"(2x,'azimuthal incident angle, beta = ',f7.2,';')") alpha * grd                          
  write (iOut,"(2x,'polarization angle of the incident wave, alphap = ',f7.2,';')") &
         alphap * grd
  if (ComputeDSCS) write (iOut,"(2x,'scattering plane, phiGS = ', f7.2,';')")       &
                          phiGS * grd  
  write (iOut,"(2x, a, 1pe10.3, a)")                                                & 
 'characteristic length of the scatterer, anorm = ', anorm, ';'   
  if (normalized) write (iOut,"(2x, a, 1pe13.4, a)")                                &
 'normalization constant, pi * anorm**2 = ', Pi * anorm * anorm, ';'
  write (iOut,*)        
end subroutine inputDSCS_EMF3DFILM




  





























