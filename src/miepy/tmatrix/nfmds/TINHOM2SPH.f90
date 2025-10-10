subroutine TINHOM2SPH
!------------------------------------------------------------------------------------
! 1. General Considerations                                                         !
! --------------------------                                                        !
! TINHOM2SPH is a routine for computing the T matrix and the scattering             !
! characteristics of an an inhomogeneous, dielectric sphere with a spherical        !
! inclusion. The program supports calculation for dielectric host spheres with      !
! REAL refractive indices.                                                          !
!                                                                                   !
! For axisymmetric inhomogeneous particles, the scattering problem decouples over   !
! the azimuthal modes m and the T matrix can be computed separately for each m.     !
!                                                                                   !
! 2. Convergence Test                                                               !
! --------------------                                                              !
! Convergence tests precede the scattering characteristics calculation. The         !
! parameters which control the T-matrix computation are the maximum expansion order !
! Nrank, and the maximum azimuthal order Mrank for the host particle.               ! 
!                                                                                   !
! For convergence tests, the incident wave is assumed to be a plane wave traveling  !
! along the Z-axis of the global coordinate system and the scattering               !
! characteristics are computed in the azimuthal plane  phi = 0°.  The convergence   !
! test over Nrank is interactive, while the convergence test over Mrank is          !
! automatically performed.                                                          !
!                                                                                   !
! Nrank is assumed to not depend on Nrank1, which is the maximum expansion order    !
! of the inclusion. The convergence tests are carried out over Nrank and Mrank,     !
! while Nrank1 is kept constant. For the expansion order test, an axisymmetric      !
! orientation of the inhomogeneous particle is considered, i.e., the incident wave  !
! propagates along the axis of symmetry of the particle. In this case, all          !
! azimuthal modes are zero except for m = - 1 and m = 1. The scattering problem is  !
! solved for Nrank and Nrank - 1 and the normalized differential scattering cross   !
! section (DSCS) will be checked at 20° increments for convergence within epsNrank  !
! tolerance. If the calculated results converge within this tolerance at 80% of     !
! the scattering angles, then  convergence is achieved.                             !
!                                                                                   ! 
! After Nrank has been determined we pass to the azimuthal order test. The program  !
! automatically sets the particle to a more general orientation, i.e., alpha =      !
! beta = 45°, and solves the scattering problem  for increasing m values until      !
! convergence of the angular scattering is achieved. The T matrix is stored for     !
! later use by other programs, and the values of Nrank and Mrank are printed to     !
! the screen and to the T-matrix information file (see "Description.txt"). These    !
! values together with the T matrix serve as INPUT PARAMETERS for other programs.   !
!                                                                                   ! 
! The above convergence tests require an estimate of Nrank. This estimate must be   !
! supplied by the user, and for this purpose, Wiscomb's truncation limit criterion  !
! [W. J. Wiscombe, Improved Mie scattering algorithms, Applied Optics, 19,          !
! 1505-1509, 1980] can be used.                                                     !
!                                                                                   !
! The truncation limit criterion proposed by Wiscombe provides the estimate         !
!                                                                                   !
!                NrankW = int(x + 4.05 * x**0.33 + 2),                              !
!                                                                                   ! 
! where x is the size parameter, x = k * r, k is the wave number and r is the       !
! radius of the host sphere.                                                        !
!                                                                                   !
! The convergence test over Nrank can be switched off by setting the logical        !
! variable DoConvTest to false. In this case, the value of Nrank must be specified  !
! in the input file.                                                                !
!                                                                                   !
! 3. Input Parameters                                                               !
! --------------------                                                              !
! The parameters specified in the input file "/INPUTFILES/InputINHOM2SPH.dat"       !
! are listed below.                                                                 !
!                                                                                   !
! - wavelength (real) - wavelength of the incident light in vacuo.                  !
!                                                                                   !
! - ind_refMed (real) - refractive index of the ambient medium (nonabsorbing        !
!   dielectric medium).                                                             !
!                                                                                   !
! - ind_refRel (REAL) - relative refractive index of the host sphere with respect   !
!   to the ambient medium.                                                          !
!                                                                                   !
! - ind_refRel1 (complex) - relative refractive index of the spherical inclusion    !
!   with respect to the host sphere.                                                !
!                                                                                   !
! - r (real) - radius of the host sphere.                                           !
!                                                                                   !
! - anorm (real) - characteristic length of the host sphere which is used to        !
!   normalize the differential scattering cross sections.                           !
!                                                                                   !
! - Nrank1 (integer) - maximum expansion order for the spherical inclusion.         !
!                                                                                   !
! - r1 (real) - radius of the spherical inclusion.                                  !
!                                                                                   !
! - z1 (real) - axial position of the spherical inclusion origin with respect to    !
!   the host sphere.                                                                !
!                                                                                   !
! - DoConvTest (logical) - if DoConvTest = t, the interactive convergence test      !
!   over Nrank is performed. An estimates of Nrank is given by Wiscombe's           !
!   truncation limit criterion. If DoConvTest = f, the value of Nrank must be       !
!   supplied in the input file.                                                     !
!                                                                                   !
! - Nrank (integer) - maximum expansion order for the host sphere. This             !
!   parameter is used if the convergence test is not performed (DoConvTest = f).    !
!                                                                                   !
!   The following table explain the significance of the variable Nrank.             !
!                                                                                   !
!          DoConvTest           Nrank             TypeConvTest                      !
!              t             console input          1 or 2                          !
!              f              file input              2                             !
!                                                                                   !  
!   The significance of the variable TypeConvTest is summarized below.              !
!                                                                                   !
!           TypeConvTest                    Significance                            !  
!               1                      convergence test over Nrank                  ! 
!               2                      convergence test over Mrank                  !
!                                                                                   !
! - epsNrank (real) - error tolerance for the expansion order test.                 !
!                                                                                   !
! - epsMrank (real) - error tolerance for the azimuthal order test.                 !
!                                                                                   !
! - FileTmat (character(80)) - name of the file to which the T matrix of the        !
!   inhomogeneous spherical particle is written.                                    !
!                                                                                   !
! - PrnProgress (logical) - if PrnProgress = t, the progress of calculation         !
!   is printed.                                                                     !
!                                                                                   !
! Note that all input parameters specifying lengths must be provided in the same    !
! units of lengths.                                                                 !
!                                                                                   !
! 4. Logical Scheme                                                                 !
! ------------------                                                                !  
! The organization of the code is as follows:                                       !
! 1. < read input data; the input file provides                                     !
!    - the maximum expansion order Nrank1 for the inclusion,                        ! 
!    - the radius of the spherical inclusion, and                                   !
!    - the axial position of the inclusion with respect to the host particle >      !  
!                                                                                   !
! 2. < T-matrix calculation >                                                       !
!                                                                                   !
!    if ( do convergence test ) then                                                !
!                                                                                   !
!       < the code computes an estimate of Nrank accordingly to                     !
!         Wiscombe's truncation criterion >                                         !
!       < the code prompts for the estimated value of Nrank >                       !
!       < the code prompts for the type of convergence test:                        !
!         1 - over Nrank, 2 - over Mrank >                                          !
!                                                                                   !
!    else if ( .not. do convergence test ) then                                     !
!                                                                                   !
!       < the input file provides the value of Nrank >                              !
!       type of convergence test = 2 (convergence test over Mrank)                  !
!                                                                                   !
!    end if                                                                         !
!    if ( do expansion order test, i.e., type of convergence test = 1 ) then        !
!       < the code computes the DSCS for Nrank and Nrank - 1 and                    !
!         write the results to the file "Output.dat" >                              !
!       < the execution is terminated >                                             !
!    end if                                                                         !
!    if ( do azimuthal order test, i.e., type of convergence test = 2 ) then        !
!       < the code computes the DSCS for increasing m values and                    !
!         write the results to the file "Output.dat" >                              !
!       < the T matrix is stored in the file FileTmat and the T-matrix information  !
!         file is created >                                                         !
!       < Mrank and Nrank are written to the screen and to the T-matrix             !
!         information file >                                                        ! 
!       < go to Step 3 >                                                            !
!    end if                                                                         !
! 3. < scattering characteristics calculation >                                     !    
!------------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none 
  integer       :: Nrank, Nrank1, TypeConvTest                     
  real(O)       :: ks, ind_refMed, ind_refRel, wavelength, r, r1, anorm, snorm,     &
                   epsNrank, epsMrank, z1 
  complex(O)    :: ind_refRel1
  character(80) :: FileTmat
  logical       :: DoConvTest, PrnProgress
! -----------------------------------------------------------------------------------
!                                Read the input file                                ! 
! -----------------------------------------------------------------------------------
  call readinputINHOM2SPH ( wavelength, ind_refMed, ind_refRel, ind_refRel1, r,     &
       anorm, Nrank1, r1, z1, DoConvTest, Nrank, epsNrank, epsMrank, FileTmat,      &
       PrnProgress, TypeConvTest, ks, snorm ) 
! -----------------------------------------------------------------------------------
!                                      Main                                         !
! -----------------------------------------------------------------------------------                   
  open (unit = iOutput, file = FileOutput, status = "replace") 
  call printinputINHOM2SPH (wavelength, anorm, ind_refMed, ind_refRel, ind_refRel1, &
       r, r1, z1, Nrank1, epsNrank, epsMrank)
  if (DoConvTest) then
    if (TypeConvTest == 1) then
      call convergence_NrankINHOM2SPH (ks, ind_refRel, ind_refRel1, r, r1, snorm,   &
           z1, Nrank, Nrank1, epsNrank, PrnProgress )
    else
      call convergence_MrankINHOM2SPH (ks, ind_refRel, ind_refRel1, r, r1, snorm,   &
           z1, Nrank, Nrank1, epsMrank, FileTmat, PrnProgress )    
    end if 
  else
    call convergence_MrankINHOM2SPH (ks, ind_refRel, ind_refRel1, r, r1, snorm, z1, &
         Nrank, Nrank1, epsMrank, FileTmat, PrnProgress )
  end if
  close (unit = iOutput)  
end subroutine TINHOM2SPH
!***********************************************************************************
subroutine readinputINHOM2SPH ( wavelength, ind_refMed, ind_refRel, ind_refRel1,    &
           r, anorm, Nrank1, r1, z1, DoConvTest, Nrank, epsNrank, epsMrank,         &
           FileTmat, PrnProgress, TypeConvTest, ks, snorm) 
  use parameters
  use derived_parameters
  implicit none 
  integer       :: Nrank, Nrank1, NrankW, TypeConvTest, ios                     
  real(O)       :: ks, ind_refMed, ind_refRel, wavelength, r, r1, anorm, xpart,     &
                   snorm, epsNrank, epsMrank, z1, x 
  complex(O)    :: ind_refRel1
  character(80) :: FileTmat, string
  logical       :: DoConvTest, PrnProgress, XFindPar
! -----------------------------------------------------------------------------------
!                       Read the input file FileInputINHOM2SPH                      ! 
! -----------------------------------------------------------------------------------    
  call DrvParameters  
  open (unit = iInputINHOM2SPH, file = FileInputINHOM2SPH, status = "old",          &
        position = "rewind")   
  wavelength  = 0.1_O * 2._O * Pi
  ind_refMed  = 1._O
  ind_refRel  = 1.2_O
  ind_refRel1 = (1.5_O,0._O)
  string     = 'OptProp'
  if (XFindPar (iInputINHOM2SPH, string)) then
    read (iInputINHOM2SPH, *, iostat = ios) wavelength
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable wavelength;')"
      stop
    end if
    read (iInputINHOM2SPH, *, iostat = ios) ind_refMed
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refMed;')"
      stop
    end if
    read (iInputINHOM2SPH, *, iostat = ios) ind_refRel
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refRel;')"
      stop
    end if 
    read (iInputINHOM2SPH, *, iostat = ios) ind_refRel1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refRel1;')"
      stop
    end if
  else
    print "(/,2x,'Group name OptProp not found;')"
    stop  
  end if  
  call check_ind_ref (ind_refRel1) 
  ks = 2._O * Pi * ind_refMed / wavelength
!     
  r = 1._O  
  anorm  = 1._O
  string = 'GeomPropHost'
  if (XFindPar (iInputINHOM2SPH, string)) then
    read (iInputINHOM2SPH, *, iostat = ios) r
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable r;')"
      stop
    end if
    read (iInputINHOM2SPH, *, iostat = ios) anorm
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable anorm;')"
      stop
    end if
  else
    print "(/,2x,'Group name GeomPropHost not found;')"
    stop  
  end if  
  call check_anorm (anorm)    
  xpart = ks * anorm
  snorm = Pi * xpart * xpart        
!
  Nrank1 = 12
  string = 'TmatIncl'
  if (XFindPar (iInputINHOM2SPH, string)) then
    read (iInputINHOM2SPH, *, iostat = ios) Nrank1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nrank1;')"
      stop
    end if
  else
    print "(/,2x,'Group name TmatIncl not found;')"
    stop  
  end if                            
!
  r1 = 0.5_O  
  z1 = 0.3_O  
  string = 'GeomPropIncl'
  if (XFindPar (iInputINHOM2SPH, string)) then
    read (iInputINHOM2SPH, *, iostat = ios) r1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable r1;')"
      stop
    end if
    read (iInputINHOM2SPH, *, iostat = ios) z1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable z1;')"
      stop
    end if
  else
    print "(/,2x,'Group name GeomPropIncl not found;')"
    stop  
  end if         
!
  DoConvTest = .true.
  string     = 'ConvTest'
  if (XFindPar (iInputINHOM2SPH, string)) then
    read (iInputINHOM2SPH, *, iostat = ios) DoConvTest
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
    print "(/,2x,'Convergence Test for an Inhomogeneous Sphere')"
    print "(  2x,'--------------------------------------------')"
  else
    print "(/,2x,'Convergence Test for an Inhomogeneous Sphere over Mrank')"
    print "(  2x,'-------------------------------------------------------')"
  end if      
!
  x = ks * r
  NrankW = int(x + 4.05_O * x**0.33_O + 2._O)
  if (.not. DoConvTest) then        
    Nrank  =  17
    string = 'NrankHost'
    if (XFindPar (iInputINHOM2SPH, string)) then
      read (iInputINHOM2SPH, *, iostat = ios) Nrank
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Nrank;')"
        stop
      end if
    else
      print "(/,2x,'Group name NrankHost not found;')"
      stop  
    end if
    print "(/,2x,'Input values:')"
    print "(  2x,'the input value of Nrank is ', i3,', while the estimated')", Nrank 
    print "(  2x,'value of Nrank from Wiscombe''s criterion is ', i3,';')", NrankW                      
  else    
    print "(/,2x,'Nrank estimate:')"                                                    
    print "(  2x, a, i3, a)",                                                       &
   'the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';'                             
    print "(/,2x,'- enter the estimated value of Nrank;')"
    call read_integer (Nrank)  
  end if
  if (Nrank < Nrank1) print "(/,2x,'Warning: the value of Nrank is too low;')"
!
  if (DoConvTest) then  
    print "(/,2x,'- enter the type of convergence test: 1 - Nrank, 2 - Mrank;')"
    call read_integerbound (TypeConvTest, 1, 2)
  else
    TypeConvTest = 2
  end if 
!
  epsNrank = 5.e-2_O
  epsMrank = 5.e-2_O 
  string   = 'Errors'
  if (XFindPar (iInputINHOM2SPH, string)) then
    read (iInputINHOM2SPH, *, iostat = ios) epsNrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsNrank;')"
      stop
    end if
    read (iInputINHOM2SPH, *, iostat = ios) epsMrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsMrank;')"
      stop
    end if
  else
    print "(/,2x,'Group name Errors not found;')"
    stop  
  end if       
!
  FileTmat = '../TMATFILES/T.dat'
  string   = 'Tmat'
  if (XFindPar (iInputINHOM2SPH, string)) then
    read (iInputINHOM2SPH, *, iostat = ios) FileTmat
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
  string      = 'PrintProgress'
  if (XFindPar (iInputINHOM2SPH, string)) then
    read (iInputINHOM2SPH, *, iostat = ios) PrnProgress
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable PrnProgress;')"
      stop
    end if
  else
    print "(/,2x,'Group name PrintProgress not found;')"
    stop  
  end if    
  close (unit = iInputINHOM2SPH)
end subroutine readinputINHOM2SPH
!***********************************************************************************
subroutine printinputINHOM2SPH (wavelength, anorm, ind_refMed, ind_refRel,          &
           ind_refRel1, r, r1, z1, Nrank1, epsNrank, epsMrank)
  use parameters
  implicit none
  integer       :: Nrank1 
  real(O)       :: wavelength, anorm, ind_refMed, ind_refRel, r, r1, z1, epsNrank,  &
                   epsMrank
  complex(O)    :: ind_refRel1 
!
  write (iOutput,"(/,2x,'Input Parameters of the Scattering Problem:',/)")
  write (iOutput,"(2x,'wavelength of the free space, wavelength = ',1pe13.4,';')")  &
         wavelength
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'refractive index of the ambient medium, ind_refMed = ', ind_refMed, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'relative refractive index of the host particle, ind_refRel = ', ind_refRel, ';'
  write (iOutput,"(2x,'relative refractive index of the inclusion with respect ')")
  write (iOutput,"(2x, a, 1pe10.3,',', 1pe10.3, a)")                                &
 'to the host particle, ind_refRel1 = (', ind_refRel1, ');' 
  write (iOutput,*)
  write (iOutput,"(2x,'radius of the host particle, r = ',1pe10.3,';')") r
  write (iOutput,"(2x,'radius of the inclusion, r1 = ',1pe10.3,';')") r1
  write (iOutput,"(2x, a, i3, a)")                                                  &
 'maximum expansion order for the inclusion: Nrank = ', Nrank1, ';'
  write (iOutput,"(2x,'axial position of the inclusion, z1 = ', 1pe10.3,';')") z1
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'characteristic length of the host particle, anorm = ', anorm, ';'  
  write (iOutput,*) 
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum expansion order tolerance, epsNrank = ', epsNrank, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum azimuthal order tolerance, epsMrank = ', epsMrank, '.'            
  write (iOutput,"(/)")                       
end subroutine printinputINHOM2SPH
!***********************************************************************************
subroutine convergence_NrankINHOM2SPH (ks, ind_ref, ind_ref1, r, r1, snorm, z1,     &
           Nrank, Nrank1, epsNrank, PrnProgress )          
  use parameters
  implicit none
  integer       :: Mrank, Nrank, Nrank1 
  real(O)       :: ks, ind_ref, r, r1, snorm, epsNrank, z1
  complex(O)    :: ind_ref1           
  logical       :: PrnProgress 
!       
  integer       :: Nteta, MM, NN, Mstart, Nmaxmax, m, Nmax, Nmax1, i, NthetaConv,   &
                   NrankAL
  real(O)       :: k, alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,  &
                   Cext, Qext 
  complex(O)    :: ind_refC     
  real(O),allocatable    :: h(:), v(:), oldh(:), oldv(:)
  complex(O),allocatable :: CC(:,:), S1(:,:), S2(:,:), S(:,:)
  complex(O),allocatable :: t(:), t1(:), t2(:), t3(:), t4(:), a(:,:), b(:,:), c(:), &
                            c1(:), ce(:) 
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = 0._O
  beta   = 0._O
  gama   = 0._O
  alfap  = Pi / 4._O
  Mstart = 1
  Mrank  = 1 
  MM = max(Nrank, Nrank1)
  NN = MM + 5
  Nmax  = Nrank
  Nmax1 = Nrank1 
  NrankAL = MM
  k = ks * ind_ref  
  ind_refC = cmplx(ind_ref,0.0,O)
  call write_TypeConvHead (2)
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)             
  allocate (CC(0:NN,0:2*NN+1))
  allocate (S1(2*NrankAL,2*NrankAL), S2(2*NrankAL,2*NrankAL),                       &
            S(2*NrankAL,2*NrankAL))
  allocate (t(2*NrankAL), t1(2*NrankAL), t2(2*NrankAL), t3(2*NrankAL),              &
            t4(2*NrankAL), a(2*NrankAL,2*NrankAL), b(2*NrankAL,2*NrankAL),          &
            c(2*NrankAL), c1(2*NrankAL))
  allocate (ce(2*Nmaxmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do
  if (PrnProgress) call write_progress (.true., 1, 7)    
! --- Nrank configuration ---
  do m = 0, 1                      
    call MatTransCmnn1 (1, k, -z1, m, NN, CC, NN)
    if (m /= 0) then
      call coefficients_fg_m (k, r1, ind_ref1, m, Nrank1, Nmax1, c)
      call MatTransAB_mnn1 (k, .true., -z1, m, CC, NN, Nmax1, Nmax, S2,             &
           NrankAL, NrankAL)
      call matrix_S1 (Nmax, Nmax1, S2, NrankAL, NrankAL, S1, NrankAL, NrankAL)
      call product_matrix_vector2 (2*Nmax1, 2*Nmax, c, S2, 2*NrankAL, 2*NrankAL)
      call product_matrices1 (2*Nmax, 2*Nmax1, 2*Nmax, S1, 2*NrankAL, 2*NrankAL,    &
           S2, 2*NrankAL, 2*NrankAL, S, 2*NrankAL, 2*NrankAL)
      if (PrnProgress) call write_progress (.false., 2, 7)
      call vector_Q_sphere1_m (1, 3, ks, r, ind_refC, m, Nrank, Nmax, t)
      call copy_vector (t, t1, 2*Nmax)          
      call product_matrix_vector1 (2*Nmax, 2*Nmax, t, S, 2*NrankAL, 2*NrankAL,      &
           b, 2*NrankAL, 2*NrankAL)
      call vector_Q_sphere1_m (1, 1, ks, r, ind_refC, m, Nrank, Nmax, t)
      call copy_vector (t, t2, 2*Nmax)  
      call sum_diagonal_elements (2*Nmax, t, b, 2*NrankAL, 2*NrankAL) 
      if (PrnProgress) call write_progress (.false., 3, 7)
      call vector_Q_sphere1_m (3, 3, ks, r, ind_refC, m, Nrank, Nmax, t)
      call copy_vector (t, t3, 2*Nmax)                  
      call product_matrix_vector1 (2*Nmax, 2*Nmax, t, S, 2*NrankAL, 2*NrankAL,      &
           a, 2*NrankAL, 2*NrankAL)
      call vector_Q_sphere1_m (3, 1, ks, r, ind_refC, m, Nrank, Nmax, t)
      call copy_vector (t, t4, 2*Nmax)  
      call sum_diagonal_elements (2*Nmax, t, a, 2*NrankAL, 2*NrankAL)
      if (PrnProgress) call write_progress (.false., 4, 7)
      call LU_SYSTEM (a, 2*NrankAL, 2*NrankAL, b, 2*NrankAL, 2*NrankAL, 2*Nmax)
      if (PrnProgress) call write_progress (.false., 5, 7)
      call minus_matrix (2*Nmax, 2*Nmax, b, 2*NrankAL, 2*NrankAL)
      call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,   &
           Nmax, c)
      call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*NrankAL, 2*NrankAL, c, c1)
      call extend_vector_positive (c1, ce, m, Mstart, Nrank, Nmax, Nmaxmax)
      call matrix_m_negativ (Nmax, Nmax, b, NrankAL, NrankAL)         
      call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,  &
           Nmax, c) 
      call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*NrankAL, 2*NrankAL, c, c1)
      call extend_vector_negative (c1, ce, m, Nrank, Nmax, Nmaxmax)
    end if
  end do                    
  call DSCS (ce, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, ks, snorm,  &
      .false.,.true., h, v)
  call CQscat (ce, Mrank, Nrank, Nmaxmax, ks, snorm, Cscat, Qscat)
  call CQext (ce, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama, alfap,    &
       ks, snorm, Cext, Qext)
  do i = 1, Nteta
    oldh(i) = h(i)
    oldv(i) = v(i)
  end do  
  call write_2ConvParam (Nrank, Mrank)
  call write_DSCS (Nteta,.false., h, v)
  call write_Effic (Qscat, Qext)
! --- (Nrank - 1) configuration ---
  m = 1 
  call matrix_Nrank_m (Nmax, S, NrankAL, NrankAL)
  call copy_vector (t1, t, 2*Nmax)
  call vector_Nrank_m (Nmax, t) 
  call product_matrix_vector1 (2*Nmax, 2*Nmax, t, S, 2*NrankAL, 2*NrankAL, b,       &
       2*NrankAL, 2*NrankAL)
  call copy_vector (t2, t, 2*Nmax)
  call vector_Nrank_m (Nmax, t)
  call sum_diagonal_elements (2*Nmax, t, b, 2*NrankAL, 2*NrankAL) 
  call copy_vector (t3, t, 2*Nmax)      
  call vector_Nrank_m (Nmax, t)
  call product_matrix_vector1 (2*Nmax, 2*Nmax, t, S, 2*NrankAL, 2*NrankAL, a,       &
       2*NrankAL, 2*NrankAL)
  call copy_vector (t4, t, 2*Nmax)      
  call vector_Nrank_m (Nmax, t)
  call sum_diagonal_elements (2*Nmax, t, a, 2*NrankAL, 2*NrankAL) 
  call matrix_Nrank_m_left (Nmax, a, NrankAL, NrankAL)
  if (PrnProgress) call write_progress (.false., 6, 7)
  call LU_SYSTEM (a, 2*NrankAL, 2*NrankAL, b, 2*NrankAL, 2*NrankAL, 2*Nmax)
  if (PrnProgress) call write_progress (.false., 7, 7)
  call minus_matrix (2*Nmax, 2*Nmax, b, 2*NrankAL, 2*NrankAL)
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,       &
       Nmax, c)
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*NrankAL, 2*NrankAL, c, c1)
  call extend_vector_positive (c1, ce, m, Mstart, Nrank, Nmax, Nmaxmax)
  call matrix_m_negativ (Nmax, Nmax, b, NrankAL, NrankAL)             
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,      &
       Nmax, c) 
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*NrankAL, 2*NrankAL, c, c1)
  call extend_vector_negative (c1, ce, m, Nrank, Nmax, Nmaxmax)     
  call DSCS (ce, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, ks, snorm,  &
      .false.,.true., h, v)
  call CQscat (ce, Mrank, Nrank, Nmaxmax, ks, snorm, Cscat, Qscat)
  call CQext (ce, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,           &
       alfap, ks, snorm, Cext, Qext)
  call delta_DSCS (Nteta, h, v, oldh, oldv, epsNrank, NthetaConv)  
  call write_2ConvParam (Nrank - 1, Mrank)
  call write_DSCS (Nteta,.false., h, v) 
  call write_Effic (Qscat, Qext)
  call write_NrankConvRes (NthetaConv, Nteta, epsNrank)
  deallocate (CC, S1, S2, S, t, t1, t2, t3, t4, a, b, c, c1, ce, h, v, oldh, oldv) 
end subroutine convergence_NrankINHOM2SPH
!***********************************************************************************
subroutine convergence_MrankINHOM2SPH (ks, ind_ref, ind_ref1, r, r1, snorm, z1,     &
           Nrank, Nrank1, epsMrank, FileTmat, PrnProgress )          
  use parameters
  implicit none
  integer       :: Mrank, Nrank, Nrank1 
  real(O)       :: ks, ind_ref, r, r1, snorm, epsMrank, z1
  complex(O)    :: ind_ref1           
  character(80) :: FileTmat
  logical       :: PrnProgress 
!       
  integer       :: Nteta, MM, NN, Mstart, Nmaxmax, m, Nmax, Nmax1, i, NthetaConv,   &
                   NrankAL
  real(O)       :: k, alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,  &
                   Cext, Qext      
  complex(O)    :: ind_refC
  real(O),allocatable    :: h(:), v(:), oldh(:), oldv(:)
  complex(O),allocatable :: CC(:,:), S1(:,:), S2(:,:), S(:,:)
  complex(O),allocatable :: t(:), a(:,:), b(:,:), c(:), c1(:), ce(:) 
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = Pi / 4._O
  beta   = Pi / 4._O
  gama   = 0._O
  alfap  = Pi / 4._O
  Mstart = 0
  Mrank  = Nrank 
  MM = max(Nrank, Nrank1)
  NN = MM + 5 
  NrankAL = MM
  k = ks * ind_ref  
  ind_refC = cmplx(ind_ref,0.0,O)  
  open (unit = iTmat, file = FileTmat, status = 'replace')
  call write_HeadFileTmat (NrankAL, NrankAL) 
  call write_TypeConvHead (3)
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)             
  allocate (CC(0:NN,0:2*NN+1))
  allocate (S1(2*NrankAL,2*NrankAL), S2(2*NrankAL,2*NrankAL),                       &
            S(2*NrankAL,2*NrankAL))
  allocate (t(2*NrankAL), a(2*NrankAL,2*NrankAL), b(2*NrankAL,2*NrankAL),           &
            c(2*NrankAL), c1(2*NrankAL))
  allocate (ce(2*Nmaxmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do  
  Mrank = - 1
  do m = Mstart, MM    
    call write_1ConvParam (m)
    Mrank = Mrank + 1        
    if (m > Nrank1) Nrank1 = m
    if (m == 0) then
      Nmax  = Nrank
      Nmax1 = Nrank1
    else
      Nmax  = Nrank - m + 1
      Nmax1 = Nrank1 - m + 1
    end if
    if (PrnProgress) call write_progress_m (.true., m, 1, 5)        
    call MatTransCmnn1 (1, k, -z1, m, NN, CC, NN)
    call coefficients_fg_m (k, r1, ind_ref1, m, Nrank1, Nmax1, c)
    call MatTransAB_mnn1 (k,.true., -z1, m, CC, NN, Nmax1, Nmax, S2,                &
         NrankAL, NrankAL)
    call matrix_S1 (Nmax, Nmax1, S2, NrankAL, NrankAL, S1, NrankAL, NrankAL)
    call product_matrix_vector2 (2*Nmax1, 2*Nmax, c, S2, 2*NrankAL, 2*NrankAL)
    call product_matrices1 (2*Nmax, 2*Nmax1, 2*Nmax, S1, 2*NrankAL, 2*NrankAL,      &
         S2, 2*NrankAL, 2*NrankAL, S, 2*NrankAL, 2*NrankAL)
    if (PrnProgress) call write_progress_m (.false., m, 2, 5)   
    call vector_Q_sphere1_m (1, 3, ks, r, ind_refC, m, Nrank, Nmax, t)              
    call product_matrix_vector1 (2*Nmax, 2*Nmax, t, S, 2*NrankAL, 2*NrankAL,        &
         b, 2*NrankAL, 2*NrankAL)
    call vector_Q_sphere1_m (1, 1, ks, r, ind_refC, m, Nrank, Nmax, t)
    call sum_diagonal_elements (2*Nmax, t, b, 2*NrankAL, 2*NrankAL) 
    if (PrnProgress) call write_progress_m (.false., m, 3, 5)   
    call vector_Q_sphere1_m (3, 3, ks, r, ind_refC, m, Nrank, Nmax, t)          
    call product_matrix_vector1 (2*Nmax, 2*Nmax, t, S, 2*NrankAL, 2*NrankAL,        &
         a, 2*NrankAL, 2*NrankAL)
    call vector_Q_sphere1_m (3, 1, ks, r, ind_refC, m, Nrank, Nmax, t)
    call sum_diagonal_elements (2*Nmax, t, a, 2*NrankAL, 2*NrankAL)
    if (PrnProgress) call write_progress_m (.false., m, 4, 5)
    call LU_SYSTEM (a, 2*NrankAL, 2*NrankAL, b, 2*NrankAL, 2*NrankAL, 2*Nmax)
    if (PrnProgress) call write_progress_m (.false., m, 5, 5)
    call minus_matrix (2*Nmax, 2*Nmax, b, 2*NrankAL, 2*NrankAL)    
    call write_FileTmat (NrankAL, NrankAL, b)       
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,     &
         Nmax, c)
    call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*NrankAL, 2*NrankAL, c, c1)
    call extend_vector_positive (c1, ce, m, Mstart, Nrank, Nmax, Nmaxmax)
    if (m /= 0) then      
      call matrix_m_negativ (Nmax, Nmax, b, NrankAL, NrankAL)         
      call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,  &
           Nmax, c) 
      call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*NrankAL, 2*NrankAL, c, c1)
      call extend_vector_negative (c1, ce, m, Nrank, Nmax, Nmaxmax)     
    end if
    call DSCS (ce, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, ks,       &
         snorm,.false.,.true., h, v)
    call delta_DSCS (Nteta, h, v, oldh, oldv, epsMrank, NthetaConv)     
    call write_DSCS (Nteta,.false., h, v)
    if (NthetaConv >= int(0.8*Nteta)) exit 
  end do
  close (unit = iTmat)
  call CQscat (ce, Mrank, Nrank, Nmaxmax, ks, snorm, Cscat, Qscat)
  call CQext (ce, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,           &
       alfap, ks, snorm, Cext, Qext)
  call write_Effic (Qscat, Qext)
  call write_MrankConvRes (NthetaConv, epsMrank)
  if (NthetaConv >= int(0.8*Nteta)) then
    print "(/,2x,'Convergence criterion for Mrank is satisfied;')"                                     
  else
    print "(/,2x,'Convergence criterion for Mrank is not satisfied;')"
  end if
  call write_InfoFileTmat (FileTmat, Mrank, Nrank, .true., .false., .false.)
  call ScatCharact (ks, FileTmat, Mrank, Nrank, .true., .false., .false.)
  print "(/,2x,'T matrix is stored in ',a50)", FileTmat
  print "(  2x,'The dimensions of the T matrix are given by:')" 
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", Nrank
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank                     
  deallocate (CC, S1, S2, S, t, a, b, c, c1, ce, h, v, oldh, oldv) 
end subroutine convergence_MrankINHOM2SPH
