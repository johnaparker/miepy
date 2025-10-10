Subroutine TMULT2SPH
!------------------------------------------------------------------------------------
! 1. General Considerations                                                         !
! --------------------------                                                        !
! TMULT2SPH is a routine for computing the T matrix and the scattering              !
! characteristics of two homogeneous, dielectric spheres .                          !  
!                                                                                   !
! For axisymmetric systems of particles, the scattering problem decouples over      !
! the azimuthal modes m and the T matrix can be computed separately for each m.     !
!                                                                                   !
! 2. Convergence Test                                                               !
! --------------------                                                              !
! In the code, the maximum expansion and azimuthal orders Nrank and Mrank of the    !
! cluster are independent on the maximum expansion and azimuthal orders             !
! NrankPart(i) and MrankPart(i), i = 1,2, of the spherical particles forming the    !
! cluster. The convergence tests are carried out over Nrank and Mrank, while        !
! NrankPart(i) and MrankPart(i), i = 1,2, are kept constant.                        !
!                                                                                   !
! For convergence tests, the incident wave is assumed to be a plane wave traveling  !
! along the Z-axis of the global coordinate system and the scattering               !
! characteristics are computed in the azimuthal plane  phi = 0°. The convergence    !
! test over Nrank is interactive, while the convergence test over Mrank is          !
! automatically performed.                                                          !
!                                                                                   !
! For the expansion order test, an axisymmetric orientation of the cluster is       !
! considered, i.e., the incident wave propagates along the axis of symmetry of      !
! the cluster. In this case, all azimuthal modes are zero except for m = - 1 and    !
! m = 1. For the convergence test over Nrank, the scattering problem is solved for  !
! Nrank and Nrank - 1. The normalized differential scattering cross section (DSCS)  !
! will be checked at 20° increments for convergence within epsNrank tolerance.      !
! If the calculated results converge within this tolerance at 80% of the            !
! scattering angles, then convergence is achieved.                                  !
!                                                                                   ! 
! After Nrank has been determined we pass to the azimuthal order test. The program  !
! automatically sets the cluster to a more general orientation, i.e., alpha = beta  !
! = 45°, and solves the scattering problem  for increasing m values until           !
! convergence of the angular scattering is achieved. The T matrix is stored for     !
! later use by other programs, and the values of Nrank and Mrank are printed to     !
! the screen and to the T-matrix information file (see "Description.txt"). These    !
! values together with the T matrix serve as INPUT PARAMETERS for other programs.   ! 
!                                                                                   !
! An estimate of Nrank is given by Wiscombe's truncation limit criterion            !
! [W. J. Wiscombe, Improved Mie scattering algorithms, Applied Optics, 19,          !
! 1505-1509, 1980], i.e.,                                                           !
!                                                                                   !
!                NrankW = int(x + 4.05 * x**0.33 + 2),                              !
!                                                                                   !
! where x is the size parameter, x = k * Rcirc, k is the wave number and Rcirc is   !
! the radius of the smallest sphere enclosing the cluster.                          !
!                                                                                   ! 
! The convergence test over Nrank can be switched off by setting the logical        !
! variable DoConvTest to false. In this case, the value of Nrank must be specified  !
! in the input file.                                                                !
!                                                                                   !
! 3. Input Parameters                                                               !
! --------------------                                                              !
! The parameters specified in the input file "/INPUTFILES/InputMULT.dat" are listed !
! below.                                                                            !
!                                                                                   !
! - wavelength (real) - wavelength of the incident light in vacuo.                  !
!                                                                                   !
! - ind_refMed (real) - refractive index of the ambient medium (nonabsorbing        !
!   dielectric medium).                                                             !
!                                                                                   !
! - anorm (real) - characteristic length of the cluster which is used to normalize  !
!   the differential scattering cross sections.                                     !
!                                                                                   !
! - Rcirc (real) - characteristic length of the cluster (usually the radius of the  !
!   smallest circumscribing sphere) which is used to compute an estimate of the     !
!   maximum  expansion order by using Wiscombe's truncation limit criterion (the    !
!   size parameter is x = k * Rcirc, where k is the wave number in the ambient      !
!   medium). This parameter must be specified if the interactive convergence tests  !
!   are performed (DoConvTest = t).                                                 !  
!                                                                                   !
! The next parameters (specified in a group statement) correspond to each            !
! particle. THE GROUP STATEMENT MUST BE REPEATED FOR BOTH PARTICLES.                !
!                                                                                   !
! - rp (real) - radius of the actual spherical particle.                            !
!                                                                                   !
! - Nrankp (integer) - maximum expansion order for the actual spherical particle.   !
!                                                                                   !
! - ind_refp (complex) - relative refractive index of the actual spherical          !
!   particle with respect to the ambient medium.                                    !
!                                                                                   !
! - zp (real) - axial position of the actual spherical particle with respect        !
!   to the coordinate system of the cluster.                                        !
!                                                                                   ! 
! - DoConvTest (logical) - if DoConvTest = t, the interactive convergence test      !
!   over Nrank is performed. An estimates of Nrank is given by Wiscombe's           !
!   truncation limit criterion. If DoConvTest = f, the values of Nrank must be      !
!   supplied in the namelist input.                                                 !
!                                                                                   !
! - Nrank (integer) - maximum expansion order for the cluster. This parameter is    !
!   used if the convergence test is not performed (DoConvTest = f).                 !
!                                                                                   !
!   The following table explain the significance of the variable Nrank.             !
!          DoConvTest           Nrank             TypeConvTest                      !
!              t             console input          1 or 2                          !
!              f              file input              2                             !  
!   where the significance of the variable TypeConvTest is                          !
!           TypeConvTest                    Significance                            !  
!               1                      convergence test over Nrank                  ! 
!               2                      convergence test over Mrank                  !
!                                                                                   !
! - epsNrank (real) - error tolerance for the expansion order test.                 !
!                                                                                   !
! - epsMrank (real) - error tolerance for the azimuthal order test.                 !
!                                                                                   !
! - FileTmat (character(80)) - name of the file to which the T matrix of the        !
!   cluster is written.                                                             !
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
! 1. < read input data; the input file provides:                                    !
!    - the maximum expansion order Nrankp,                                          ! 
!    - the radius of each spherical particle, and                                   !                               
!    - the axial position of each spherical particle >                              !
!                                                                                   !   
! 2. < T-matrix calculation >                                                       !
!                                                                                   ! 
!    if ( do convergence test ) then                                                !
!                                                                                   !
!      < the code computes an estimate of Nrank accordingly to                      !
!        Wiscombe's truncation criterion >                                          !
!      < the code prompts for the estimated value of Nrank >                        !
!      < the code prompts for the type of convergence test:                         !
!        1 - over Nrank, 2 - over Mrank >                                           !
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
  implicit none 
  integer       :: Nrank1, Nrank2, Nrank, TypeConvTest                                   
  real(O)       :: k, ind_refMed, wavelength, anorm, snorm, epsNrank, epsMrank, r1, &
                   r2, z1, z2, Rcirc
  complex(O)    :: ind_ref1, ind_ref2
  logical       :: DoConvTest, PrnProgress
  character(80) :: FileTmat
! -----------------------------------------------------------------------------------
!                               Read the input file                                 ! 
! -----------------------------------------------------------------------------------
  call readinputMULT2SPH (wavelength, ind_refMed, anorm, Rcirc, r1, r2,             &
       z1, z2, Nrank1, Nrank2, ind_ref1, ind_ref2, DoConvTest, Nrank, epsNrank,     &
       epsMrank, FileTmat, TypeConvTest, PrnProgress, k, snorm)  
! -----------------------------------------------------------------------------------
!                                      Main                                         !
! -----------------------------------------------------------------------------------                   
  open (unit = iOutput,file = FileOutput,status = "replace") 
  call printinputMULT2SPH (wavelength, ind_refMed, r1, r2, z1, z2, ind_ref1,        &
       ind_ref2, Nrank1, Nrank2, epsNrank, epsMrank, anorm, Rcirc)
  if (DoConvTest) then
    if (TypeConvTest == 1) then
      call convergence_NrankMULT2SPH (k, snorm, r1, r2, ind_ref1, ind_ref2, z1, z2, &
           Nrank1, Nrank2, Nrank, epsNrank, PrnProgress )
    else 
      call convergence_MrankMULT2SPH (k, snorm, r1, r2, ind_ref1, ind_ref2, z1, z2, &
           Nrank1, Nrank2, Nrank, epsMrank, FileTmat, PrnProgress )  
    end if   
  else
    call convergence_MrankMULT2SPH (k, snorm, r1, r2, ind_ref1, ind_ref2, z1, z2,   &
         Nrank1, Nrank2, Nrank, epsMrank, FileTmat, PrnProgress )
  end if
  close (unit = iOutput)
end subroutine TMULT2SPH
!***********************************************************************************
subroutine readinputMULT2SPH (wavelength, ind_refMed, anorm, Rcirc, r1, r2,         &
           z1, z2, Nrank1, Nrank2, ind_ref1, ind_ref2, DoConvTest, Nrank, epsNrank, &
           epsMrank, FileTmat, TypeConvTest, PrnProgress, k, snorm)   
  use parameters
  use derived_parameters
  implicit none 
  integer       :: Nrankp, Nrank1, Nrank2, Nrank, TypeConvTest, ios, NrankW                                    
  real(O)       :: k, ind_refMed, wavelength, anorm, x, snorm, epsNrank, epsMrank,  &
                   rp, zp, r1, r2, z1, z2, Rcirc, xR
  complex(O)    :: ind_refp, ind_ref1, ind_ref2
  logical       :: DoConvTest, PrnProgress, XFindPar
  character(80) :: FileTmat, string
! -----------------------------------------------------------------------------------
!                        Read the input file FileInputMULT2SPH                      ! 
! -----------------------------------------------------------------------------------    
  call DrvParameters  
  open (unit = iInputMULT2SPH, file = FileInputMULT2SPH, status = "old",            &
        position = "rewind")   
  wavelength = 0.1_O * 2._O * Pi                                                                    
  ind_refMed = 1._O
  string     = 'OptProp'
  if (XFindPar (iInputMULT2SPH, string)) then
    read (iInputMULT2SPH, *, iostat = ios) wavelength
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable wavelength;')"
      stop
    end if
    read (iInputMULT2SPH, *, iostat = ios) ind_refMed
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refMed;')"
      stop
    end if         
  else
    print "(/,2x,'Group name OptProp not found;')"
    stop  
  end if  
  k = 2._O * Pi * ind_refMed / wavelength         
!
  anorm  = 1._O     
  Rcirc  = 1._O
  string = 'GenProp'
  if (XFindPar (iInputMULT2SPH, string)) then
    read (iInputMULT2SPH, *, iostat = ios) anorm
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable anorm;')"
      stop
    end if
    read (iInputMULT2SPH, *, iostat = ios) Rcirc
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Rcirc;')"
      stop
    end if         
  else
    print "(/,2x,'Group name GenProp not found;')"
    stop  
  end if    
  call check_anorm (anorm)     
  x = k * anorm
  snorm = Pi * x * x  
!
  rp = 0.2_O
  zp = 0._O             
  Nrankp   = 9
  ind_refp = (1.5_O,0._O) 
  string   = 'PartProp'
  if (XFindPar (iInputMULT2SPH, string)) then
    read (iInputMULT2SPH, *, iostat = ios) rp
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable rp;')"
      print "(  2x,'for the first particle;')"
      stop
    end if
    read (iInputMULT2SPH, *, iostat = ios) zp
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable zp;')"
      print "(  2x,'for the first particle;')"
      stop
    end if
    read (iInputMULT2SPH, *, iostat = ios) Nrankp
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nrankp;')"
      print "(  2x,'for the first particle;')"
      stop
    end if
    read (iInputMULT2SPH, *, iostat = ios) ind_refp
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refp;')"
      print "(  2x,'for the first particle;')"
      stop
    end if
  else
    print "(/,2x,'Group name PartProp not found;')"
    stop  
  end if        
  r1 = rp
  z1 = zp
  Nrank1   = Nrankp
  ind_ref1 = ind_refp
  call check_ind_ref1 (1, ind_ref1)
!
  rp = 0.2_O
  zp = 0._O             
  Nrankp   = 9
  ind_refp = (1.5_O,0._O) 
  string   = 'PartProp'
  if (XFindPar (iInputMULT2SPH, string)) then
    read (iInputMULT2SPH, *, iostat = ios) rp
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable rp;')"
      print "(  2x,'for the second particle;')"
      stop
    end if
    read (iInputMULT2SPH, *, iostat = ios) zp
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable zp;')"
      print "(  2x,'for the second particle;')"
      stop
    end if
    read (iInputMULT2SPH, *, iostat = ios) Nrankp
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nrankp;')"
      print "(  2x,'for the second particle;')"
      stop
    end if
    read (iInputMULT2SPH, *, iostat = ios) ind_refp
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refp;')"
      print "(  2x,'for the second particle;')"
      stop
    end if
  else
    print "(/,2x,'Group name PartProp not found;')"
    stop  
  end if     
  r2 = rp
  z2 = zp
  Nrank2   = Nrankp
  ind_ref2 = ind_refp  
  call check_ind_ref1 (2, ind_ref2)
!
  DoConvTest = .true.
  string     = 'ConvTest'
  if (XFindPar (iInputMULT2SPH, string)) then
    read (iInputMULT2SPH, *, iostat = ios) DoConvTest
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
    print "(/,2x,'Convergence Test for an Ensemble of Two Spheres')"
    print "(  2x,'-----------------------------------------------')"
  else
    print "(/,2x,'Convergence Test for an Ensemble of Two Spheres over Mrank')"
    print "(  2x,'----------------------------------------------------------')"
  end if   
!
  xR = k * Rcirc 
  NrankW = int(xR + 4.05_O * xR**0.33_O + 2._O)
  if (.not. DoConvTest) then        
    Nrank  =  17
    string = 'NrankCluster'
    if (XFindPar (iInputMULT2SPH, string)) then
      read (iInputMULT2SPH, *, iostat = ios) Nrank
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Nrank;')"
        stop
      end if          
    else
      print "(/,2x,'Group name NrankCluster not found;')"
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
  if (Nrank < Nrank1 .or. Nrank < Nrank2)                                           &
      print "(/,2x,'Warning: the value of Nrank is too low;')"
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
  if (XFindPar (iInputMULT2SPH, string)) then
     read (iInputMULT2SPH, *, iostat = ios) epsNrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsNrank;')"
      stop
    end if 
    read (iInputMULT2SPH, *, iostat = ios) epsMrank
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
  if (XFindPar (iInputMULT2SPH, string)) then
    read (iInputMULT2SPH, *, iostat = ios) FileTmat
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
  if (XFindPar (iInputMULT2SPH, string)) then
    read (iInputMULT2SPH, *, iostat = ios) PrnProgress
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable PrnProgress;')"
      stop
    end if          
  else
    print "(/,2x,'Group name PrintProgress not found;')"
    stop  
  end if    
  close (unit = iInputMULT2SPH)
end subroutine readinputMULT2SPH
!***********************************************************************************
subroutine printinputMULT2SPH (wavelength, ind_refMed, r1, r2, z1, z2, ind_ref1,    &
           ind_ref2, Nrank1, Nrank2, epsNrank, epsMrank, anorm, Rcirc)
  use parameters
  implicit none  
  integer       :: Nrank1, Nrank2, Nrank
  real(O)       :: wavelength, ind_refMed, r1, r2, z1, z2, epsNrank, epsMrank, x,   &
                   anorm, Rcirc, k
  complex(O)    :: ind_ref1, ind_ref2  
!
  write (iOutput,"(/,2x,'Input Parameters of the Scattering Problem:',/)")
  write (iOutput,"(2x,'wavelength of the free space, wavelength = ',1pe13.4,';')")  &
         wavelength
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'refractive index of the ambient medium, ind_refMed = ', ind_refMed, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'characteristic length of the cluster, anorm = ', anorm, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'radius of the circumscribing sphere, Rcirc = ', Rcirc, ';'
  write (iOutput,*)
  write (iOutput,"(2x,'particle 1: ')")
  write (iOutput,"(2x,'radius: r = ',1pe10.3,';')") r1
  write (iOutput,"(2x,'axial position: z = ',1pe10.3,';')") z1
  write (iOutput,"(2x, a, 1pe10.3, ',', 1pe10.3, a)")                               &
 'relative refractive index of the particle, ind_refRel = (', ind_ref1, ');' 
  write (iOutput,"(2x,'maximum expansion order: Nrank = ',i3,';')") Nrank1
  k = 2._O * Pi * ind_refMed / wavelength  
  x = k * r1
  Nrank = int(x + 4.05_O * x**0.33_O + 2._O)
  write (iOutput,"(2x, a, i3, a)")                                                  &
 'the estimated value of Nrank from Wiscombe criterion is ', Nrank, ';' 
  write (iOutput,*)
  write (iOutput,"(2x,'particle 2: ')")
  write (iOutput,"(2x,'radius: r = ',1pe10.3,';')") r2
  write (iOutput,"(2x,'axial position: z = ',1pe10.3,';')") z2
  write (iOutput,"(2x, a, 1pe10.3, ',', 1pe10.3, a)")                               &
 'relative refractive index of the particle, ind_refRel = (', ind_ref2, ');' 
  write (iOutput,"(2x,'maximum expansion order: Nrank = ',i3,';')") Nrank2
  x = k * r2
  Nrank = int(x + 4.05_O * x**0.33_O + 2._O)
  write (iOutput,"(2x, a, i3, a)")                                                  &
 'the estimated value of Nrank from Wiscombe criterion is ', Nrank, ';' 
  write (iOutput,*)
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum expansion order tolerance, epsNrank = ', epsNrank, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum azimuthal order tolerance, epsMrank = ', epsMrank, '.'            
  write (iOutput,"(/)")                       
end subroutine printinputMULT2SPH
!***********************************************************************************
subroutine convergence_NrankMULT2SPH (k, snorm, r1, r2, ind_ref1, ind_ref2, z1, z2, &
           Nrank1, Nrank2, Nrank, epsNrank, PrnProgress )
  use parameters
  implicit none
  integer       :: Nrank1, Nrank2, Nrank  
  real(O)       :: k, snorm, r1, r2, z1, z2, epsNrank
  complex(O)    :: ind_ref1, ind_ref2
  logical       :: PrnProgress 
!       
  integer       :: m, Mstart, Mrank, MM, NN, Nmaxim, Nrankrank, Nmax, Nmax1, Nmax2, &
                   Nmaxmax, Nteta, i, NthetaConv, NrankAL  
  real(O)       :: alfap, tetaGI, phiGI, phiGS, alfa, beta, gama, Cscat, Qscat,     &
                   Cext, Qext, z1z2
  logical       :: type1, type2                             
  real(O),allocatable    :: h(:), v(:), oldh(:), oldv(:)
  complex(O),allocatable :: CC1(:,:), CC2(:,:), CC3(:,:), aa(:,:), bb(:,:), a(:,:), &
                            b(:,:), c(:,:), t1(:), t2(:), c1(:), c2(:), ce(:)
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
  MM = Nrank
  NN = Nrank + 5
  Nmax   = Nrank
  Nmax1  = Nrank1
  Nmax2  = Nrank2
  Nmaxim = Nrank  + Mrank * (2 * Nrank - Mrank + 1) 
  Nrankrank = Nrank1 + Nrank2 
  Nmaxmax = Nmax1  + Nmax2 
  NrankAL = max(Nrank, Nrank1, Nrank2)
  call write_TypeConvHead (2)
  allocate (CC1(0:NN,0:2*NN+1),  CC2(0:NN,0:2*NN+1),  CC3(0:NN,0:2*NN+1))   
  allocate (aa(2*Nrankrank,2*Nrankrank), bb(2*Nrankrank,2*NrankAL))        
  allocate (a(2*NrankAL,2*NrankAL), b(2*NrankAL,2*NrankAL), c(2*NrankAL,2*NrankAL))
  allocate (t1(2*Nrank1), t2(2*Nrank2))
  allocate (c1(2*NrankAL), c2(2*NrankAL), ce(2*Nmaxim))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do 
  z1z2 = abs(z1 - z2)
  if (z1 > z2) then
    type1 = .false.
    type2 = .true.
  else
    type1 = .true.
    type2 = .false.
  end if
  if (PrnProgress) call write_progress (.true., 1, 8)
! --- Nrank configuration ---  
  do m = 0, 1    
    call MatTransCmnn1 (3, k, z1z2, m, NN, CC3,  NN)
    call MatTransCmnn1 (1, k, z1,   m, NN, CC1,  NN)          
    call MatTransCmnn1 (1, k, z2,   m, NN, CC2,  NN)      
    if (m /= 0) then
      call coefficients_fg_m (k, r1, ind_ref1, m, Nrank1, Nmax1, t1)    
      call coefficients_fg_m (k, r2, ind_ref2, m, Nrank2, Nmax2, t2)        
      call identity_matrix (2*Nmaxmax, aa, 2*Nrankrank, 2*Nrankrank)    
      call MatTransAB_mnn1 (k,type1, z1z2, m, CC3, NN, Nmax1, Nmax2, a,             &
           NrankAL, NrankAL)        
      call product_matrix_vector1 (2*Nmax1, 2*Nmax2, t1, a, 2*NrankAL, 2*NrankAL,   &
           b, 2*NrankAL, 2*NrankAL) 
      call extend_matrix3 (Nmax1, Nmax2, 0, 2*Nmax1, b, NrankAL, NrankAL, aa,       &
           Nrankrank, Nrankrank)
      if (PrnProgress) call write_progress (.false., 2, 8)                                    
      call MatTransAB_mnn1 (k,type2, z1z2, m, CC3, NN, Nmax2, Nmax1, a,             &
           NrankAL, NrankAL)          
      call product_matrix_vector1 (2*Nmax2, 2*Nmax1, t2, a, 2*NrankAL, 2*NrankAL,   &
           b, 2*NrankAL, 2*NrankAL)
      call extend_matrix3 (Nmax2, Nmax1, 2*Nmax1, 0, b, NrankAL, NrankAL, aa,       &
           Nrankrank, Nrankrank)
      if (PrnProgress) call write_progress (.false., 3, 8)                      
      call MatTransAB_mnn1 (k,.false., z1, m, CC1, NN, Nmax1, Nmax, a,              &
           NrankAL, NrankAL)           
      call product_matrix_vector1 (2*Nmax1, 2*Nmax, t1, a, 2*NrankAL, 2*NrankAL,    &
           b, 2*NrankAL, 2*NrankAL)
      call extend_matrix4 (1, Nmax1, Nmax, Nmaxmax, 0, b, NrankAL, NrankAL, bb,     &
           Nrankrank, NrankAL) 
      if (PrnProgress) call write_progress (.false., 4, 8)
      call MatTransAB_mnn1 (k,.false., z2, m, CC2, NN, Nmax2, Nmax, a,              &
           NrankAL, NrankAL)      
      call product_matrix_vector1 (2*Nmax2, 2*Nmax, t2, a, 2*NrankAL, 2*NrankAL,    &
           b, 2*NrankAL, 2*NrankAL)
      call extend_matrix4 (2, Nmax2, Nmax, Nmaxmax, 2*Nmax1, b, NrankAL, NrankAL,   &
           bb, Nrankrank, NrankAL)  
      if (PrnProgress) call write_progress (.false., 5, 8)
      call LU_SYSTEM_DIRECT (aa, 2*Nrankrank, 2*Nrankrank, bb, 2*Nrankrank,         &
           2*NrankAL, 2*Nmaxmax, 2*Nmax)
      if (PrnProgress) call write_progress (.false., 6, 8)                   
      call extract_matrix2 (Nmax1, Nmax, 0, b, NrankAL, NrankAL, bb, Nrankrank,     &
           NrankAL)
      call MatTransAB_mnn1 (k,.true., z1, m, CC1, NN, Nmax, Nmax1, a,               &
           NrankAL, NrankAL)
      call product_matrices2 (1, 2*Nmax, 2*Nmax1, 2*Nmax, a, 2*NrankAL, 2*NrankAL,  &
           b, 2*NrankAL, 2*NrankAL, c, 2*NrankAL, 2*NrankAL)    
      if (PrnProgress) call write_progress (.false., 7, 8)
      call extract_matrix2 (Nmax2, Nmax, 2*Nmax1, b, NrankAL, NrankAL,              &
           bb, Nrankrank, NrankAL)
      call MatTransAB_mnn1 (k,.true., z2, m, CC2, NN, Nmax, Nmax2, a,               &
           NrankAL, NrankAL)
      call product_matrices2 (2, 2*Nmax, 2*Nmax2, 2*Nmax, a, 2*NrankAL, 2*NrankAL,  &
           b, 2*NrankAL, 2*NrankAL, c, 2*NrankAL, 2*NrankAL)                
      if (PrnProgress) call write_progress (.false., 8, 8)
      call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,   &
           Nmax, c1)
      call product_matrix_vector (2*Nmax, 2*Nmax, c, 2*NrankAL, 2*NrankAL, c1, c2)
      call extend_vector_positive (c2, ce, m, Mstart, Nrank, Nmax, Nmaxim)        
      call matrix_m_negativ (Nmax, Nmax, c, NrankAL, NrankAL)         
      call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,  &
           Nmax, c1)        
      call product_matrix_vector (2*Nmax, 2*Nmax, c, 2*NrankAL, 2*NrankAL, c1, c2)
      call extend_vector_negative (c2, ce, m, Nrank, Nmax, Nmaxim)      
    end if
  end do 
  call DSCS (ce, Mrank, Nrank, Nmaxim, Nteta, phiGS, alfa, beta, gama, k, snorm,    &
      .false.,.true., h, v)    
  call CQscat (ce, Mrank, Nrank, Nmaxim, k, snorm, Cscat, Qscat)
  call CQext (ce, Mrank, Nrank, Nmaxim, tetaGI, phiGI, alfa, beta, gama,            &
       alfap, k, snorm, Cext, Qext)
  do i = 1, Nteta
    oldh(i) = h(i)
    oldv(i) = v(i)
  end do  
  call write_2ConvParam (Nrank, Mrank)
  call write_DSCS (Nteta,.false., h, v)
  call write_Effic (Qscat, Qext)   
! --- (Nrank - 1) configuration ---
  m = 1
  if (PrnProgress) call write_progress_low
  call matrix_m_negativ (Nmax, Nmax, c, NrankAL, NrankAL)
  call matrix_Nrank_m (Nmax, c, NrankAL, NrankAL)
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,       &
       Nmax, c1)
  call product_matrix_vector (2*Nmax, 2*Nmax, c, 2*NrankAL, 2*NrankAL, c1, c2)
  call extend_vector_positive (c2, ce, m, Mstart, Nrank, Nmax, Nmaxim)
  call matrix_m_negativ (Nmax, Nmax, c, NrankAL, NrankAL)             
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap,-m, Nrank,       &
       Nmax, c1)        
  call product_matrix_vector (2*Nmax, 2*Nmax, c, 2*NrankAL, 2*NrankAL, c1, c2)
  call extend_vector_negative (c2, ce, m, Nrank, Nmax, Nmaxim)
  call DSCS (ce, Mrank, Nrank, Nmaxim, Nteta, phiGS, alfa, beta, gama, k, snorm,    &
      .false.,.true., h, v)
  call CQscat (ce, Mrank, Nrank, Nmaxim, k, snorm, Cscat, Qscat)
  call CQext (ce, Mrank, Nrank, Nmaxim, tetaGI, phiGI, alfa, beta, gama,            &
       alfap, k, snorm, Cext, Qext)
  call delta_DSCS (Nteta, h, v, oldh, oldv, epsNrank, NthetaConv)  
  call write_2ConvParam (Nrank - 1, Mrank)
  call write_DSCS (Nteta,.false., h, v) 
  call write_Effic (Qscat, Qext)
  call write_NrankConvRes (NthetaConv, Nteta, epsNrank)
  deallocate (CC1, CC2, CC3, aa, bb, a, b, c, c1, c2, ce, t1, t2, h, v, oldh, oldv) 
end subroutine convergence_NrankMULT2SPH
!***********************************************************************************
subroutine convergence_MrankMULT2SPH (k, snorm, r1, r2, ind_ref1, ind_ref2, z1, z2, &
           Nrank1, Nrank2, Nrank, epsMrank, FileTmat, PrnProgress )
  use parameters
  implicit none
  integer       :: Nrank1, Nrank2, Nrank  
  real(O)       :: k, snorm, r1, r2, z1, z2, epsMrank
  complex(O)    :: ind_ref1, ind_ref2
  character(80) :: FileTmat
  logical       :: PrnProgress 
!       
  integer       :: m, Mstart, Mrank, MM, NN, Nmaxim, Nrankrank, Nmax, Nmax1, Nmax2, &
                   Nmaxmax, Nteta, i, NthetaConv, NrankAL  
  real(O)       :: alfap, tetaGI, phiGI, phiGS, alfa, beta, gama, Cscat, Qscat,     &
                   Cext, Qext, z1z2 
  logical       :: type1, type2
  real(O),allocatable    :: h(:), v(:), oldh(:), oldv(:)
  complex(O),allocatable :: CC1(:,:), CC2(:,:), CC3(:,:), aa(:,:), bb(:,:), a(:,:), &
                            b(:,:), c(:,:), t1(:), t2(:), c1(:), c2(:), ce(:)   
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
  MM = Nrank
  NN = Nrank  + 5
  Nmaxim  = Nrank  + Mrank * (2 * Nrank - Mrank + 1) 
  Nrankrank = Nrank1 + Nrank2  
  NrankAL = max(Nrank,Nrank1,Nrank2)
  call write_TypeConvHead (3)
  open (unit = iTmat, file = FileTmat, status = 'replace')
  call  write_HeadFileTmat (NrankAL, NrankAL) 
  allocate (CC1(0:NN,0:2*NN+1), CC2(0:NN,0:2*NN+1), CC3(0:NN,0:2*NN+1))   
  allocate (aa(2*Nrankrank,2*Nrankrank), bb(2*Nrankrank,2*NrankAL))        
  allocate (a(2*NrankAL,2*NrankAL), b(2*NrankAL,2*NrankAL), c(2*NrankAL,2*NrankAL))
  allocate (t1(2*Nrank1), t2(2*Nrank2))
  allocate (c1(2*NrankAL), c2(2*NrankAL), ce(2*Nmaxim))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do 
  z1z2 = abs(z1 - z2)
  if (z1 > z2) then
    type1 = .false.
    type2 = .true.
  else
    type1 = .true.
    type2 = .false.
  end if
  Mrank = - 1
  do m = Mstart, MM    
    call write_1ConvParam (m)
    Mrank = Mrank + 1       
    if (m > Nrank1) Nrank1 = m
    if (m > Nrank2) Nrank2 = m
    if (m == 0) then
      Nmax  = Nrank
      Nmax1 = Nrank1
      Nmax2 = Nrank2
    else
      Nmax  = Nrank - m + 1
      Nmax1 = Nrank1 - m + 1
      Nmax2 = Nrank2 - m + 1
    end if
    Nmaxmax = Nmax1 + Nmax2
    if (PrnProgress) call write_progress_m (.true., m, 1, 8)        
    call MatTransCmnn1 (3, k, z1z2, m, NN, CC3, NN)
    call MatTransCmnn1 (1, k, z1, m, NN, CC1, NN)     
    call MatTransCmnn1 (1, k, z2, m, NN, CC2, NN) 
    call coefficients_fg_m (k, r1, ind_ref1, m, Nrank1, Nmax1, t1)  
    call coefficients_fg_m (k, r2, ind_ref2, m, Nrank2, Nmax2, t2)  
    call identity_matrix (2*Nmaxmax, aa, 2*Nrankrank, 2*Nrankrank)     
    call MatTransAB_mnn1 (k, type1, z1z2, m, CC3, NN, Nmax1, Nmax2, a,              &
         NrankAL, NrankAL)  
    call product_matrix_vector1 (2*Nmax1, 2*Nmax2, t1, a, 2*NrankAL, 2*NrankAL,     &
         b, 2*NrankAL, 2*NrankAL)   
    call extend_matrix3 (Nmax1, Nmax2, 0, 2*Nmax1, b, NrankAL, NrankAL, aa,         &
         Nrankrank, Nrankrank)   
    if (PrnProgress) call write_progress_m (.false., m, 2, 8)         
    call MatTransAB_mnn1 (k, type2, z1z2, m, CC3, NN, Nmax2, Nmax1, a,              &
         NrankAL, NrankAL)
    call product_matrix_vector1 (2*Nmax2, 2*Nmax1, t2, a, 2*NrankAL, 2*NrankAL,     &
         b, 2*NrankAL, 2*NrankAL)
    call extend_matrix3 (Nmax2, Nmax1, 2*Nmax1, 0, b, NrankAL, NrankAL, aa,         &
         Nrankrank, Nrankrank)      
    if (PrnProgress) call write_progress_m (.false., m, 3, 8)   
    call MatTransAB_mnn1 (k,.false., z1, m, CC1, NN, Nmax1, Nmax, a,                &
         NrankAL, NrankAL)
    call product_matrix_vector1 (2*Nmax1, 2*Nmax, t1, a, 2*NrankAL, 2*NrankAL,      &
         b, 2*NrankAL, 2*NrankAL)
    call extend_matrix4 (1, Nmax1, Nmax, Nmaxmax, 0, b, NrankAL, NrankAL, bb,       &
         Nrankrank, NrankAL) 
    if (PrnProgress) call write_progress_m (.false., m, 4, 8)   
    call MatTransAB_mnn1 (k,.false., z2, m, CC2, NN, Nmax2, Nmax, a,                &  
         NrankAL, NrankAL)
    call product_matrix_vector1 (2*Nmax2, 2*Nmax, t2, a, 2*NrankAL, 2*NrankAL,      &
         b, 2*NrankAL, 2*NrankAL)
    call extend_matrix4 (2, Nmax2, Nmax, Nmaxmax, 2*Nmax1, b, NrankAL, NrankAL,     &
         bb, Nrankrank, NrankAL)             
    if (PrnProgress) call write_progress_m (.false., m, 5, 8)   
    call LU_SYSTEM_DIRECT (aa, 2*Nrankrank, 2*Nrankrank, bb, 2*Nrankrank,           &
         2*NrankAL, 2*Nmaxmax,2*Nmax)  
    if (PrnProgress) call write_progress_m (.false., m, 6, 8)   
    call extract_matrix2 (Nmax1, Nmax, 0, b, NrankAL, NrankAL, bb, Nrankrank,       &
         NrankAL)
    call MatTransAB_mnn1 (k,.true., z1, m, CC1, NN, Nmax, Nmax1, a,                 &
         NrankAL, NrankAL)
    call product_matrices2 (1, 2*Nmax, 2*Nmax1, 2*Nmax, a, 2*NrankAL, 2*NrankAL,    &
         b, 2*NrankAL, 2*NrankAL, c, 2*NrankAL, 2*NrankAL)    
    if (PrnProgress) call write_progress_m (.false., m, 7, 8)   
    call extract_matrix2 (Nmax2, Nmax, 2*Nmax1, b, NrankAL, NrankAL, bb,            &
         Nrankrank, NrankAL)
    call MatTransAB_mnn1 (k,.true., z2, m, CC2, NN, Nmax, Nmax2, a,                 &
         NrankAL, NrankAL)
    call product_matrices2 (2, 2*Nmax, 2*Nmax2, 2*Nmax, a, 2*NrankAL, 2*NrankAL,    &
         b, 2*NrankAL, 2*NrankAL, c, 2*NrankAL, 2*NrankAL)                      
    if (PrnProgress) call write_progress_m (.false., m, 8, 8)   
    call write_FileTmat (NrankAL, NrankAL, c)
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,     &
         Nmax, c1)
    call product_matrix_vector (2*Nmax, 2*Nmax, c, 2*NrankAL, 2*NrankAL, c1, c2)
    call extend_vector_positive (c2, ce, m, Mstart, Nrank, Nmax, Nmaxim)
    if (m /= 0) then      
      call matrix_m_negativ (Nmax, Nmax, c, NrankAL, NrankAL)         
      call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap,-m, Nrank,   &
           Nmax, c1)        
      call product_matrix_vector (2*Nmax, 2*Nmax, c, 2*NrankAL, 2*NrankAL, c1, c2)
      call extend_vector_negative (c2, ce, m, Nrank, Nmax, Nmaxim)      
    end if
    call DSCS (ce, Mrank, Nrank, Nmaxim, Nteta, phiGS, alfa, beta, gama, k, snorm,  &
        .false.,.true., h, v)
    call delta_DSCS (Nteta, h, v, oldh, oldv, epsMrank, NthetaConv)     
    call write_DSCS (Nteta,.false., h, v)
    if (NthetaConv >= int(0.8*Nteta)) exit 
  end do
  close (unit = iTmat)
  call CQscat (ce, Mrank, Nrank, Nmaxim, k, snorm, Cscat, Qscat)
  call CQext (ce, Mrank, Nrank, Nmaxim, tetaGI, phiGI, alfa, beta, gama,            &
       alfap, k, snorm, Cext, Qext)
  call write_Effic (Qscat, Qext)
  call write_MrankConvRes (NthetaConv, epsMrank)
  if (NthetaConv >= int(0.8*Nteta)) then
    print "(/,2x,'Convergence criterion for Mrank is satisfied;')"
  else
    print "(/,2x,'Convergence criterion for Mrank is not satisfied;')"
  end if
  call write_InfoFileTmat (FileTmat, Mrank, Nrank, .true., .false., .false.)
  call ScatCharact (k, FileTmat, Mrank, Nrank, .true., .false., .false.)
  print "(/,2x,'The T matrix is stored in ',a50)", FileTmat
  print "(  2x,'The dimensions of the T matrix are given by:')"
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", Nrank 
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank                     
  deallocate (CC1, CC2, CC3, aa, bb, a, b, c, c1, c2, ce, t1, t2, h, v, oldh, oldv) 
end subroutine convergence_MrankMULT2SPH

       
