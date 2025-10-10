subroutine TINHOMSPH
!------------------------------------------------------------------------------------
! 1. General Considerations                                                         !
! --------------------------                                                        !
! TINHOMSPH is a code for computing the T matrix and the scattering characteristics !
! of an an inhomogeneous, dielectric sphere with an arbitrarily shaped inclusion.   !
! The program supports calculation for dielectric host spheres with REAL refractive !
! index. The main feature of this routine is that the T matrix of the inclusion is  !
! provided as input parameter in the file FileTmat. The inclusion can be a          !
! homogeneous, axisymmetric or nonaxisymmetric particle, a composite or a layered   !
! scatterer and an aggregate. For computing the T matrix of the inclusion, the      !
! refractive index of the ambient medium (the parameter ind_refMed in the           !
! corresponding T-matrix program) must be identical with the relative refractive    !
! index of host particle (the parameter indrefRel in the present T-matrix program). ! 
! Note that the T matrix of a spherical inclusion must be computed with the TAXSYM  !
! routine and not with the TSPHERE routine (the T matrix of the inclusion must be   !
! specified as a two-dimensional array and the TSPHERE routine provides an          !
! one-dimensional array).                                                           ! 
!                                                                                   !
! 2. Convergence Test                                                               !
! --------------------                                                              !
! Convergence tests precede the scattering characteristics calculation. The         !
! parameters which control the T-matrix computation are the maximum expansion order !
! Nrank and the maximum azimuthal order Mrank for the host particle.                !
!                                                                                   !
! For convergence tests, the incident wave is assumed to be a plane wave traveling  !
! along the Z-axis of the global coordinate system and the scattering               !
! characteristics are computed in the azimuthal plane  phi = 0°. The Euler          !
! orientation angles of the host sphere are alpha = beta = 45° and gamma = 0°, and  !
! interactive convergence tests over Nrank and Mrank are performed.                 !
!                                                                                   !  
! Nrank and Mrank are assumed to not depend on Nrank1 and Mrank1, which are the     !
! maximum expansion and azimuthal orders of the inclusion. The convergence tests    !
! are carried out over Nrank and Mrank, while Nrank1 and Mrank1 are kept constant.  !
! For the convergence test over the expansion order, the scattering problem is      !
! solved for the pairs (Nrank, Mrank) and (Nrank - 1, Mrank), while for the         !
! azimuthal order test, the cases (Nrank, Mrank) and (Nrank, Mrank - 1) are         !
! considered. The normalized differential scattering cross section will be checked  !
! at 20° increments for convergence within epsX (epsNint,epsNrank or epsMrank)      !
! tolerance. If the calculated results are converged within this tolerance at       !
! 80% of the scattering angles, then convergence is achieved. The T matrix is       !
! stored for later use by other programs, and the values of Nrank and Mrank are     !
! printed to the screen and to the T-matrix information file (see                   !
! "Description.txt"). These values together with the T matrix serve as INPUT        !
! PARAMETERS for other programs.                                                    !
!                                                                                   !
! 3. Estimates of Nrank and Mrank                                                   !
! --------------------------------                                                  !
! The above convergence tests require estimates of Nrank and Mrank. These estimates !
! must be supplied by the user, and for this purpose, Wiscombe's truncation limit   !
! criterion [W. J. Wiscombe, Improved Mie scattering algorithms, Applied Optics,    !
! 19, 1505-1509, 1980] can be used.                                                 !
!                                                                                   !
! The truncation limit criterion proposed by Wiscombe provides the estimate         !
!                                                                                   ! 
!                NrankW = int(x + 4.05 * x**0.33 + 2),                              !
!                                                                                   !
! where x is the size parameter, x = k * r, k is the wave number and r is the       !
! radius of the smallest circumscribing sphere. For the maximum azimuthal order     !
! test, we may use the conservative estimate Mrank = Nrank - 2,..,Nrank.            ! 
!                                                                                   !
! 4. Strategies for Performing Convergence Tests                                    !
! -----------------------------------------------                                   !
! The following strategy for performing the convergence test can be employed:       !
! 1. choose a value of Nrank close to the value predicted by Wiscombe's truncation  !
!    limit criterion (or smaller) and set Mrank = Nrank - 2 (or smaller);           !
! 2. perform the convergence tests over Nrank and Mrank;                            ! 
! 3. if convergence is not achieved, set Nrank = Nrank + 1, maintain the relation   !
!    Mrank = Nrank - 2 and go to Step 2.                                            !
! For particles which are too extreme in size, the errors of the extinction and     !
! scattering cross sections (instead of the differential scattering cross section)  !
! can be analyzed.                                                                  ! 
!                                                                                   !
! 5. Additional Comments                                                            !
! -----------------------                                                           !
! The convergence tests can be switched off by setting the logical variable         !
! DoConvTest to false. In this case, the values of Nrank and Mrank must be          !
! specified in the input file.                                                      !
!                                                                                   !
! 6. Input Parameters                                                               !
! --------------------                                                              !
! The parameters specified in the input file "/INPUTFILES/InputINHOMSPH.dat"        !
! are listed below.                                                                 !
!                                                                                   !
! - wavelength (real) - wavelength of the incident light in vacuo.                  !
!                                                                                   !
! - ind_refMed (real) - refractive index of the ambient medium (nonabsorbing        !
!   dielectric medium).                                                             !
!                                                                                   !
! - ind_refRel (REAL) - relative refractive index of the host sphere with           ! 
!   respect to the ambient medium.                                                  !
!                                                                                   !
! - r (real) - radius of the host sphere.                                           ! 
!                                                                                   !
! - anorm (real) - characteristic length of the host sphere which is used to        !
!   normalize the differential scattering cross sections.                           !
!                                                                                   !
! - FileTmat (character(80)) - name of the file containing the T matrix of the      !
!   inclusion.                                                                      !
!                                                                                   !
! - axsym1 (logical) - if axsym1 = t, the inclusion is a rotationally symmetric     !
!   particle (axisymmetric particle).                                               !
!                                                                                   !
! - chiral1 (logical) - if chiral1 = t, the inclusion is an optical active particle  !
!   (chiral particle).                                                              !
!                                                                                   !
! - Nrank1, Mrank1 (integer variables) - the maximum expansion and azimuthal        !
!   orders of the inclusion.                                                        !
!                                                                                   !   
! - x1, y1, z1 (real variables) - Cartesian coordinates specifying the position     !
!   of the inclusion with respect to the coordinate system of the host sphere.      !
!                                                                                   ! 
! - alpha1, beta1, gamma1 (real variables) - Euler angles specifying the            !
!   orientation of the coordinate system of the inclusion with respect to the       !
!   coordinate system of the host sphere.                                           !
!                                                                                   !
! - DoConvTest (logical) - if DoConvTest = t, the interactive convergence tests     !
!   over Nrank and Mrank are invoked. An estimates of Nrank is given by             !
!   Wiscombe's truncation limit criterion. If DoConvTest = f, the values of Nrank   !
!   and Mrank must be supplied in the input file.                                   !
!                                                                                   !
! - ExtThetaDom (logical) - if ExtThetaDom = t the DSCS is computed for             !
!   scattering angles ranging from 0° to 180° in the azimuthal plane phiGS = 0°,    !
!   and from 180° to 0° in the azimuthal plane phiGS = 180°. The total number of    !
!   scattering angles is 10. If ExtThetaDom = f the DSCS is computed for            !
!   scattering angles ranging from 0 to 180° in the azimuthal plane phiGS = 0.      !
!                                                                                   !
! - Nrank (integer) - maximum expansion order for the host sphere. This             !
!   parameter is used if the convergence tests are not performed (DoConvTest = f).  !
!                                                                                   !
! - Mrank (integer) - maximum azimuthal order for the host sphere. This             !
!   parameter is used if the convergence test are not performed (DoConvTest = f).   !
!                                                                                   !
!   The following table explain the significance of the variables Nrank and Mrank.  !
!                                                                                   !
!                  DoConvTest       Nrank and Mrank                                 !
!                     t             console input                                   !
!                     f              file input                                     ! 
!                                                                                   !
! - epsNrank (real) - error tolerance for the expansion order test.                 !
!                                                                                   !
! - epsMrank (real) - error tolerance for the azimuthal order test.                 !
!                                                                                   !
! - FileTmatG (character(80)) - name of the file to which the T matrix of the       !
!   inhomogeneous spherical particle is written.                                    !
!                                                                                   !
! - PrnProgress (logical) - if PrnProgress = t, the progress of calculation         !
!   is printed.                                                                     !
!                                                                                   !
! Note that all input parameters specifying lengths must be provided in the same    !
! units of lengths.                                                                 !
!                                                                                   !
! 7. Logical Scheme                                                                 !
! ------------------                                                                !  
! The organization of the code is as follows:                                       !
! 1. < read input data; the input file provides:                                    !
!    - the name of the file containing the T matrix of the inclusion,               !
!    - the maximum expansion and azimuthal orders Nrank1 and Mrank1,                ! 
!    - the Cartesian coordinates specifying the position of the inclusion, and      !
!    - the Euler angles specifying the orientation of the inclusion >               !   
!                                                                                   !
! 2. < T-matrix calculation >                                                       !
!                                                                                   !
!    if ( do convergence test ) then                                                !
!                                                                                   !
!      < the code computes an estimate of Nrank accordingly to                      !
!        Wiscombe's truncation criterion >                                          !
!      < the code prompts for the estimated values of Nrank and Mrank >             !
!                                                                                   !
!    else if ( .not. do convergence test ) then                                     !
!                                                                                   !
!       < the input file provides the values of Nrank and Mrank >                   !
!                                                                                   !
!    end if                                                                         !
!    if ( do convergence test ) then                                                !
!       < the code computes the DSCS for (Nrank,Mrank), (Nrank - 1,Mrank) and       ! 
!         (Nrank, Mrank - 1) and write the results to the file "Output.dat" >       !
!       < the T matrix is stored in the file FileTmatG and the T-matrix             !
!         information file is created >                                             !
!       < Mrank and Nrank are written to the screen and to the T-matrix             !
!         information file >                                                        !
!       < go to Step 3 >                                                            !
!    end if                                                                         !
!    if ( .not. do convergence test ) then                                          !
!       < the T matrix is computed and stored in the file FileTmatG >               !
!       < the T-matrix information file is created >                                !
!       < Mrank and Nrank are written to the screen and to the T-matrix             !
!         information file >                                                        !
!       < go to Step 3 >                                                            !
!    end if                                                                         !
! 3. < scattering characteristics calculation >                                     !   
!------------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none 
  integer       :: Mrank, Nrank, Nrank1, Mrank1                    
  real(O)       :: ks, ind_refMed, ind_refRel, wavelength, r, anorm, snorm,         &
                   epsNrank, epsMrank, x1, y1, z1, alpha1, beta1, gamma1 
  character(80) :: FileTmat, FileTmatG
  logical       :: axsym1, chiral1, DoConvTest, ExtThetaDom, PrnProgress 
! -----------------------------------------------------------------------------------
!                              Read the input file                                  ! 
! -----------------------------------------------------------------------------------  
  call readinputINHOMSPH ( wavelength, ind_refMed, ind_refRel, r, anorm,             &
       FileTmat, axsym1, chiral1, Nrank1, Mrank1, x1, y1, z1, alpha1, beta1,         &
       gamma1, DoConvTest, ExtThetaDom, Nrank, Mrank, epsNrank, epsMrank,            &
       FileTmatG, PrnProgress, ks, snorm ) 
! -----------------------------------------------------------------------------------
!                                      Main                                         !
! -----------------------------------------------------------------------------------
  if (DoConvTest) then    
    open (unit = iOutput, file = FileOutput, status = "replace") 
    call printinputINHOMSPH (r, FileTmat, axsym1, x1, y1, z1, alpha1, beta1, gamma1,&
         Mrank1, Nrank1, wavelength, anorm, ind_refMed, ind_refRel, epsNrank,       &
         epsMrank)
    call convergence_Nrank_MrankINHOMSPH (r, ks, ind_refRel, snorm, Mrank, Nrank,   &
         x1, y1, z1, alpha1, beta1, gamma1, Mrank1, Nrank1, FileTmat, axsym1,       &
         chiral1, epsNrank, epsMrank, ExtThetaDom, FileTmatG, PrnProgress)   
    close (unit = iOutput)
  else
    call TMatrix_Nrank_MrankINHOMSPH (r, ks, ind_refRel, Mrank, Nrank, x1, y1, z1,  &
         alpha1, beta1, gamma1, Mrank1, Nrank1, FileTmat, axsym1, chiral1,          &
         FileTmatG, PrnProgress)
  end if
end subroutine TINHOMSPH
!***********************************************************************************
subroutine readinputINHOMSPH ( wavelength, ind_refMed, ind_refRel, r, anorm,        &
           FileTmat, axsym1, chiral1, Nrank1, Mrank1, x1, y1, z1, alpha1, beta1,    &
           gamma1, DoConvTest, ExtThetaDom, Nrank, Mrank, epsNrank, epsMrank,       &
           FileTmatG, PrnProgress, ks, snorm )     
  use parameters
  use derived_parameters
  implicit none 
  integer       :: Mrank, Nrank, NrankW, Mrank1, Nrank1, ios                     
  real(O)       :: ks, ind_refMed, ind_refRel, wavelength, r, anorm, xpart, snorm,  &
                   epsNrank, epsMrank, x1, y1, z1, alpha1, beta1, gamma1, x, grd 
  character(80) :: FileTmat, FileTmatG, string
  logical       :: axsym1, chiral1, DoConvTest, ExtThetaDom, PrnProgress, XFindPar 
! -----------------------------------------------------------------------------------
!                        Read the input file FileInputINHOMSPH                      ! 
! ----------------------------------------------------------------------------------- 
  call DrvParameters 
  open (unit = iInputINHOMSPH, file = FileInputINHOMSPH, status = "old",            &
        position = "rewind")   
  wavelength = 0.1_O * 2._O * Pi
  ind_refMed = 1._O                                                                 
  ind_refRel = 1.2_O
  string     = 'OptProp'
  if (XFindPar (iInputINHOMSPH, string)) then
    read (iInputINHOMSPH, *, iostat = ios) wavelength
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable wavelength;')"
      stop
    end if
    read (iInputINHOMSPH, *, iostat = ios) ind_refMed
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refMed;')"
      stop
    end if
    read (iInputINHOMSPH, *, iostat = ios) ind_refRel
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refRel;')"
      stop
    end if      
  else
    print "(/,2x,'Group name OptProp not found;')"
    stop  
  end if  
  ks = 2._O * Pi * ind_refMed / wavelength
!       
  r = 1._O
  anorm  = 1._O
  string = 'GeomPropHost'
  if (XFindPar (iInputINHOMSPH, string)) then
    read (iInputINHOMSPH, *, iostat = ios) r
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable r;')"
      stop
    end if
    read (iInputINHOMSPH, *, iostat = ios) anorm
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
  FileTmat = '../TMATFILES/T.dat'
  axsym1   = .true.
  chiral1  = .false. 
  Nrank1   = 6
  Mrank1   = 4 
  string   = 'TmatIncl'
  if (XFindPar (iInputINHOMSPH, string)) then
    read (iInputINHOMSPH, *, iostat = ios) FileTmat
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable FileTmat;')"
      stop
    end if
    read (iInputINHOMSPH, *, iostat = ios) axsym1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable axsym1;')"
      stop
    end if 
    read (iInputINHOMSPH, *, iostat = ios) chiral1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable chiral1;')"
      stop
    end if
    read (iInputINHOMSPH, *, iostat = ios) Nrank1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nrank1;')"
      stop
    end if
    read (iInputINHOMSPH, *, iostat = ios) Mrank1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Mrank1;')"
      stop
    end if
  else
    print "(/,2x,'Group name TmatIncl not found;')"
    stop  
  end if  
  call check_MrankNrank (Mrank1, Nrank1) 
!
  x1 = 0.1_O
  y1 = 0.1_O
  z1 = 0.1_O
  alpha1 = 45._O
  beta1  = 45._O
  gamma1 = 0._O
  string = 'GeomPropIncl'
  if (XFindPar (iInputINHOMSPH, string)) then
    read (iInputINHOMSPH, *, iostat = ios) x1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable x1;')"
      stop
    end if
    read (iInputINHOMSPH, *, iostat = ios) y1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable y1;')"
      stop
    end if 
    read (iInputINHOMSPH, *, iostat = ios) z1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable z1;')"
      stop
    end if
    read (iInputINHOMSPH, *, iostat = ios) alpha1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable alpha1;')"
      stop
    end if
    read (iInputINHOMSPH, *, iostat = ios) beta1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable beta1;')"
      stop
    end if
    read (iInputINHOMSPH, *, iostat = ios) gamma1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable gamma1;')"
      stop
    end if
  else
    print "(/,2x,'Group name GeomPropIncl not found;')"
    stop  
  end if    
  grd    = Pi / 180._O 
  alpha1 = alpha1 * grd
  beta1  = beta1 * grd
  gamma1 = gamma1 * grd
!
  DoConvTest  = .true.
  ExtThetaDom = .true.
  string      = 'ConvTest'
  if (XFindPar (iInputINHOMSPH, string)) then
    read (iInputINHOMSPH, *, iostat = ios) DoConvTest
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable DoConvTest;')"
      stop
    end if
    read (iInputINHOMSPH, *, iostat = ios) ExtThetaDom
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ExtThetaDom;')"
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
    print "(/,2x,'T-Matrix Computation for an Inhomogeneous Sphere')"
    print "(  2x,'------------------------------------------------')"
  end if
!  
  x = ks * r
  NrankW = int(x + 4.05_O * x**0.33_O + 2._O)
  if (.not. DoConvTest) then
    Nrank  =  16
    Mrank  =   8
    string = 'NrankMrankHost'
    if (XFindPar (iInputINHOMSPH, string)) then
      read (iInputINHOMSPH, *, iostat = ios) Nrank
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Nrank;')"
        stop
      end if
      read (iInputINHOMSPH, *, iostat = ios) Mrank
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Mrank;')"
        stop
      end if         
    else
      print "(/,2x,'Group name NrankMrankHost not found;')"
      stop  
    end if    
    print "(/,2x,'Input values:')"
    print "(  2x, a, i3, a, i3, a)",                                                &
   'the input values of Nrank and Mrank are ', Nrank, ' and ', Mrank, ', while'
    print "(  2x, a, i3, a)",                                                       &
   'the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';'                           
  else       
    print "(/,2x,'Nrank estimate:')"                                                    
    print "(  2x, a, i3, a)",                                                       &
   'the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';'
    print "(/,2x,'- enter the estimated values of Nrank and Mrank, where')"
    print "(  2x,'  Mrank = Nrank - 2,...,Nrank;')"     
    call read_integer2 (Nrank, Mrank)           
  end if  
  call check_MrankNrank (Mrank, Nrank)
  if (Nrank < Nrank1) print "(/,2x,'warning: Nrank is too low;')"
  if (Mrank < Mrank1) print "(/,2x,'warning: Mrank is too low;')" 
!
  epsNrank = 5.e-2_O
  epsMrank = 5.e-2_O 
  string   = 'Errors'
  if (XFindPar (iInputINHOMSPH, string)) then
    read (iInputINHOMSPH, *, iostat = ios) epsNrank 
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsNrank;')"
      stop
    end if
    read (iInputINHOMSPH, *, iostat = ios) epsMrank 
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsMrank;')"
      stop
    end if         
  else
    print "(/,2x,'Group name Errors not found;')"
    stop  
  end if        
!
  FileTmatG = '../TMATFILES/TG.dat'
  string    = 'Tmat'
  if (XFindPar (iInputINHOMSPH, string)) then
    read (iInputINHOMSPH, *, iostat = ios) FileTmatG 
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable FileTmatG;')"
      stop
    end if            
  else
    print "(/,2x,'Group name Tmat not found;')"
    stop  
  end if   
!
  PrnProgress = .true.
  string      = 'PrintProgress'
  if (XFindPar (iInputINHOMSPH, string)) then
    read (iInputINHOMSPH, *, iostat = ios) PrnProgress 
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable PrnProgress;')"
      stop
    end if            
  else
    print "(/,2x,'Group name PrintProgress not found;')"
    stop  
  end if        
  close (unit = iInputINHOMSPH)  
end subroutine readinputINHOMSPH
!***********************************************************************************
subroutine printinputINHOMSPH (r, FileTmat, axsym1, x1, y1, z1, alfa1, beta1, gama1,&
           Mrank1, Nrank1, wavelength, anorm, ind_refMed, ind_refRel, epsNrank,     &
           epsMrank)
  use parameters
  implicit none
  integer       :: Mrank1, Nrank1, LenString
  real(O)       :: wavelength, anorm, ind_refMed, ind_refRel, r, x1, y1, z1, alfa1, &
                   beta1, gama1, epsNrank, epsMrank
  character(80) :: FileTmat, FileTmatWrite
  logical       :: axsym1
!
  write (iOutput,"(/,2x,'Input Parameters of the Scattering Problem:',/)")
  write (iOutput,"(2x,'wavelength of the free space, wavelength = ',1pe13.4,';')")  &
         wavelength
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'refractive index of the ambient medium, ind_refMed = ', ind_refMed, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'relative refractive index of the host particle, ind_refRel = ', ind_refRel, ';'
  write (iOutput,*)
  write (iOutput,"(2x,'radius of the host particle, r = ',1pe10.3,';')") r  
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'characteristic length of the host particle, anorm = ', anorm, ';'
  write (iOutput,*)
  FileTmatWrite = FileTmat(14:LenString(FileTmat))
  write (iOutput,"(2x, a, a)")                                                      &
 'name of the file containing the T matrix of the inclusion = ', FileTmatWrite
  if (axsym1) then
    write (iOutput,"(2x,'axisymmetric inclusion;')")     
  else
    write (iOutput,"(2x,'nonaxisymmetric inclusion;')")          
  end if
  write (iOutput,"(2x,'maximum expansion order, Nrank =',i3,';')") Nrank1
  write (iOutput,"(2x,'maximum azimuthal order, Mrank =',i3,';')") Mrank1
  write (iOutput,"(2x, a, 1pe10.3, a, 1pe10.3, a, 1pe10.3, a)")                     &
 'Cartesian coordinates of the inclusion: x = ', x1, ', y = ', y1, ', z = ', z1, ';'             
  write (iOutput,"(2x, a, 1pe10.3, a, 1pe10.3, a, 1pe10.3, a)")                     &
 'Euler angles: alpha = ', alfa1 * 180._O / Pi, ', beta = ', beta1 * 180._O / Pi,   &
 ', gamma = ', gama1 * 180._O / Pi, ';'  
  write (iOutput,*)
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum expansion order tolerance, epsNrank = ', epsNrank, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum azimuthal order tolerance, epsMrank = ', epsMrank, '.'          
  write (iOutput,"(/)")                       
end subroutine printinputINHOMSPH
!***********************************************************************************
subroutine convergence_Nrank_MrankINHOMSPH (Rb, ks, ind_ref, snorm, Mrank, Nrank,   &
           x1, y1, z1, alfa1, beta1, gama1, Mrank1, Nrank1, FileTmat, axsym1,       &
           chiral1, epsNrank, epsMrank, ExtThetaDom, FileTmatG, PrnProgress)
  use parameters
  implicit none
  integer       :: Mrank, Nrank, Mrank1, Nrank1
  real(O)       :: ks, ind_ref, Rb, snorm, x1, y1, z1, alfa1, beta1, gama1,         &
                   epsNrank, epsMrank  
  character(80) :: FileTmat, FileTmatG
  logical       :: axsym1, chiral1, ExtThetaDom, PrnProgress 
!
  integer       :: Nmax, Nmax1, Nteta, i, NthetaConvN, NthetaConvM, NmaxAL,         &
                   nt1g, mt1g
  real(O)       :: k, alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,  &
                   Cext, Qext, r1, sumang
  complex(O)    :: ind_refC
  real(O),allocatable    :: h(:), v(:), oldh(:), oldv(:), oldh0(:), oldv0(:)
  complex(O),allocatable :: a(:,:), b(:,:), s(:,:), c(:), c1(:), t1(:,:), t(:),     &
                            t11(:), t13(:), t31(:), t33(:)     
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
  Nmax1  = Nrank1 + Mrank1 * (2 * Nrank1 - Mrank1 + 1) 
  NmaxAL = max(Nmax, Nmax1)
  k = ks * ind_ref  
  ind_refC = cmplx(ind_ref,0.0,O)  
  call write_TypeConvHead (4)
  allocate (a(2*NmaxAL,2*NmaxAL), b(2*NmaxAL,2*NmaxAL), s(2*NmaxAL,2*NmaxAL),       &
            c(2*NmaxAL), c1(2*NmaxAL))
  allocate (t(2*NmaxAL), t11(2*NmaxAL), t13(2*NmaxAL), t31(2*NmaxAL), t33(2*NmaxAL))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta), oldh0(Nteta), oldv0(Nteta))
  open (unit = iTmat, file = FileTmat, status = "old", position = "rewind")
  if (.not. axsym1) then
    call read_HeadFileTmat (nt1g, mt1g)
    call check_dimensionMat (nt1g, mt1g, Nmax1)       
    allocate (t1(2*nt1g,2*mt1g))
    call read_FileTmat (nt1g, mt1g, t1) 
  else
    nt1g = Nmax1
    mt1g = Nmax1
    allocate (t1(2*nt1g,2*mt1g))
    call read_Tmatrix (chiral1, Nrank1, Mrank1, Nmax1, t1, nt1g, mt1g)            
  end if
  close (unit = iTmat)
  open  (unit = iTmat, file = FileTmatG, status = 'replace')
  call  write_HeadFileTmat (NmaxAL, NmaxAL)  
  if (PrnProgress) call  write_progress (.true., 1, 6)            
  r1 = sqrt(x1**2 + y1**2 + z1**2)
  sumang = abs(alfa1) + abs(beta1) + abs(gama1)
  if (r1 /=0 .and. sumang /= 0) then
    call MatTransRot_TR_mn_m1n1 (1, k, x1, y1, z1, alfa1, beta1, gama1, Mrank,      &
         Nrank, Nmax, Mrank1, Nrank1, Nmax1, s, NmaxAL, NmaxAL)
  else if (r1 /=0 .and. sumang == 0) then 
    call MatTransAB_mn_m1n1 (1, k, x1, y1, z1, Mrank, Nrank, Nmax, Mrank1,          &
         Nrank1, Nmax1, s, NmaxAL, NmaxAL)
  else if (r1 == 0 .and. sumang /= 0) then    
    call MatRot_R_mn_m1n1 (alfa1, beta1, gama1, Mrank, Nrank, Nmax, Mrank1,         &
         Nrank1, Nmax1, s, NmaxAL, NmaxAL, .true.)
  else if (r1 == 0 .and. sumang == 0) then
    call identity_transformation (Mrank, Nrank, Nmax, Mrank1, Nrank1, Nmax1, s,     &
         NmaxAL, NmaxAL)    
  end if  
  call product_matrices (2*Nmax, 2*Nmax1, 2*Nmax1, s, 2*NmaxAL, 2*NmaxAL, t1,       &
       2*nt1g, 2*mt1g)
  if (PrnProgress) call write_progress (.false., 2, 6)
  if (r1 /=0 .and. sumang /= 0) then
    call MatRotTrans_RT_mn_m1n1 (1, k, x1, y1, z1, alfa1, beta1, gama1, Mrank1,     &
         Nrank1, Nmax1, Mrank, Nrank, Nmax, a, NmaxAL, NmaxAL)
  else if (r1 /=0 .and. sumang == 0) then 
    call MatTransAB_mn_m1n1 (1, k, -x1, -y1, -z1, Mrank1, Nrank1, Nmax1, Mrank,     &
         Nrank, Nmax, a, NmaxAL, NmaxAL)
  else if (r1 == 0 .and. sumang /= 0) then    
    call MatRot_R_mn_m1n1 (-gama1, -beta1, -alfa1, Mrank1, Nrank1, Nmax1, Mrank,    &
         Nrank, Nmax, a, NmaxAL, NmaxAL, .true.)
  else if (r1 == 0 .and. sumang == 0) then
    call identity_transformation (Mrank1, Nrank1, Nmax1, Mrank, Nrank, Nmax,        &
         a, NmaxAL, NmaxAL)    
  end if   
  call product_matrices (2*Nmax, 2*Nmax1, 2*Nmax, s, 2*NmaxAL, 2*NmaxAL, a,         &
       2*NmaxAL, 2*NmaxAL) 
  if (PrnProgress) call write_progress (.false., 3, 6)      
  call vector_Q_sphere (3, 3, ks, Rb, ind_refC, Mrank, Nrank, Nmax, t)
  call copy_vector (t, t33, 2*Nmax)     
  call product_matrix_vector1 (2*Nmax, 2*Nmax, t, s, 2*NmaxAL, 2*NmaxAL,            &
       a, 2*NmaxAL, 2*NmaxAL)
  call vector_Q_sphere (3, 1, ks, Rb, ind_refC, Mrank, Nrank, Nmax, t)
  call copy_vector (t, t31, 2*Nmax)     
  call sum_diagonal_elements (2*Nmax, t, a, 2*NmaxAL, 2*NmaxAL) 
  if (PrnProgress) call write_progress (.false., 4, 6)
  call vector_Q_sphere (1, 3, ks, Rb, ind_refC, Mrank, Nrank, Nmax, t)
  call copy_vector (t, t13, 2*Nmax)     
  call product_matrix_vector1 (2*Nmax, 2*Nmax, t, s, 2*NmaxAL, 2*NmaxAL,            &
       b, 2*NmaxAL, 2*NmaxAL)
  call vector_Q_sphere (1, 1, ks, Rb, ind_refC, Mrank, Nrank, Nmax, t)
  call copy_vector (t, t11, 2*Nmax)     
  call sum_diagonal_elements (2*Nmax, t, b, 2*NmaxAL, 2*NmaxAL)  
  if (PrnProgress) call write_progress (.false., 5, 6)
  call LU_SYSTEM (a, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL, 2*Nmax)
  if (PrnProgress) call write_progress (.false., 6, 6)
  call minus_matrix (2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL)
  call write_FileTmat (NmaxAL, NmaxAL, b)       
  call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrank, Nrank,     &
       Nmax, c)    
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL, c, c1)
  call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfa, beta, gama, ks, snorm,     &
       ExtThetaDom,.true., h, v)
  call CQscat (c1, Mrank, Nrank, Nmax, ks, snorm, Cscat, Qscat)
  call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfa, beta, gama,              &
       alfap, ks, snorm, Cext, Qext)       
  call write_2ConvParam (Nrank, Mrank)  
  call write_DSCS (Nteta, ExtThetaDom, h, v)
  call write_Effic (Qscat, Qext)
  do i = 1, Nteta
    oldh(i)  = h(i)
    oldv(i)  = v(i)
    oldh0(i) = h(i)
    oldv0(i) = v(i)
  end do 
  close (unit = iTmat) 
! --- (Nrank - 1) configuration ---
  if (PrnProgress) call write_progress_low
  call copy_vector (t33, t, 2*Nmax)
  call vector_Nrank_1 (Mrank, Nrank, Nmax, t)   
  call product_matrix_vector1 (2*Nmax, 2*Nmax, t, s, 2*NmaxAL, 2*NmaxAL,            &
       a, 2*NmaxAL, 2*NmaxAL)
  call copy_vector (t31, t, 2*Nmax)
  call vector_Nrank_1 (Mrank, Nrank, Nmax, t)   
  call sum_diagonal_elements (2*Nmax, t, a, 2*NmaxAL, 2*NmaxAL)
  call matrix_Nrank_1_left (Mrank, Nrank, Nmax, a, NmaxAL, NmaxAL)  
  call copy_vector (t13, t, 2*Nmax)
  call vector_Nrank_1 (Mrank, Nrank, Nmax, t)   
  call product_matrix_vector1 (2*Nmax, 2*Nmax, t, s, 2*NmaxAL, 2*NmaxAL,            &
       b, 2*NmaxAL, 2*NmaxAL)
  call copy_vector (t11, t, 2*Nmax)
  call vector_Nrank_1 (Mrank, Nrank, Nmax, t)   
  call sum_diagonal_elements (2*Nmax, t, b, 2*NmaxAL, 2*NmaxAL)  
  call matrix_Nrank_1_right (Mrank, Nrank, Nmax, b, NmaxAL, NmaxAL)  
  call LU_SYSTEM (a, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL, 2*Nmax)
  call minus_matrix (2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL) 
  call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrank, Nrank,     &
       Nmax, c)    
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL, c, c1)
  call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfa, beta, gama, ks, snorm,     &
       ExtThetaDom,.true., h, v)
  call CQscat (c1, Mrank, Nrank, Nmax, ks, snorm, Cscat, Qscat)
  call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfa, beta, gama,              &
       alfap, ks, snorm, Cext, Qext)
  call delta_DSCS (Nteta, h, v, oldh, oldv, epsNrank, NthetaConvN)
  call write_2ConvParam (Nrank - 1, Mrank)
  call write_DSCS (Nteta, ExtThetaDom, h, v)
  call write_Effic (Qscat, Qext)
  call write_NrankConvRes (NthetaConvN, Nteta, epsNrank)
! --- (Mrank - 1) configuration ---
  call copy_vector (t33, t, 2*Nmax)
  call vector_Mrank_1 (Mrank, Nrank, Nmax, t)   
  call product_matrix_vector1 (2*Nmax, 2*Nmax, t, s, 2*NmaxAL, 2*NmaxAL,            &
       a, 2*NmaxAL, 2*NmaxAL)
  call copy_vector (t31, t, 2*Nmax)
  call vector_Mrank_1 (Mrank, Nrank, Nmax, t)   
  call sum_diagonal_elements (2*Nmax, t, a, 2*NmaxAL, 2*NmaxAL)
  call matrix_Mrank_1_left (Mrank, Nrank, Nmax, a, NmaxAL, NmaxAL) 
  call copy_vector (t13, t, 2*Nmax)
  call vector_Mrank_1 (Mrank, Nrank, Nmax, t)   
  call product_matrix_vector1 (2*Nmax, 2*Nmax, t, s, 2*NmaxAL, 2*NmaxAL,            &
       b, 2*NmaxAL, 2*NmaxAL)
  call copy_vector (t11, t, 2*Nmax)
  call vector_Mrank_1 (Mrank, Nrank, Nmax, t)   
  call sum_diagonal_elements (2*Nmax, t, b, 2*NmaxAL, 2*NmaxAL)  
  call matrix_Mrank_1_right (Mrank, Nrank, Nmax, b, NmaxAL, NmaxAL) 
  call LU_SYSTEM (a, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL, 2*Nmax)
  call minus_matrix (2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL)
  call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrank, Nrank,     &
       Nmax, c)    
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL, c, c1)
  call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfa, beta, gama, ks, snorm,     &
       ExtThetaDom,.true., h, v)
  call CQscat (c1, Mrank, Nrank, Nmax, ks, snorm, Cscat, Qscat)
  call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfa, beta, gama,              &
       alfap, ks, snorm, Cext, Qext)
  call delta_DSCS (Nteta, h, v, oldh0, oldv0, epsMrank, NthetaConvM)
  call write_2ConvParam (Nrank, Mrank - 1)
  call write_DSCS (Nteta, ExtThetaDom, h, v)
  call write_Effic (Qscat, Qext)
  call write_MrankConvRes (NthetaConvM, epsMrank)
  if (NthetaConvN >= int(0.8*Nteta) .and. NthetaConvM >= int(0.8*Nteta)) then
    print "(/,2x,'Convergence criteria for Nrank and Mrank are satisfied;')"                       
  else
    print "(/,2x,'Convergence criteria for Nrank and Mrank are not satisfied;')"
  end if 
  call write_InfoFileTmat (FileTmatG, Mrank, Nrank, .false., .false., .false.)
  call ScatCharact (ks, FileTmatG, Mrank, Nrank, .false., .false., .false.)  
  print "(/,2x,'T matrix is stored in ',a50)", FileTmatG
  print "(  2x,'The dimensions of the T matrix are given by:')"
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", Nrank                  
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank              
  deallocate (a, b, s, c, c1, t1, t, t11, t13, t31, t33, h, v, oldh, oldv,          &
              oldh0, oldv0)     
end subroutine convergence_Nrank_MrankINHOMSPH
!***********************************************************************************
subroutine TMatrix_Nrank_MrankINHOMSPH (Rb, ks, ind_ref, Mrank, Nrank, x1, y1, z1,  &
           alfa1, beta1, gama1, Mrank1, Nrank1, FileTmat, axsym1, chiral1,          &
           FileTmatG, PrnProgress)
  use parameters
  implicit none
  integer       :: Mrank, Nrank, Mrank1, Nrank1
  real(O)       :: ks, ind_ref, Rb, x1, y1, z1, alfa1, beta1, gama1                    
  character(80) :: FileTmat, FileTmatG
  logical       :: axsym1, chiral1, PrnProgress 
!
  integer       :: Nmax, Nmax1, NmaxAL, nt1g, mt1g
  real(O)       :: k, r1, sumang
  complex(O)    :: ind_refC
  complex(O),allocatable :: a(:,:), b(:,:), s(:,:), t1(:,:), t(:)                            
!  
  Nmax   = Nrank + Mrank * (2 * Nrank - Mrank + 1) 
  Nmax1  = Nrank1 + Mrank1 * (2 * Nrank1 - Mrank1 + 1) 
  NmaxAL = max(Nmax, Nmax1)
  k = ks * ind_ref  
  ind_refC = cmplx(ind_ref,0.0,O)  
  allocate (a(2*NmaxAL,2*NmaxAL), b(2*NmaxAL,2*NmaxAL), s(2*NmaxAL,2*NmaxAL))
  allocate (t(2*NmaxAL))
  open (unit = iTmat, file = FileTmat, status = "old", position = "rewind")
  if (.not. axsym1) then
    call read_HeadFileTmat (nt1g, mt1g)
    call check_dimensionMat (nt1g, mt1g, Nmax1)       
    allocate (t1(2*nt1g,2*mt1g))
    call read_FileTmat (nt1g, mt1g, t1) 
  else
    nt1g = Nmax1
    mt1g = Nmax1
    allocate (t1(2*nt1g,2*mt1g))
    call read_Tmatrix (chiral1, Nrank1, Mrank1, Nmax1, t1, nt1g, mt1g)            
  end if
  close (unit = iTmat)
  open  (unit = iTmat, file = FileTmatG, status = 'replace')
  call  write_HeadFileTmat (NmaxAL, NmaxAL)  
  if (PrnProgress) call  write_progress (.true., 1, 6)            
  r1 = sqrt(x1**2 + y1**2 + z1**2)
  sumang = abs(alfa1) + abs(beta1) + abs(gama1)
  if (r1 /=0 .and. sumang /= 0) then
    call MatTransRot_TR_mn_m1n1 (1, k, x1, y1, z1, alfa1, beta1, gama1, Mrank,      &
         Nrank, Nmax, Mrank1, Nrank1, Nmax1, s, NmaxAL, NmaxAL)
  else if (r1 /=0 .and. sumang == 0) then 
    call MatTransAB_mn_m1n1 (1, k, x1, y1, z1, Mrank, Nrank, Nmax, Mrank1,          &
         Nrank1, Nmax1, s, NmaxAL, NmaxAL)
  else if (r1 == 0 .and. sumang /= 0) then    
    call MatRot_R_mn_m1n1 (alfa1, beta1, gama1, Mrank, Nrank, Nmax, Mrank1,         &
         Nrank1, Nmax1, s, NmaxAL, NmaxAL, .true.)
  else if (r1 == 0 .and. sumang == 0) then
    call identity_transformation (Mrank, Nrank, Nmax, Mrank1, Nrank1, Nmax1, s,     &
         NmaxAL, NmaxAL)    
  end if  
  call product_matrices (2*Nmax, 2*Nmax1, 2*Nmax1, s, 2*NmaxAL, 2*NmaxAL, t1,       &
       2*nt1g, 2*mt1g)
  if (PrnProgress) call write_progress (.false., 2, 6)
  if (r1 /=0 .and. sumang /= 0) then
    call MatRotTrans_RT_mn_m1n1 (1, k, x1, y1, z1, alfa1, beta1, gama1, Mrank1,     &
         Nrank1, Nmax1, Mrank, Nrank, Nmax, a, NmaxAL, NmaxAL)
  else if (r1 /=0 .and. sumang == 0) then 
    call MatTransAB_mn_m1n1 (1, k, -x1, -y1, -z1, Mrank1, Nrank1, Nmax1, Mrank,     &
         Nrank, Nmax, a, NmaxAL, NmaxAL)
  else if (r1 == 0 .and. sumang /= 0) then    
    call MatRot_R_mn_m1n1 (-gama1, -beta1, -alfa1, Mrank1, Nrank1, Nmax1, Mrank,    &
         Nrank, Nmax, a, NmaxAL, NmaxAL, .true.)
  else if (r1 == 0 .and. sumang == 0) then
    call identity_transformation (Mrank1, Nrank1, Nmax1, Mrank, Nrank, Nmax,        &
         a, NmaxAL, NmaxAL)    
  end if   
  call product_matrices (2*Nmax, 2*Nmax1, 2*Nmax, s, 2*NmaxAL, 2*NmaxAL, a,         &
       2*NmaxAL, 2*NmaxAL) 
  if (PrnProgress) call write_progress (.false., 3, 6)      
  call vector_Q_sphere (3, 3, ks, Rb, ind_refC, Mrank, Nrank, Nmax, t)
  call product_matrix_vector1 (2*Nmax, 2*Nmax, t, s, 2*NmaxAL, 2*NmaxAL,            &
       a, 2*NmaxAL, 2*NmaxAL)
  call vector_Q_sphere (3, 1, ks, Rb, ind_refC, Mrank, Nrank, Nmax, t)
  call sum_diagonal_elements (2*Nmax, t, a, 2*NmaxAL, 2*NmaxAL) 
  if (PrnProgress) call write_progress (.false., 4, 6)
  call vector_Q_sphere (1, 3, ks, Rb, ind_refC, Mrank, Nrank, Nmax, t)
  call product_matrix_vector1 (2*Nmax, 2*Nmax, t, s, 2*NmaxAL, 2*NmaxAL,            &
       b, 2*NmaxAL, 2*NmaxAL)
  call vector_Q_sphere (1, 1, ks, Rb, ind_refC, Mrank, Nrank, Nmax, t)
  call sum_diagonal_elements (2*Nmax, t, b, 2*NmaxAL, 2*NmaxAL)  
  if (PrnProgress) call write_progress (.false., 5, 6)
  call LU_SYSTEM (a, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL, 2*Nmax)
  if (PrnProgress) call write_progress (.false., 6, 6)
  call minus_matrix (2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL)
  call write_FileTmat (NmaxAL, NmaxAL, b)         
  close (unit = iTmat) 
  call write_InfoFileTmat (FileTmatG, Mrank, Nrank, .false., .false., .false.)
  call ScatCharact (ks, FileTmatG, Mrank, Nrank, .false., .false., .false.) 
  print "(/,2x,'T matrix is stored in ',a50)", FileTmatG
  print "(  2x,'The dimensions of the T matrix are given by:')"         
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", Nrank          
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank              
  deallocate (a, b, s, t1, t)   
end subroutine TMatrix_Nrank_MrankINHOMSPH

