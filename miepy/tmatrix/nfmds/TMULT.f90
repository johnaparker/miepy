subroutine TMULT
!------------------------------------------------------------------------------------
! 1. General Considerations                                                         !
! --------------------------                                                        !
! TMULT is a routine for computing the T matrix and the scattering characteristics  !
! of clusters of arbitrarily shaped particles. As for inhomogeneous scatterers, the !
! individual T matrices of the particles are input parameters of the code (the      !
! names of the files containing the T matrices are supplied in the input file), and !
! they may correspond to homogeneous, axisymmetric or nonaxisymmetric particles and !
! inhomogeneous, composite or layered particles. The refractive indices of          !
! the ambient medium in the present T-matrix program and the programs computing     !
! individual T matrices must coincide. Note that the T matrix of a spherical        !
! particle must be computed with the routine TAXSYM and not with the program        !
! TSPHERE (because the T matrix of the particle must be specified as a              !
! two-dimensional array and the routine TSPHERE provides an one-dimensional array). ! 
!                                                                                   !
! The T matrix of a system of two particles 1 and 2 can be expressed in terms of    !
! the T matrices of each particle as                                                !
!                                                                                   !
!            T =  TR1(O,O1) * T1 * [I - TR3(O1,O2) * T2 * TR3(O2,O1) * T1 ]^(-1)    !
!               * [TR3(O1,O2) * T2 * TR1(O2,O) + TR1(O1,O) ]                        !
!               + TR1(O,O2) * T3 * [I - TR3(O2,O1) * T1 * TR3(O1,O2) * T2 ]^(-1)    ! 
!               * [TR3(O2,O1) * T1 * TR1(O1,O) + TR1(O2,O) ],                       !
!                                                                                   !
! where O is the origin of the cluster, O1 and O2 are the origins of the particles  !
! 1 and 2, respectively, T1 and T2 are the individual transition matrices, and      !
! TR(O1,O2) is the matrix transforming the vector spherical wave functions from     !
! the origin O1 to the origin O2 (through rotations and translations).              !
!                                                                                   !
! 2. Convergence Test                                                               !
! --------------------                                                              !
! In the code, the maximum expansion and azimuthal orders Nrank and Mrank of the    !
! cluster are independent on the maxium expansion and azimuthal orders              !
! NrankPart(i) and MrankPart(i), i = 1,...,Npart, of all particles forming the      !
! cluster. The convergence tests are carried out over Nrank and Mrank, while        !
! NrankPart(i) and MrankPart(i), i = 1,...,Npart, are kept constant.                !
!                                                                                   !
! For convergence tests, the incident wave is assumed to be a plane wave traveling  !
! along the Z-axis of the global coordinate system and the scattering               !
! characteristics are computed in the azimuthal plane  phi = 0°. The Euler          !
! orientation angles of the cluster are alphaC = betaC = 45° and gammaC = 0°, and   !
! interactive convergence tests over Nrank and Mrank are performed.                 !
!                                                                                   !  
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
! An estimate of Nrank is given by Wiscombe's truncation limit criterion            !
! [W. J. Wiscombe, Improved Mie scattering algorithms, Applied Optics, 19,          !
! 1505-1509, 1980], i.e.,                                                           !
!                                                                                   !
!                NrankW = int(x + 4.05 * x**0.33 + 2),                              !
!                                                                                   !
! where x is the size parameter, x = k * Rcirc, k is the wave number and Rcirc is   !
! the radius of the smallest sphere enclosing the cluster. For the maximum azimuthal! 
! order test, we may use the conservative estimate Mrank = Nrank - 2,..,Nrank.      !  
!                                                                                   ! 
! The convergence tests can be switched off by setting the logical variable         !
! DoConvTest to false. In this case, the values of Nrank and Mrank must be          !
! specified in the input file.                                                      !
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
! - Npart (integer) - number of particles.                                          !
!                                                                                   !
! - anorm (real) - characteristic length of the cluster which is used to normalize  !
!   the differential scattering cross sections.                                     !
!                                                                                   !
! - Rcirc (real) - characteristic length of the cluster (usually the radius of the  !
!   smallest circumscribing sphere) which is used to compute an estimate of the     !
!   maximum expansion order by using Wiscombe's truncation limit criterion (the size!
!   parameter is x = k * Rcirc, where k is the wave number in the ambient medium).  !
!   This parameter must be specified if the interactive convergence tests are       !
!   performed (DoConvTest = t).                                                     ! 
!                                                                                   !
! The next parameters (specified in two group statements) correspond to each        !
! particle. THE GROUP STATEMENTS MUST BE REPEATED FOR ALL Npart PARTICLES.          !
!                                                                                   !  
! - FileTmatPart (character(80)) - name of the file containing the T matrix of      !
!   the actual particle.                                                            !
!                                                                                   !
! - axsymPart (logical) - if axsymPart = t, the actual scatterer is a rotationally  !
!   symmetric particle (axisymmetric particle).                                     !
!                                                                                   !
! - chiralPart (logical) - if chiralPart = t, the actual scatterer is an optical    !
!   active particle (chiral particle).                                              !
!                                                                                   !
! - NrankPart, MrankPart (integer variables) - the maximum expansion and azimuthal  !
!   orders of the actual particle.                                                  !
!                                                                                   !   
! - xPart, yPart, zPart (real variables) - Cartesian coordinates specifying the     !
!   position of the actual particle in the cluster coordinate system.               !
!                                                                                   ! 
! - alphaPart, betaPart, gammaPart (real variables) - Euler angles specifying the   !
!   orientation of the coordinate system of the particle with respect to the        !
!   coordinate system of the cluster.                                               !
!                                                                                   !
! - DoConvTest (logical) - if DoConvTest = t, the interactive convergence tests     !
!   over Nrank and Mrank are performed. An estimates of Nrank is given by           !
!   Wiscombe's truncation limit criterion. If DoConvTest = f, the values of Nrank   !
!   and Mrank must be supplied in the input file.                                   !
!                                                                                   !
! - ExtThetaDom (logical) - if ExtThetaDom = t the DSCS is computed for             !
!   scattering angles ranging from 0 to 180° in the azimuthal plane phiGS = 0° and  !
!   from 180° to 0° in the azimuthal plane phiGS = 180°. The total number of        !
!   scattering angles is 10. If ExtThetaDom = f the DSCS is computed for            !
!   scattering angles ranging from 0° to 180° in the azimuthal plane phiGS = 0°.    !
!                                                                                   !
! - Nrank (integer) - maximum expansion order for the cluster. This parameter is    !
!   used if the convergence tests are not performed (DoConvTest = f).               !
!                                                                                   !
! - Mrank (integer) - maximum azimuthal order for the cluster. This parameter is    !
!   used if the convergence tests are not performed (DoConvTest = f).               !
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
!    - the name of the file containing the T matrix of each particle,               !
!    - the maximum expansion and azimuthal orders NrankPart and MrankPart,          ! 
!    - the Cartesian coordinates specifying the position of each particle, and      !
!    - the Euler angles specifying the orientation of each particle >               !
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
!       < the T matrix is stored in the file FileTmat and the T-matrix information  !
!         file is created >                                                         !
!       < Mrank and Nrank are written to the screen and to the T-matrix             !
!         information file >                                                        !
!       < go to Step 3 >                                                            !
!    end if                                                                         !
!    if ( .not. do convergence test ) then                                          !
!       < the T matrix is stored in the file FileTmat and the T-matrix information  !
!         file is created >                                                         !
!       < Mrank and Nrank are written to the screen and to the T-matrix             !
!         information file >                                                        !  
!       < go to Step 3 >                                                            !
!    end if                                                                         !
! 3. < scattering characteristics calculation >                                     !    
!------------------------------------------------------------------------------------
  use parameters
  use allocation, only: Mrankp, Nrankp, xp, yp, zp, alphap, betap, gammap,          &
                        FileTmatp, axsymp, chiralp  
  implicit none                                                                         
  integer       :: Npart, Mrank, Nrank                                   
  real(O)       :: k, ind_refMED, wavelength, anorm, snorm, epsNrank, epsMrank, Rcirc
  character(80) :: FileTmat
  logical       :: DoConvTest, ExtThetaDom, PrnProgress  
! -----------------------------------------------------------------------------------
!                               Read the input file                                 ! 
! -----------------------------------------------------------------------------------
  call readinputMULT ( wavelength, ind_refMed, Npart, anorm, Rcirc, DoConvTest,     &
       ExtThetaDom, Nrank, Mrank, epsNrank, epsMrank, FileTmat, PrnProgress,        &
       k, snorm ) 
! -----------------------------------------------------------------------------------
!                                      Main                                         !
! ----------------------------------------------------------------------------------- 
  if (DoConvTest) then            
    open (unit = iOutput,file = FileOutput,status = "replace")
    call printinputMULT (wavelength, ind_refMed, Npart, FileTmatp, axsymp, xp, yp,  &
         zp, alphap, betap, gammap, Mrankp, Nrankp, epsNrank, epsMrank, anorm, Rcirc)
    call convergence_Nrank_MrankMULT (k, snorm, Npart, FileTmatp, axsymp, chiralp,  &
         xp, yp, zp, alphap, betap, gammap, Mrankp, Nrankp, Mrank, Nrank, epsNrank, &
         epsMrank, ExtThetaDom, FileTmat, PrnProgress)
    close (unit = iOutput)
  else
    call TMatrix_Nrank_MrankMULT (k, Npart, FileTmatp, axsymp, chiralp, xp, yp, zp, &
         alphap, betap, gammap, Mrankp, Nrankp, Mrank, Nrank, FileTmat, PrnProgress)
  end if  
  deallocate (FileTmatp, axsymp, chiralp, Mrankp, Nrankp, xp, yp, zp, alphap,       &
              betap, gammap)  
end subroutine TMULT
!***********************************************************************************
subroutine readinputMULT ( wavelength, ind_refMed, Npart, anorm, Rcirc, DoConvTest, &
           ExtThetaDom, Nrank, Mrank, epsNrank, epsMrank, FileTmat, PrnProgress,    &
           k, snorm )    
  use parameters
  use derived_parameters
  use allocation, only: Mrankp, Nrankp, xp, yp, zp, alphap, betap, gammap,          &
                        FileTmatp, axsymp, chiralp
  implicit none                                                                         
  integer       :: Npart, MrankPart, NrankPart, ipart, Mrank, Nrank, NrankW, ios                                    
  real(O)       :: k, ind_refMED, wavelength, anorm, x, snorm, epsNrank, epsMrank,  &
                   xPart, yPart, zPart, alphaPart, betaPart, gammaPart, Rcirc, xR,  &
                   grd		   
  character(80) :: FileTmatPart, FileTmat, string
  logical       :: axsymPart, chiralPart, DoConvTest, ExtThetaDom, PrnProgress,     &
                   XFindPar  
! -----------------------------------------------------------------------------------
!                          Read the input file FileInputMULT                        ! 
! ----------------------------------------------------------------------------------- 
  call DrvParameters 
  open (unit = iInputMULT, file = FileInputMULT, status = "old", position = "rewind")   
  wavelength = 0.1_O * 2._O * Pi                                                                    
  ind_refMed = 1.5_O
  string     = 'OptProp'
  if (XFindPar (iInputMULT, string)) then
    read (iInputMULT, *, iostat = ios) wavelength
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable wavelength;')"
      stop
    end if
    read (iInputMULT, *, iostat = ios) ind_refMed
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
  Npart  = 2
  anorm  = 1._O  
  Rcirc  = 1._O 
  string = 'GenProp'
  if (XFindPar (iInputMULT, string)) then
    read (iInputMULT, *, iostat = ios) Npart
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Npart')"
      stop
    end if
    read (iInputMULT, *, iostat = ios) anorm
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable anorm;')"
      stop
    end if
    read (iInputMULT, *, iostat = ios) Rcirc
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
  allocate (FileTmatp(Npart), axsymp(Npart), chiralp(Npart), Mrankp(Npart),         &
            Nrankp(Npart), xp(Npart), yp(Npart), zp(Npart), alphap(Npart),          &
            betap(Npart), gammap(Npart))
  grd = Pi / 180._O	    
  do ipart = 1, Npart
    FileTmatPart = '../TMATFILES/T.dat'
    axsymPart  = .true.
    chiralPart = .false.     
    NrankPart  = 6
    MrankPart  = 3
    string     = 'TmatPart'
    if (XFindPar (iInputMULT, string)) then
      read (iInputMULT, *, iostat = ios) FileTmatPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable FileTmatPart')"
        print "(  2x,'for the particle ',i3,';')", ipart
        stop
      end if
      read (iInputMULT, *, iostat = ios) axsymPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable axsymPart;')"
        print "(  2x,'for the particle ',i3,';')", ipart
        stop
      end if
      read (iInputMULT, *, iostat = ios) chiralPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable chiralPart;')"
        print "(  2x,'for the particle ',i3,';')", ipart
        stop
      end if
      read (iInputMULT, *, iostat = ios) NrankPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable NrankPart;')"
        print "(  2x,'for the particle ',i3,';')", ipart 
        stop
      end if
      read (iInputMULT, *, iostat = ios) MrankPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable MrankPart;')"
        print "(  2x,'for the particle ',i3,';')", ipart 
        stop
      end if
    else
      print "(/,2x,'Group name TmatPart not found;')"
      stop  
    end if      
    call check_MrankNrank (MrankPart, NrankPart)    
    FileTmatp(ipart) = FileTmatPart
    axsymp(ipart)  = axsymPart
    chiralp(ipart) = chiralPart     
    Nrankp(ipart)  = NrankPart  
    Mrankp(ipart)  = MrankPart      
!
    xPart = 0.1_O
    yPart = 0.1_O
    zPart = 0.1_O
    alphaPart = 45._O
    betaPart  = 45._O
    gammaPart = 0._O 
    string    = 'GeomPartProp'
    if (XFindPar (iInputMULT, string)) then
      read (iInputMULT, *, iostat = ios) xPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable xPart')"
        print "(  2x,'for the particle ',i3,';')", ipart
        stop
      end if
      read (iInputMULT, *, iostat = ios) yPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable yPart;')"
        print "(  2x,'for the particle ',i3,';')", ipart
        stop
      end if
      read (iInputMULT, *, iostat = ios) zPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable zPart;')"
        print "(  2x,'for the particle ',i3,';')", ipart
        stop
      end if
      read (iInputMULT, *, iostat = ios) alphaPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable alphaPart;')"
        print "(  2x,'for the particle ',i3,';')", ipart 
        stop
      end if
      read (iInputMULT, *, iostat = ios) betaPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable betaPart;')"
        print "(  2x,'for the particle ',i3,';')", ipart 
        stop
      end if
      read (iInputMULT, *, iostat = ios) gammaPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable gammaPart;')"
        print "(  2x,'for the particle ',i3,';')", ipart 
        stop
      end if
    else
      print "(/,2x,'Group name GeomPartProp not found;')"
      stop  
    end if            
    xp(ipart) = xPart
    yp(ipart) = yPart
    zp(ipart) = zPart
    alphap(ipart) = alphaPart * grd
    betap(ipart)  = betaPart  * grd
    gammap(ipart) = gammaPart * grd
  end do
!
  DoConvTest  = .true.
  ExtThetaDom = .true.
  string      = 'ConvTest'
  if (XFindPar (iInputMULT, string)) then
    read (iInputMULT, *, iostat = ios) DoConvTest
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable DoConvTest')"
      stop
    end if
    read (iInputMULT, *, iostat = ios) ExtThetaDom
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
    print "(/,2x,'Convergence Test for a Cluster of Particles')"
    print "(  2x,'-------------------------------------------')"
  else
    print "(/,2x,'T-Matrix Computation for a Cluster of Particles')"
    print "(  2x,'-----------------------------------------------')"
  end if
!
  xR = k * Rcirc
  NrankW = int(xR + 4.05_O * xR**0.33_O + 2._O)
  if (.not. DoConvTest) then
    Nrank  =  16
    Mrank  =   8    
    string = 'NrankMrankCluster'
    if (XFindPar (iInputMULT, string)) then
      read (iInputMULT, *, iostat = ios) Nrank
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Nrank')"
        stop
      end if
      read (iInputMULT, *, iostat = ios) Mrank
      if (ios /= 0) then
       print "(/,2x,'Error by reading the input variable Mrank;')"
       stop
      end if            
    else
      print "(/,2x,'Group name NrankMrankCluster not found;')"
      stop  
    end if     
    print "(/,2x,'Input values:')"
    print "(  2x, a, i4, a, i4, a)",                                                &
   'the input values of Nrank and Mrank are ', Nrank, ' and', Mrank, ', while' 
    print "(  2x, a, i3, a)",                                                       &
   'the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';'                 
  else   
    print "(/,2x,'Nrank estimate:')"                                                    
    print "(  2x, a, i3, a)",                                                       &
   'the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';'
    print "(/,2x,'- enter the estimated values of Nrank and Mrank, where')"     
    print "(  2x,'  it is recommended to set Mrank = Nrank - 2,...,Nrank;')"        
    call read_integer2 (Nrank, Mrank)           
  end if 
  call check_MrankNrank (Mrank, Nrank)    
!
  epsNrank = 5.e-2_O
  epsMrank = 5.e-2_O
  string   = 'Errors'
  if (XFindPar (iInputMULT, string)) then
    read (iInputMULT, *, iostat = ios) epsNrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsNrank;')"
      stop
    end if
    read (iInputMULT, *, iostat = ios) epsMrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsMrank;')"
      stop
    end if            
  else
    print "(/,2x,'Group name Errors not found;')"
    stop  
  end if          
!
  FileTmat = '../TMATFILES/TG.dat'
  string   = 'Tmat'
  if (XFindPar (iInputMULT, string)) then
    read (iInputMULT, *, iostat = ios) FileTmat
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
  if (XFindPar (iInputMULT, string)) then
    read (iInputMULT, *, iostat = ios) PrnProgress
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable PrnProgress;')"
      stop
    end if             
  else
    print "(/,2x,'Group name PrintProgress not found;')"
    stop  
  end if        
  close (unit = iInputMULT)
end subroutine readinputMULT
!***********************************************************************************
subroutine printinputMULT (wavelength, ind_refMed, Npart, FileTmatp, axsymp, xp, yp,&
           zp, alfap, betap, gamap, Mrankp, Nrankp, epsNrank, epsMrank, anorm, Rcirc)
  use parameters
  implicit none  
  integer       :: Npart, i, Mrankp(Npart), Nrankp(Npart), LenString
  real(O)       :: wavelength, ind_refMed, xp(Npart), yp(Npart), zp(Npart),         &
                   alfap(Npart), betap(Npart), gamap(Npart), epsNrank, epsMrank,    &
                   anorm, Rcirc
  character(80) :: FileTmatp(Npart), FileTmatWrite              
  logical       :: axsymp(Npart)
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
  do i = 1, Npart    
    write (iOutput,"(2x,'particle: ',i3)") i
           FileTmatWrite = FileTmatp(i)(14:LenString(FileTmatp(i)))
    write (iOutput,"(2x,'name of the file containing the T matrix = ',a)")          &
           FileTmatWrite 
    if (axsymp(i)) then
      write (iOutput,"(2x,'axisymmetric particle;')")              
    else
      write (iOutput,"(2x,'nonaxisymmetric particle;')")             
    end if
    write (iOutput,"(2x,'maximum expansion order, Nrank = ',i3,';')") Nrankp(i)
    write (iOutput,"(2x,'maximum azimuthal order, Mrank = ',i3,';')") Mrankp(i)
    write (iOutput,"(2x, a, 1pe10.3, a, 1pe10.3, a, 1pe10.3, a)")                   &
   'Cartesian coordinates of the particle: x = ', xp(i), ', y = ', yp(i),           &
   ', z = ', zp(i), ';'          
    write (iOutput,"(2x, a, 1pe10.3, a, 1pe10.3, a, 1pe10.3, a)")                   &
   'Euler angles: alpha = ', alfap(i) * 180._O / Pi, ', beta = ',                   &
    betap(i) * 180._O / Pi, ', gamma = ', gamap(i) * 180._O / Pi, ';'  
    write (iOutput,*)
  end do
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum expansion order tolerance, epsNrank = ', epsNrank, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum azimuthal order tolerance, epsMrank = ', epsMrank, '.'            
  write (iOutput,"(/)")                       
end subroutine printinputMULT
!***********************************************************************************
subroutine convergence_Nrank_MrankMULT (k, snorm, Npart, FileTmatp, axsymp, chiralp,&
           xp, yp, zp, alfap, betap, gamap, Mrankp, Nrankp, Mrank, Nrank, epsNrank, &
           epsMrank, ExtThetaDom, FileTmat, PrnProgress)
  use parameters
  implicit none
  integer       :: Npart, Mrankp(Npart), Nrankp(Npart), Mrank, Nrank                   
  real(O)       :: k, snorm, xp(Npart), yp(Npart), zp(Npart), alfap(Npart),         &
                   betap(Npart), gamap(Npart), epsNrank, epsMrank
  character(80) :: FileTmatp(Npart), FileTmat         
  logical       :: axsymp(Npart), chiralp(Npart), ExtThetaDom, PrnProgress
!      
  integer       :: Nmax, Nmaxmax, NmaxAL, Nteta, ipart, jpart, NthetaConvN,         &
                   NthetaConvM, Nl, Nc, Mrankpl, Nrankpl, Nmaxpl, Mrankpl1,         &
                   Nrankpl1, Nmaxpl1, i, ntpg, mtpg  
  real(O)       :: alfa, beta, gama, alfaPol, x, y, z, tetaGI, phiGI, phiGS,        &
                   alfa1, beta1, gama1, x1, y1, z1, alfaC, betaC, gamaC, Cscat,     &
                   Qscat, Cext, Qext
  integer,allocatable    :: Nmaxp(:)
  real(O),allocatable    :: h(:), v(:), oldh(:), oldv(:), oldh0(:), oldv0(:)
  complex(O),allocatable :: aa(:,:), bb(:,:), a(:,:), b(:,:), t(:,:), c(:), c1(:)   
!
  tetaGI  = 0._O
  phiGI   = 0._O
  phiGS   = 0._O
  Nteta   = 10
  alfaC   = Pi / 4._O
  betaC   = Pi / 4._O
  gamaC   = 0._O
  alfaPol = Pi / 4._O
  Nmax = Nrank + Mrank * (2 * Nrank - Mrank + 1) 
  allocate (Nmaxp(Npart))
  NmaxAL  = Nmax
  Nmaxmax = 0
  do ipart = 1, Npart
    Nmaxp(ipart) = Nrankp(ipart) + Mrankp(ipart) * (2 * Nrankp(ipart) -             &
                   Mrankp(ipart) + 1)
    Nmaxmax      = Nmaxmax + Nmaxp(ipart)
    if (Nmaxp(ipart) > NmaxAL) NmaxAL = Nmaxp(ipart)
    if (Nmaxp(ipart) > Nmax) then
      print "(/,2x, a,i3,  a)",                                                     &
     'Warning: the input values of Nrank and Mrank are too low for particle ',ipart,';'        
    end if
  end do   
  call write_TypeConvHead (4) 
  allocate (aa(2*Nmaxmax,2*Nmaxmax), bb(2*Nmaxmax,2*NmaxAL))        
  allocate (a(2*NmaxAL,2*NmaxAL), b(2*NmaxAL,2*NmaxAL), c(2*NmaxAL), c1(2*NmaxAL))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta), oldh0(Nteta),             &
            oldv0(Nteta))       
  call identity_matrix (2*Nmaxmax, aa, 2*Nmaxmax, 2*Nmaxmax)
  if (PrnProgress) call write_progress (.true., 1, 2*Npart+2)            
  Nl = 0  
  do ipart = 1, Npart
    Mrankpl = Mrankp(ipart)
    Nrankpl = Nrankp(ipart)
    Nmaxpl  = Nmaxp(ipart)
    open (unit = iTmat, file = FileTmatp(ipart), status = "old", position = "rewind")
    if (.not. axsymp(ipart)) then
      call read_HeadFileTmat (ntpg, mtpg)
      call check_dimensionMat (ntpg, mtpg, Nmaxpl)            
      allocate (t(2*ntpg,2*mtpg))
      call read_FileTmat (ntpg, mtpg, t) 
    else
      ntpg = Nmaxpl
      mtpg = Nmaxpl
      allocate (t(2*ntpg,2*mtpg))
      call read_Tmatrix (chiralp(ipart), Nrankpl, Mrankpl, Nmaxpl, t, ntpg, mtpg)         
    end if
    close (unit = iTmat)                                                                                                       
    x = xp(ipart)
    y = yp(ipart)
    z = zp(ipart)
    alfa = alfap(ipart)
    beta = betap(ipart)
    gama = gamap(ipart)
    Nc   = 0      
    do jpart = 1, Npart            
      if (jpart /= ipart) then
        Mrankpl1 = Mrankp(jpart)
        Nrankpl1 = Nrankp(jpart)
        Nmaxpl1  = Nmaxp(jpart)
        x1 = xp(jpart)
        y1 = yp(jpart)
        z1 = zp(jpart)
        alfa1 = alfap(jpart)
        beta1 = betap(jpart)
        gama1 = gamap(jpart)            
        call MatRotTransRot_RTR_mn_m1n1 (3, k, x, y, z, alfa, beta, gama, Mrankpl,  &
             Nrankpl, Nmaxpl, x1, y1, z1, alfa1, beta1, gama1, Mrankpl1, Nrankpl1,  &
             Nmaxpl1, a, NmaxAL, NmaxAL)        
        call product_matrices1 (2*Nmaxpl, 2*Nmaxpl, 2*Nmaxpl1, t, 2*ntpg, 2*mtpg,   &
             a, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL)                      
        call extend_matrix1 (ipart, jpart, Npart, Nmaxp, Nl, Nc, b, NmaxAL, NmaxAL, &
             aa, Nmaxmax, Nmaxmax)                                                      
      end if
      Nc = Nc + 2 * Nmaxp(jpart)      
    end do                 
    call MatRotTrans_RT_mn_m1n1 (1, k, x, y, z, alfa, beta, gama, Mrankpl, Nrankpl, &
         Nmaxpl, Mrank, Nrank, Nmax, a, NmaxAL, NmaxAL)
    call product_matrices1 (2*Nmaxpl, 2*Nmaxpl, 2*Nmax, t, 2*ntpg, 2*mtpg,          &
         a, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL)              
    call extend_matrix2 (ipart, Npart, Nmaxp, Nmax, Nmaxmax, Nl, b, NmaxAL, NmaxAL, &
         bb, Nmaxmax, NmaxAL)            
    Nl = Nl + 2 * Nmaxp(ipart)
    if (PrnProgress) call write_progress (.false., ipart+1, 2*Npart+2)
    deallocate (t)
  end do
  call LU_SYSTEM_DIRECT (aa, 2*Nmaxmax, 2*Nmaxmax, bb, 2*Nmaxmax, 2*NmaxAL,         &
       2*Nmaxmax, 2*Nmax)
  if (PrnProgress) call write_progress (.false., Npart+2, 2*Npart+2)
  open (unit = iTmat, file = FileTmat, status = 'replace')
  call  write_HeadFileTmat (NmaxAL, NmaxAL) 
  allocate (t(2*NmaxAL,2*NmaxAL))
  Nl = 0
  do ipart = 1, Npart
    Mrankpl = Mrankp(ipart)
    Nrankpl = Nrankp(ipart)
    Nmaxpl  = Nmaxp(ipart)          
    x = xp(ipart)
    y = yp(ipart)
    z = zp(ipart)
    alfa = alfap(ipart)
    beta = betap(ipart)
    gama = gamap(ipart)
    call extract_matrix1 (ipart, Npart, Nmaxp, Nmax, Nl, b, NmaxAL, NmaxAL,         &
         bb, Nmaxmax, NmaxAL)       
    call MatTransRot_TR_mn_m1n1 (1, k, x, y, z, alfa, beta, gama, Mrank, Nrank,     &
         Nmax, Mrankpl, Nrankpl, Nmaxpl, t, NmaxAL, NmaxAL)    
    call product_matrices2 (ipart, 2*Nmax, 2*Nmaxpl, 2*Nmax, t, 2*NmaxAL, 2*NmaxAL, &
         b, 2*NmaxAL, 2*NmaxAL, a, 2*NmaxAL, 2*NmaxAL)                      
    Nl = Nl + 2 * Nmaxpl
    if (PrnProgress) call write_progress (.false., Npart+2+ipart, 2*Npart+2)
  end do
  deallocate (t)  
  call copy_matrix (2*Nmax, 2*Nmax, a, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL)
  call write_FileTmat (NmaxAL, NmaxAL, a)
  call PWcoefficients_ab (tetaGI, phiGI, alfaC, betaC, gamaC, alfaPol, Mrank,       &
       Nrank, Nmax, c)    
  call product_matrix_vector (2*Nmax, 2*Nmax, a, 2*NmaxAL, 2*NmaxAL, c, c1)
  call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfaC, betaC, gamaC, k, snorm,   & 
       ExtThetaDom,.true., h, v)
  call CQscat (c1, Mrank, Nrank, Nmax, k, snorm, Cscat, Qscat)
  call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfaC, betaC, gamaC,           &
       alfaPol, k, snorm, Cext, Qext)        
  call write_2ConvParam (Nrank, Mrank)  
  call write_DSCS (Nteta,ExtThetaDom, h, v)        
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
  call matrix_Nrank_1_right (Mrank, Nrank, Nmax, a, NmaxAL, NmaxAL)  
  call PWcoefficients_ab (tetaGI, phiGI, alfaC, betaC, gamaC, alfaPol, Mrank,       &
       Nrank, Nmax, c)    
  call product_matrix_vector (2*Nmax, 2*Nmax, a, 2*NmaxAL, 2*NmaxAL, c, c1)
  call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfaC, betaC, gamaC, k, snorm,   &
       ExtThetaDom,.true., h, v)
  call delta_DSCS (Nteta, h, v, oldh, oldv, epsNrank, NthetaConvN)
  call CQscat (c1, Mrank, Nrank, Nmax, k, snorm, Cscat, Qscat)
  call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfaC, betaC, gamaC,           &
       alfaPol, k, snorm, Cext, Qext) 
  call write_2ConvParam (Nrank - 1, Mrank)
  call write_DSCS (Nteta, ExtThetaDom, h, v) 
  call write_Effic (Qscat, Qext)
  call write_NrankConvRes (NthetaConvN, Nteta, epsNrank)
! --- (Mrank - 1) configuration ---  
  call matrix_Mrank_1_right (Mrank, Nrank, Nmax, b, NmaxAL, NmaxAL)  
  call PWcoefficients_ab (tetaGI, phiGI, alfaC, betaC, gamaC, alfaPol, Mrank,       &
       Nrank, Nmax, c)    
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL, c, c1)
  call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfaC, betaC, gamaC, k, snorm,   &
       ExtThetaDom,.true., h, v)
  call delta_DSCS (Nteta, h, v, oldh0, oldv0, epsMrank, NthetaConvM)
  call CQscat (c1, Mrank, Nrank, Nmax, k, snorm, Cscat, Qscat)
  call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfaC, betaC, gamaC,           &
       alfaPol, k, snorm, Cext, Qext) 
  call write_2ConvParam (Nrank, Mrank - 1)
  call write_DSCS (Nteta, ExtThetaDom, h, v) 
  call write_Effic (Qscat, Qext)
  call write_MrankConvRes (NthetaConvM, epsMrank)
  if (NthetaConvN >= int(0.8*Nteta) .and. NthetaConvM >= int(0.8*Nteta)) then
    print "(/,2x,'Convergence criteria for Nrank and Mrank are satisfied;')"                                    
  else
    print "(/,2x,'Convergence criteria for Nrank and Mrank are not satisfied;')"
  end if  
  call write_InfoFileTmat (FileTmat, Mrank, Nrank, .false., .false., .false.)
  call ScatCharact (k, FileTmat, Mrank, Nrank, .false., .false., .false.)
  print "(/,2x,'T matrix is stored in ',a50)", FileTmat
  print "(  2x,'The dimensions of the T matrix are given by:')"           
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", Nrank
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank  
  deallocate (aa, bb, a, b, c, c1, h, v, oldh, oldv, oldh0, oldv0, Nmaxp)      
end subroutine convergence_Nrank_MrankMULT
!***********************************************************************************
subroutine TMatrix_Nrank_MrankMULT (k, Npart, FileTmatp, axsymp, chiralp, xp, yp,   &
           zp, alfap, betap, gamap, Mrankp, Nrankp, Mrank, Nrank, FileTmat,         &
           PrnProgress)
  use parameters
  implicit none
  integer       :: Npart, Mrankp(Npart), Nrankp(Npart), Mrank, Nrank                   
  real(O)       :: k, xp(Npart), yp(Npart), zp(Npart), alfap(Npart), betap(Npart),  &
                   gamap(Npart)
  character(80) :: FileTmatp(Npart), FileTmat           
  logical       :: axsymp(Npart), chiralp(Npart), PrnProgress
!      
  integer       :: Nmax, Nmaxmax, NmaxAL, ipart, jpart, Nl, Nc, Mrankpl, Nrankpl,   &
                   Nmaxpl, Mrankpl1, Nrankpl1, Nmaxpl1, ntpg, mtpg  
  real(O)       :: alfa, beta, gama, x, y, z, alfa1, beta1, gama1, x1, y1, z1                              
  integer,allocatable    :: Nmaxp(:)
  complex(O),allocatable :: aa(:,:), bb(:,:), a(:,:), b(:,:), t(:,:)   
!  
  Nmax = Nrank + Mrank * (2 * Nrank - Mrank + 1) 
  allocate (Nmaxp(Npart))
  NmaxAL  = Nmax
  Nmaxmax = 0
  do ipart = 1, Npart
    Nmaxp(ipart) = Nrankp(ipart) + Mrankp(ipart) * (2 * Nrankp(ipart) -             &
                   Mrankp(ipart) + 1)
    Nmaxmax      = Nmaxmax + Nmaxp(ipart)
    if (Nmaxp(ipart) > NmaxAL) NmaxAL = Nmaxp(ipart)
    if (Nmaxp(ipart) > Nmax) then
      print "(/,2x, a, i3, a)",                                                     &
     'Warning: the input values of Nrank and Mrank are too low for particle ',ipart,';'        
    end if
  end do   
  allocate (aa(2*Nmaxmax,2*Nmaxmax), bb(2*Nmaxmax,2*NmaxAL))        
  allocate (a(2*NmaxAL,2*NmaxAL), b(2*NmaxAL,2*NmaxAL))         
  call identity_matrix (2*Nmaxmax, aa, 2*Nmaxmax, 2*Nmaxmax)
  if (PrnProgress) call write_progress (.true., 1, 2*Npart+2)            
  Nl = 0  
  do ipart = 1, Npart
    Mrankpl = Mrankp(ipart)
    Nrankpl = Nrankp(ipart)
    Nmaxpl  = Nmaxp(ipart)
    open (unit = iTmat, file = FileTmatp(ipart), status = "old", position = "rewind")
    if (.not. axsymp(ipart)) then
      call read_HeadFileTmat (ntpg, mtpg)
      call check_dimensionMat (ntpg, mtpg, Nmaxpl)            
      allocate (t(2*ntpg,2*mtpg))
      call read_FileTmat (ntpg, mtpg, t) 
    else
      ntpg = Nmaxpl
      mtpg = Nmaxpl
      allocate (t(2*ntpg,2*mtpg))
      call read_Tmatrix (chiralp(ipart), Nrankpl, Mrankpl, Nmaxpl, t, ntpg, mtpg)         
    end if
    close (unit = iTmat)                                                                                                       
    x = xp(ipart)
    y = yp(ipart)
    z = zp(ipart)
    alfa = alfap(ipart)
    beta = betap(ipart)
    gama = gamap(ipart)
    Nc   = 0    
    do jpart = 1, Npart     
      if (jpart /= ipart) then
        Mrankpl1 = Mrankp(jpart)
        Nrankpl1 = Nrankp(jpart)
        Nmaxpl1  = Nmaxp(jpart)
        x1 = xp(jpart)
        y1 = yp(jpart)
        z1 = zp(jpart)
        alfa1 = alfap(jpart)
        beta1 = betap(jpart)
        gama1 = gamap(jpart)            
        call MatRotTransRot_RTR_mn_m1n1 (3, k, x, y, z, alfa, beta, gama, Mrankpl,  &
             Nrankpl, Nmaxpl, x1, y1, z1, alfa1, beta1, gama1, Mrankpl1, Nrankpl1,  &
             Nmaxpl1, a, NmaxAL, NmaxAL)        
        call product_matrices1 (2*Nmaxpl, 2*Nmaxpl, 2*Nmaxpl1, t, 2*ntpg, 2*mtpg,   &
             a, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL)                      
        call extend_matrix1 (ipart, jpart, Npart, Nmaxp, Nl, Nc, b, NmaxAL, NmaxAL, &
             aa, Nmaxmax, Nmaxmax)                                                      
      end if
      Nc = Nc + 2 * Nmaxp(jpart)      
    end do              
    call MatRotTrans_RT_mn_m1n1 (1, k, x, y, z, alfa, beta, gama, Mrankpl, Nrankpl, &
         Nmaxpl, Mrank, Nrank, Nmax, a, NmaxAL, NmaxAL)
    call product_matrices1 (2*Nmaxpl, 2*Nmaxpl, 2*Nmax, t, 2*ntpg, 2*mtpg,          &
         a, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL)              
    call extend_matrix2 (ipart, Npart, Nmaxp, Nmax, Nmaxmax, Nl, b, NmaxAL, NmaxAL, &
         bb, Nmaxmax, NmaxAL)            
    Nl = Nl + 2 * Nmaxp(ipart)
    if (PrnProgress) call write_progress (.false., ipart+1, 2*Npart+2)
    deallocate (t)
  end do
  call LU_SYSTEM_DIRECT (aa, 2*Nmaxmax, 2*Nmaxmax, bb, 2*Nmaxmax, 2*NmaxAL,         &
       2*Nmaxmax, 2*Nmax)
  if (PrnProgress) call write_progress (.false., Npart+2, 2*Npart+2)
  open (unit = iTmat, file = FileTmat, status = 'replace')
  call  write_HeadFileTmat (NmaxAL, NmaxAL) 
  allocate (t(2*NmaxAL,2*NmaxAL))
  Nl = 0
  do ipart = 1, Npart
    Mrankpl = Mrankp(ipart)
    Nrankpl = Nrankp(ipart)
    Nmaxpl  = Nmaxp(ipart)        
    x = xp(ipart)
    y = yp(ipart)
    z = zp(ipart)
    alfa = alfap(ipart)
    beta = betap(ipart)
    gama = gamap(ipart)
    call extract_matrix1 (ipart, Npart, Nmaxp, Nmax, Nl, b, NmaxAL, NmaxAL,         &
         bb, Nmaxmax, NmaxAL)       
    call MatTransRot_TR_mn_m1n1 (1, k, x, y, z, alfa, beta, gama, Mrank, Nrank,     &
         Nmax, Mrankpl, Nrankpl, Nmaxpl, t, NmaxAL, NmaxAL)    
    call product_matrices2 (ipart, 2*Nmax, 2*Nmaxpl, 2*Nmax, t, 2*NmaxAL, 2*NmaxAL, &
         b, 2*NmaxAL, 2*NmaxAL, a, 2*NmaxAL, 2*NmaxAL)                      
    Nl = Nl + 2 * Nmaxpl
    if (PrnProgress) call write_progress (.false., Npart+2+ipart, 2*Npart+2)
  end do  
  call write_FileTmat (NmaxAL, NmaxAL, a)  
  close (unit = iTmat)   
  call write_InfoFileTmat (FileTmat, Mrank, Nrank, .false., .false., .false.) 
  call ScatCharact (k, FileTmat, Mrank, Nrank, .false., .false., .false.)
  print "(/,2x,'T matrix is stored in ',a50)", FileTmat
  print "(  2x,'The dimensions of the T matrix are given by:')"           
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", Nrank
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank  
  deallocate (aa, bb, a, b, t, Nmaxp)   
end subroutine TMatrix_Nrank_MrankMULT

