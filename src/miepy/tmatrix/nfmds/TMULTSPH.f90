subroutine TMULTSPH
!------------------------------------------------------------------------------------
! 1. General Considerations                                                         !
! --------------------------                                                        !
! TMULTSPH is a routine for computing the T matrix and the scattering               !
! characteristics of clusters of homogeneous, dielectric spheres. The particles in  !
! the cluster can have the same radius or different radii.                          !
!                                                                                   !
! 2. Convergence Test                                                               !
! --------------------                                                              !
! In the code, the maximum expansion and azimuthal orders Nrank and Mrank of the    !
! cluster are independent on the maxium expansion and azimuthal orders              !
! NrankPart(i) and MrankPart(i), i = 1,...,Npart, of all spherical particles        !
! forming the cluster. The convergence tests are carried out over Nrank and Mrank,  !
! while NrankPart(i) and MrankPart(i), i = 1,...,Npart, are kept constant.          !
!                                                                                   !
! For convergence tests, the incident wave is assumed to be a plane wave traveling  !
! along the Z-axis of the global coordinate system and the scattering               !
! characteristics are computed in the azimuthal plane phi = 0°. The Euler           !
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
! the radius of the smallest sphere enclosing the cluster. For the azimuthal order  !
! test, Mrank can be chosen as Nrank - 2,..,Nrank. The same prescriptions hold for  !
! NrankPart(i) and MrankPart(i) for all i = 1,2,...,Npart.                          !
!                                                                                   ! 
! The convergence tests can be switched off by setting the logical variable         !
! DoConvTest to false. In this case, the values of Nrank and Mrank must be          !
! specified in the input file.                                                      !
!                                                                                   !
! 3. Input Parameters                                                               !
! --------------------                                                              !
! The parameters specified in the input file "/INPUTFILES/InputMULTSPH.dat" are     !
! listed below.                                                                     !
!                                                                                   !
! - wavelength (real) - wavelength of the incident light in vacuo.                  !
!                                                                                   !
! - ind_refMed (real) - refractive index of the ambient medium (nonabsorbing        !
!   dielectric medium).                                                             !
!                                                                                   !
! - Npart (integer) - number of spherical particles.                                !
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
! - identical (logical) - if identical = true, the spherical particles have the     !
!   same radius. Otherwise, the radii are different.                                !
!                                                                                   !
! The next parameters (specified in two group statements) correspond to each        !
! spherical particle. THE GROUP STATEMENTS MUST BE REPEATED FOR ALL Npart           !
! PARTICLES. FOR IDENTICAL PARTICLES, USES AFTER THE TmatPart KEYWORD, A BLANK      !
! LINE.                                                                             ! 
!                                                                                   !
! - rPart (real) - radius of the actual spherical particle.                         !
!                                                                                   !
! - NrankPart, MrankPart (integer variables) - maximum expansion and azimuthal      !
!   orders for the actual spherical particle.                                       !
!                                                                                   !
! - ind_refPart (complex) - relative refractive index of the actual spherical       !
!   particle with respect to the ambient medium.                                    !
!                                                                                   !
! - xPart, yPart, zPart (real variables) - Cartesian coordinates specifying the     !
!   position of the actual spherical particle with respect to the coordinate        !
!   system of the cluster.                                                          !
!                                                                                   !
! - DoConvTest (logical) - if DoConvTest = t, the interactive convergence tests     !
!   over Nrank and Mrank are performed. An estimates of Nrank is given by           !
!   Wiscombe's truncation limit criterion. If DoConvTest = f, the values of Nrank   !
!   and Mrank must be supplied in the input file.                                   !
!                                                                                   !
! - ExtThetaDom (logical) - if ExtThetaDom = t the DSCS is computed for             !
!   scattering angles ranging from 0° to 180° in the azimuthal plane phiGS = 0°,    !
!   and from 180° to 0° in the azimuthal plane phiGS = 180°. The total number of    !
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
!    - the radius of each spherical particle,                                       !
!    - the maximum expansion and azimuthal orders NrankPart and MrankPart, and      ! 
!    - the Cartesian coordinates specifying the position of each spherical          !
!      particle >                                                                   ! 
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
!        (Nrank, Mrank - 1) and write the results to the file "Output.dat" >        !
!       < the T matrix is stored in the file FileTmat and the T-matrix information  !
!         file is created >                                                         !
!       < Mrank and Nrank are written to the screen and to the T-matrix             !
!         information file >                                                        ! 
!       < go to Step 3 >                                                            !
!    end if                                                                         !
!    if ( .not. convergence test ) then                                             !
!       < the T matrix is stored in the file FileTmat and the T-matrix information  !
!         file is created >                                                         !
!       < Mrank and Nrank are written to the screen and to the T-matrix             !
!         information file >                                                        ! 
!       < go to Step 3 >                                                            !
!    end if                                                                         !
! 3. < scattering characteristics calculation >                                     !   
!------------------------------------------------------------------------------------
  use parameters
  use allocation, only: Mrankp, Nrankp, rp, xp, yp, zp, ind_refp
  implicit none 
  integer       :: Npart, Mrank, Nrank
  real(O)       :: k, ind_refMED, wavelength, anorm, snorm, epsNrank, epsMrank, Rcirc
  character(80) :: FileTmat
  logical       :: identical, DoConvTest, ExtThetaDom, PrnProgress  
! -----------------------------------------------------------------------------------
!                      Read the input file FileInputMULTSPH                         ! 
! -----------------------------------------------------------------------------------
  call readinputMULTSPH ( wavelength, ind_refMed, Npart, anorm, Rcirc,              &
       identical, DoConvTest, ExtThetaDom, Nrank, Mrank, epsNrank, epsMrank,        &
       FileTmat, PrnProgress, k, snorm)  
! -----------------------------------------------------------------------------------
!                                      Main                                         !
! ----------------------------------------------------------------------------------- 
  if (DoConvTest) then     
    open (unit = iOutput, file = FileOutput, status = "replace")
    call printinputMULTSPH (wavelength, ind_refMed, Npart, identical, rp, xp, yp,   &
         zp, ind_refp, Mrankp, Nrankp, epsNrank, epsMrank, anorm, Rcirc)
    call convergence_Nrank_MrankMULTSPH (k, snorm, Npart, identical, rp, xp, yp,    &
         zp, ind_refp, Mrankp, Nrankp, Mrank, Nrank, epsNrank, epsMrank,            &
         ExtThetaDom, FileTmat, PrnProgress)
    close (unit = iOutput)
  else
    call TMatrix_Nrank_MrankMULTSPH (k, Npart, identical, rp, xp, yp, zp, ind_refp, &
         Mrankp, Nrankp, Mrank, Nrank, FileTmat, PrnProgress)
  end if
  deallocate (Mrankp, Nrankp, rp, xp, yp, zp, ind_refp)  
end subroutine TMULTSPH
!***********************************************************************************
subroutine readinputMULTSPH ( wavelength, ind_refMed, Npart, anorm, Rcirc,          &
           identical, DoConvTest, ExtThetaDom, Nrank, Mrank, epsNrank, epsMrank,    &
           FileTmat, PrnProgress, k, snorm)   
  use parameters
  use derived_parameters
  use allocation, only: Mrankp, Nrankp, rp, xp, yp, zp, ind_refp
  implicit none 
  integer       :: Npart, MrankPart, NrankPart, ipart, Mrank, Nrank, ios, NrankW
  real(O)       :: k, ind_refMED, wavelength, anorm, x, snorm, epsNrank, epsMrank,  &
                   rPart, xPart, yPart, zPart, Rcirc, xR
  complex(O)    :: ind_refPart
  character(80) :: FileTmat, string
  logical       :: identical, DoConvTest, ExtThetaDom, PrnProgress, XFindPar
! -----------------------------------------------------------------------------------
!                         Read the input file FileInputMULTSPH                      ! 
! ----------------------------------------------------------------------------------- 
  call DrvParameters 
  open (unit = iInputMULTSPH, file = FileInputMULTSPH, status = "old",              &
        position = "rewind")   
  wavelength = 0.1_O * 2._O * Pi                                                                    
  ind_refMed = 1.5_O
  string     = 'OptProp'
  if (XFindPar (iInputMULTSPH, string)) then
    read (iInputMULTSPH, *, iostat = ios) wavelength
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable wavelength;')"
      stop
    end if
    read (iInputMULTSPH, *, iostat = ios) ind_refMed
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
  Npart = 2
  anorm = 1._O   
  Rcirc = 1._O
  identical = .false.
  string    = 'GenProp'
  if (XFindPar (iInputMULTSPH, string)) then
    read (iInputMULTSPH, *, iostat = ios) Npart
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Npart')"
      stop
    end if
    read (iInputMULTSPH, *, iostat = ios) anorm
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable anorm;')"
      stop
    end if
    read (iInputMULTSPH, *, iostat = ios) Rcirc
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Rcirc;')"
      stop
    end if 
    read (iInputMULTSPH, *, iostat = ios) identical
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable identical;')"
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
  allocate (Mrankp(Npart), Nrankp(Npart), rp(Npart), xp(Npart), yp(Npart),          &
            zp(Npart), ind_refp(Npart)) 
  do ipart = 1, Npart    
    rPart = 1.0_O    
    NrankPart   = 6
    MrankPart   = 4
    ind_refPart = (1.5_O,0._O) 
    string     = 'TmatPart'
    if (.not. identical .or. ipart == 1) then
      if (XFindPar (iInputMULTSPH, string)) then
        read (iInputMULTSPH, *, iostat = ios) rPart
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable rPart')"
          print "(  2x,'for the particle ',i3,';')", ipart
          stop
        end if      
        read (iInputMULTSPH, *, iostat = ios) NrankPart
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable NrankPart;')"
          print "(  2x,'for the particle ',i3,';')", ipart 
          stop
        end if
        read (iInputMULTSPH, *, iostat = ios) MrankPart
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable MrankPart;')"
          print "(  2x,'for the particle ',i3,';')", ipart 
          stop
        end if
        read (iInputMULTSPH, *, iostat = ios) ind_refPart
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable ind_refPart;')"
          print "(  2x,'for the particle ',i3,';')", ipart 
          stop
        end if
      else
        print "(/,2x,'Group name TmatPart not found;')"
        stop  
      end if
      call check_MrankNrank (MrankPart, NrankPart)  
      rp(ipart) = rPart    
      Nrankp(ipart)   = NrankPart
      Mrankp(ipart)   = MrankPart
      ind_refp(ipart) = ind_refPart
    else
      rp(ipart) = rp(ipart-1)     
      Nrankp(ipart)   = Nrankp(ipart-1)
      Mrankp(ipart)   = Mrankp(ipart-1)
      ind_refp(ipart) = ind_refp(ipart-1)
    end if
!
    xPart = 0.1_O
    yPart = 0.1_O
    zPart = 0.1_O  
    string    = 'GeomPartProp'
    if (XFindPar (iInputMULTSPH, string)) then
      read (iInputMULTSPH, *, iostat = ios) xPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable xPart')"
        print "(  2x,'for the particle ',i3,';')", ipart
        stop
      end if
      read (iInputMULTSPH, *, iostat = ios) yPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable yPart;')"
        print "(  2x,'for the particle ',i3,';')", ipart
        stop
      end if
      read (iInputMULTSPH, *, iostat = ios) zPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable zPart;')"
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
  end do  
!
  DoConvTest  = .true.
  ExtThetaDom = .true.
  string      = 'ConvTest'
  if (XFindPar (iInputMULTSPH, string)) then
    read (iInputMULTSPH, *, iostat = ios) DoConvTest
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable DoConvTest')"
      stop
    end if
    read (iInputMULTSPH, *, iostat = ios) ExtThetaDom
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
    print "(/,2x,'Convergence Test for a Cluster of Spheres')"
    print "(  2x,'-----------------------------------------')"
  else
    print "(/,2x,'T-Matrix Computation for a Cluster of Spheres')"
    print "(  2x,'---------------------------------------------')"
  end if
!
  xR = k * Rcirc
  NrankW = int(xR + 4.05_O * xR**0.33_O + 2._O)
  if (.not. DoConvTest) then
    Nrank  =  16
    Mrank  =   8
    string = 'NrankMrankCluster'
    if (XFindPar (iInputMULTSPH, string)) then
      read (iInputMULTSPH, *, iostat = ios) Nrank
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Nrank')"
        stop
      end if
      read (iInputMULTSPH, *, iostat = ios) Mrank
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
  if (XFindPar (iInputMULTSPH, string)) then
    read (iInputMULTSPH, *, iostat = ios) epsNrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsNrank;')"
      stop
    end if
    read (iInputMULTSPH, *, iostat = ios) epsMrank
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
  if (XFindPar (iInputMULTSPH, string)) then
    read (iInputMULTSPH, *, iostat = ios) FileTmat
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
  if (XFindPar (iInputMULTSPH, string)) then
    read (iInputMULTSPH, *, iostat = ios) PrnProgress
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable PrnProgress;')"
      stop
    end if             
  else
    print "(/,2x,'Group name PrintProgress not found;')"
    stop  
  end if         
  close (unit = iInputMULTSPH)  
end subroutine readinputMULTSPH
!***********************************************************************************
subroutine printinputMULTSPH (wavelength, ind_refMed, Npart, identical, rp, xp, yp, &
           zp, ind_refp, Mrankp, Nrankp, epsNrank, epsMrank, anorm, Rcirc)
  use parameters
  implicit none  
  integer       :: Npart, i, Mrankp(Npart), Nrankp(Npart)
  real(O)       :: wavelength, ind_refMed, xp(Npart), yp(Npart), zp(Npart),         &
                   epsNrank, epsMrank, rp(Npart), anorm, Rcirc
  complex(O)    :: ind_refp(Npart)
  logical       :: identical  
!
  write (iOutput,"(/,2x,'Input Parameters of the Scattering Problem:',/)")
  write (iOutput,"(2x,'wavelength of the free space, wavelength = ',1pe13.4,';')")  &
         wavelength
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'refractive index of the ambient medium, ind_refMed = ', ind_refMed, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'characteristic length of the cluster, anorm = ', anorm, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'radius of the circumscribing sphere, Rcirc =', Rcirc, ';'
  write (iOutput,*)
  if (identical) then
    write (iOutput,"(2x,'identical particles;')")
    write (iOutput,"(2x,'radius: r = ',1pe10.3,';')") rp(1)
    write (iOutput,"(2x, a, 1pe10.3,',', 1pe10.3, a)")                              &
   'relative refractive index, ind_refRel = (', ind_refp(1), ');'  
    write (iOutput,"(2x,'number of particles: Npart = ',i3,';')") Npart 
    write (iOutput,"(2x,'maximum expansion order, Nrank = ',i3,';')") Nrankp(1)
    write (iOutput,"(2x,'maximum azimuthal order, Mrank = ',i3,';')") Mrankp(1)      
    do i = 1, Npart        
      write (iOutput,"(2x, a, i3, a, 1pe10.3, a, 1pe10.3, a, 1pe10.3, a)")          &
     'Cartesian coordinates of the particle ', i, ': x = ', xp(i), ', y = ', yp(i), &
     ', z = ', zp(i), ';'                   
    end do
    write (iOutput,*)       
  else  
    write (iOutput,"(2x,'number of particles: Npart = ',i3,';')") Npart
    write (iOutput,*)
    do i = 1, Npart  
      write (iOutput,"(2x,'particle: ',i3)") i    
      write (iOutput,"(2x,'radius, r = ',1pe10.3,';')") rp(i)
      write (iOutput,"(2x, a, 1pe10.3, ',', 1pe10.3, a)")                           &
     'relative refractive index, ind_refRel = (', ind_refp(i), ');'  
      write (iOutput,"(2x,'maximum expansion order, Nrank = ',i3,';')") Nrankp(i)
      write (iOutput,"(2x,'maximum azimuthal order, Mrank = ',i3,';')") Mrankp(i)       
      write (iOutput,"(2x, a, i3, a, 1pe10.3, a, 1pe10.3, a, 1pe10.3, a)")          &
     'Cartesian coordinates of the particle ', i, ': x = ', xp(i), ', y = ', yp(i), &
     ', z = ', zp(i), ';'                             
      write (iOutput,*) 
    end do
  end if  
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum expansion order tolerance, epsNrank = ', epsNrank, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum azimuthal order tolerance, epsMrank = ', epsMrank, '.'            
  write (iOutput,"(/)")                       
end subroutine printinputMULTSPH
!***********************************************************************************
subroutine convergence_Nrank_MrankMULTSPH (k, snorm, Npart, identical, rp, xp, yp,  &
           zp, ind_refp, Mrankp, Nrankp, Mrank, Nrank, epsNrank, epsMrank,          &
           ExtThetaDom, FileTmat, PrnProgress)
  use parameters
  implicit none
  integer       :: Npart, Mrankp(Npart), Nrankp(Npart), Mrank, Nrank
  real(O)       :: k, snorm, rp(Npart), xp(Npart), yp(Npart), zp(Npart), epsNrank,  &
                   epsMrank
  complex(O)    :: ind_refp(Npart)
  character(80) :: FileTmat
  logical       :: identical, ExtThetaDom, PrnProgress
!       
  integer       :: Nmax, Nmaxmax, Nteta, ipart, jpart, NthetaConvN, NthetaConvM,    &
                   Nl, Nc, Mrankpl, Nrankpl, Nmaxpl, Mrankpl1, Nrankpl1, Nmaxpl1,   &
                   i, NmaxAL  
  real(O)       :: r, alfaPol, x, y, z, tetaGI, phiGI, phiGS,                       &
                   x1, y1, z1, dx, dy, dz, alfaC, betaC, gamaC, Cscat, Qscat,       &
                   Cext, Qext
  complex(O)    :: ind_ref
  integer,allocatable    :: Nmaxp(:)
  real(O),allocatable    :: h(:), v(:), oldh(:), oldv(:), oldh0(:), oldv0(:)
  complex(O),allocatable :: aa(:,:), bb(:,:), a(:,:), b(:,:), t(:), c(:), c1(:)   
!
  tetaGI  = 0._O
  phiGI   = 0._O
  phiGS   = 0._O
  Nteta   = 10
  alfaC   = Pi / 4._O
  betaC   = Pi / 4._O
  gamaC   = 0._O
  alfaPol = Pi / 4._O
  Nmax    = Nrank + Mrank * (2 * Nrank - Mrank + 1) 
  NmaxAL  = Nmax
  allocate (Nmaxp(Npart))
  if (.not. identical) then
    Nmaxmax = 0
    do ipart = 1, Npart
      Nmaxp(ipart) = Nrankp(ipart) + Mrankp(ipart) * (2 * Nrankp(ipart) -           &
                     Mrankp(ipart) + 1)
      Nmaxmax      = Nmaxmax + Nmaxp(ipart)      
      if (Nmaxp(ipart) > NmaxAL) NmaxAL = Nmaxp(ipart)
      if (Nmaxp(ipart) > Nmax) then
        print "(/,2x, a, i3, a)",                                                   &
       'Warning: the input values of Nrank and Mrank are too low for particle ',ipart, ';'        
      end if 
    end do
  else 
    Nmaxp(1) = Nrankp(1) + Mrankp(1) * (2 * Nrankp(1) - Mrankp(1) + 1)
    Nmaxmax  = Nmaxp(1) * Npart
    if (Nmaxp(1) > NmaxAL) then
      NmaxAL = Nmaxp(1)     
      print "(/,2x,'Warning: the input values of Nrank and Mrank are too low;')"
    end if  
  end if     
  call write_TypeConvHead (4)
  open (unit = iTmat, file = FileTmat, status = 'replace')
  call  write_HeadFileTmat (NmaxAL, NmaxAL) 
  allocate (aa(2*Nmaxmax,2*Nmaxmax), bb(2*Nmaxmax,2*NmaxAL))        
  allocate (a(2*NmaxAL,2*NmaxAL), b(2*NmaxAL,2*NmaxAL), c(2*NmaxAL), c1(2*NmaxAL))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta), oldh0(Nteta),             &
            oldv0(Nteta))       
  call identity_matrix (2*Nmaxmax, aa, 2*Nmaxmax, 2*Nmaxmax)
  if (PrnProgress) call write_progress (.true., 1, 2*Npart+2)
  if (.not. identical) then              
    Nl = 0  
    do ipart = 1, Npart   
      Mrankpl = Mrankp(ipart)
      Nrankpl = Nrankp(ipart)
      Nmaxpl  = Nmaxp(ipart)
      allocate (t(2*Nmaxpl))
      r       = rp(ipart)
      ind_ref = ind_refp(ipart)
      call coefficients_fg (k, r, ind_ref, Mrankpl, Nrankpl, Nmaxpl, t)         
      x  = xp(ipart)
      y  = yp(ipart)
      z  = zp(ipart)    
      Nc = 0    
      do jpart = 1, Npart           
        if (jpart /= ipart) then
          Mrankpl1 = Mrankp(jpart)
          Nrankpl1 = Nrankp(jpart)
          Nmaxpl1  = Nmaxp(jpart)
          x1 = xp(jpart)
          y1 = yp(jpart)
          z1 = zp(jpart)
          dx = x1 - x
          dy = y1 - y
          dz = z1 - z       
          call MatTransAB_mn_m1n1 (3, k, dx, dy, dz, Mrankpl, Nrankpl, Nmaxpl,      &
               Mrankpl1, Nrankpl1, Nmaxpl1, a, NmaxAL, NmaxAL)
          call product_matrix_vector1 (2*Nmaxpl, 2*Nmaxpl1, t, a, 2*NmaxAL,         &
               2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL)                         
          call extend_matrix1 (ipart, jpart, Npart, Nmaxp, Nl, Nc, b, NmaxAL,       &
               NmaxAL, aa, Nmaxmax, Nmaxmax)                                                    
        end if
        Nc = Nc + 2 * Nmaxp(jpart)      
      end do        
      call MatTransAB_mn_m1n1 (1, k, -x, -y, -z, Mrankpl, Nrankpl, Nmaxpl, Mrank,   &
           Nrank, Nmax, a, NmaxAL, NmaxAL)
      call product_matrix_vector1 (2*Nmaxpl, 2*Nmax, t, a, 2*NmaxAL, 2*NmaxAL,      &
           b, 2*NmaxAL, 2*NmaxAL)    
      call extend_matrix2 (ipart, Npart, Nmaxp, Nmax, Nmaxmax, Nl, b, NmaxAL,       &
           NmaxAL, bb, Nmaxmax, NmaxAL)          
      Nl = Nl + 2 * Nmaxp(ipart)          
      if (PrnProgress) call write_progress (.false., ipart+1, 2*Npart+2)     
      deallocate (t)
    end do
  else     
    Mrankpl = Mrankp(1)
    Nrankpl = Nrankp(1)
    Nmaxpl  = Nmaxp(1)
    allocate (t(2*Nmaxpl))
    r       = rp(1)
    ind_ref = ind_refp(1)   
    call coefficients_fg (k, r, ind_ref, Mrankpl, Nrankpl, Nmaxpl, t)
    Nl = 0  
    do ipart = 1, Npart                 
      x  = xp(ipart)
      y  = yp(ipart)
      z  = zp(ipart)    
      Nc = 0    
      do jpart = 1, Npart           
        if (jpart > ipart) then          
          x1 = xp(jpart)
          y1 = yp(jpart)
          z1 = zp(jpart)
          dx = x1 - x
          dy = y1 - y
          dz = z1 - z
          call MatTransAB_mn_m1n1 (3, k, dx, dy, dz, Mrankpl, Nrankpl, Nmaxpl,      &
               Mrankpl, Nrankpl, Nmaxpl, a, NmaxAL, NmaxAL)
          call product_matrix_vector1 (2*Nmaxpl, 2*Nmaxpl, t, a, 2*NmaxAL,          &
               2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL)
          call extend_matrix3 (Nmaxpl, Nmaxpl, Nl, Nc, b, NmaxAL, NmaxAL,           &
               aa, Nmaxmax, Nmaxmax)
          call matrix_inverse (Mrankpl, Nrankpl, Nmaxpl, Mrankpl, Nrankpl, Nmaxpl,  &
               a, NmaxAL, NmaxAL)
          call product_matrix_vector1 (2*Nmaxpl, 2*Nmaxpl, t, a, 2*NmaxAL,          &
               2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL)
          call extend_matrix3 (Nmaxpl, Nmaxpl, Nc, Nl, b, NmaxAL, NmaxAL,           &
               aa, Nmaxmax, Nmaxmax)                                                                    
        end if
        Nc = Nc + 2 * Nmaxpl      
      end do        
      call MatTransAB_mn_m1n1 (1, k, -x, -y, -z, Mrankpl, Nrankpl, Nmaxpl, Mrank,   &
           Nrank, Nmax, a, NmaxAL, NmaxAL)
      call product_matrix_vector1 (2*Nmaxpl, 2*Nmax, t, a, 2*NmaxAL, 2*NmaxAL,      &
           b, 2*NmaxAL, 2*NmaxAL)
      call extend_matrix4 (ipart, Nmaxpl, Nmax, Nmaxmax, Nl, b, NmaxAL, NmaxAL,     &
           bb, Nmaxmax, NmaxAL)                  
      Nl = Nl + 2 * Nmaxpl        
      if (PrnProgress) call write_progress (.false., ipart+1, 2*Npart+2)      
    end do
    deallocate (t)
  end if        
  call LU_SYSTEM_DIRECT (aa, 2*Nmaxmax, 2*Nmaxmax, bb, 2*Nmaxmax, 2*NmaxAL,         &
       2*Nmaxmax, 2*Nmax)
  if (PrnProgress) call write_progress (.false., Npart+2, 2*Npart+2)         
  deallocate (aa)
  allocate (aa(2*NmaxAL,2*NmaxAL))
  Nl = 0
  Mrankpl = Mrankp(1)
  Nrankpl = Nrankp(1)
  Nmaxpl  = Nmaxp(1)
  do ipart = 1, Npart    
    if (.not. identical) then
      Mrankpl = Mrankp(ipart)
      Nrankpl = Nrankp(ipart)
      Nmaxpl  = Nmaxp(ipart)      
    end if
    x = xp(ipart)
    y = yp(ipart)
    z = zp(ipart)    
    if (.not. identical) then
      call extract_matrix1 (ipart, Npart, Nmaxp, Nmax, Nl, b, NmaxAL, NmaxAL,       &
           bb, Nmaxmax, NmaxAL)
    else 
      call extract_matrix2 (Nmaxpl, Nmax, Nl, b, NmaxAL, NmaxAL, bb, Nmaxmax, NmaxAL)
    end if
    call MatTransAB_mn_m1n1 (1, k, x, y, z, Mrank, Nrank, Nmax, Mrankpl, Nrankpl,   &
         Nmaxpl, aa, NmaxAL, NmaxAL)                              
    call product_matrices2 (ipart, 2*Nmax, 2*Nmaxpl, 2*Nmax, aa, 2*NmaxAL,          &
         2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL, a, 2*NmaxAL, 2*NmaxAL)                            
    Nl = Nl + 2 * Nmaxpl
    if (PrnProgress) call write_progress (.false., Npart+2+ipart, 2*Npart+2)
  end do   
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
  print "(/,2x,'The T matrix is stored in ',a50)", FileTmat
  print "(  2x,'The dimensions of the T matrix are given by:')"
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", Nrank
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank  
  deallocate (aa, bb, a, b, c, c1, h, v, oldh, oldv, oldh0, oldv0, Nmaxp)       
end subroutine convergence_Nrank_MrankMULTSPH
!***********************************************************************************
subroutine TMatrix_Nrank_MrankMULTSPH (k, Npart, identical, rp, xp, yp, zp,         &
           ind_refp, Mrankp, Nrankp, Mrank, Nrank, FileTmat, PrnProgress)
  use parameters
  implicit none
  integer       :: Npart, Mrankp(Npart), Nrankp(Npart), Mrank, Nrank
  real(O)       :: k, rp(Npart), xp(Npart), yp(Npart), zp(Npart)                   
  complex(O)    :: ind_refp(Npart)
  character(80) :: FileTmat
  logical       :: identical, PrnProgress
!       
  integer       :: Nmax, Nmaxmax, ipart, jpart, Nl, Nc, Mrankpl, Nrankpl, Nmaxpl,   &
                   Mrankpl1, Nrankpl1, Nmaxpl1, NmaxAL  
  real(O)       :: r, x, y, z, x1, y1, z1, dx, dy, dz                             
  complex(O)    :: ind_ref
  integer,allocatable    :: Nmaxp(:)
  complex(O),allocatable :: aa(:,:), bb(:,:), a(:,:), b(:,:), t(:)   
!  
  Nmax    = Nrank + Mrank * (2 * Nrank - Mrank + 1) 
  NmaxAL  = Nmax
  allocate (Nmaxp(Npart))
  if (.not. identical) then
    Nmaxmax = 0
    do ipart = 1, Npart
      Nmaxp(ipart) = Nrankp(ipart) + Mrankp(ipart) * (2 * Nrankp(ipart) -           &
                     Mrankp(ipart) + 1)
      Nmaxmax      = Nmaxmax + Nmaxp(ipart)      
      if (Nmaxp(ipart) > NmaxAL) NmaxAL = Nmaxp(ipart)
      if (Nmaxp(ipart) > Nmax) then
        print "(/,2x, a, i3, a)",                                                   &
       'Warning: the input values of Nrank and Mrank are too low for particle ',ipart,';'        
      end if 
    end do
  else 
    Nmaxp(1) = Nrankp(1) + Mrankp(1) * (2 * Nrankp(1) - Mrankp(1) + 1)
    Nmaxmax  = Nmaxp(1) * Npart
    if (Nmaxp(1) > NmaxAL) then
      NmaxAL = Nmaxp(1)     
      print "(/,2x,'Warning: the input values of Nrank and Mrank are too low;')"
    end if  
  end if     
  open (unit = iTmat, file = FileTmat, status = 'replace')
  call  write_HeadFileTmat (NmaxAL, NmaxAL) 
  allocate (aa(2*Nmaxmax,2*Nmaxmax), bb(2*Nmaxmax,2*NmaxAL))        
  allocate (a(2*NmaxAL,2*NmaxAL), b(2*NmaxAL,2*NmaxAL))
  call identity_matrix (2*Nmaxmax, aa, 2*Nmaxmax, 2*Nmaxmax)
  if (PrnProgress) call write_progress (.true., 1, 2*Npart+2)
  if (.not. identical) then              
    Nl = 0  
    do ipart = 1, Npart   
      Mrankpl = Mrankp(ipart)
      Nrankpl = Nrankp(ipart)
      Nmaxpl  = Nmaxp(ipart)
      allocate (t(2*Nmaxpl))
      r       = rp(ipart)
      ind_ref = ind_refp(ipart)
      call coefficients_fg (k, r, ind_ref, Mrankpl, Nrankpl, Nmaxpl, t)         
      x  = xp(ipart)
      y  = yp(ipart)
      z  = zp(ipart)    
      Nc = 0    
      do jpart = 1, Npart           
        if (jpart /= ipart) then
          Mrankpl1 = Mrankp(jpart)
          Nrankpl1 = Nrankp(jpart)
          Nmaxpl1  = Nmaxp(jpart)
          x1 = xp(jpart)
          y1 = yp(jpart)
          z1 = zp(jpart)
          dx = x1 - x
          dy = y1 - y
          dz = z1 - z       
          call MatTransAB_mn_m1n1 (3, k, dx, dy, dz, Mrankpl, Nrankpl, Nmaxpl,      &
               Mrankpl1, Nrankpl1, Nmaxpl1, a, NmaxAL, NmaxAL)
          call product_matrix_vector1 (2*Nmaxpl, 2*Nmaxpl1, t, a, 2*NmaxAL,         &
               2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL)                         
          call extend_matrix1 (ipart, jpart, Npart, Nmaxp, Nl, Nc, b, NmaxAL,       &
               NmaxAL, aa, Nmaxmax, Nmaxmax)                                                    
        end if
        Nc = Nc + 2 * Nmaxp(jpart)      
      end do        
      call MatTransAB_mn_m1n1 (1, k, -x, -y, -z, Mrankpl, Nrankpl, Nmaxpl, Mrank,   &
           Nrank, Nmax, a, NmaxAL, NmaxAL)
      call product_matrix_vector1 (2*Nmaxpl, 2*Nmax, t, a, 2*NmaxAL, 2*NmaxAL,      &
           b, 2*NmaxAL, 2*NmaxAL)    
      call extend_matrix2 (ipart, Npart, Nmaxp, Nmax, Nmaxmax, Nl, b, NmaxAL,       &
           NmaxAL, bb, Nmaxmax, NmaxAL)          
      Nl = Nl + 2 * Nmaxp(ipart)          
      if (PrnProgress) call write_progress (.false., ipart+1, 2*Npart+2)     
      deallocate (t)
    end do
  else     
    Mrankpl = Mrankp(1)
    Nrankpl = Nrankp(1)
    Nmaxpl  = Nmaxp(1)
    allocate (t(2*Nmaxpl))
    r       = rp(1)
    ind_ref = ind_refp(1)   
    call coefficients_fg (k, r, ind_ref, Mrankpl, Nrankpl, Nmaxpl, t)
    Nl = 0  
    do ipart = 1, Npart                 
      x  = xp(ipart)
      y  = yp(ipart)
      z  = zp(ipart)    
      Nc = 0    
      do jpart = 1, Npart           
        if (jpart > ipart) then          
          x1 = xp(jpart)
          y1 = yp(jpart)
          z1 = zp(jpart)
          dx = x1 - x
          dy = y1 - y
          dz = z1 - z
          call MatTransAB_mn_m1n1 (3, k, dx, dy, dz, Mrankpl, Nrankpl, Nmaxpl,      &
               Mrankpl, Nrankpl, Nmaxpl, a, NmaxAL, NmaxAL)
          call product_matrix_vector1 (2*Nmaxpl, 2*Nmaxpl, t, a, 2*NmaxAL,          &
               2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL)
          call extend_matrix3 (Nmaxpl, Nmaxpl, Nl, Nc, b, NmaxAL, NmaxAL,           &
               aa, Nmaxmax, Nmaxmax)
          call matrix_inverse (Mrankpl, Nrankpl, Nmaxpl, Mrankpl, Nrankpl, Nmaxpl,  &
               a, NmaxAL, NmaxAL)
          call product_matrix_vector1 (2*Nmaxpl, 2*Nmaxpl, t, a, 2*NmaxAL,          &
               2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL)
          call extend_matrix3 (Nmaxpl, Nmaxpl, Nc, Nl, b, NmaxAL, NmaxAL,           &
               aa, Nmaxmax, Nmaxmax)                                                                    
        end if
        Nc = Nc + 2 * Nmaxpl      
      end do        
      call MatTransAB_mn_m1n1 (1, k, -x, -y, -z, Mrankpl, Nrankpl, Nmaxpl, Mrank,   &
           Nrank, Nmax, a, NmaxAL, NmaxAL)
      call product_matrix_vector1 (2*Nmaxpl, 2*Nmax, t, a, 2*NmaxAL, 2*NmaxAL,      &
           b, 2*NmaxAL, 2*NmaxAL)
      call extend_matrix4 (ipart, Nmaxpl, Nmax, Nmaxmax, Nl, b, NmaxAL, NmaxAL,     &
           bb, Nmaxmax, NmaxAL)                  
      Nl = Nl + 2 * Nmaxpl        
      if (PrnProgress) call write_progress (.false., ipart+1, 2*Npart+2)      
    end do
    deallocate (t)
  end if        
  call LU_SYSTEM_DIRECT (aa, 2*Nmaxmax, 2*Nmaxmax, bb, 2*Nmaxmax, 2*NmaxAL,         &
       2*Nmaxmax, 2*Nmax)
  if (PrnProgress) call write_progress (.false., Npart+2, 2*Npart+2)         
  deallocate (aa)
  allocate (aa(2*NmaxAL,2*NmaxAL))
  Nl = 0
  Mrankpl = Mrankp(1)
  Nrankpl = Nrankp(1)
  Nmaxpl  = Nmaxp(1)
  do ipart = 1, Npart    
    if (.not. identical) then
      Mrankpl = Mrankp(ipart)
      Nrankpl = Nrankp(ipart)
      Nmaxpl  = Nmaxp(ipart)      
    end if
    x = xp(ipart)
    y = yp(ipart)
    z = zp(ipart)    
    if (.not. identical) then
      call extract_matrix1 (ipart, Npart, Nmaxp, Nmax, Nl, b, NmaxAL, NmaxAL,       &
           bb, Nmaxmax, NmaxAL)
    else 
      call extract_matrix2 (Nmaxpl, Nmax, Nl, b, NmaxAL, NmaxAL, bb, Nmaxmax, NmaxAL)
    end if
    call MatTransAB_mn_m1n1 (1, k, x, y, z, Mrank, Nrank, Nmax, Mrankpl, Nrankpl,   &
         Nmaxpl, aa, NmaxAL, NmaxAL)                              
    call product_matrices2 (ipart, 2*Nmax, 2*Nmaxpl, 2*Nmax, aa, 2*NmaxAL,          &
         2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL, a, 2*NmaxAL, 2*NmaxAL)                            
    Nl = Nl + 2 * Nmaxpl
    if (PrnProgress) call write_progress (.false., Npart+2+ipart, 2*Npart+2)
  end do   
  call write_FileTmat (NmaxAL, NmaxAL, a)  
  close (unit = iTmat)  
  call write_InfoFileTmat (FileTmat, Mrank, Nrank, .false., .false., .false.)
  call ScatCharact (k, FileTmat, Mrank, Nrank, .false., .false., .false.) 
  print "(/,2x,'The T matrix is stored in ',a50)", FileTmat
  print "(  2x,'The dimensions of the T matrix are given by:')"
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", Nrank
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank  
  deallocate (aa, bb, a, b, Nmaxp)      
end subroutine TMatrix_Nrank_MrankMULTSPH

