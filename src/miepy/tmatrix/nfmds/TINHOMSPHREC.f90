subroutine TINHOMSPHREC
!------------------------------------------------------------------------------------
! 1. General Considerations                                                         !
! --------------------------                                                        !
! TINHOMSPHREC is a code for computing the T matrix and the scattering              !
! characteristics of an inhomogeneous, dielectric sphere with multiple spherical    !
! inclusions. The program supports calculation for dielectric host spheres with     !
! real refractive index and uses a recursive algorithm for computing the T matrix.  !
! The inclusions are IDENTICALLY and randomly distributed inside the host sphere.   !
! The routine "SphConfig" from the file "Random.f90" generates a random             !
! distribution of particles using the sequential addition method. The user can      !
! implement other particle distributions by modifying this routine. However, the    !
! list of parameters must be maintained.                                            !
!                                                                                   !
! The T matrix of an inhomogeneous sphere can be expressed in terms of the T        !
! matrices of the host sphere and the inclusions as                                 !
!                                                                                   !
!            T =  { T_host - Q13 * T_incl * Q31^(-1) }                              !
!               * { I - Q33 * T_incl * * Q31^(-1) }^(-1),                           !
!                                                                                   ! 
! where T_host is the transition matrix of the host sphere and T_incl is the        !
! transition matrix of the inclusions. The recursive algorithm is used to compute   !
! the T matrix of the inclusions. Essentially, the T matrix of a system of Npart    !
! particles is computed by using the T matrices of the new added q particles and    !
! the T matrix of the previous system of Npart - q particles. Let us consider Ncs   !
! concentric spheres with radii Rcirc(i), i = 1,2,..,Ncs, in increasing order.      !
! Inside the sphere of radius Rcirc(1) there are N(1) particles, and in the         !
! spherical shell bounded by Rcirc(i-1) and Rcirc(i) there are N(i) particles. At   !
! the first iteration step we compute the T matrix T(1) of all particles situated   !
! inside the sphere of radius Rcirc(1). At the iteration step i, we compute the     ! 
! transition matrix T(i) of all particles situated inside the sphere of radius      ! 
! Rcirc(i), that is, the system T-matrix of the N(1) + N(2) + ...+ N(i) particles.  !
! The system T-matrix T(i) is computed by using the system T-matrix T(i-1) and the  !
! individual transition matrices of all N(i) particles situated in the spherical    !
! shell between Rcirc(i-1) and Rcirc(i). The procedure is repeated until all        !
! particles are exhausted, and it is apparent that at each iteration step, only a   !
! N(i)+1-scatterer problem needs to be solved, where N(i) usually is much smaller   !
! than Npart. Thus, we can keep the dimension of the problem manageable even when   !
! Npart is very large. The geometric constraint of the T-matrix method requires     !
! that the N(i) particles completely reside inside the spherical shell between      !
! Rcirc(i-1) and Rcirc(i). However, numerical simulations certify the accuracy      !
! of the recursive scheme for small size parameters and large particle numbers.     !
!                                                                                   !
! 2. Convergence Test                                                               !
! --------------------                                                              !
! Convergence tests precede the scattering characteristics calculation. The         !
! parameters which control the T-matrix computation are the maximum expansion order !
! Nrank, and the maximum azimuthal order Mrank for the host particle.               ! 
!                                                                                   !
! For convergence tests, the incident wave is assumed to be a plane wave traveling  !
! along the Z-axis of the global coordinate system and the scattering               !
! characteristics are computed in the azimuthal plane  phi = 0° . The Euler         !
! orientation angles of the host sphere are alpha = beta = 45° and gamma = 0°, and  !
! convergence tests over Nrank and Mrank are performed.                             !
!                                                                                   !
! The input file generates a sequence of Ncs concentric spheres with radii          !  
! Rcirc(i), i = 1,2,...,Ncs, in increasing order. The largest sphere is the host    !
! sphere, i.e., the radius of the host sphere is Rcirc(Ncs). For each auxiliary     !
! sphere i, the maximum expansion order NrankCirc(i) and the number of azimuthal    !
! modes MrankCirc(i) are specified in the input file. The maximum expansion         !  
! order Nrank and the number of azimuthal modes Mrank for the inclusion are also    !
! supplied in the input file. The convergence tests are carried out over            !
! NrankCirc(Ncs) and MrankCirc(Ncs) (corresponding to the host sphere), while       !
! Nrank and Mrank are kept constant. For the convergence test over the expansion    !
! order, the scattering problem is solved for the pairs                             !
! (NrankCirc(Ncs),MrankCirc(Ncs)) and (NrankCirc(Ncs) - 1,MrankCirc(Ncs)), while    !
! for the azimuthal order test, the cases (NrankCirc(Ncs),MrankCirc(Ncs)) and       !
! (NrankCirc(Ncs),MrankCirc(Ncs) - 1) are considered. The normalized differential   !
! scattering cross section will be checked at 20° increments for convergence        !
! within epsX (epsNint,epsNrank or epsMrank) tolerance. If the calculated results   !
! are converged within this tolerance at 80% of the scattering angles, then         !
! convergence is achieved. The T matrix is stored for later use by other programs,  !
! and the values of NrankCirc(Ncs) and MrankCirc(Ncs)  are printed to the screen    !
! and to the T-matrix information file (see "Description.txt"). These values        !
! together with the T matrix serve as INPUT PARAMETERS for other programs.          !
!                                                                                   !
! The values of NrankCirc(i), i = 1,2,...,Ncs, and Nrank can be computed with       !
! Wiscombe's truncation limit criterion [W. J. Wiscombe, Improved Mie scattering    !
! algorithms, Applied Optics, 19, 1505-1509, 1980]. MrankCirc(i) can be chosen as   !
! MrankCirc(i) = NrankCirc(i) - 2,...,NrankCirc(i) and Mrank as Mrank = Nrank.      !
! In general, these prescriptions guarantee convergence. The main problem which     !
! has to be solved is the selection of the concentric spheres. The radii must be    !
! chosen such that the numbers of particles in each spherical shell do not vary     ! 
! significantly. Because the total computer time can be expressed as                !
!                                                                                   !
!                      CPU_total = Ncs * CPU_shell,                                 !
!                                                                                   !
! where CPU_shell is the computer time for a specific shell calculation, and        !
!                                                                                   !    
!                      CPU_shell ~ Npart / Ncs,                                     !
!                                                                                   !
! an optimum value of Ncs should be found by repeated simulations.                  !
!                                                                                   !
! The convergence tests can be switched off by setting the logical variable         !
! DoConvTest to false.                                                              !  
!                                                                                   !
! 3. Input Parameters                                                               !
! --------------------                                                              !
! The parameters specified in the input file "/INPUTFILES/InputINHOMSPH.dat"        !
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
! - ind_refRel1 (complex) - relative refractive index of the spherical inclusions   !
!   with respect to the host sphere. The imaginary part of the relative refractive  !
!   index must be zero for nonabsorbing particles and positive for absorbing        !
!   particles.                                                                      !
!                                                                                   !
! - Npart (integer) - number of spherical inclusions.                               !
!                                                                                   !
! - anorm (real) - characteristic length of the host sphere which is used to        !
!   normalize the differential scattering cross sections.                           !
!                                                                                   !
! - Ncs (integer) - number of concentric spheres.                                   !
!                                                                                   !
! - Ntry (integer) - maximum number of calls of the sequential addition method      !
!   routine for generating the desired random distribution of spheres.              !
!                                                                                   !
! - DoConvTest (logical) - if DoConvTest = t, the convergence tests over            !
!   NrankCirc(Ncs) and MrankCirc(Ncs) are performed.                                !
!                                                                                   !
! - ExtThetaDom (logical) - if ExtThetaDom = t the DSCS is computed for             !
!   scattering angles ranging from 0° to 180° in the azimuthal plane phiGS = 0°,    !
!   and from 180° to 0° in the azimuthal plane phiGS = 180°. The total number of    !
!   scattering angles is 10. If ExtThetaDom = f the DSCS is computed for            !
!   scattering angles ranging from 0° to 180° in the azimuthal plane phiGS = 0°.    !
!                                                                                   !
! - r (real) - radius of the spherical inclusions.                                  !
!                                                                                   !
! - Nrank, Mrank (integer variables) - maximum expansion and azimuthal orders for   !
!   the spherical inclusions.                                                       !
!                                                                                   !
! The next parameters (specified in a group statement) correspond to each            !
! auxiliary sphere. THE GROUP STATEMENT MUST BE REPEATED FOR ALL Ncs SPHERES.       !
!                                                                                   !
! - Rcirc (real) - radius of the actual auxiliary sphere. Rcs must be specified     !
!   in increasing order. Note that Rcirc(Ncs) is the radius of the host sphere.     ! 
!                                                                                   !
! - NrankCirc, MrankCirc (integer variables) - maximum expansion and azimuthal      !
!   orders corresponding to the actual auxiliary sphere.                            ! 
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
!------------------------------------------------------------------------------------
  use parameters
  use allocation, only: Mrankcs, Nrankcs, Rcs
  implicit none 
  integer       :: Npart, Ncs, Mrank, Nrank, Ntry                                    
  real(O)       :: ks, ind_refMED, wavelength, ind_refRel, anorm, snorm, epsNrank,  &
                   epsMrank, r
  complex(O)    :: ind_refRel1
  character(80) :: FileTmat
  logical       :: DoConvTest, ExtThetaDom, PrnProgress   
  integer,allocatable    :: Npartcs(:)
  real(O),allocatable    :: xp(:,:), yp(:,:), zp(:,:) 
! -----------------------------------------------------------------------------------
!                                 Read the input file                               ! 
! -----------------------------------------------------------------------------------   
  call readinputINHOMSPHREC ( wavelength, ind_refMed, ind_refRel, ind_refRel1,      &
       Npart, anorm, Ncs, Ntry, DoConvTest, ExtThetaDom, r, Nrank, Mrank, epsNrank, &
       epsMrank, FileTmat, PrnProgress, ks, snorm)
! -----------------------------------------------------------------------------------
!                                      Main                                         !
! -----------------------------------------------------------------------------------                  
  allocate (Npartcs(Ncs), xp(Ncs,Npart), yp(Ncs,Npart), zp(Ncs,Npart))
  call SphConfig (Ntry, Npart, Ncs, Npartcs, r, Rcs, xp, yp, zp)    
  if (DoConvTest) then
    open (unit = iOutput, file = FileOutput, status = "replace")    
    call printinputINHOMSPHREC (wavelength, ind_refMed, ind_refRel, ind_refRel1, r, &
         Npart, Ncs, Npartcs, Rcs, Mrankcs, Nrankcs, Mrank, Nrank, epsNrank,        &
         epsMrank, anorm)
    call conv_Nrank_MrankINHOMSPHREC (ks, snorm, r, ind_refRel, ind_refRel1, Npart, &
         Ncs, Npartcs, Rcs, xp, yp, zp, Mrankcs, Nrankcs, Mrank, Nrank, epsNrank,   &
         epsMrank, ExtThetaDom, FileTmat, PrnProgress)
    close (unit = iOutput)
  else
    call TMatrix_Nrank_MrankINHOMSPHREC (ks, r, ind_refRel, ind_refRel1, Npart, Ncs,&
         Npartcs, Rcs, xp, yp, zp, Mrankcs, Nrankcs, Mrank, Nrank, FileTmat,        &
         PrnProgress)
  end if
  deallocate (Rcs, Mrankcs, Nrankcs, Npartcs, xp, yp, zp)  
end subroutine TINHOMSPHREC
!***********************************************************************************
subroutine readinputINHOMSPHREC ( wavelength, ind_refMed, ind_refRel, ind_refRel1,  &
           Npart, anorm, Ncs, Ntry, DoConvTest, ExtThetaDom, r, Nrank, Mrank,       &
           epsNrank, epsMrank, FileTmat, PrnProgress, ks, snorm)     
  use parameters
  use derived_parameters
  use allocation, only: Mrankcs, Nrankcs, Rcs
  implicit none 
  integer       :: Npart, Ncs, MrankCirc, NrankCirc, Mrank, Nrank, NrankW, ics,     &
                   Ntry, ios                                    
  real(O)       :: ks, ind_refMED, wavelength, ind_refRel, anorm, x, snorm,         &
                   epsNrank, epsMrank, r, Rcirc, xR
  complex(O)    :: ind_refRel1
  character(80) :: FileTmat, string
  logical       :: DoConvTest, ExtThetaDom, PrnProgress, XFindPar
! -----------------------------------------------------------------------------------
!                        Read the input file FileInputINHOMSPHREC                   ! 
! ----------------------------------------------------------------------------------- 
  call DrvParameters 
  open (unit = iInputINHOMSPHREC, file = FileInputINHOMSPHREC, status = "old",      &
        position = "rewind")   
  wavelength  = 0.1_O * 2._O * Pi                                                                   
  ind_refMED  = 1._O
  ind_refRel  = 1.2_O
  ind_refRel1 = (1.5_O,0._O)
  string     = 'OptProp'
  if (XFindPar (iInputINHOMSPHREC, string)) then
    read (iInputINHOMSPHREC, *, iostat = ios) wavelength
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable wavelength;')"
      stop
    end if
    read (iInputINHOMSPHREC, *, iostat = ios) ind_refMed
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refMed;')"
      stop
    end if
    read (iInputINHOMSPHREC, *, iostat = ios) ind_refRel
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refRel;')"
      stop
    end if 
    read (iInputINHOMSPHREC, *, iostat = ios) ind_refRel1
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
  Npart  = 2
  anorm  = 1._O 
  Ncs    = 2  
  Ntry   = 1000
  string = 'GenProp'
  if (XFindPar (iInputINHOMSPHREC, string)) then
    read (iInputINHOMSPHREC, *, iostat = ios) Npart
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Npart;')"
      stop
    end if
    read (iInputINHOMSPHREC, *, iostat = ios) anorm
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable anorm;')"
      stop
    end if
    read (iInputINHOMSPHREC, *, iostat = ios) Ncs
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Ncs;')"
      stop
    end if 
    read (iInputINHOMSPHREC, *, iostat = ios) Ntry
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Ntry;')"
      stop
    end if
  else
    print "(/,2x,'Group name GenProp not found;')"
    stop  
  end if   
  call check_anorm (anorm) 
  x  = ks * anorm
  snorm  = Pi * x * x  
!
  DoConvTest  = .true.
  ExtThetaDom = .true.
  string      = 'ConvTest'
  if (XFindPar (iInputINHOMSPHREC, string)) then
    read (iInputINHOMSPHREC, *, iostat = ios) DoConvTest
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable DoConvTest;')"
      stop
    end if
    read (iInputINHOMSPHREC, *, iostat = ios) ExtThetaDom
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ExtThetaDom;')"
      stop
    end if
  else
    print "(/,2x,'Group name ConvTest not found;')"
    stop  
  end if      
!
  r = 0.2
  Nrank  = 3
  Mrank  = 3
  string = 'TmatPart'
  if (XFindPar (iInputINHOMSPHREC, string)) then
    read (iInputINHOMSPHREC, *, iostat = ios) r
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable r;')"
      stop
    end if
    read (iInputINHOMSPHREC, *, iostat = ios) Nrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nrank;')"
      stop
    end if
    read (iInputINHOMSPHREC, *, iostat = ios) Mrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Mrank;')"
      stop
    end if
  else
    print "(/,2x,'Group name TmatPart not found;')"
    stop  
  end if  
!  
  if (DoConvTest) then    
    print "(/,2x,'Convergence Test for a Sphere with Spherical Inclusions')"
    print "(  2x,'-------------------------------------------------------')"                 
  else
    print "(/,2x,'T-Matrix Computation for a Sphere with Spherical Inclusions')"
    print "(  2x,'-----------------------------------------------------------')"
  end if    
!
  allocate (Rcs(Ncs), Mrankcs(Ncs), Nrankcs(Ncs))
  do ics = 1, Ncs
    Rcirc = 0.5 
    NrankCirc = 7
    MrankCirc = 5   
    string    = 'RecProp'
    if (XFindPar (iInputINHOMSPHREC, string)) then
      read (iInputINHOMSPHREC, *, iostat = ios) Rcirc
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Rcirc;')"
        print "(  2x,'for the circumscribing sphere ',i3,';')", ics
        stop
      end if
      read (iInputINHOMSPHREC, *, iostat = ios) NrankCirc
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable NrankCirc;')"
        print "(  2x,'for the circumscribing sphere ',i3,';')", ics
        stop
      end if
      read (iInputINHOMSPHREC, *, iostat = ios) MrankCirc
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable MrankCirc;')"
        print "(  2x,'for the circumscribing sphere ',i3,';')", ics
        stop
      end if
    else
      print "(/,2x,'Group name RecProp not found;')"
      stop  
    end if
    call check_MrankNrank (MrankCirc, NrankCirc)
    Rcs(ics) = Rcirc
    Nrankcs(ics) = NrankCirc
    Mrankcs(ics) = MrankCirc
  end do
  if (Ncs > 1) call check_circum_radii (Ncs, Rcs)   
  xR = ks * Rcs(Ncs)
  NrankW = int(xR + 4.05_O * xR**0.33_O + 2._O)
  print "(/,2x,'Input values:')"
  print "(  2x, a, i3, a, i3, a)",                                                  &
 'the input values of Nrank and Mrank for the host sphere are ', Nrankcs(Ncs),      &
 ' and ', Mrankcs(Ncs), ','
  print "(  2x, a, i3, a)",                                                         &
 'while the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';'       
!
  epsNrank = 5.e-2_O
  epsMrank = 5.e-2_O
  string   = 'Errors'
  if (XFindPar (iInputINHOMSPHREC, string)) then
    read (iInputINHOMSPHREC, *, iostat = ios) epsNrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsNrank;')"
      stop
    end if
    read (iInputINHOMSPHREC, *, iostat = ios) epsMrank
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
  if (XFindPar (iInputINHOMSPHREC, string)) then
    read (iInputINHOMSPHREC, *, iostat = ios) FileTmat
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
  if (XFindPar (iInputINHOMSPHREC, string)) then
    read (iInputINHOMSPHREC, *, iostat = ios) PrnProgress
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable PrnProgress;')"
      stop
    end if
  else
    print "(/,2x,'Group name PrintProgress not found;')"
    stop  
  end if        
  close (unit = iInputINHOMSPHREC) 
end subroutine readinputINHOMSPHREC
!***********************************************************************************
subroutine printinputINHOMSPHREC (wavelength, ind_refMed, ind_ref, ind_ref1, r,     &
           Npart, Ncs, Npartcs, Rcs, Mrankcs, Nrankcs, Mrank, Nrank, epsNrank,      &
           epsMrank, anorm)
  use parameters
  implicit none  
  integer       :: Npart, Ncs, Npartcs(Ncs), Mrankcs(Ncs), Nrankcs(Ncs), Mrank,     &
                   Nrank, ics
  real(O)       :: wavelength, ind_refMed, ind_ref, anorm, r, Rcs(Ncs), epsNrank,   &
                   epsMrank
  complex(O)    :: ind_ref1  
!
  write (iOutput,"(/,2x,'Input Parameters of the Scattering Problem:',/)")
  write (iOutput,"(2x,'wavelength of the free space, wavelength = ',1pe13.4,';')")  &
         wavelength
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'refractive index of the ambient medium, ind_refMed = ', ind_refMed, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'relative refractive index of the host particle, ind_refRel = ', ind_ref, ';'
  write (iOutput,"(2x,'relative refractive index of the inclusion with respect ')")
  write (iOutput,"(2x, a, 1pe10.3,',',1pe10.3, a)")                                 &
 'to the host particle, ind_refRel1 = (', ind_ref1, ');' 
  write (iOutput,*)
  write (iOutput,"(2x,'radius of the host particle, R = ',1pe10.3,';')") Rcs(Ncs)
  write (iOutput,"(2x,'radius of each spherical inclusion, r = ',1pe10.3,';')") r
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'characteristic length of the particle, anorm = ', anorm, ';'
  write (iOutput,"(2x,'number of spherical inclusions: Npart = ',i3,';')") Npart
  write (iOutput,"(2x, a, i3, a)")                                                  &
 'maximum expansion order for each inclusion, Nrank = ', Nrank, ';'
  write (iOutput,"(2x, a, i3, a)")                                                  &
 'maximum azimuthal order for each inclusion, Mrank = ', Mrank, ';'
  write (iOutput,"(2x,'number of auxiliary surfaces, Ncs = ',i3,';')") Ncs
  write (iOutput,*)
  do ics = 1, Ncs
    write (iOutput,"(2x, a, i3, a, i3, a)")                                         &
   'number of inclusions contained in the region ', ics, ', Npart = ',              &
    Npartcs(ics), ';'
    write (iOutput,"(2x, a, i3, a, 1pe10.3, a)")                                    &
   'radius of the region ', ics, ', Rcs = ', Rcs(ics), ';'
    write (iOutput,"(2x, a, i3, a, i3, a)")                                         &
   'maximum expansion order for region ', ics, ', Nrank = ', Nrankcs(ics), ';'
    write (iOutput,"(2x, a, i3, a, i3, a)")                                         &
   'maximum azimuthal order for region ', ics, ', Mrank = ', Mrankcs(ics), ';'
    write (iOutput,*)
  end do
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum expansion order tolerance, epsNrank = ', epsNrank, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum azimuthal order tolerance, epsMrank = ', epsMrank, '.'            
  write (iOutput,"(/)")                               
end subroutine printinputINHOMSPHREC
!***********************************************************************************
subroutine conv_Nrank_MrankINHOMSPHREC (ks, snorm, r, ind_ref, ind_ref1, Npart,     &
           Ncs, Npartcs, Rcs, xp, yp, zp, Mrankcs, Nrankcs, Mrank, Nrank, epsNrank, &
           epsMrank, ExtThetaDom, FileTmat, PrnProgress)
  use parameters
  implicit none
  integer       :: Npart, Ncs, Npartcs(Ncs), Mrankcs(Ncs), Nrankcs(Ncs), Mrank,     &
                   Nrank
  real(O)       :: ks, snorm, ind_ref, r, Rcs(Ncs), xp(Ncs,Npart), yp(Ncs,Npart),   &
                   zp(Ncs,Npart), epsNrank,epsMrank
  complex(O)    :: ind_ref1
  character(80) :: FileTmat
  logical       :: ExtThetaDom, PrnProgress
!       
  integer       :: Nmax, Nmaxmax, Nteta, ipart, jpart, NthetaConvN, NthetaConvM,    &
                   Nl, Nc, i, NmaxAL, ics  
  real(O)       :: k, alfap, x, y, z, tetaGI, phiGI, phiGS, x1, y1, z1, dx, dy, dz, &
                   alfa, beta, gama, Cscat, Qscat, Cext, Qext
  complex(O)    :: ind_refC
  integer,allocatable    :: Nmaxcs(:)
  real(O),allocatable    :: h(:), v(:), oldh(:), oldv(:), oldh0(:), oldv0(:)
  complex(O),allocatable :: aa(:,:), bb(:,:), a(:,:), b(:,:), c(:,:), tt(:,:),      &
                            t(:), t11(:), t13(:), t31(:), t33(:), cc(:), cc1(:)   
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = Pi / 4._O
  beta   = Pi / 4._O
  gama   = 0._O
  alfap  = Pi / 4._O
  k = ks * ind_ref
  ind_refC = cmplx(ind_ref,0.0,O)  
  Nmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  allocate (Nmaxcs(Ncs))
  do ics = 1, Ncs
    Nmaxcs(ics) = Nrankcs(ics) + Mrankcs(ics) * (2 * Nrankcs(ics) -                 &
                  Mrankcs(ics) + 1)
  end do 
  call write_TypeConvHead (4)
  open (unit = iTmat, file = FileTmat, status = 'replace')
  call  write_HeadFileTmat (Nmaxcs(Ncs), Nmaxcs(Ncs)) 
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta), oldh0(Nteta),             &
            oldv0(Nteta))
  allocate (tt(2*Nmaxcs(1),2*Nmaxcs(1)), t(2*Nmax))      
  call coefficients_fg (k, r, ind_ref1, Mrank, Nrank, Nmax, t)
  if (PrnProgress)                                                                  &
      print "(/,2x,'progress of main calculation for the inclusion T matrix;')"
  do ics = 1, Ncs      
    if (PrnProgress)                                                                &
        print "(2x,'progress of main calculation for region',i3,':')", ics 
    if (ics == 1) then
      Nmaxmax = Nmax * Npartcs(ics)
    else
      Nmaxmax = Nmax * Npartcs(ics) + Nmaxcs(ics-1)
    end if
    NmaxAL = Nmaxcs(ics)
    if (ics > 1) NmaxAL = max(Nmaxcs(ics),Nmaxcs(ics-1))
    if (Nmax > NmaxAL) then
      NmaxAL = Nmax 
      print "(/,2x,'Warning: the input values of Nrank and Mrank for the')"
      print "(  2x,'auxiliary sphere ',i2,' are too low;')", ics
    end if 
      allocate (aa(2*Nmaxmax,2*Nmaxmax), bb(2*Nmaxmax,2*NmaxAL))        
    allocate (a(2*NmaxAL,2*NmaxAL), b(2*NmaxAL,2*NmaxAL), c(2*NmaxAL,2*NmaxAL))
    call identity_matrix (2*Nmaxmax, aa, 2*Nmaxmax,2*Nmaxmax)    
    Nl = 0  
    do ipart = 1, Npartcs(ics)              
      x  = xp(ics,ipart)
      y  = yp(ics,ipart)
      z  = zp(ics,ipart)    
      Nc = 0    
      do jpart = 1, Npartcs(ics)            
        if (jpart > ipart) then          
          x1 = xp(ics,jpart)
          y1 = yp(ics,jpart)
          z1 = zp(ics,jpart)
          dx = x1 - x
          dy = y1 - y
          dz = z1 - z
          call MatTransAB_mn_m1n1 (3, k, dx, dy, dz, Mrank, Nrank, Nmax, Mrank,     &
               Nrank, Nmax, a, NmaxAL, NmaxAL)
          call product_matrix_vector1 (2*Nmax, 2*Nmax, t, a, 2*NmaxAL, 2*NmaxAL,    &
               b, 2*NmaxAL, 2*NmaxAL)
          call extend_matrix3 (Nmax, Nmax, Nl, Nc, b, NmaxAL, NmaxAL, aa, Nmaxmax,  &
               Nmaxmax)
          call matrix_inverse (Mrank, Nrank, Nmax, Mrank, Nrank, Nmax, a, NmaxAL,   &
               NmaxAL)
          call product_matrix_vector1 (2*Nmax, 2*Nmax, t, a, 2*NmaxAL, 2*NmaxAL,    &
               b, 2*NmaxAL, 2*NmaxAL)
          call extend_matrix3 (Nmax, Nmax, Nc, Nl, b, NmaxAL, NmaxAL, aa, Nmaxmax,  &
               Nmaxmax)                                                                 
        end if
        Nc = Nc + 2 * Nmax      
      end do      
      if (ics > 1) then
        call MatTransAB_mn_m1n1 (3, k, -x, -y, -z, Mrank, Nrank, Nmax,              &
             Mrankcs(ics-1), Nrankcs(ics-1), Nmaxcs(ics-1), a, NmaxAL, NmaxAL)
        call product_matrix_vector1 (2*Nmax, 2*Nmaxcs(ics-1), t, a, 2*NmaxAL,       &
             2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL)
        call extend_matrix3 (Nmax, Nmaxcs(ics-1), Nl, Nc, b, NmaxAL, NmaxAL, aa,    &
             Nmaxmax, Nmaxmax)      
        call MatTransAB_mn_m1n1 (3, k, x, y, z, Mrankcs(ics-1), Nrankcs(ics-1),     &
             Nmaxcs(ics-1), Mrank, Nrank, Nmax, a, NmaxAL, NmaxAL)
        call product_matrices1 (2*Nmaxcs(ics-1), 2*Nmaxcs(ics-1), 2*Nmax, tt,       &
             2*Nmaxcs(ics-1), 2*Nmaxcs(ics-1), a, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL,  &
             2*NmaxAL)       
        call extend_matrix3 (Nmaxcs(ics-1), Nmax, Nc, Nl, b, NmaxAL, NmaxAL, aa,    &
             Nmaxmax, Nmaxmax)      
      end if  
      call MatTransAB_mn_m1n1 (1, k, -x, -y, -z, Mrank, Nrank, Nmax, Mrankcs(ics),  &
           Nrankcs(ics), Nmaxcs(ics), a, NmaxAL, NmaxAL)
      call product_matrix_vector1 (2*Nmax, 2*Nmaxcs(ics), t, a, 2*NmaxAL, 2*NmaxAL, &
           b, 2*NmaxAL, 2*NmaxAL)
      call extend_matrix4 (ipart, Nmax, Nmaxcs(ics), Nmaxmax, Nl, b, NmaxAL,        &
           NmaxAL, bb, Nmaxmax, NmaxAL)                  
      Nl = Nl + 2 * Nmax
      if (PrnProgress) call write_progress (.false., ipart, 2*Npartcs(ics)+1)      
    end do
    if (ics > 1) then
      call identity_transformation (Mrankcs(ics-1), Nrankcs(ics-1), Nmaxcs(ics-1),  &
           Mrankcs(ics), Nrankcs(ics), Nmaxcs(ics), a, NmaxAL, NmaxAL)
      call product_matrices1 (2*Nmaxcs(ics-1), 2*Nmaxcs(ics-1), 2*Nmaxcs(ics), tt,  &
           2*Nmaxcs(ics-1), 2*Nmaxcs(ics-1), a, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL,    &
           2*NmaxAL)
      call extend_matrix4 (2, Nmaxcs(ics-1), Nmaxcs(ics), Nmaxmax, Nl, b, NmaxAL,   &
           NmaxAL, bb, Nmaxmax, NmaxAL)             
    end if
    call LU_SYSTEM_DIRECT (aa, 2*Nmaxmax, 2*Nmaxmax, bb, 2*Nmaxmax, 2*NmaxAL,       &
         2*Nmaxmax, 2*Nmaxcs(ics)) 
    if (PrnProgress) call write_progress (.false., Npartcs(ics)+1, 2*Npartcs(ics)+1)              
    Nl = 0  
    do ipart = 1, Npartcs(ics)    
      x = xp(ics,ipart)
      y = yp(ics,ipart)
      z = zp(ics,ipart)         
      call extract_matrix2 (Nmax, Nmaxcs(ics), Nl, b, NmaxAL, NmaxAL, bb, Nmaxmax,  &
           NmaxAL)
      call MatTransAB_mn_m1n1 (1, k, x, y, z, Mrankcs(ics), Nrankcs(ics),           &
           Nmaxcs(ics), Mrank, Nrank, Nmax, a, NmaxAL, NmaxAL)                          
      call product_matrices2 (ipart, 2*Nmaxcs(ics), 2*Nmax, 2*Nmaxcs(ics), a,       &
           2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL, c, 2*NmaxAL, 2*NmaxAL)                        
      Nl = Nl + 2 * Nmax
      if (PrnProgress) call write_progress (.false., Npartcs(ics)+1+ipart,          &
                            2*Npartcs(ics)+1)
    end do   
    if (ics > 1) then
      call extract_matrix2 (Nmaxcs(ics-1), Nmaxcs(ics), Nl, b, NmaxAL, NmaxAL,      &
           bb, Nmaxmax, NmaxAL)
      call identity_transformation (Mrankcs(ics), Nrankcs(ics), Nmaxcs(ics),        &
           Mrankcs(ics-1), Nrankcs(ics-1), Nmaxcs(ics-1), a, NmaxAL, NmaxAL)
      call product_matrices2 (2, 2*Nmaxcs(ics), 2*Nmaxcs(ics-1), 2*Nmaxcs(ics),     &
           a, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL, c, 2*NmaxAL, 2*NmaxAL)
      deallocate (tt)
      allocate (tt(2*Nmaxcs(ics),2*Nmaxcs(ics)))
    end if
    call copy_matrix (2*Nmaxcs(ics), 2*Nmaxcs(ics), c, 2*NmaxAL, 2*NmaxAL, tt,      &
         2*NmaxAL, 2*NmaxAL)
    deallocate (aa, bb, a, b, c)
  end do 
  deallocate (t)
  allocate (t(2*Nmaxcs(Ncs)), t11(2*Nmaxcs(Ncs)), t13(2*Nmaxcs(Ncs)),               &
            t31(2*Nmaxcs(Ncs)), t33(2*Nmaxcs(Ncs)))
  allocate (a(2*Nmaxcs(Ncs),2*Nmaxcs(Ncs)), b(2*Nmaxcs(Ncs),2*Nmaxcs(Ncs)),         &
            cc(2*Nmaxcs(Ncs)), cc1(2*Nmaxcs(Ncs)))  
  if (PrnProgress)                                                                  &
      print "(/,2x,'progress of main calculation for the global T matrix;')"                    
  call vector_Q_sphere (3, 3, ks, Rcs(Ncs), ind_refC, Mrankcs(Ncs), Nrankcs(Ncs),   &
       Nmaxcs(Ncs), t)
  call copy_vector (t, t33, 2*Nmaxcs(Ncs))      
  call product_matrix_vector1 (2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs), t, tt, 2*Nmaxcs(Ncs),  &
       2*Nmaxcs(Ncs), a, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs))
  call vector_Q_sphere (3, 1, ks, Rcs(Ncs), ind_refC, Mrankcs(Ncs), Nrankcs(Ncs),   &
       Nmaxcs(Ncs), t)
  call copy_vector (t, t31, 2*Nmaxcs(Ncs))      
  call sum_diagonal_elements (2*Nmaxcs(Ncs), t, a, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs)) 
  if (PrnProgress) call write_progress (.false., 1, 3)
  call vector_Q_sphere (1, 3, ks, Rcs(Ncs), ind_refC, Mrankcs(Ncs), Nrankcs(Ncs),   &
       Nmaxcs(Ncs), t)
  call copy_vector (t, t13, 2*Nmaxcs(Ncs))      
  call product_matrix_vector1 (2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs), t, tt, 2*Nmaxcs(Ncs),  &
       2*Nmaxcs(Ncs), b, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs))
  call vector_Q_sphere (1, 1, ks, Rcs(Ncs), ind_refC, Mrankcs(Ncs), Nrankcs(Ncs),   &
       Nmaxcs(Ncs), t)
  call copy_vector (t, t11, 2*Nmaxcs(Ncs))      
  call sum_diagonal_elements (2*Nmaxcs(Ncs), t, b, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs))  
  if (PrnProgress) call write_progress (.false., 2, 3)
  call LU_SYSTEM (a, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs), b, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs), &
       2*Nmaxcs(Ncs))
  if (PrnProgress) call write_progress (.false., 3, 3)
  call minus_matrix (2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs), b, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs))
  call write_FileTmat (Nmaxcs(Ncs), Nmaxcs(Ncs), b)     
  call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrankcs(Ncs),     &
       Nrankcs(Ncs), Nmaxcs(Ncs), cc)    
  call product_matrix_vector (2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs), b, 2*Nmaxcs(Ncs),       &
       2*Nmaxcs(Ncs), cc, cc1)
  call DSCS (cc1, Mrankcs(Ncs), Nrankcs(Ncs), Nmaxcs(Ncs), Nteta, phiGS, alfa,      &
       beta, gama, ks, snorm, ExtThetaDom,.true., h, v)
  call CQscat (cc1, Mrankcs(Ncs), Nrankcs(Ncs), Nmaxcs(Ncs), ks, snorm, Cscat, Qscat)
  call CQext (cc1, Mrankcs(Ncs), Nrankcs(Ncs), Nmaxcs(Ncs), tetaGI, phiGI, alfa,    &
       beta, gama, alfap, ks, snorm, Cext, Qext)
  call write_2ConvParam (Nrankcs(Ncs), Mrankcs(Ncs))              
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
  call copy_vector (t33, t, 2*Nmaxcs(Ncs))
  call vector_Nrank_1 (Mrankcs(Ncs), Nrankcs(Ncs), Nmaxcs(Ncs), t)      
  call product_matrix_vector1 (2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs), t, tt, 2*Nmaxcs(Ncs),  &
       2*Nmaxcs(Ncs), a, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs))
  call copy_vector (t31, t, 2*Nmaxcs(Ncs))
  call vector_Nrank_1 (Mrankcs(Ncs), Nrankcs(Ncs), Nmaxcs(Ncs), t)      
  call sum_diagonal_elements (2*Nmaxcs(Ncs), t, a, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs))
  call matrix_Nrank_1_left (Mrankcs(Ncs), Nrankcs(Ncs), Nmaxcs(Ncs), a, Nmaxcs(Ncs),&
       Nmaxcs(Ncs))  
  call copy_vector (t13, t, 2*Nmaxcs(Ncs))
  call vector_Nrank_1 (Mrankcs(Ncs), Nrankcs(Ncs), Nmaxcs(Ncs), t)      
  call product_matrix_vector1 (2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs), t, tt, 2*Nmaxcs(Ncs),  &
       2*Nmaxcs(Ncs), b, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs))
  call copy_vector (t11, t, 2*Nmaxcs(Ncs))
  call vector_Nrank_1 (Mrankcs(Ncs), Nrankcs(Ncs), Nmaxcs(Ncs), t)      
  call sum_diagonal_elements (2*Nmaxcs(Ncs), t, b, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs))  
  call matrix_Nrank_1_right (Mrankcs(Ncs), Nrankcs(Ncs), Nmaxcs(Ncs), b,            &
       Nmaxcs(Ncs), Nmaxcs(Ncs))  
  call LU_SYSTEM (a, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs), b, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs), &
       2*Nmaxcs(Ncs))
  call minus_matrix (2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs), b, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs)) 
  call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrankcs(Ncs),     &
       Nrankcs(Ncs), Nmaxcs(Ncs), cc)    
  call product_matrix_vector (2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs), b, 2*Nmaxcs(Ncs),       &
       2*Nmaxcs(Ncs), cc, cc1)
  call DSCS (cc1, Mrankcs(Ncs), Nrankcs(Ncs), Nmaxcs(Ncs), Nteta, phiGS, alfa,      &
       beta, gama, ks, snorm, ExtThetaDom,.true., h, v)
  call CQscat (cc1, Mrankcs(Ncs), Nrankcs(Ncs), Nmaxcs(Ncs), ks, snorm, Cscat, Qscat)
  call CQext (cc1, Mrankcs(Ncs), Nrankcs(Ncs), Nmaxcs(Ncs), tetaGI, phiGI, alfa,    &
       beta, gama, alfap, ks, snorm, Cext, Qext)
  call delta_DSCS (Nteta, h, v, oldh, oldv, epsNrank, NthetaConvN)
  call write_2ConvParam (Nrankcs(Ncs) - 1, Mrankcs(Ncs))
  call write_DSCS (Nteta, ExtThetaDom, h, v)
  call write_Effic (Qscat, Qext)
  call write_NrankConvRes (NthetaConvN, Nteta, epsNrank)
! --- (Mrank - 1) configuration ---
  call copy_vector (t33, t, 2*Nmaxcs(Ncs))
  call vector_Mrank_1 (Mrankcs(Ncs), Nrankcs(Ncs), Nmaxcs(Ncs), t)      
  call product_matrix_vector1 (2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs), t, tt, 2*Nmaxcs(Ncs),  &
       2*Nmaxcs(Ncs), a, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs))
  call copy_vector (t31, t, 2*Nmaxcs(Ncs))
  call vector_Mrank_1 (Mrankcs(Ncs), Nrankcs(Ncs), Nmaxcs(Ncs), t)      
  call sum_diagonal_elements (2*Nmaxcs(Ncs), t, a, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs))
  call matrix_Mrank_1_left (Mrankcs(Ncs), Nrankcs(Ncs), Nmaxcs(Ncs), a, Nmaxcs(Ncs),&
       Nmaxcs(Ncs)) 
  call copy_vector (t13, t, 2*Nmaxcs(Ncs))
  call vector_Mrank_1 (Mrankcs(Ncs), Nrankcs(Ncs), Nmaxcs(Ncs), t)      
  call product_matrix_vector1 (2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs), t, tt, 2*Nmaxcs(Ncs),  &
       2*Nmaxcs(Ncs), b, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs))
  call copy_vector (t11, t, 2*Nmaxcs(Ncs))
  call vector_Mrank_1 (Mrankcs(Ncs), Nrankcs(Ncs), Nmaxcs(Ncs), t)      
  call sum_diagonal_elements (2*Nmaxcs(Ncs), t, b, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs))  
  call matrix_Mrank_1_right (Mrankcs(Ncs), Nrankcs(Ncs), Nmaxcs(Ncs), b,            &
       Nmaxcs(Ncs), Nmaxcs(Ncs)) 
  call LU_SYSTEM (a, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs), b, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs), &
       2*Nmaxcs(Ncs))
  call minus_matrix (2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs), b, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs))
  call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrankcs(Ncs),     &
       Nrankcs(Ncs), Nmaxcs(Ncs), cc)    
  call product_matrix_vector (2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs), b, 2*Nmaxcs(Ncs),       &
       2*Nmaxcs(Ncs), cc, cc1)
  call DSCS (cc1, Mrankcs(Ncs), Nrankcs(Ncs), Nmaxcs(Ncs), Nteta, phiGS, alfa,      &
       beta, gama, ks, snorm, ExtThetaDom,.true., h, v)
  call CQscat (cc1, Mrankcs(Ncs), Nrankcs(Ncs), Nmaxcs(Ncs), ks, snorm, Cscat, Qscat)
  call CQext (cc1, Mrankcs(Ncs), Nrankcs(Ncs), Nmaxcs(Ncs), tetaGI, phiGI, alfa,    &
       beta, gama, alfap, ks, snorm, Cext, Qext)
  call delta_DSCS (Nteta, h, v, oldh0, oldv0, epsMrank, NthetaConvM)
  call write_2ConvParam (Nrankcs(Ncs), Mrankcs(Ncs) - 1)
  call write_DSCS (Nteta, ExtThetaDom, h, v)
  call write_Effic (Qscat, Qext)
  call write_MrankConvRes (NthetaConvM, epsMrank)
  if (NthetaConvN >= int(0.8*Nteta) .and. NthetaConvM >= int(0.8*Nteta)) then
    print "(/,2x,'Convergence criteria for Nrank and Mrank are satisfied;')"                               
  else
    print "(/,2x,'Convergence criteria for Nrank and Mrank are not satisfied;')"
  end if 
  call write_InfoFileTmat (FileTmat, Mrankcs(Ncs), Nrankcs(Ncs), .false., .false.,  &
      .false.)
  call ScatCharact (ks, FileTmat, Mrankcs(Ncs), Nrankcs(Ncs), .false., .false.,     &
      .false.)          
  print "(/,2x,'T matrix is stored in ',a50)", FileTmat
  print "(  2x,'The dimensions of the T matrix are given by:')"                          
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", Nrankcs(Ncs)
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrankcs(Ncs)                       
  deallocate (a, b, tt, t, t11, t13, t31, t33, cc, cc1, h, v, oldh, oldv, oldh0,    &
              oldv0, Nmaxcs)    
end subroutine conv_Nrank_MrankINHOMSPHREC
!***********************************************************************************
subroutine TMatrix_Nrank_MrankINHOMSPHREC (ks, r, ind_ref, ind_ref1, Npart, Ncs,    &
           Npartcs, Rcs, xp, yp, zp, Mrankcs, Nrankcs, Mrank, Nrank, FileTmat,      &
           PrnProgress)
  use parameters
  implicit none
  integer       :: Npart, Ncs, Npartcs(Ncs), Mrankcs(Ncs), Nrankcs(Ncs), Mrank,     &
                   Nrank
  real(O)       :: ks, ind_ref, r, Rcs(Ncs), xp(Ncs,Npart), yp(Ncs,Npart),          &
                   zp(Ncs,Npart)
  complex(O)    :: ind_ref1
  character(80) :: FileTmat
  logical       :: PrnProgress
!       
  integer       :: Nmax, Nmaxmax, ipart, jpart, Nl, Nc, NmaxAL, ics  
  real(O)       :: k, x, y, z, x1, y1, z1, dx, dy, dz                   
  complex(O)    :: ind_refC
  integer,allocatable    :: Nmaxcs(:)
  complex(O),allocatable :: aa(:,:), bb(:,:), a(:,:), b(:,:), c(:,:), tt(:,:), t(:)
!  
  k = ks * ind_ref
  ind_refC = cmplx(ind_ref,0.0,O)  
  Nmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  allocate (Nmaxcs(Ncs))
  do ics = 1, Ncs
    Nmaxcs(ics) = Nrankcs(ics) + Mrankcs(ics) * (2 * Nrankcs(ics) -                 &
                  Mrankcs(ics) + 1)
  end do 
  open (unit = iTmat, file = FileTmat, status = 'replace')
  call  write_HeadFileTmat (Nmaxcs(Ncs), Nmaxcs(Ncs)) 
  allocate (tt(2*Nmaxcs(1),2*Nmaxcs(1)), t(2*Nmax))      
  call coefficients_fg (k, r, ind_ref1, Mrank, Nrank, Nmax, t)
  if (PrnProgress)                                                                  &
      print "(/,2x,'progress of main calculation for the inclusion T matrix;')"
  do ics = 1, Ncs      
    if (PrnProgress)                                                                &
        print "(2x,'progress of main calculation for region',i3,':')", ics 
    if (ics == 1) then
      Nmaxmax = Nmax * Npartcs(ics)
    else
      Nmaxmax = Nmax * Npartcs(ics) + Nmaxcs(ics-1)
    end if
    NmaxAL = Nmaxcs(ics)
    if (ics > 1) NmaxAL = max(Nmaxcs(ics),Nmaxcs(ics-1))
    if (Nmax > NmaxAL) then
      NmaxAL = Nmax 
      print "(/,2x,'Warning: the input values of Nrank and Mrank for the')"
      print "(  2x,'auxiliary sphere ',i2,' are too low;')", ics
    end if 
    allocate (aa(2*Nmaxmax,2*Nmaxmax), bb(2*Nmaxmax,2*NmaxAL))        
    allocate (a(2*NmaxAL,2*NmaxAL), b(2*NmaxAL,2*NmaxAL), c(2*NmaxAL,2*NmaxAL))
    call identity_matrix (2*Nmaxmax, aa, 2*Nmaxmax,2*Nmaxmax)    
    Nl = 0  
    do ipart = 1, Npartcs(ics)              
      x  = xp(ics,ipart)
      y  = yp(ics,ipart)
      z  = zp(ics,ipart)    
      Nc = 0    
      do jpart = 1, Npartcs(ics)            
        if (jpart > ipart) then          
          x1 = xp(ics,jpart)
          y1 = yp(ics,jpart)
          z1 = zp(ics,jpart)
          dx = x1 - x
          dy = y1 - y
          dz = z1 - z
          call MatTransAB_mn_m1n1 (3, k, dx, dy, dz, Mrank, Nrank, Nmax, Mrank,     &
               Nrank, Nmax, a, NmaxAL, NmaxAL)
          call product_matrix_vector1 (2*Nmax, 2*Nmax, t, a, 2*NmaxAL, 2*NmaxAL,    &
               b, 2*NmaxAL, 2*NmaxAL)
          call extend_matrix3 (Nmax, Nmax, Nl, Nc, b, NmaxAL, NmaxAL, aa, Nmaxmax,  &
               Nmaxmax)
          call matrix_inverse (Mrank, Nrank, Nmax, Mrank, Nrank, Nmax, a, NmaxAL,   &
               NmaxAL)
          call product_matrix_vector1 (2*Nmax, 2*Nmax, t, a, 2*NmaxAL, 2*NmaxAL,    &
               b, 2*NmaxAL, 2*NmaxAL)
          call extend_matrix3 (Nmax, Nmax, Nc, Nl, b, NmaxAL, NmaxAL, aa, Nmaxmax,  &
               Nmaxmax)                                                                 
        end if
        Nc = Nc + 2 * Nmax      
      end do      
      if (ics > 1) then
        call MatTransAB_mn_m1n1 (3, k, -x, -y, -z, Mrank, Nrank, Nmax,              &
             Mrankcs(ics-1), Nrankcs(ics-1), Nmaxcs(ics-1), a, NmaxAL, NmaxAL)
        call product_matrix_vector1 (2*Nmax, 2*Nmaxcs(ics-1), t, a, 2*NmaxAL,       &
             2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL)
        call extend_matrix3 (Nmax, Nmaxcs(ics-1), Nl, Nc, b, NmaxAL, NmaxAL, aa,    &
             Nmaxmax, Nmaxmax)      
        call MatTransAB_mn_m1n1 (3, k, x, y, z, Mrankcs(ics-1), Nrankcs(ics-1),     &
             Nmaxcs(ics-1), Mrank, Nrank, Nmax, a, NmaxAL, NmaxAL)
        call product_matrices1 (2*Nmaxcs(ics-1), 2*Nmaxcs(ics-1), 2*Nmax, tt,       &
             2*Nmaxcs(ics-1), 2*Nmaxcs(ics-1), a, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL,  &
             2*NmaxAL)       
        call extend_matrix3 (Nmaxcs(ics-1), Nmax, Nc, Nl, b, NmaxAL, NmaxAL, aa,    &
             Nmaxmax, Nmaxmax)      
      end if  
      call MatTransAB_mn_m1n1 (1, k, -x, -y, -z, Mrank, Nrank, Nmax, Mrankcs(ics),  &
           Nrankcs(ics), Nmaxcs(ics), a, NmaxAL, NmaxAL)
      call product_matrix_vector1 (2*Nmax, 2*Nmaxcs(ics), t, a, 2*NmaxAL, 2*NmaxAL, &
           b, 2*NmaxAL, 2*NmaxAL)
      call extend_matrix4 (ipart, Nmax, Nmaxcs(ics), Nmaxmax, Nl, b, NmaxAL,        &
           NmaxAL, bb, Nmaxmax, NmaxAL)                  
      Nl = Nl + 2 * Nmax
      if (PrnProgress) call write_progress (.false., ipart, 2*Npartcs(ics)+1)      
    end do
    if (ics > 1) then
      call identity_transformation (Mrankcs(ics-1), Nrankcs(ics-1), Nmaxcs(ics-1),  &
           Mrankcs(ics), Nrankcs(ics), Nmaxcs(ics), a, NmaxAL, NmaxAL)
      call product_matrices1 (2*Nmaxcs(ics-1), 2*Nmaxcs(ics-1), 2*Nmaxcs(ics), tt,  &
           2*Nmaxcs(ics-1), 2*Nmaxcs(ics-1), a, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL,    &
           2*NmaxAL)
      call extend_matrix4 (2, Nmaxcs(ics-1), Nmaxcs(ics), Nmaxmax, Nl, b, NmaxAL,   &
           NmaxAL, bb, Nmaxmax, NmaxAL)             
    end if
    call LU_SYSTEM_DIRECT (aa, 2*Nmaxmax, 2*Nmaxmax, bb, 2*Nmaxmax, 2*NmaxAL,       &
         2*Nmaxmax, 2*Nmaxcs(ics)) 
    if (PrnProgress) call write_progress (.false., Npartcs(ics)+1, 2*Npartcs(ics)+1)              
    Nl = 0  
    do ipart = 1, Npartcs(ics)    
      x = xp(ics,ipart)
      y = yp(ics,ipart)
      z = zp(ics,ipart)         
      call extract_matrix2 (Nmax, Nmaxcs(ics), Nl, b, NmaxAL, NmaxAL, bb, Nmaxmax,  &
           NmaxAL)
      call MatTransAB_mn_m1n1 (1, k, x, y, z, Mrankcs(ics), Nrankcs(ics),           &
           Nmaxcs(ics), Mrank, Nrank, Nmax, a, NmaxAL, NmaxAL)                          
      call product_matrices2 (ipart, 2*Nmaxcs(ics), 2*Nmax, 2*Nmaxcs(ics), a,       &
           2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL, c, 2*NmaxAL, 2*NmaxAL)                        
      Nl = Nl + 2 * Nmax
      if (PrnProgress) call write_progress (.false., Npartcs(ics)+1+ipart,          &
                            2*Npartcs(ics)+1)
    end do   
    if (ics > 1) then
      call extract_matrix2 (Nmaxcs(ics-1), Nmaxcs(ics), Nl, b, NmaxAL, NmaxAL,      &
           bb, Nmaxmax, NmaxAL)
      call identity_transformation (Mrankcs(ics), Nrankcs(ics), Nmaxcs(ics),        &
           Mrankcs(ics-1), Nrankcs(ics-1), Nmaxcs(ics-1), a, NmaxAL, NmaxAL)
      call product_matrices2 (2, 2*Nmaxcs(ics), 2*Nmaxcs(ics-1), 2*Nmaxcs(ics),     &
           a, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL, c, 2*NmaxAL, 2*NmaxAL)
      deallocate (tt)
      allocate (tt(2*Nmaxcs(ics),2*Nmaxcs(ics)))
    end if
    call copy_matrix (2*Nmaxcs(ics), 2*Nmaxcs(ics), c, 2*NmaxAL, 2*NmaxAL, tt,      &
         2*NmaxAL, 2*NmaxAL)
    deallocate (aa, bb, a, b, c)
  end do 
  deallocate (t)
  allocate (t(2*Nmaxcs(Ncs)))
  allocate (a(2*Nmaxcs(Ncs),2*Nmaxcs(Ncs)), b(2*Nmaxcs(Ncs),2*Nmaxcs(Ncs)))  
  if (PrnProgress)                                                                  &
      print "(/,2x,'progress of main calculation for the global T matrix;')"                    
  call vector_Q_sphere (3, 3, ks, Rcs(Ncs), ind_refC, Mrankcs(Ncs), Nrankcs(Ncs),   &
       Nmaxcs(Ncs), t)
  call product_matrix_vector1 (2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs), t, tt, 2*Nmaxcs(Ncs),  &
       2*Nmaxcs(Ncs), a, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs))
  call vector_Q_sphere (3, 1, ks, Rcs(Ncs), ind_refC, Mrankcs(Ncs), Nrankcs(Ncs),   &
       Nmaxcs(Ncs), t)
  call sum_diagonal_elements (2*Nmaxcs(Ncs), t, a, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs)) 
  if (PrnProgress) call write_progress (.false., 1, 3)
  call vector_Q_sphere (1, 3, ks, Rcs(Ncs), ind_refC, Mrankcs(Ncs), Nrankcs(Ncs),   &
       Nmaxcs(Ncs), t)
  call product_matrix_vector1 (2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs), t, tt, 2*Nmaxcs(Ncs),  &
       2*Nmaxcs(Ncs), b, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs))
  call vector_Q_sphere (1, 1, ks, Rcs(Ncs), ind_refC, Mrankcs(Ncs), Nrankcs(Ncs),   &
       Nmaxcs(Ncs), t)
  call sum_diagonal_elements (2*Nmaxcs(Ncs), t, b, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs))  
  if (PrnProgress) call write_progress (.false., 2, 3)
  call LU_SYSTEM (a, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs), b, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs), &
       2*Nmaxcs(Ncs))
  if (PrnProgress) call write_progress (.false., 3, 3)
  call minus_matrix (2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs), b, 2*Nmaxcs(Ncs), 2*Nmaxcs(Ncs))
  call write_FileTmat (Nmaxcs(Ncs), Nmaxcs(Ncs), b)     
  close (unit = iTmat) 
  call write_InfoFileTmat (FileTmat, Mrankcs(Ncs), Nrankcs(Ncs), .false., .false.,  &
      .false.)
  call ScatCharact (ks, FileTmat, Mrankcs(Ncs), Nrankcs(Ncs), .false., .false.,     &
      .false.)   
  print "(/,2x,'T matrix is stored in ',a50)", FileTmat
  print "(  2x,'The dimensions of the T matrix are given by:')"                          
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", Nrankcs(Ncs)
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrankcs(Ncs)                       
  deallocate (a, b, tt, t, Nmaxcs)      
end subroutine TMatrix_Nrank_MrankINHOMSPHREC


