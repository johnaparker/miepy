subroutine TSPHERE
!------------------------------------------------------------------------------------
! 1. General Considerations                                                         !
! --------------------------                                                        !
! TSPHERE is a routine for computing the T vector and the scattering                !
! characteristics of concentrically layered spheres. The spherical particle         !
! consists of Npart layers, each characterized by arbitrary but constant values     !
! of electric permittivity and magnetic permeability. Each layer is bounded by two  ! 
! surfaces. The enclosing (exterior) surface will be referred to as the layer       !
! surface. The first layer surface is the largest surface and corresponds to the    !
! host sphere.                                                                      !
!                                                                                   !
! The T matrix of a concentrically layered sphere is diagonal and is completely     !
! described by the diagonal elements corresponding to the azimuthal mode m = 1.     !
! Therefore, for concentrically layered spheres, the T matrix simplifies to a       !
! T vector.                                                                         !
!                                                                                   !
! The domain of analysis is divided into homogeneous layers. The electric and       !
! magnetic fields of a specific layer i are approximated by using localized vector  !
! spherical wave functions, and Nrankp(i) will be referred to as the maximum        !
! expansion order of the layer i.                                                   !
!                                                                                   !
! In contrast to the program for computing the T matrix of a layered particle,      !
! the present code assumes that the maximum expansion orders of all layers are      !
! identically, i.e., Nrankp(1) = Nrankp(2) = ...= Nrankp(Npart) = Nrank.            !
!                                                                                   !
! 2. Convergence Test                                                               !
! --------------------                                                              !
! For convergence tests, the incident wave is assumed to be a vector plane wave     !
! traveling along the Z-axis of the global coordinate system and the scattering     !
! characteristics  are computed in the azimuthal plane  phi = 0° . Note that only   !
! two azimuthal modes m = 1 and m = - 1 are required to compute the scattering      !
! characteristics.                                                                  !
!                                                                                   !
! The convergence test is carried out over the maximum expansion order Nrank.       !
! The scattering problem is solved for Nrank and Nrank - 1 and the normalized       !
! differential scattering cross section (DSCS) will be checked at 20° increments    !
! for convergence within epsNrank tolerance. If the calculated results converge     !
! within this tolerance at 80% of the scattering angles, then convergence is        !
! achieved.                                                                         ! 
!                                                                                   !
! An estimate of Nrank is given by Wiscombe's truncation limit criterion            !
! [W. J. Wiscombe, Improved Mie scattering algorithms, Applied Optics, 19,          !
! 1505-1509, 1980], i.e.,                                                           !
!                                                                                   !                     
!                NrankW = int(xp + 4.05 * xp**0.33 + 2),                            !
!                                                                                   !
! where xp is the size parameter of the layer, xp = k * r(1), k is the wave         !
! number and r(1) is the radius of the host sphere.                                 !
!                                                                                   !
! The convergence tests over Nrank can be switched off by setting the logical       !
! variable DoConvTest to false. In this case, the value of Nrank must be specified  !
! in the input file.                                                                !
!                                                                                   !
! 3. Input Parameters                                                               !
! --------------------                                                              !
! The parameters specified in the input file "/INPUTFILES/InputSPHERE.dat" are      !
! listed below.                                                                     !
!                                                                                   !
! - wavelength (real) - wavelength of the incident light in vacuo.                  !
!                                                                                   !
! - ind_refMed (real) - refractive index of the ambient medium (nonabsorbing        !
!   dielectric medium).                                                             !
!                                                                                   !
! - Npart (integer) - number of spherical layers.                                   !
!                                                                                   !
! - anorm (real) - characteristic length of the layered sphere which is used        !
!   to normalize the differential scattering cross sections.                        !
!                                                                                   !
! The next parameters (specified in a sequence of group statements) characterize     !
! each layer.                                                                       !
! - ind_refPartRel (complex) - relative refractive index of the actual spherical    !
!   layer with respect to the ambient medium. The imaginary part of the relative    !
!   refractive index must be zero for nonabsorbing particles and positive for       !
!   absorbing particles.                                                            ! 
!                                                                                   !
! - rPart (real) - radius of the actual spherical layer.                            !
!                                                                                   !
! - DoConvTest (logical) - if DoConvTest = t, the interactive convergence test      !
!   over Nrank is invoked. An estimate of Nrank for all layers is provided by       !
!   Wiscombe's truncation limit criterion. If DoConvTest = f, the values of Nrank   !
!   must be supplied in the input file.                                             !
!                                                                                   !
! - Nrank (integer) - maximum expansion order for the host sphere. This parameter   !
!   is used if the convergence test is not performed (DoConvTest = f)               !  
!                                                                                   !
! - epsNrank (real) - error tolerance for the expansion order test.                 !
!                                                                                   !
! - FileTmat (character(80)) - name of the file to which the T matrix is written.   !
!                                                                                   !
! - PrnProgress (logical) - if PrnProgress = t, the progress of calculation         !
!   is printed.                                                                     !
!                                                                                   !
! Note that all input parameters specifying lengths must be provided in the same    !
! units of lengths.                                                                 ! 
!------------------------------------------------------------------------------------
  use parameters
  use allocation, only: rp, ind_refp
  implicit none 
  integer       :: Npart, Nrank                                   
  real(O)       :: k, ind_refMed, wavelength, anorm, snorm, epsNrank
  character(80) :: FileTmat
  logical       :: DoConvTest, PrnProgress  
! -----------------------------------------------------------------------------------
!                                Read the input file                                ! 
! -----------------------------------------------------------------------------------
  call readinputSPHERE ( wavelength, ind_refMed, Npart, anorm, DoConvTest, Nrank,   &
       epsNrank, FileTmat, PrnProgress, k, snorm)  
! -----------------------------------------------------------------------------------
!                                      Main                                         !
! ----------------------------------------------------------------------------------- 
  if (DoConvTest) then            
    open (unit = iOutput, file = FileOutput, status = "replace") 
    call printinputSPHERE (wavelength, ind_refMed, Npart, rp, ind_refp, epsNrank,   &
         anorm)
    call convergence_NrankSPHERE (k, Npart, ind_refp, rp, snorm, Nrank, epsNrank,   &
         FileTmat, PrnProgress)      
    close (unit = iOutput) 
  else
    call Tvector_NrankSPHERE (k, Npart, ind_refp, rp, Nrank, FileTmat, PrnProgress) 
  end if
  deallocate (rp,ind_refp)
end subroutine TSPHERE
!************************************************************************************
subroutine readinputSPHERE ( wavelength, ind_refMed, Npart, anorm, DoConvTest,      &
           Nrank, epsNrank, FileTmat, PrnProgress, k, snorm)    
  use parameters
  use derived_parameters
  use allocation, only: rp, ind_refp
  implicit none 
  integer       :: Npart, Nrank, NrankW, ipart, ios                                   
  real(O)       :: k, ind_refMed, rPart, wavelength, anorm, xpart, snorm,           &
                   epsNrank, x
  complex(O)    :: ind_refPartRel
  character(80) :: FileTmat, string
  logical       :: DoConvTest, PrnProgress, XFindPar
! -----------------------------------------------------------------------------------
!                          Read the input file FileInputSPHERE                      ! 
! -----------------------------------------------------------------------------------    
  call DrvParameters
  open (unit = iInputSPHERE, file = FileInputSPHERE, status = "old",                &
        position = "rewind")   
  wavelength = 0.1_O * 2._O * Pi  
  ind_refMed = 1._O
  string     = 'OptProp'
  if (XFindPar (iInputSPHERE, string)) then
    read (iInputSPHERE, *, iostat = ios) wavelength
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable wavelength;')"
      stop
    end if
    read (iInputSPHERE, *, iostat = ios) ind_refMed
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
  Npart  = 1
  anorm  = 1._O
  string = 'GeomProp'
  if (XFindPar (iInputSPHERE, string)) then
    read (iInputSPHERE, *, iostat = ios) Npart
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Npart;')"
      stop
    end if
    read (iInputSPHERE, *, iostat = ios) anorm
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
  allocate (ind_refp(Npart), rp(Npart))
  do ipart = 1, Npart
    ind_refPartRel = (1.5_O,0._O)
    string = 'OptPartProp'
    if (XFindPar (iInputSPHERE, string)) then
      read (iInputSPHERE, *, iostat = ios) ind_refPartRel
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable ind_refPartRel;')"
        print "(  2x,'for the particle ',i3,';')", ipart
        stop
      end if             
    else
      print "(/,2x,'Group name OptPartProp not found;')"
      stop  
    end if
    call check_ind_ref1 (ipart, ind_refPartRel)
    ind_refp(ipart) = ind_refPartRel
    rPart  = 1._O
    string = 'GeomPartProp'
    if (XFindPar (iInputSPHERE, string)) then
      read (iInputSPHERE, *, iostat = ios) rPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable rPart;')"
        print "(  2x,'for the particle ',i3,';')", ipart
        stop
      end if             
    else
      print "(/,2x,'Group name GeomPartProp not found;')"
      stop  
    end if
    rp(ipart) = rPart
  end do 
!
  DoConvTest = .true.
  string     = 'ConvTest'
  if (XFindPar (iInputSPHERE, string)) then
    read (iInputSPHERE, *, iostat = ios) DoConvTest
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
    print "(/,2x,'Convergence Test for a Spherical Layered Particle')"
    print "(  2x,'-------------------------------------------------')"
  else
    print "(/,2x,'T-vector calculation for a Spherical Layered Particle')"
    print "(  2x,'-----------------------------------------------------')"
  end if
!
  x = k * rp(1)
  NrankW = int(x + 4.05_O * x**0.33_O + 2._O)
  if (.not. DoConvTest) then
    Nrank  = 20
    string = 'NrankSph'
    if (XFindPar (iInputSPHERE, string)) then
      read (iInputSPHERE, *, iostat = ios) Nrank
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Nrank;')"
        stop
      end if          
    else
      print "(/,2x,'Group name NrankSph not found;')"
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
!
  epsNrank = 5.e-2_O
  string   = 'Errors'
  if (XFindPar (iInputSPHERE, string)) then
    read (iInputSPHERE, *, iostat = ios) epsNrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsNrank;')"
      stop
    end if          
  else
    print "(/,2x,'Group name Errors not found;')"
    stop  
  end if  
!
  FileTmat = 'T.dat'
  string   = 'Tmat'
  if (XFindPar (iInputSPHERE, string)) then
    read (iInputSPHERE, *, iostat = ios) FileTmat
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
  string = 'PrintProgress'
  if (XFindPar (iInputSPHERE, string)) then
    read (iInputSPHERE, *, iostat = ios) PrnProgress
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable PrnProgress;')"
      stop
    end if          
  else
    print "(/,2x,'Group name PrintProgress not found;')"
    stop  
  end if
  close (unit = iInputSPHERE)
end subroutine readinputSPHERE
!***********************************************************************************
subroutine printinputSPHERE (wavelength, ind_refMed, Npart, r, ind_ref, epsNrank,   &
           anorm)
  use parameters
  implicit none  
  integer       :: Npart,i
  real(O)       :: wavelength, ind_refMed, r(Npart), epsNrank, anorm                   
  complex(O)    :: ind_ref(Npart)  
!
  write (iOutput,"(/,2x,'Input Parameters of the Scattering Problem:',/)")
  write (iOutput,"(2x,'wavelength of the free space, wavelength = ',1pe13.4,';')")  &
         wavelength
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'refractive index of the ambient medium, ind_refMed = ', ind_refMed, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'characteristic length of the particle, anorm = ', anorm, ';'
  write (iOutput,*) 
  write (iOutput,"(2x,'number of layers, Nlayer = ',i3,';')") Npart 
  do i = 1, Npart      
    write (iOutput,"(2x,'radius of the layer ',i3,': r = ',1pe10.3,';')") i,r(i)
    write (iOutput,"(2x, a, i3, a, 1pe10.3, ',', 1pe10.3, a)")                      &
   'relative refractive index of the layer ', i, ', ind_refRel = (',                &
    ind_ref(i), ');'           
  end do
  write (iOutput,*)     
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum expansion order tolerance, epsNrank = ', epsNrank, '.'              
  write (iOutput,"(/)")  
end subroutine printinputSPHERE
! **********************************************************************************
subroutine convergence_NrankSPHERE (k, Npart, ind_ref, r, snorm, Nrank, epsNrank,   &
           FileTmat, PrnProgress)   
  use parameters               
  implicit none 
  integer       :: Npart, Nrank      
  real(O)       :: k, r(Npart), snorm, epsNrank
  complex(O)    :: ind_ref(Npart)
  character(80) :: FileTmat
  logical       :: PrnProgress
!  
  integer       :: Mstart, Nteta, Mrank, Nmax, Nmaxmax, i, m, NthetaConv
  real(O)       :: alfa, beta, gama, alfap, phiGS, tetaGI, phiGI, Cscat, Qscat,     &
                   Cext, Qext
  real(O),allocatable    :: h(:), v(:), oldh(:), oldv(:)
  complex(O),allocatable :: t(:), c(:), c1(:), cc(:)             
!
  phiGS  = 0._O
  Nteta  = 10
  alfa   = 0._O        
  beta   = 0._O
  gama   = 0._O
  alfap  = Pi / 4._O
  tetaGI = 0._O
  phiGI  = 0._O
  Mstart = 1
  Mrank  = 1
  Nmax   = Nrank 
  m = 1 
  open (unit = iTmat, file = FileTmat, status = 'replace')
  call write_HeadFileTmatVct (Nrank) 
  call write_TypeConvHead (2) 
  allocate (t(2*Nrank), c(2*Nrank), c1(2*Nrank))
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  allocate (cc(2*Nmaxmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))
  if (PrnProgress) call write_progress (.true., 1, 5)             
! --- Nrank configuration ---
  call coefficients_fg_m_coated (k, Npart, r, ind_ref, m, Nrank, Nmax, t) 
  if (PrnProgress) call write_progress (.false., 2, 5)  
  call write_FileTmatVct (Nrank, t)
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,       &
       Nmax, c)             
  call product_vector_vector (2*Nmax, t, c, c1)
  call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)         
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap,-m, Nrank,       &
       Nmax, c)
  call product_vector_vector (2*Nmax, t, c, c1)
  call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax) 
  if (PrnProgress) call write_progress (.false., 3, 5)                                             
  call DSCS (cc, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm,   &
      .false.,.true., h, v)
  call CQscat (cc, Mrank, Nrank, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,           &
       alfap, k, snorm, Cext, Qext)   
  do i = 1, Nteta
    oldh(i) = h(i)
    oldv(i) = v(i)
  end do  
  write (iOutput,"(17x,'Nrank = ',i3,/)") Nrank
  call  write_DSCS (Nteta,.false., h, v)
  call  write_Effic (Qscat, Qext)
  close (unit = iTmat)  
! --- (Nrank - 1) configuration ---  ...  
  call vector_Nrank_m (Nmax, t)
  if (PrnProgress) call write_progress (.false., 4, 5)
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,       &
       Nmax, c)             
  call product_vector_vector (2*Nmax, t, c, c1)
  call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)         
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,      &
       Nmax, c)
  call product_vector_vector (2*Nmax, t, c, c1)
  call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)                                                  
  if (PrnProgress) call write_progress (.false., 5, 5)
  call DSCS (cc, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm,   &
      .false.,.true., h, v)
  call CQscat (cc, Mrank, Nrank, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,           &
       alfap, k, snorm, Cext, Qext) 
  call delta_DSCS (Nteta, h, v, oldh, oldv, epsNrank, NthetaConv)  
  write (iOutput,"(17x,'Nrank = ',i3,' ---',/)") Nrank - 1
  call write_DSCS (Nteta,.false., h, v)
  call write_Effic (Qscat, Qext)
  call write_MrankConvRes (NthetaConv, epsNrank)
  if (NthetaConv >= int(0.8*Nteta)) then
    print "(/,2x,'Convergence criterion for Nrank is satisfied;')"    
  else
    print "(/,2x,'Convergence criterion for Nrank is not satisfied;')"
  end if
  Mrank = Nrank 
  call write_InfoFileTmat (FileTmat, Mrank, Nrank, .true., .true., .false.) 
  call ScatCharact (k, FileTmat, Mrank, Nrank, .true., .true., .false.)  
  print "(/,2x,'The T vector is stored in ',a50)", FileTmat
  print "(  2x,'The dimension of the T vector is given by:')"
  print "(  2x,'- maximum expansion order, Nrank = ',i3,';')", Nrank                      
  deallocate (t, c, c1, cc, h, v, oldh, oldv)   
end subroutine convergence_NrankSPHERE 
! **********************************************************************************
subroutine Tvector_NrankSPHERE (k, Npart, ind_ref, r, Nrank, FileTmat, PrnProgress)   
  use parameters                    
  implicit none 
  integer       :: Npart, Nrank    
  real(O)       :: k, r(Npart)
  complex(O)    :: ind_ref(Npart)
  character(80) :: FileTmat
  logical       :: PrnProgress
!  
  integer       :: Mrank, m
  complex(O),allocatable :: t(:)
!  
  m = 1 
  open (unit = iTmat, file = FileTmat, status = 'replace')
  call write_HeadFileTmatVct (Nrank) 
  allocate (t(2*Nrank))  
  if (PrnProgress) call write_progress (.true., 1, 2)             
  call coefficients_fg_m_coated (k, Npart, r, ind_ref, m, Nrank, Nrank, t) 
  if (PrnProgress) call write_progress (.false., 2, 2)    
  call write_FileTmatVct (Nrank, t)  
  close (unit = iTmat)  
  Mrank = Nrank 
  call write_InfoFileTmat (FileTmat, Mrank, Nrank, .true., .true., .false.) 
  call ScatCharact (k, FileTmat, Mrank, Nrank, .true., .true., .false.)  
  print "(/,2x,'The T vector is stored in ',a50)", FileTmat
  print "(  2x,'The dimension of the T vector is given by:')"
  print "(  2x,'- maximum expansion order, Nrank = ',i3,';')", Nrank                      
  deallocate (t)      
end subroutine Tvector_NrankSPHERE     
