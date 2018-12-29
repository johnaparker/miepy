! **********************************************************************************
! *           SCATTERING CHARACTERISTICS: RANDOMLY ORIENTED PARTICLES              *
! *   --------------------------------------------------------------------------   *
! *   Partial list of subroutines:                                                 *
! *     DiffScatCrossSectAvrg,      ScatteringMatrix,           HovTest,           *     
! *     ExtendedScatteringMatrix,   AsymmetryParameter,         SijSCpqEff,        *
! *     AvCscatCextEff,             tnScatEff,                  tnScatAxsymEff,    *
! *     tnExtEff,                   tnExtAxsymEff,              tnExtOffEff,       *
! *     RandomExtinctionMatrix,     AvCscatCextV,               MNinfiniteMatrix,  * 
! *     IdentifyMatrixElement,      IdentifyMatrixElementAxsym, IdentifyIndex,     *
! *     Hv,                         dimensionTV,                formTV,            *  
! *     SijSCpqInt,                 Smatrix                                        *
! **********************************************************************************
subroutine DiffScatCrossSectAvrg (MirorSym, SS, Ntheta, anorm, epol_beta,           &
           epol_alpha, phiGS, normalized, Cscat, Cext, Qscat, Qext,                 &
           CscatV, CextV, QscatV, QextV)
! -----------------------------------------------------------------------------------
! The routine computes the average differential scattering cross section at a set   !
! of Ntheta scattering angles in the azimuthal plane phiGS. The incident            !
! polarization state is described by the polarization vector                        !
!                 epol = epol_beta * e_beta + epol_alpha * e_alpha.                 !
! The average differential scattering cross section for parallel and perpendicular  !
! polarizations are given by                                                        !
!           <h> = <|FS_theta|**2>                                                   !    
!               = <|S(1,1)|**2> * |e_beta|**2 + <|S(1,2)|**2> * |e_alpha|**2        !
!                 + 2*Re{ <S(1,1)*conjg(S(1,2))> * e_beta * conjg(e_alpha) }        !
! and                                                                               !                                             
!           <v> = <|FS_phi|**2>                                                     !
!               = <|S(2,1)|**2> * |e_beta|**2 + <|S(2,2)|**2> * |e_alpha|**2        !                   
!                 + 2*Re{ <S(2,1)*conjg(S(2,2))> * e_beta * conjg(e_alpha) }        !
! where                                                                             !
!           e_beta  =   epol_beta * cos(phiGS) + epol_alpha * sin(phiGS),           !
!           e_alpha = - epol_beta * sin(phiGS) + epol_alpha * cos(phiGS).           !
! For isotropic and mirror-symmetric media                                          !
!              Re{ <S(1,1)*conjg(S(1,2))> } = 0                                     !    
! and                                                                               !
!              Re{ <S(2,1)*conjg(S(2,2))> } = 0.                                    !
! -----------------------------------------------------------------------------------
  use parameters
  implicit none 
  integer       :: Ntheta, i
  real(O)       :: anorm, phiGS, Cscat, Cext, Qscat, Qext, CscatV, CextV, QscatV,   &
                   QextV, eb2, ea2, fact, h, v, theta, norm     
  complex(O)    :: SS(10,Ntheta), epol_beta, epol_alpha, e_beta, e_alpha, eba,      &
                   S11S12eba, S21S22eba 
  logical       :: MirorSym, normalized
!  
  write (iDSCS,"(/,2x,'Results:',/)")
  write (iDSCS,"(  2x,'Cross Sections and Efficiencies:')")      
  if (.not.MirorSym) then
    write (iDSCS,"(2x, a, 1pe13.4, 5x, a, 1pe13.4)")                                &
          '<Cscat>_I = ', Cscat,  '<Qscat>_I = ', Qscat 
    write (iDSCS,"(2x, a, 1pe13.4, 5x, a, 1pe13.4)")                                &
          '<Cscat>_V = ', CscatV, '<Qscat>_V = ', QscatV 
    write (iDSCS,*)	    
    write (iDSCS,"(2x, a, 1pe13.4, 5x, a, 1pe13.4)")                                &
          '<Cext>_I  = ', Cext,   '<Qext>_I  = ', Qext
    write (iDSCS,"(2x, a, 1pe13.4, 5x, a, 1pe13.4)")                                &
          '<Cext>_V  = ', CextV,  '<Qext>_V  = ', QextV	    
    write (iDSCS,*)           	    	                          
  else   
    write (iDSCS,"(2x, a, 1pe13.4, 5x, a, 1pe13.4)")                                &
          '<Cscat> = ', Cscat,  '<Qscat> = ', Qscat
    write (iDSCS,"(2x, a, 1pe13.4, 5x, a, 1pe13.4)")                                &
          '<Cext>  = ', Cext,   '<Qext>  = ', Qext       
  end if 
  write (iDSCS,*)        
  write (iDSCS,"(2x,'Differential Scattering Cross Section:')")      
  if (normalized) then
    write (iDSCS,"(2x, a, /, 2x, a, 9x, a, 8x, a, /)")                              &
   'normalized average  differential scattering cross section', 'theta',            &
   'parallel', 'perpendicular'                  
  else
    write (iDSCS,"(8x, a, /, 2x, a, 9x, a, 8x, a, /)")                              &
   'average differential scattering cross section', 'theta',                        &
   'parallel', 'perpendicular'                            
  end if  
  norm       =   sqrt ( abs(epol_beta)**2 + abs(epol_alpha)**2 )  
  epol_beta  =   epol_beta  / norm
  epol_alpha =   epol_alpha / norm       
  e_beta     =   epol_beta * cos(phiGS) + epol_alpha * sin(phiGS)  
  e_alpha    = - epol_beta * sin(phiGS) + epol_alpha * cos(phiGS)    
  eb2  = abs(e_beta )**2
  ea2  = abs(e_alpha)**2
  eba  = e_beta * conjg(e_alpha)               
  fact = Pi * anorm * anorm
  fact = 1._O / fact
  do i = 1, Ntheta    
    theta = real(i - 1,O) * 180._O / real(Ntheta - 1,O)
    if (.not.MirorSym) then
      S11S12eba = SS(2,i) * eba
      h = real(SS(1,i),O) * eb2 + real(SS( 5,i),O) * ea2 +                          &
          2._O * real(S11S12eba,O)	    
      S21S22eba = SS(9,i) * eba	    
      v = real(SS(8,i),O) * eb2 + real(SS(10,i),O) * ea2 +                          &
          2._O * real(S21S22eba,O)          
    else
      h = real(SS(1,i),O) * eb2 + real(SS( 5,i),O) * ea2             
      v = real(SS(8,i),O) * eb2 + real(SS(10,i),O) * ea2        
    end if               
    if (normalized) then
      h = h * fact
      v = v * fact                                          
    end if  
    write (iDSCS,"(1x,f6.2,5x,1pe13.4,5x,1pe13.4)") theta, h, v
  end do
end subroutine DiffScatCrossSectAvrg 
! **********************************************************************************
subroutine ScatteringMatrix (MirorSym, theta, Ntheta, ZE, Z)
!-----------------------------------------------------------------------------------
! The routine computes the scattering matrix Z at a specified scattering angle     !
! theta. ZE is the extended scattering matrix and Ntheta is the number of zenith   !
! angles (at which the extended scattering matrix is computed) in the azimuthal    !
! plane phi = 0°.                                                                  !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: Ntheta
  real(O)    :: theta, ZE(4,4,Ntheta), Z(4,4)
  logical    :: MirorSym   
!
  integer    :: i
  real(O)    :: Z11, Z12, Z13, Z14, Z22, Z23, Z24, Z33, Z34, Z44    
  real(O),allocatable :: thetaV(:), ZV(:) 
! 
  allocate (thetaV(Ntheta), ZV(Ntheta))
  do i = 1, Ntheta
    thetaV(i) = real(i - 1,O) * Pi / real(Ntheta - 1,O)            
  end do  
  do i = 1, Ntheta
    ZV(i) = ZE(1,1,i)
  end do      
  call Interp (Ntheta, thetaV, ZV, theta, Z11)                          
  do i = 1, Ntheta
    ZV(i) = ZE(1,2,i)
  end do
  call Interp (Ntheta, thetaV, ZV, theta, Z12)  
!  
  do i = 1, Ntheta
    ZV(i) = ZE(2,2,i)
  end do   
  call Interp (Ntheta, thetaV, ZV, theta, Z22)      
!  
  do i = 1, Ntheta
    ZV(i) = ZE(3,3,i)
  end do
  call Interp (Ntheta, thetaV, ZV, theta, Z33)    
  do i = 1, Ntheta
    ZV(i) = ZE(3,4,i)
  end do
  call Interp (Ntheta, thetaV, ZV, theta, Z34)    
!  
  do i = 1, Ntheta
    ZV(i) = ZE(4,4,i)
  end do
  call Interp (Ntheta, thetaV, ZV, theta, Z44)  
!
  if (.not.MirorSym) then
    do i = 1, Ntheta
      ZV(i) = ZE(1,3,i)
    end do    
    call Interp (Ntheta, thetaV, ZV, theta, Z13)    
    do i = 1, Ntheta
      ZV(i) = ZE(1,4,i)
    end do
    call Interp (Ntheta, thetaV, ZV, theta, Z14)  
!
    do i = 1, Ntheta
      ZV(i) = ZE(2,3,i)
    end do
    call Interp (Ntheta, thetaV, ZV, theta, Z23)
    do i = 1, Ntheta
      ZV(i) = ZE(2,4,i)
    end do
    call Interp (Ntheta, thetaV, ZV, theta, Z24)         
  end if   
  Z(1,1) = Z11    
  Z(1,2) = Z12
  if (.not.MirorSym) then
    Z(1,3) = Z13
    Z(1,4) = Z14              
  else
    Z(1,3) = 0._O
    Z(1,4) = 0._O            
  end if     
!  
  Z(2,1) = Z12
  Z(2,2) = Z22
  if (.not.MirorSym) then
    Z(2,3) = Z23
    Z(2,4) = Z24              
  else
    Z(2,3) = 0._O
    Z(2,4) = 0._O            
  end if    
! 
  if (.not.MirorSym) then
    Z(3,1) = - Z13
    Z(3,2) = - Z23              
  else
    Z(3,1) = 0._O
    Z(3,2) = 0._O            
  end if   
  Z(3,3) = Z33
  Z(3,4) = Z34
! 
  if (.not.MirorSym) then
    Z(4,1) = Z14
    Z(4,2) = Z24              
  else
    Z(4,1) = 0._O
    Z(4,2) = 0._O            
  end if  
  Z(4,3) = - Z34
  Z(4,4) =   Z44
  deallocate (thetaV, ZV)
end subroutine ScatteringMatrix
! **********************************************************************************
subroutine HovTest (Z, FailHovTest, Nfail, FailMessage)
! ----------------------------------------------------------------------------------
! The routine performs the Hovenier and van der Mee test of the phase matrix.      !
! ----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer        :: Nfail
  real(O)        :: Z(4,4)
  logical        :: FailHovTest
  character(256) :: FailMessage(10)
!
  integer :: ifail
  real(O) :: tol, tmp 
!  
  ifail = 0
  tol   = 1.e-3_O
  if (Z(1,1) < 0._O) then
    FailHovTest = .true.      
    ifail = ifail + 1
    FailMessage (ifail) = 'the condition Z11 > 0 is not satisfied;'
  end if
  if (abs(Z(1,2)) > Z(1,1)) then
    FailHovTest = .true.
    ifail = ifail + 1
    FailMessage (ifail) = 'the condition |Z12| < Z11 is not satisfied;'
  end if
  if (abs(Z(2,2)) > Z(1,1)) then
    FailHovTest = .true.
    ifail = ifail + 1
    FailMessage (ifail) = 'the condition |Z22| < Z11 is not satisfied;'
  end if
  if (abs(Z(3,3)) > Z(1,1)) then
    FailHovTest = .true.
    ifail = ifail + 1
    FailMessage (ifail) = 'the condition |Z33| < Z11 is not satisfied;'
  end if
  if (abs(Z(3,4)) > Z(1,1)) then
    FailHovTest = .true.
    ifail = ifail + 1
    FailMessage (ifail) = 'the condition |Z34| < Z11 is not satisfied;'
  end if
  if (abs(Z(4,4)) > Z(1,1)) then
    FailHovTest = .true.
    ifail = ifail + 1
    FailMessage (ifail) = 'the condition |Z44| < Z11 is not satisfied;'
  end if   
  tmp = (Z(3,3) + Z(4,4))**2 + 4._O * Z(3,4)**2 - (Z(1,1) + Z(2,2))**2 +            &
         4._O * Z(1,2)**2 
  if (tmp > tol) then
    FailHovTest = .true.
    ifail = ifail + 1
    FailMessage (ifail) =                                                           &
   'the condition (Z33 + Z44)**2 + 4*Z34**2 < (Z11 + Z22)**2 - 4*Z12**2  is not satisfied;'
  end if
  tmp = abs(Z(3,3) - Z(4,4)) - Z(1,1) + Z(2,2)
  if (tmp > tol) then
    FailHovTest = .true.
    ifail = ifail + 1
    FailMessage (ifail) =                                                           &
   'the condition |Z33 - Z44| < Z11 - Z22 is not satisfied;'
  end if
  tmp = abs(Z(2,2) - Z(1,2)) - Z(1,1) + Z(1,2)
  if (tmp > tol) then
    FailHovTest = .true.
    ifail = ifail + 1
    FailMessage (ifail) =                                                           &
   'the condition |Z22 - Z12| < Z11 - Z12 is not satisfied;'
  end if
  tmp = abs(Z(2,2) + Z(1,2)) - Z(1,1) - Z(1,2)
  if (tmp > tol) then
    FailHovTest = .true.
    ifail = ifail + 1
    FailMessage (ifail) =                                                           &
   'the condition |Z22 + Z12| < Z11 + Z12 is not satisfied;'
  end if
  Nfail = ifail
end subroutine HovTest
! **********************************************************************************
subroutine ExtendedScatteringMatrix (MirorSym, Ntheta, SS, ZE)
!-----------------------------------------------------------------------------------
! The routine computes the scattering matrix ZE for all scattering angles in the   !
! azimuthal plane phi = 0°. SS is the average matrix <SpqSp1q1*> and Ntheta is the !
! number of zenith angles in the azimuthale plane phi = 0°.                        !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: Ntheta, i
  real(O)    :: ZE(4,4,Ntheta)
  complex(O) :: SS(10,Ntheta)
  complex(O) :: S11SC11, S11SC12, S11SC21, S11SC22, S12SC12, S12SC21, S12SC22,      &
                S21SC21, S21SC22, S22SC22, S21SC11, S21SC12, S22SC11, S22SC12,      &
                S22SC21   
  logical    :: MirorSym		
!
  do i = 1, Ntheta        
    S11SC11 = SS(1,i)                    
    S11SC22 = SS(4,i)        
    S12SC12 = SS(5,i)        
    S12SC21 = SS(6,i)            
    S21SC21 = SS(8,i)          
    S22SC22 = SS(10,i)
    if (.not.MirorSym) then
      S11SC12 = SS(2,i)
      S11SC21 = SS(3,i)
      S12SC22 = SS(7,i)
      S21SC22 = SS(9,i)                
    end if 
!         
    S21SC12 = conjg(S12SC21)
    S22SC11 = conjg(S11SC22)             
    if (.not.MirorSym) then
      S21SC11 = conjg(S11SC21)     
      S22SC12 = conjg(S12SC22) 
      S22SC21 = conjg(S21SC22)    
    end if     
!
    ZE(1,1,i) =   0.5_O * real(S11SC11 + S12SC12 + S21SC21 + S22SC22,O)     
    ZE(1,2,i) =   0.5_O * real(S11SC11 - S12SC12 + S21SC21 - S22SC22,O)    
    if (.not.MirorSym) then
      ZE(1,3,i) = - real (S11SC12 + S22SC21,O)
      ZE(1,4,i) = - aimag(S11SC12 - S22SC21)            
    else
      ZE(1,3,i) =   0._O
      ZE(1,4,i) =   0._O        
    end if         
!    
    ZE(2,1,i) =   ZE(1,2,i)
    ZE(2,2,i) =   0.5_O * real(S11SC11 - S12SC12 - S21SC21 + S22SC22,O)      
    if (.not.MirorSym) then
      ZE(2,3,i) = - real (S11SC12 - S22SC21,O)  
      ZE(2,4,i) = - aimag(S11SC12 + S22SC21)    
    else
      ZE(2,3,i) =   0._O
      ZE(2,4,i) =   0._O        
    end if     
!        
    if (.not.MirorSym) then
      ZE(3,1,i) = - real (S11SC21 + S22SC12,O)  
      ZE(3,2,i) = - real (S11SC21 - S22SC12,O)    
    else
      ZE(3,1,i) =   0._O
      ZE(3,2,i) =   0._O        
    end if     
    ZE(3,3,i) =   real (S11SC22 + S12SC21,O)  
    ZE(3,4,i) =   aimag(S11SC22 + S21SC12)        
!        
    if (.not.MirorSym) then
      ZE(4,1,i) = - aimag(S21SC11 + S22SC12)  
      ZE(4,2,i) = - aimag(S21SC11 - S22SC12)      
    else
      ZE(4,1,i) =   0._O
      ZE(4,2,i) =   0._O        
    end if     
    ZE(4,3,i) =   - ZE(3,4,i)
    ZE(4,4,i) =   real (S22SC11 - S12SC21,O) 
  end do
end subroutine ExtendedScatteringMatrix
! **********************************************************************************
subroutine AsymmetryParameter (MirorSym, Ntheta, ZE, Cscat, AsymPar, AsymParV)
!------------------------------------------------------------------------------------     
! The routine computes the asymmetry parameter for a randomly oriented particle.    ! 
! The relations are                                                                 !
!                        2*Pi     |Pi                                               ! 
!     <cos THETA>    = ---------- |  <F11(theta)> cos(theta)sin(theta) d(theta)     !
!                        <Cscat>  |0                                                !
!  and                                                                              !
!                        2*Pi     |Pi                                               !
!     <cos THETA >_V = ---------- |  <F14(theta)> cos(theta)sin(theta) d(theta)     !
!                        <Cscat>  |0                                                ! 
!------------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: Ntheta
  real(O)    :: ZE(4,4,Ntheta), Cscat, AsymPar, AsymParV
  logical    :: MirorSym 
!
  integer    :: i
  real(O)    :: atheta, btheta, theta, sum, sumV
  real(O),allocatable :: wtheta(:), xtheta(:)
!
  allocate (wtheta(Ntheta), xtheta(Ntheta))   
  atheta = 0._O
  btheta = Pi
  call Simpson (atheta, btheta, Ntheta, xtheta, wtheta)
  sum  = 0._O 
  sumV = 0._O
  do i = 1, Ntheta    
    theta = xtheta(i)                         
    sum  = sum  + wtheta(i) * sin(theta) * cos(theta) * ZE(1,1,i)
    if (.not.MirorSym) sumV = sumV + wtheta(i) * sin(theta) * cos(theta) * ZE(1,4,i)                           
  end do                     
  AsymPar = 2._O * Pi * sum / Cscat
  if (.not.MirorSym) then
    AsymParV = 2._O * Pi * sumV / Cscat  
  else
    AsymParV = 0._O  
  end if     
  deallocate (wtheta, xtheta)     
end subroutine AsymmetryParameter
! **********************************************************************************
! *               ANALYTICAL ORIENTATION AVERAGING PROCEDURE                       *
! **********************************************************************************
subroutine SijSCpqEff (MirorSym, wavenumber, axsym, chiral, FileTmat, Nrank, Mrank, &
           Nrankeff, Mrankeff, Ntheta, S, PrnProgress)
!-----------------------------------------------------------------------------------
! The routine computes the matrix <Spq Sp1q1*> by using an analytical averaging    !
! procedure.                                                                       !      
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer       :: Nrank, Mrank, Nrankeff, Mrankeff, Ntheta
  real(O)       :: wavenumber 
  complex(O)    :: S(10,Ntheta)
  character(80) :: FileTmat
  logical       :: MirorSym, axsym, chiral, PrnProgress
!      
  integer       :: i, j, n, nt, m, m1tl, m1l, u, dmp, m1, mp, m1p, n1inf, n1sup, n1,&
                   k, kt, nTV, m1t, mt, dimensionTV, IdentifyIndex, ntg, mtg,       &
                   dm, dmmin, dmmax, mpmin, mpmax, mmin, mmax, Nmax, Nmaxeff 
  real(O)       :: factr, m11t, Hv, factur, factmr
  complex(O)    :: T11, T12, T21, T22, fact, AXX, AXY, AYX, AYY, BXX, BXY, BYX, BYY,&
                   CXX, CXY, CYX, CYY, DXX, DXY, DYX, DYY, MTh, MPh, MThC, MPhC,    &
                   m1c, m1tc
  real(O),allocatable    :: C(:)
  complex(O),allocatable :: T(:,:), TV(:), M2inf(:,:), M3inf(:,:), SA(:,:,:,:,:),   &
                            SB(:,:,:,:,:), SC(:,:,:,:,:), SD(:,:,:,:,:),            &
                            A(:,:,:), B(:,:,:), Al(:,:), Bl(:,:)
!
  Nmax    = Nrank + Mrank * (2 * Nrank - Mrank + 1)   
  Nmaxeff = Nrankeff + Mrankeff * (2 * Nrankeff - Mrankeff + 1)   
  if (.not. axsym) then    
    open (unit = iTmat, file = FileTmat, status = 'unknown')       
    call read_HeadFileTmat (ntg, mtg)
    call check_dimensionMat (ntg, mtg, Nmax)              
    allocate (T(2*ntg,2*mtg))
    call read_FileTmat (ntg, mtg, T) 
    close (unit = iTmat)
  else 
    nTV = dimensionTV (Mrank, Nrank, chiral)
    allocate (TV(nTV))
    call formTV (FileTmat, chiral, Mrank, Nrank, TV, nTV)    
  end if                  
  allocate (M2inf(Ntheta,Nmaxeff), M3inf(Ntheta,Nmaxeff)) 
  call MNinfiniteMatrix (Mrankeff, Nrankeff, Nmaxeff, Ntheta, M2inf, M3inf)
  allocate (SA(Nrankeff,Nrankeff,-Mrankeff:Mrankeff,2,2),                           &
            SB(Nrankeff,Nrankeff,-Mrankeff:Mrankeff,2,2),                           &
            SC(Nrankeff,Nrankeff,-Mrankeff:Mrankeff,2,2),                           &
            SD(Nrankeff,Nrankeff,-Mrankeff:Mrankeff,2,2))
  allocate (A(Nrankeff,-Mrankeff:Mrankeff,2), B(Nrankeff,-Mrankeff:Mrankeff,2),     &
            Al(2,Nrankeff), Bl(2,Nrankeff))
  do n = 1, Nrankeff
    do nt = 1, Nrankeff
      do m = - Mrankeff, Mrankeff          
        do m1l = 1, 2
          do m1tl = 1, 2
            SA(n,nt,m,m1l,m1tl) = zero
            SB(n,nt,m,m1l,m1tl) = zero
            SC(n,nt,m,m1l,m1tl) = zero
            SD(n,nt,m,m1l,m1tl) = zero                 
          end do  
        end do    
      end do
    end do
  end do
  if (.not. axsym) then
    do u = 0, 2*Nrankeff
      if (PrnProgress) call write_progress (.false., u+1, 2*Nrankeff+3)
      factur = sqrt(real(2 * u + 1,O)) / wavenumber
      do dmp = - u, u
        do n = 1, Nrankeff
          do m = - Mrankeff, Mrankeff
            do m1l = 1, 2
              A(n,m,m1l) = zero
              B(n,m,m1l) = zero
            end do
          end do
        end do              
        do n = 1, Nrankeff
          allocate (C(0:n+u))
          n1inf = max(1, abs(n - u))                      
          n1sup = min(Nrankeff, n + u)
          do m1l = 1, 2
            do n1 = n1inf, n1sup                          
              Al(m1l,n1) = zero
              Bl(m1l,n1) = zero 
            end do
          end do                                                 
          mpmax = min(n, Mrankeff, Mrankeff + dmp)
          mpmin = min(n, Mrankeff, Mrankeff - dmp)
          do mp = - mpmin, mpmax 
            m1p = mp - dmp  
            call CouplingCoef (dmp, u, - mp, n, C)
            n1inf  = max(1, abs(n - u), abs(m1p))                         
            n1sup  = min(Nrankeff, n + u)
            factmr = Hv(mp) * Hv(m1p)                                       
            do n1 = n1inf, n1sup                                                                            
              call IdentifyMatrixElement (mp, n, m1p, n1, Nrank, Nmax, T, ntg,  &
                   mtg, T11, T12, T21, T22)                                  
              factr = factmr * C(n1) 
              do m1l = 1, 2
                m1 = -3 + 2 * m1l
                Al(m1l,n1) = Al(m1l,n1) + factr * (m1 * T11 + T12)
                Bl(m1l,n1) = Bl(m1l,n1) + factr * (m1 * T21 + T22)
              end do              
            end do
          end do                           
          do m1l = 1, 2
            m1 = - 3 + 2 * m1l
            mmax = min(n, Mrankeff, u + m1)
            mmin = min(n, Mrankeff, u - m1)
            do m = - mmin, mmax                     
              call CouplingCoef (m1 - m, u, m, n, C)
              n1inf  = max(1, abs(n - u))                     
              n1sup  = min(Nrankeff, n + u)                                     
              factmr = Hv(m) * Hv(m1)
              do n1 = n1inf, n1sup
                factr = factmr * factur / sqrt(real(2 * n1 + 1,O))                                      
                fact  = (- 1._O)**(n + n1) * im**(n1 - 1) * C(n1) * factr                                        
                A(n,m,m1l) = A(n,m,m1l) + fact * Al(m1l,n1)
                B(n,m,m1l) = B(n,m,m1l) + fact * Bl(m1l,n1)
              end do
            end do
          end do
          deallocate (C)
        end do                                                                  
        do m1l = 1, 2
          m1 = - 3 + 2 * m1l
          do m1tl = 1, 2
            m1t = - 3 + 2 * m1tl
            dmmax = min(u, Mrankeff - m1, Mrankeff - m1t)
            dmmin = min(u, Mrankeff + m1, Mrankeff + m1t)
            do dm = - dmmin, dmmax
              m  = dm + m1
              mt = dm + m1t
              do n = max(1,abs(m)), Nrankeff
                do nt = max(1,abs(mt)), Nrankeff
                  SA(n,nt,m,m1l,m1tl) = SA(n,nt,m,m1l,m1tl) +                       &
                                        A(n,m,m1l) * conjg(A(nt,mt,m1tl))
                  SB(n,nt,m,m1l,m1tl) = SB(n,nt,m,m1l,m1tl) +                       &
                                        A(n,m,m1l) * conjg(B(nt,mt,m1tl))
                  SC(n,nt,m,m1l,m1tl) = SC(n,nt,m,m1l,m1tl) +                       &
                                        B(n,m,m1l) * conjg(A(nt,mt,m1tl))
                  SD(n,nt,m,m1l,m1tl) = SD(n,nt,m,m1l,m1tl) +                       &
                                        B(n,m,m1l) * conjg(B(nt,mt,m1tl))           
                end do  
              end do
            end do
          end do
        end do                                                    
      end do
    end do
  else 
    do u = 0, 2*Nrankeff
      if (PrnProgress) call write_progress (.false., u+1, 2*Nrankeff+3)
      factur = sqrt(real(2 * u + 1,O)) / wavenumber
      do n = 1, Nrankeff
        do m = - Mrankeff, Mrankeff
          do m1l = 1, 2
            A(n,m,m1l) = zero
            B(n,m,m1l) = zero
          end do
        end do
      end do                  
      do n = 1, Nrankeff
        allocate (C(0:n+u))
        n1inf = max(1, abs(n - u))                        
        n1sup = min(Nrankeff, n + u)
        do m1l = 1, 2
          do n1 = n1inf, n1sup
            Al(m1l,n1) = zero
            Bl(m1l,n1) = zero 
          end do
        end do                                                    
        do mp = - min(n,Mrankeff), min(n,Mrankeff)              
          call CouplingCoef (0, u, - mp, n, C)
          n1inf = max(1, abs(n - u), abs(mp))                     
          n1sup = min(Nrankeff, n + u)                                  
          do n1 = n1inf, n1sup                
            call IdentifyMatrixElementAxsym (mp, n, n1, Nrank, chiral, TV, nTV,     &
                 T11, T12, T21, T22)                
            do m1l = 1, 2
              m1 = - 3 + 2 * m1l
              Al(m1l,n1) = Al(m1l,n1) + C(n1) * (m1 * T11 + T12)
              Bl(m1l,n1) = Bl(m1l,n1) + C(n1) * (m1 * T21 + T22)
            end do
          end do
        end do    
        do m1l = 1, 2
          m1 = - 3 + 2 * m1l
          mmax = min(n, Mrankeff, u + m1)
          mmin = min(n, Mrankeff, u - m1)
          do m = - mmin, mmax   
            call CouplingCoef (m1 - m, u, m, n, C)
            n1inf  = max(1, abs(n - u))                       
            n1sup  = min(Nrankeff, n + u)
            factmr = Hv(m) * Hv(m1)                                 
            do n1 = n1inf, n1sup
              factr = factmr * factur / sqrt(real(2 * n1 + 1,O))                                    
              fact  = (- 1._O)**(n + n1) * im**(n1 - 1) * C(n1) * factr                                        
              A(n,m,m1l) = A(n,m,m1l) + fact * Al(m1l,n1)
              B(n,m,m1l) = B(n,m,m1l) + fact * Bl(m1l,n1)
            end do                  
          end do
        end do
        deallocate (C)
      end do      
      do m1l = 1, 2
        m1 = - 3 + 2 * m1l
        do m1tl = 1, 2
          m1t = - 3 + 2 * m1tl
          dmmax = min(u, Mrankeff - m1, Mrankeff - m1t)
          dmmin = min(u, Mrankeff + m1, Mrankeff + m1t)
          do dm = - dmmin, dmmax
            m  = dm + m1
            mt = dm + m1t
            do n = max(1,abs(m)), Nrankeff
              do nt = max(1,abs(mt)), Nrankeff
                SA(n,nt,m,m1l,m1tl) = SA(n,nt,m,m1l,m1tl) +                         &
                                      A(n,m,m1l) * conjg(A(nt,mt,m1tl))
                SB(n,nt,m,m1l,m1tl) = SB(n,nt,m,m1l,m1tl) +                         &
                                      A(n,m,m1l) * conjg(B(nt,mt,m1tl))
                SC(n,nt,m,m1l,m1tl) = SC(n,nt,m,m1l,m1tl) +                         &
                                      B(n,m,m1l) * conjg(A(nt,mt,m1tl))
                SD(n,nt,m,m1l,m1tl) = SD(n,nt,m,m1l,m1tl) +                         &
                                      B(n,m,m1l) * conjg(B(nt,mt,m1tl))                                       
              end do  
            end do
          end do
        end do
      end do
    end do                  
  end if    
  do j = 1, 10        
    do i = 1, Ntheta
      S(j,i) = zero       
    end do
  end do
  do m1l = 1, 2
    if (PrnProgress) call write_progress (.false., m1l+2*Nrankeff+1, 2*Nrankeff+3) 
    m1 = - 3 + 2 * m1l
    do m1tl = 1, 2
      m1t   = - 3 + 2 * m1tl 
      m1c   = im * real(m1,O)
      m1tc  = im * real(m1t,O)
      m11t  = real(m1 * m1t,O)   
      dmmax = min(Mrankeff - m1, Mrankeff - m1t)
      dmmin = min(Mrankeff + m1, Mrankeff + m1t)
      do dm = - dmmin, dmmax
        m  = dm + m1
        mt = dm + m1t
        do n = max(1,abs(m)), Nrankeff
          do nt = max(1,abs(mt)), Nrankeff                                        
            AXX  =   SA(n,nt,m,m1l,m1tl)
            AXY  =   m1tc * SA(n,nt,m,m1l,m1tl)
            AYX  = - m1c  * SA(n,nt,m,m1l,m1tl)
            AYY  =   m11t * SA(n,nt,m,m1l,m1tl)           
            BXX  =   im * SB(n,nt,m,m1l,m1tl)
            BXY  =   im * m1tc * SB(n,nt,m,m1l,m1tl)
            BYX  = - im * m1c  * SB(n,nt,m,m1l,m1tl)
            BYY  =   im * m11t * SB(n,nt,m,m1l,m1tl)              
            CXX  =   im * SC(n,nt,m,m1l,m1tl)
            CXY  =   im * m1tc * SC(n,nt,m,m1l,m1tl)
            CYX  = - im * m1c  * SC(n,nt,m,m1l,m1tl)
            CYY  =   im * m11t * SC(n,nt,m,m1l,m1tl)              
            DXX  =   SD(n,nt,m,m1l,m1tl)
            DXY  =   m1tc * SD(n,nt,m,m1l,m1tl)
            DYX  = - m1c  * SD(n,nt,m,m1l,m1tl)
            DYY  =   m11t * SD(n,nt,m,m1l,m1tl)
            k    = IdentifyIndex (m, n, Nrankeff)
            kt   = IdentifyIndex (mt, nt, Nrankeff) 
            do i = 1, Ntheta
              MTh  = M2inf(i,k)                                
              MPh  = M3inf(i,k)
              MThC = conjg(M2inf(i,kt))                            
              MPhC = conjg(M3inf(i,kt))
              S(1,i)  = S(1,i)  + AXX * MTh * MThC + BXX * MTh * MPhC -             &
                                  CXX * MPh * MThC + DXX * MPh * MPhC 
              S(4,i)  = S(4,i)  + AXY * MTh * MPhC - BXY * MTh * MThC -             &
                                  CXY * MPh * MPhC - DXY * MPh * MThC
              S(5,i)  = S(5,i)  + AYY * MTh * MThC + BYY * MTh * MPhC -             &
                                  CYY * MPh * MThC + DYY * MPh * MPhC				
              S(6,i)  = S(6,i)  + AYX * MTh * MPhC - BYX * MTh * MThC -             &
                                  CYX * MPh * MPhC - DYX * MPh * MThC			
              S(8,i)  = S(8,i)  + AXX * MPh * MPhC - BXX * MPh * MThC +             &
                                  CXX * MTh * MPhC + DXX * MTh * MThC
              S(10,i) = S(10,i) + AYY * MPh * MPhC - BYY * MPh * MThC +             & 
                                  CYY * MTh * MPhC + DYY * MTh * MThC
              if (.not.MirorSym) then 				 				 				 
                S(2,i) = S(2,i) + AXY * MTh * MThC + BXY * MTh * MPhC -             &
                                  CXY * MPh * MThC + DXY * MPh * MPhC
                S(3,i) = S(3,i) + AXX * MTh * MPhC - BXX * MTh * MThC -             &
                                  CXX * MPh * MPhC - DXX * MPh * MThC
                S(7,i) = S(7,i) + AYY * MTh * MPhC - BYY * MTh * MThC -             &
                                  CYY * MPh * MPhC - DYY * MPh * MThC				    				 
                S(9,i) = S(9,i) + AXY * MPh * MPhC - BXY * MPh * MThC +             &
                                  CXY * MTh * MPhC + DXY * MTh * MThC	 
              end if				  
            end do                  
          end do
        end do          
      end do
    end do
  end do                      
  deallocate (M2inf, M3inf, SA, SB, SC, SD, A, B)
  if (.not.axsym) then
    deallocate (T)
  else 
    deallocate (TV)
  end if  
end subroutine SijSCpqEff
! **********************************************************************************
subroutine AvCscatCextEff (ReducedOrder, wavenumber, snorm, FileTmat, axsym, chiral,&
           delta, Nrank, Mrank, Nrankeff, Mrankeff, Cext, Cscat, Qext, Qscat)
!------------------------------------------------------------------------------------     
! The routine computes the average optical cross sections and efficiencies.         !
! The optical cross sections are:                                                   !
!                        2*Pi                                                       !
!            <Cext> = - ------ Re {SUM (2n+1) * (t11_n + t22_n)                     !
!                        k**2                                                       !                         
! and                                                                               !
!                         2*Pi                                                      !
!            <Cscat> =   ------ SUM (2n+1) * (t'11_n + t'22_n).                     !
!                         k**2                                                      !  
! The efficiencies are:                                                             ! 												 !
! 			  <Cext>         <Cext> * k**2		                    !
!            <Qext>  = ------------- = -----------------                            !
! 			Pi*anorm**2	     snorm                                  !
! and 			                                                            !
! 			 <Cscat>         <Cscat> * k**2	                            !
!            <Qscat> = ------------- = -----------------.                           !
! 			Pi*anorm**2	     snorm                                  !
!------------------------------------------------------------------------------------
  use parameters                                            
  implicit none
  integer       :: Nrank, Mrank, Nrankeff, Mrankeff
  real(O)       :: wavenumber, snorm, delta, Cext, Cscat, Qext, Qscat       
  character(80) :: FileTmat
  logical       :: ReducedOrder, axsym, chiral
!      
  integer       :: Nmax, nTV, n, dimensionTV, ntg, mtg
  real(O)       :: wavenumber2, Cexteff, Cscateff, dCext, dCscat 
  complex(O)    :: sum1, sum2, t11, t22
  logical       :: more, changeMrank
  complex(O),allocatable :: T(:,:), TV(:)
!
  Nmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)   
  if (.not. axsym) then                                           
    open (unit = iTmat, file = FileTmat, status = 'unknown')       
    call read_HeadFileTmat (ntg, mtg)
    call check_dimensionMat (ntg, mtg, Nmax)              
    allocate (T(2*ntg,2*mtg))
    call read_FileTmat (ntg, mtg, T) 
    close (unit = iTmat)
  else 
    nTV = dimensionTV (Mrank, Nrank, chiral)
    allocate (TV(nTV))
    call formTV (FileTmat, chiral, Mrank, Nrank, TV, nTV)  
  end if          
  Nrankeff = Nrank
  Mrankeff = Mrank
  sum1 = zero
  sum2 = zero
  do n = 1, Nrank
    if (.not. axsym) then
      call tnExtEff (n, Nrank, Mrankeff, Nmax, T, ntg, mtg, t11, t22)
    else 
      call tnExtAxsymEff (n, Nrank, Mrankeff, chiral, TV, nTV, t11, t22)
    end if
    sum1 = sum1 + (t11 + t22)
    if (.not. axsym) then
      call tnScatEff (n, Nrank, Mrankeff, Nrankeff, Nmax, T, ntg, mtg, t11, t22)
    else 
      call tnScatAxsymEff (n, Nrank, Mrankeff, Nrankeff, chiral, TV, nTV, t11, t22)
    end if
    sum2 = sum2 + (t11 + t22)
  end do                
  wavenumber2 = wavenumber * wavenumber
  Cext  = - 2._O * Pi * real(sum1,O) / wavenumber2
  Cscat =   2._O * Pi * real(sum2,O) / wavenumber2  
  if (ReducedOrder) then   
    Cexteff  = Cext
    Cscateff = Cscat
    dCext    = abs(Cext  - Cexteff)
    dCscat   = abs(Cscat - Cscateff)
    more     = (dCext < delta * abs(Cext)) .and. (dCscat < delta * abs(Cscat))
    changeMrank = .false.
    do while (more .and. Nrankeff > 1)
      Nrankeff = Nrankeff - 1
      if (Nrankeff < Mrank) then
        Mrankeff    =  Mrankeff - 1
        changeMrank = .true.
      end if
      sum1 = zero
      sum2 = zero
      do n = 1, Nrankeff
        if (.not. axsym) then
          call tnExtEff (n, Nrank, Mrankeff, Nmax, T, ntg, mtg, t11, t22)
        else 
          call tnExtAxsymEff (n, Nrank, Mrankeff, chiral, TV, nTV, t11, t22)
        end if
        sum1 = sum1 + (t11 + t22)
        if (.not. axsym) then
          call tnScatEff (n, Nrank, Mrankeff, Nrankeff, Nmax, T, ntg, mtg,          &
               t11, t22)
        else 
          call tnScatAxsymEff (n, Nrank, Mrankeff, Nrankeff, chiral, TV, nTV,       &
               t11, t22)
        end if
        sum2 = sum2 + (t11 + t22)
      end do                
      Cexteff  = - 2._O * Pi * real(sum1,O) / wavenumber2
      Cscateff =   2._O * Pi * real(sum2,O) / wavenumber2
      dCext  = abs(Cext  - Cexteff)
      dCscat = abs(Cscat - Cscateff)
      more   = (dCext < delta * abs(Cext)) .and. (dCscat < delta * abs(Cscat))
    end do  
    Nrankeff = Nrankeff + 1
    if (changeMrank) Mrankeff = Mrankeff + 1
  end if
  Qscat =  Cscat * wavenumber2 / snorm
  Qext  =  Cext  * wavenumber2 / snorm
  if (.not. axsym) then
    deallocate (T)
  else 
    deallocate (TV)
  end if                
end subroutine AvCscatCextEff
! **********************************************************************************
subroutine tnScatEff (n, Nrank, Mrankeff, Nrankeff, Nmax, T, ntg, mtg, t11, t22)
  use parameters
  implicit none                
  integer    :: n, Nrank, Mrankeff, Nrankeff, Nmax, ntg, mtg
  complex(O) :: T(2*ntg,2*mtg), t11, t22
!      
  integer    :: mp, m1, n1inf, n1
  real(O)    :: sum1, sum2
  complex(O) :: TG11, TG12, TG21, TG22
!
  sum1 = 0._O
  sum2 = 0._O      
  do mp = - Mrankeff, Mrankeff
    if (mp <= n) then
      do m1 = - Mrankeff, Mrankeff
        n1inf = max(1, abs(m1))
        do n1 = n1inf, Nrankeff
          call IdentifyMatrixElement (m1, n1, mp, n, Nrank, Nmax, T, ntg, mtg,      &
               TG11, TG12, TG21, TG22)
          sum1 = sum1 + abs(TG11)**2 + abs(TG21)**2   
          sum2 = sum2 + abs(TG12)**2 + abs(TG22)**2
        end do
      end do
    end if
  end do  
  t11 = cmplx(sum1,0.0,O)
  t22 = cmplx(sum2,0.0,O)
end subroutine tnScatEff
! **********************************************************************************
subroutine tnScatAxsymEff (n, Nrank, Mrankeff, Nrankeff, chiral, TV, nTV, t11, t22)
  use parameters
  implicit none                
  integer    :: n, Nrank, Mrankeff, Nrankeff, nTV
  complex(O) :: TV(nTV), t11, t22
  logical    :: chiral
!      
  integer    :: mp, n1inf, n1
  real(O)    :: sum1, sum2
  complex(O) :: TG11, TG12, TG21, TG22
!
  sum1 = 0._O
  sum2 = 0._O      
  do mp = - Mrankeff, Mrankeff
    if (mp <= n) then     
      n1inf = max(1, abs(mp))
      do n1 = n1inf, Nrankeff       
        call IdentifyMatrixElementAxsym (mp, n1, n, Nrank, chiral, TV, nTV,         &
             TG11, TG12, TG21, TG22)
        sum1 = sum1 + abs(TG11)**2 + abs(TG21)**2   
        sum2 = sum2 + abs(TG12)**2 + abs(TG22)**2
      end do      
    end if
  end do        
  t11 = cmplx(sum1,0.0,O)
  t22 = cmplx(sum2,0.0,O)
end subroutine tnScatAxsymEff
! **********************************************************************************
subroutine tnExtEff (n, Nrank, Mrankeff, Nmax, T, ntg, mtg, t11, t22)
  use parameters
  implicit none                
  integer    :: n, Nrank, Mrankeff, Nmax, ntg, mtg
  complex(O) :: T(2*ntg,2*mtg), t11, t22
!      vmr_main_input.dat
  integer    :: mp
  complex(O) :: sum1, sum2, TG11, TG12, TG21, TG22
!
  sum1 = zero
  sum2 = zero       
  do mp = - Mrankeff, Mrankeff
    if (abs(mp) <= n) then                        
      call IdentifyMatrixElement (mp, n, mp, n, Nrank, Nmax, T, ntg, mtg,           &
           TG11, TG12, TG21, TG22)       
      sum1 = sum1 + TG11            
      sum2 = sum2 + TG22            
    end if
  end do
  t11 = sum1
  t22 = sum2
end subroutine tnExtEff 
! **********************************************************************************
subroutine tnExtAxsymEff (n, Nrank, Mrankeff, chiral, TV, nTV, t11, t22)
  use parameters
  implicit none                
  integer    :: n, Nrank, Mrankeff, nTV
  complex(O) :: TV(nTV), t11, t22
  logical    :: chiral
!      
  integer    :: mp
  complex(O) :: sum1, sum2, TG11, TG12, TG21, TG22
!
  sum1 = zero
  sum2 = zero      
  do mp = - Mrankeff, Mrankeff
    if (abs(mp) <= n) then      
      call IdentifyMatrixElementAxsym (mp, n, n, Nrank, chiral, TV, nTV,            &
           TG11, TG12, TG21, TG22)        
      sum1 = sum1 + TG11      
      sum2 = sum2 + TG22            
    end if
  end do
  t11 = sum1
  t22 = sum2
end subroutine tnExtAxsymEff
! **********************************************************************************
subroutine RandomExtinctionMatrix (wavenumber, FileTmat, MirorSym, axsym, chiral,   &
           Nrank, Mrank, KeAv)
!------------------------------------------------------------------------------------     
! The routine computes the average extinction matrix. The elements of the           !
! extinction matrix are given by                                                    !
!                       2*Pi                                                        !
!            <Kii> = - ------ Re {SUM (2n+1) * (t11_n + t22_n)                      !
!                       k**2                                                        !
! for i = 1,2,3,4 and                                                               ! 
!                                2*Pi                                               !
!            <K14> = <K41>   =  ------ Re {SUM (2n+1) * (t12_n + t21_n),            !
!                                k**2                                               !
!                                                                                   !
!                                2*Pi                                               !
!            <K23> = - <K32> =  ------ Im {SUM (2n+1) * (t12_n + t21_n).            !
!                                k**2                                               !      
!------------------------------------------------------------------------------------
  use parameters                                            
  implicit none
  integer       :: Nrank, Mrank
  real(O)       :: wavenumber, KeAv(4,4)
  character(80) :: FileTmat
  logical       :: MirorSym, axsym, chiral
!      
  integer       :: Nmax, nTV, n, dimensionTV, ntg, mtg, i, j
  real(O)       :: wavenumber2, KII, K14, K23
  complex(O)    :: sum1, sum2, t11, t22, t12, t21  
  complex(O),allocatable :: T(:,:), TV(:)
!
  Nmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)   
  if (.not. axsym) then                                           
    open (unit = iTmat, file = FileTmat, status = 'unknown')       
    call read_HeadFileTmat (ntg, mtg)
    call check_dimensionMat (ntg, mtg, Nmax)              
    allocate (T(2*ntg,2*mtg))
    call read_FileTmat (ntg, mtg, T) 
    close (unit = iTmat)
  else 
    nTV = dimensionTV (Mrank, Nrank, chiral)
    allocate (TV(nTV))
    call formTV (FileTmat, chiral, Mrank, Nrank, TV, nTV)  
  end if            
  sum1 = zero
  sum2 = zero
  do n = 1, Nrank  
    if (.not. axsym) then
      call tnExtEff (n, Nrank, Mrank, Nmax, T, ntg, mtg, t11, t22)
    else 
      call tnExtAxsymEff (n, Nrank, Mrank, chiral, TV, nTV, t11, t22)
    end if
    if (.not.MirorSym .and. .not.axsym) then
      call tnExtOffEff (n, Nrank, Mrank, Nmax, T, ntg, mtg, t12, t21)        
    end if                
    sum1 = sum1 + (t11 + t22)
    sum2 = sum2 + (t12 + t21)        
  end do                
  wavenumber2 = wavenumber * wavenumber
  KII = - 2._O * Pi * real (sum1,O) / wavenumber2
  K14 =   2._O * Pi * real (sum2,O) / wavenumber2
  K23 =   2._O * Pi * aimag(sum2)   / wavenumber2    
  do i = 1, 4
    do j = 1, 4
      KeAv(i,j) = 0._O
    end do
  end do
  do i = 1, 4
    KeAv(i,i) = KII
  end do
  if (.not. MirorSym) then  
    KeAv(1,4) =   K14
    KeAv(4,1) =   K14      
    KeAv(2,3) =   K23
    KeAv(3,2) = - K23    
  end if                                        
  if (.not. axsym) then
    deallocate (T)
  else 
    deallocate (TV)
  end if                
end subroutine RandomExtinctionMatrix
! **********************************************************************************
subroutine tnExtOffEff (n, Nrank, Mrankeff, Nmax, T, ntg, mtg, t12, t21)
  use parameters
  implicit none                
  integer    :: n, Nrank, Mrankeff, Nmax, ntg, mtg
  complex(O) :: T(2*ntg,2*mtg), t12, t21
!      
  integer    :: mp
  complex(O) :: sum1, sum2, TG11, TG12, TG21, TG22
!
  sum1 = zero
  sum2 = zero       
  do mp = - Mrankeff, Mrankeff
    if (abs(mp) <= n) then                        
      call IdentifyMatrixElement (mp, n, mp, n, Nrank, Nmax, T, ntg, mtg,           &
           TG11, TG12, TG21, TG22)       
      sum1 = sum1 + TG12            
      sum2 = sum2 + TG21            
    end if
  end do
  t12 = sum1
  t21 = sum2
end subroutine tnExtOffEff 
! **********************************************************************************
subroutine AvCscatCextV (wavenumber, snorm, Ntheta, ZE, KeAv, CscatV, CextV,         &
                         QscatV, QextV)
!------------------------------------------------------------------------------------     
! The routine computes the V-component of the average scattering cross section,     !
!                                                                                   !
!                            |Pi                                                    !
!           <Cscat>_V = 2*Pi |  <F14(theta)> sin(theta) d(theta).                   !
!                            |0                                                     !
!                                                                                   !
! and the V-component of the average extinction cross section,                      !
!                                                                                   !
!                           <Cext>_V = <K14>                                        !        
!------------------------------------------------------------------------------------		     
  use parameters
  implicit none
  integer    :: Ntheta
  real(O)    :: wavenumber, snorm, ZE(4,4,Ntheta), KeAv(4,4), CscatV, QscatV,        &
                CextV, QextV
  
!
  integer    :: i
  real(O)    :: wavenumber2, atheta, btheta, theta, sum
  real(O),allocatable :: wtheta(:), xtheta(:)
!
  wavenumber2 = wavenumber * wavenumber
  allocate (wtheta(Ntheta), xtheta(Ntheta))
  atheta = 0._O
  btheta = Pi
  call Simpson (atheta, btheta, Ntheta, xtheta, wtheta)
  sum = 0._O 
  do i = 1, Ntheta    
    theta = xtheta(i)                         
    sum  = sum + wtheta(i) * sin(theta) * ZE(1,4,i)                           
  end do                   
  CscatV = 2._O * Pi * sum  
  CextV  = KeAv(1,4)   
  QscatV = CscatV * wavenumber2 / snorm      
  QextV  = CextV  * wavenumber2 / snorm          
  deallocate (wtheta, xtheta)  
end subroutine AvCscatCextV
! **********************************************************************************
subroutine MNinfiniteMatrix (Mrank, Nrank, Nmax, Ntheta, M2inf, M3inf)
!-----------------------------------------------------------------------------------
! The routine computes the theta- and phi- components of the vector spherical wave !
! functions Mmn at infinity. The scattering direction is in the azimuthal plane    !
! phi = 0° and therefore the azimuthal dependence exp(j * m * phi) is omitted.     !
! Note that N2inf = - j * M3inf and N3inf =  j * M2inf.                            !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax, Ntheta
  complex(O) :: M2inf(Ntheta,Nmax), M3inf(Ntheta,Nmax)
!
  integer    :: i, k, m, N0, n, l, ml
  real(O)    :: theta, nm, mlr 
  complex(O) :: factc, factp, factt
  real(O),allocatable :: Pnm(:), dPnm(:), pinm(:), taunm(:) 
!
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))
  do i = 1, Ntheta    
    theta = real(i - 1,O) * Pi / real(Ntheta - 1,O)          
    do m = 0, Mrank
      call leg_normalized (theta, m, Nrank, Pnm, dPnm, pinm, taunm)
      if (m == 0) then
        do k = 1, Nrank
          n  = k
          nm = real(2 * n * (n + 1),O)          
          nm = 1._O / sqrt(nm)        
          factc  = (- im)**(n + 1) * nm           
          factt  = factc * taunm(n)            
          M2inf(i,k) =   zero
          M3inf(i,k) = - factt          
        end do
      else
        N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
        ml = m
        do l = 1, 2
          mlr  = real(ml,O)              
          do k = 1, Nrank - m + 1
            n  = m + k - 1
            nm = real(2 * n * (n + 1),O)                
            nm = 1._O / sqrt(nm)              
            factc  = (- im)**(n + 1) * nm
            factp  = factc * mlr * pinm(n)
            factt  = factc * taunm(n)            
            M2inf(i,N0+k) =   im * factp
            M3inf(i,N0+k) = - factt           
          end do          
          N0 =   N0 + Nrank - m + 1
          ml = - ml
        end do
      end if                
    end do
  end do    
  deallocate (Pnm, dPnm, pinm, taunm)      
end subroutine MNinfiniteMatrix   
! **********************************************************************************
subroutine  IdentifyMatrixElement (m, n, m1, n1, Nrank, Nmax, T, ntg, mtg, T11, T12,&
            T21, T22)
  use parameters
  implicit none
  integer    :: m, n, m1, n1, Nrank, Nmax, ntg, mtg, k, k1, N0
  complex(O) :: T(2*ntg,2*mtg), T11, T12, T21, T22 
!
  if (m == 0) then      
    k = n
  else
    N0  = Nrank + (abs(m) - 1) * (2 * Nrank - abs(m) + 2)
    if (m > 0) then        
      k  = N0 + n - abs(m) + 1
    else
      N0 = N0 + Nrank - abs(m) + 1
      k  = N0 + n - abs(m) + 1
    end if
  end if      
  if (m1 == 0) then      
    k1 = n1
  else
    N0 = Nrank + (abs(m1) - 1) * (2 * Nrank - abs(m1) + 2)
    if (m1 > 0) then       
      k1 = N0 + n1 - abs(m1) + 1
    else
      N0 = N0 + Nrank - abs(m1) + 1
      k1 = N0 + n1 - abs(m1) + 1
    end if
  end if
  T11 = T(k,k1)
  T12 = T(k,k1+Nmax)
  T21 = T(k+Nmax,k1)
  T22 = T(k+Nmax,k1+Nmax)
end subroutine IdentifyMatrixElement 
! **********************************************************************************
subroutine IdentifyMatrixElementAxsym (m, n, n1, Nrank, chiral, TV, nTV, T11, T12,  &
           T21, T22)
  use parameters
  implicit none
  integer     :: m, n, n1, Nrank, nTV, Nc, Nmax, i, j, ma, tn, tm
  complex(O)  :: TV(nTV), T11, T12, T21, T22      
  logical     :: chiral
!
  if (m == 0) then        
    Nc = 0
    Nmax = Nrank
    i  = n
    j  = n1
    T11  = TV(Nc+i+2*(j-1)*Nmax)        
    T12  = TV(Nc+i+2*(j+Nmax-1)*Nmax)   
    T21  = TV(Nc+i+Nmax+2*(j-1)*Nmax)     
    T22  = TV(Nc+i+Nmax+2*(j+Nmax-1)*Nmax) 
  else        
    ma = abs(m)
    tn = Nrank * Nrank
    tm = (ma - 2) * (ma - 1)
    if (.not. chiral) then
      Nc = 8 * tn + 4 * (ma - 2) * tn - 4 * tm * Nrank + 2 * tm * (2 * ma - 3) / 3
    else 
      Nc = 12 * tn + 8 * (ma - 2) * tn - 8 * tm * Nrank + 4 * tm * (2 * ma - 3) / 3
    end if              
    Nmax = Nrank - ma + 1       
    i = n - ma + 1
    j = n1 - ma + 1
    if (m > 0) then
      T11 = TV(Nc+i+2*(j-1)*Nmax)              
      T12 = TV(Nc+i+2*(j+Nmax-1)*Nmax)        
      T21 = TV(Nc+i+Nmax+2*(j-1)*Nmax)        
      T22 = TV(Nc+i+Nmax+2*(j+Nmax-1)*Nmax)   
    else
      if (.not. chiral) then
        T11 =   TV(Nc+i+2*(j-1)*Nmax)         
        T12 = - TV(Nc+i+2*(j+Nmax-1)*Nmax)     
        T21 = - TV(Nc+i+Nmax+2*(j-1)*Nmax)    
        T22 =   TV(Nc+i+Nmax+2*(j+Nmax-1)*Nmax) 
      else 
        Nc  = Nc + 4 * Nmax * Nmax
        T11 = TV(Nc+i+2*(j-1)*Nmax)            
        T12 = TV(Nc+i+2*(j+Nmax-1)*Nmax)       
        T21 = TV(Nc+i+Nmax+2*(j-1)*Nmax)       
        T22 = TV(Nc+i+Nmax+2*(j+Nmax-1)*Nmax)  
      end if  
    end if        
  end if
end subroutine IdentifyMatrixElementAxsym  
! **********************************************************************************
function  IdentifyIndex (m, n, Nrank) result (k)
  implicit none
  integer  :: m, n, Nrank, k, N0
!
  if (m == 0) then      
    k = n
  else
    N0  = Nrank + (abs(m) - 1) * (2 * Nrank - abs(m) + 2)
    if (m > 0) then        
      k  = N0 + n - abs(m) + 1
    else
      N0 = N0 + Nrank - abs(m) + 1
      k  = N0 + n - abs(m) + 1
    end if
  end if            
end function IdentifyIndex                      
! **********************************************************************************
function Hv (m) result (a)
  use parameters
  implicit none
  integer  :: m
  real(O)  :: a
!
  if (m >= 0 ) then
    a = 1._O
  else 
    a = (- 1._O)**m  
  end if
end function Hv
! **********************************************************************************
function dimensionTV (Mrank, Nrank, chiral) result (nTV)
  implicit none
  integer   :: nTV, m, Mrank, Nrank, Nmax
  logical   :: chiral
!
  nTV = 0
  do m = 0, Mrank
    if (m == 0) then
      Nmax = Nrank
    else
      Nmax = Nrank - m + 1
    end if
    nTV = nTV + 4 * Nmax * Nmax
    if (chiral .and. m /= 0) then
      nTV = nTV + 4 * Nmax * Nmax
    end if
  end do      
end function dimensionTV
! **********************************************************************************
subroutine formTV (FileTmat, chiral, Mrank, Nrank, TV, nTV)
  use parameters
  implicit none
  integer       :: Mrank, Nrank, nTV
  complex(O)    :: TV(nTV)
  character(80) :: FileTmat
  logical       :: chiral
!      
  integer       :: Nc, m, Nmax, i, j, ntl, mtl
  complex(O),allocatable :: TL(:,:) 
!
  open (unit = iTmat, file = FileTmat, status = 'unknown')
  call read_HeadFileTmat (ntl, mtl)
  call check_dimensionMat (ntl, mtl, Nrank)             
  allocate (TL(2*ntl, 2*mtl))
  Nc = 0
  do m = 0, Mrank
    if (m == 0) then
      Nmax = Nrank
    else
      Nmax = Nrank - m + 1
    end if    
    call read_FileTmat (ntl, mtl, TL)
    do i = 1, 2*Nmax
      do j = 1, 2*Nmax
        TV(Nc+i+2*(j-1)*Nmax) = TL(i,j)
      end do
    end do
    Nc = Nc + 4 * Nmax * Nmax
    if (chiral .and. m /= 0) then      
      call read_FileTmat (ntl, mtl, TL)
      do i = 1, 2*Nmax
        do j = 1, 2*Nmax
          TV(Nc+i+2*(j-1)*Nmax) = TL(i,j)
        end do
      end do
      Nc = Nc + 4 * Nmax * Nmax
    end if
  end do
  close (unit = iTmat)
  deallocate (TL) 
end subroutine formTV 
! **********************************************************************************
! *               NUMERICAL ORIENTATION AVERAGING PROCEDURE                        *
! **********************************************************************************
subroutine SijSCpqInt (MirorSym, wavenumber, snorm, FileTmat, axsym, chiral, Nrank, &
           Mrank, Nalpha, Nbeta, Ngamma, Ntheta, UseSimpson, CextAv, CscatAv,       &
           QextAv, QscatAv, S, PrnProgress)
!------------------------------------------------------------------------------------
! The routine computes the matrix <Spq Sp1q1*> by using a numerical averaging       !
! procedure.                                                                        !      
!------------------------------------------------------------------------------------
  use parameters
  implicit none
  integer       :: Nrank, Mrank, Nalpha, Nbeta, Ngamma, Ntheta
  real(O)       :: wavenumber, snorm, CextAv, CscatAv, QextAv, QscatAv       
  complex(O)    :: S(10,Ntheta)
  character(80) :: FileTmat
  logical       :: MirorSym, axsym, chiral, UseSimpson, PrnProgress
!      
  integer       :: nTV, i, j, ialpha, ibeta, igamma, m, Nc, Mstart, Nmaxl,          &
                   dimensionTV, ntg, mtg, Nmax
  real(O)       :: aalpha, balpha, abeta, bbeta, agamma, bgamma, thetaGI, phiGI,    &
                   alphapX, alphapY, phiAZIMUT, norm, alpha, beta, gamma, Cext,     &
                   Cscat, Qext, Qscat, fact 
  complex(O)    :: FthetaX, FphiX, FthetaY, FphiY
  logical       :: ExtThetaDom                                 
  real(O),allocatable    :: walpha(:), xalpha(:), wbeta(:), xbeta(:), wgamma(:),    &
                            xgamma(:)      
  complex(O),allocatable :: T(:,:), TV(:), cX(:), cY(:), cX1(:), cY1(:), ccX(:),    &
                            ccY(:), FtX(:), FpX(:), FtY(:), FpY(:)
!
  Nmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  if (.not. axsym) then          
    open (unit = iTmat, file = FileTmat, status = 'unknown')       
    call read_HeadFileTmat (ntg, mtg)
    call check_dimensionMat (ntg, mtg, Nmax)              
    allocate (T(2*ntg,2*mtg))
    call read_FileTmat (ntg, mtg, T) 
    close (unit = iTmat)                                          
  else 
    nTV = dimensionTV (Mrank, Nrank, chiral)
    allocate (TV(nTV))
    call formTV (FileTmat, chiral, Mrank, Nrank, TV, nTV)  
  end if
  norm = 2._O * Pi
  norm = 1._O / norm
  allocate (walpha(Nalpha), xalpha(Nalpha))  
  aalpha = 0._O
  balpha = 2._O * Pi    
  call Simpson (aalpha, balpha, Nalpha, xalpha, walpha) 
  do ialpha = 1, Nalpha
    walpha(ialpha) = walpha(ialpha) * norm
  end do         
  allocate (wbeta(Nbeta), xbeta(Nbeta))
  if (.not. UseSimpson) then
    abeta = 0._O
    bbeta = Pi   
    call Gauss_Legendre (abeta, bbeta, Nbeta, wbeta, xbeta)   
  else
    abeta = - 1._O
    bbeta =   1._O 
    call Simpson (abeta, bbeta, Nbeta, xbeta, wbeta)     
  end if 
  do ibeta = 1, Nbeta      
    if (.not. UseSimpson) then 
      beta = xbeta(ibeta)
      wbeta(ibeta) = 0.5_O * wbeta(ibeta) * sin(beta)
    else
      wbeta(ibeta) = 0.5_O * wbeta(ibeta)
    end if
  end do           
  if (.not. axsym) then
    allocate (wgamma(Ngamma), xgamma(Ngamma))
    agamma = 0._O
    bgamma = 2._O * Pi           
    call Simpson (agamma, bgamma, Ngamma, xgamma, wgamma)
    do igamma = 1, Ngamma 
      wgamma(igamma) = wgamma(igamma) * norm
    end do
  end if
  thetaGI = 0._O
  phiGi   = 0._O
  ExtThetaDom = .false.   
  alphapX   = 0._O
  alphapY   = Pi / 2._O      
  phiAZIMUT = 0._O      
  allocate (ccX(2*Nmax), ccY(2*Nmax))
  allocate (FtX(Ntheta), FpX(Ntheta), FtY(Ntheta), FpY(Ntheta))
  do j = 1, 10
    do i = 1, Ntheta
      S(j,i) = zero
    end do
  end do
  CscatAv = 0._O
  CextAv  = 0._O                                  
  if (axsym) then
    allocate (cX(2*Nrank), cY(2*Nrank), cX1(2*Nrank), cY1(2*Nrank))     
    Mstart = 0
    do ialpha = 1, Nalpha
      alpha = xalpha(ialpha)
      if (PrnProgress) call write_progress (.false., ialpha, Nalpha)
      do ibeta = 1, Nbeta          
        if (.not. UseSimpson) then 
          beta = xbeta(ibeta)
        else
          beta = acos(xbeta(ibeta))
        end if
        gamma = 0._O    
        Nc = 0
        do m = Mstart, Mrank
          if (m == 0) then
            Nmaxl = Nrank
          else
            Nmaxl = Nrank - m + 1
          end if        
          call PWcoefficients_ab_m (thetaGI, phiGI, alpha, beta, gamma,             &
               alphapX, m, Nrank, Nmaxl, cX)
          call PWcoefficients_ab_m (thetaGI, phiGI, alpha, beta, gamma,             &
               alphapY, m, Nrank, Nmaxl, cY)          
          call product_vector_vector_AVplus (Nc, Nmaxl, TV, nTV, cX, cX1)
          call product_vector_vector_AVplus (Nc, Nmaxl, TV, nTV, cY, cY1)             
          call extend_vector_positive (cX1, ccX, m, Mstart, Nrank, Nmaxl, Nmax)
          call extend_vector_positive (cY1, ccY, m, Mstart, Nrank, Nmaxl, Nmax)                    
          if (m /= 0) then
            call PWcoefficients_ab_m (thetaGI, phiGI, alpha, beta, gamma,           &
                 alphapX, - m, Nrank, Nmaxl, cX)
            call PWcoefficients_ab_m (thetaGI, phiGI, alpha, beta, gamma,           &
                 alphapY, - m, Nrank, Nmaxl, cY)                        
            if (.not. chiral) then
              call product_vector_vector_AVminus (Nc, Nmaxl, TV, nTV, cX, cX1)
              call product_vector_vector_AVminus (Nc, Nmaxl, TV, nTV, cY, cY1)              
            else 
              Nc = Nc + 4 * Nmaxl * Nmaxl
              call product_vector_vector_AVplus (Nc, Nmaxl, TV, nTV, cX, cX1)
              call product_vector_vector_AVplus (Nc, Nmaxl, TV, nTV, cY, cY1)                         
            end if  
            call extend_vector_negative (cX1, ccX, m, Nrank, Nmaxl, Nmax)
            call extend_vector_negative (cY1, ccY, m, Nrank, Nmaxl, Nmax)                   
          end if    
          Nc = Nc + 4 * Nmaxl * Nmaxl           
        end do
        call F_azimuthal (ccX, Mrank, Nrank, Nmax, Ntheta, phiAZIMUT,               &
             alpha, beta, gamma, wavenumber, ExtThetaDom, FtX, FpX)
        call F_azimuthal (ccY, Mrank, Nrank, Nmax, Ntheta, phiAZIMUT,               &
             alpha, beta, gamma, wavenumber, ExtThetaDom, FtY, FpY)                       
        fact = walpha(ialpha) * wbeta(ibeta)                  
        call Smatrix  (MirorSym, Ntheta, fact, FtX, FpX, FtY, FpY, S)
        call F_azimuthal_tetaGS (ccX, Mrank, Nrank, Nmax, thetaGI, phiGI, alpha,    &
             beta, gamma, wavenumber, FthetaX, FphiX)
        call F_azimuthal_tetaGS (ccY, Mrank, Nrank, Nmax, thetaGI, phiGI, alpha,    &
             beta, gamma, wavenumber, FthetaY, FphiY)                     
        call CQscat (ccX, Mrank, Nrank, Nmax, wavenumber, snorm, Cscat, Qscat)
        call CQext  (ccX, Mrank, Nrank, Nmax, thetaGI, phiGI, alpha, beta, gamma,   &
             alphapX, wavenumber, snorm, Cext, Qext)         
        CscatAv = CscatAv + fact * Cscat
        CextAv  = CextAv + fact * Cext
      end do
    end do 
    deallocate (cX, cY, cX1, cY1)       
  else 
    allocate (cX(2*Nmax), cY(2*Nmax))   
    do ialpha = 1, Nalpha
      alpha = xalpha(ialpha)
      if (PrnProgress) call write_progress (.false., ialpha, Nalpha)
      do ibeta = 1, Nbeta                                       
        if (.not. UseSimpson) then 
          beta = xbeta(ibeta)
        else
          beta = acos(xbeta(ibeta))
        end if       
        do igamma = 1, Ngamma 
          gamma = xgamma(igamma)
          call PWcoefficients_ab (thetaGI, phiGI, alpha, beta, gamma,               &
               alphapX, Mrank, Nrank, Nmax, cX)
          call PWcoefficients_ab (thetaGI, phiGI, alpha, beta, gamma,               &
               alphapY, Mrank, Nrank, Nmax, cY)           
          call product_matrix_vector (2*Nmax, 2*Nmax, T, 2*ntg, 2*mtg, cX, ccX)   
          call product_matrix_vector (2*Nmax, 2*Nmax, T, 2*ntg, 2*mtg, cY, ccY)                
          call F_azimuthal (ccX, Mrank, Nrank, Nmax, Ntheta, phiAZIMUT,             &
               alpha, beta, gamma, wavenumber, ExtThetaDom, FtX, FpX)
          call F_azimuthal (ccY, Mrank, Nrank, Nmax, Ntheta, phiAZIMUT,             &
               alpha, beta, gamma, wavenumber, ExtThetaDom, FtY, FpY)                                               
          fact = walpha(ialpha) * wbeta(ibeta) * wgamma(igamma)                              
          call Smatrix  (MirorSym, Ntheta, fact, FtX, FpX, FtY, FpY, S)
          call F_azimuthal_tetaGS (ccX, Mrank, Nrank, Nmax, thetaGI, phiGI, alpha,  &
               beta, gamma, wavenumber, FthetaX, FphiX)
          call F_azimuthal_tetaGS (ccY, Mrank, Nrank, Nmax, thetaGI, phiGI, alpha,  &
               beta, gamma, wavenumber, FthetaY, FphiY)                             
          call CQscat (ccX, Mrank, Nrank, Nmax, wavenumber, snorm, Cscat, Qscat)
          call CQext  (ccX, Mrank, Nrank, Nmax, thetaGI, phiGI, alpha, beta, gamma, &
               alphapX, wavenumber, snorm, Cext, Qext)          
          CscatAv = CscatAv + fact * Cscat
          CextAv  = CextAv + fact * Cext
        end do        
      end do
    end do
    deallocate (cX, cY)         
  end if        
  norm    = wavenumber * wavenumber / snorm
  QscatAv = CscatAv * norm
  QextAv  = CextAv  * norm                                    
  deallocate (ccX, ccY, FtX, FpX, FtY, FpY)
  deallocate (walpha, xalpha, wbeta, xbeta)
  if (.not. axsym) then
    deallocate (T, wgamma, xgamma)
  else 
    deallocate (TV)
  end if
end subroutine SijSCpqInt
! **********************************************************************************
subroutine Smatrix (MirorSym, Ntheta, fact, FtX, FpX, FtY, FpY, S)
  use parameters
  implicit none
  integer    :: Ntheta, i
  real(O)    :: fact
  complex(O) :: FtX(Ntheta), FpX(Ntheta), FtY(Ntheta), FpY(Ntheta), S(10,Ntheta) 
  logical    :: MirorSym  
!
  do i = 1, Ntheta
    S(1,i)  = S(1,i)  + fact * FtX(i) * conjg(FtX(i))        
    S(4,i)  = S(4,i)  + fact * FtX(i) * conjg(FpY(i))             
    S(5,i)  = S(5,i)  + fact * FtY(i) * conjg(FtY(i))        
    S(6,i)  = S(6,i)  + fact * FtY(i) * conjg(FpX(i))                
    S(8,i)  = S(8,i)  + fact * FpX(i) * conjg(FpX(i))           
    S(10,i) = S(10,i) + fact * FpY(i) * conjg(FpY(i))
    if (.not.MirorSym) then
      S(2,i) = S(2,i) + fact * FtX(i) * conjg(FtY(i))
      S(3,i) = S(3,i) + fact * FtX(i) * conjg(FpX(i))     
      S(7,i) = S(7,i) + fact * FtY(i) * conjg(FpY(i))
      S(9,i) = S(9,i) + fact * FpX(i) * conjg(FpY(i))            
    end if    
  end do
end subroutine Smatrix
