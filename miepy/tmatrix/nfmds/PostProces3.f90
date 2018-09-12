! **********************************************************************************
! *         ANALYTICAL SIZE AVERAGING PROCEDURE FOR SPHERICAL PARTICLES            *
! *    -------------------------------------------------------------------------   *
! *    Partial list of subroutines:                                                *
! *      DiffScatCrossSectAvrgSPH,         ScatteringMatrixSPH,                    *
! *      ExtendedScatteringMatrixSPH,      SijSCpqSPH,                             *
! *      AVSCTCOEF,                        AvCscatCextSPH,                         *
! *      AVSCTCOEF1,                       GeometryPars,                           *
! *      PDF                                                                       *
! **********************************************************************************
subroutine DiffScatCrossSectAvrgSPH (anorm, alphap, Ntheta, SS, normalized, Cext,   &
           Cscat, Qext, Qscat)
!-----------------------------------------------------------------------------------
! The routine computes the average differential scattering cross section in the    !
! azimuthal plane phi = 0. The incident wave propagates along the Z-axis of the    !
! global coordinate system.                                                        !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none 
  integer       :: Ntheta, i
  real(O)       :: anorm, alphap, Cext, Cscat, Qext, Qscat, fact, ca, sa, ca2, sa2, &
                   theta, h, v   
  complex(O)    :: SS(3,Ntheta)
  logical       :: normalized
!
  write (iDSCS,"(/,2x,'Results:',/)")
  write (iDSCS,"(2x,'Cross Sections and Efficiencies:')") 
  write (iDSCS,"(2x,'average scattering cross section = ',1pe13.4,';')") Cscat
  write (iDSCS,"(2x,'average extinction cross section = ',1pe13.4,';')") Cext 
  write (iDSCS,*)    
  write (iDSCS,"(2x,'average scattering efficiency    = ',1pe13.4,';')") Qscat
  write (iDSCS,"(2x,'average extinction efficiency    = ',1pe13.4,';')") Qext 
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
  ca   = cos(alphap)
  sa   = sin(alphap) 
  ca2  = ca * ca
  sa2  = sa * sa
  fact = Pi * anorm * anorm
  fact = 1._O / fact
  do i = 1, Ntheta    
    theta = real(i - 1,O) * 180._O / real(Ntheta - 1,O)       
    h = real(SS(1,i),O) * ca2
    v = real(SS(3,i),O) * sa2  
    if (normalized) then
      h = h * fact
      v = v * fact                                          
    end if  
    write (iDSCS,"(1x,f6.2,5x,1pe13.4,5x,1pe13.4)") theta, h, v
  end do
end subroutine DiffScatCrossSectAvrgSPH 
! **********************************************************************************
subroutine ScatteringMatrixSPH (theta, Ntheta, ZE, Z)
!-----------------------------------------------------------------------------------
! The routine computes the scattering matrix Z at the scattering angle theta. ZE is!
! the extended scattering matrix and Ntheta is the number of scattering angles (at !
! which the extended scattering matrix is computed) in the azimuthale plane        !
! phi = 0.                                                                         !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: Ntheta
  real(O)    :: theta, ZE(4,Ntheta), Z(4,4)
! 
  integer    :: i
  real(O)    :: Z11, Z12, Z33, Z34
  real(O),allocatable :: thetaV(:), ZV(:) 
! 
  allocate (thetaV(Ntheta), ZV(Ntheta))
  do i = 1, Ntheta
    thetaV(i) = real(i - 1,O) * Pi / real(Ntheta - 1,O)        
  end do
  do i = 1, Ntheta
    ZV(i) = ZE(1,i)
  end do
  call Interp (Ntheta, thetaV, ZV, theta, Z11)
  do i = 1, Ntheta
    ZV(i) = ZE(2,i)
  end do
  call Interp (Ntheta, thetaV, ZV, theta, Z12)
  do i = 1, Ntheta
    ZV(i) = ZE(3,i)
  end do
  call Interp (Ntheta, thetaV, ZV, theta, Z33)
  do i = 1, Ntheta
    ZV(i) = ZE(4,i)
  end do
  call Interp (Ntheta, thetaV, ZV, theta, Z34)
  Z(1,1) =   Z11     
  Z(1,2) =   Z12
  Z(1,3) =   0._O
  Z(1,4) =   0._O 
  Z(2,1) =   Z12
  Z(2,2) =   Z11  
  Z(2,3) =   0._O
  Z(2,4) =   0._O
  Z(3,1) =   0._O
  Z(3,2) =   0._O
  Z(3,3) =   Z33  
  Z(3,4) =   Z34        
  Z(4,1) =   0._O
  Z(4,2) =   0._O
  Z(4,3) = - Z34 
  Z(4,4) =   Z33 
  deallocate (thetaV, ZV)
end subroutine ScatteringMatrixSPH
! **********************************************************************************
subroutine ExtendedScatteringMatrixSPH (Ntheta, SS, ZE)
!-----------------------------------------------------------------------------------
! The routine computes the scattering matrix ZE, for all scattering angles. SS is  !
! the average size matrix and Ntheta is the number of scattering angles.           !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: Ntheta, i
  real(O)    :: ZE(4,Ntheta)
  complex(O) :: SS(3,Ntheta)
  complex(O) :: S11SC11, S11SC22, S22SC22 
!
  do i = 1, Ntheta     
    S11SC11 = SS(1,i) 
    S11SC22 = SS(2,i) 
    S22SC22 = SS(3,i)        
    ZE(1,i) = 0.5_O * real(S11SC11 + S22SC22,O)     
    ZE(2,i) = 0.5_O * real(S11SC11 - S22SC22,O)   
    ZE(3,i) = real (S11SC22,O)  
    ZE(4,i) = aimag(S11SC22)        
  end do
end subroutine ExtendedScatteringMatrixSPH
! **********************************************************************************
subroutine SijSCpqSPH (wavenumber, ind_ref, amin, amax, TypeDist, Npar, par, Cnst,  &
           Nrank, Nint, Ntheta, S, PrnProgress)
!-----------------------------------------------------------------------------------     
! The routine computes the size average matrix <Spq Sp1q1*> by using an analytical !
! averaging procedure.                                                             ! 
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer      :: TypeDist, Npar, Nrank, Nint, Ntheta
  real(O)      :: wavenumber, amin, amax, Cnst, par(Npar)
  complex(O)   :: ind_ref, S(3,Ntheta)  
  logical      :: PrnProgress
!      
  integer       :: i, j, n, nt, ml, m, mtl, mt, k, kt, Nmax, IdentifyIndex
  real(O)       :: wavenumber2, invwno
  complex(O)    :: fct, Ptt, Ptp, Ppt, Ppp, MTh, MPh, MThC, MPhC  
  complex(O),allocatable :: M2inf(:,:), M3inf(:,:), A(:,:), B(:,:), C(:,:)
!
  Nmax = 3 * Nrank
  allocate (M2inf(Ntheta,Nmax), M3inf(Ntheta,Nmax)) 
  call MNinfiniteMatrix (1, Nrank, Nmax, Ntheta, M2inf, M3inf)
  allocate (A(Nrank,Nrank), B(Nrank,Nrank), C(Nrank,Nrank))
  call AVSCTCOEF (wavenumber, ind_ref, amin, amax, TypeDist, Npar, par, Cnst,       &
       Nrank, Nint, Nrank, A, B, C)
  wavenumber2 = wavenumber * wavenumber
  invwno      = 1._O / wavenumber2
  do j = 1,3        
    do i = 1, Ntheta
      S(j,i) = zero       
    end do
  end do
  do n = 1, Nrank
    if (PrnProgress) call write_progress (.false., n, Nrank)
    do nt = 1, Nrank
      do i = 1, Ntheta
        Ptt = zero
        Ptp = zero
        Ppt = zero
        Ppp = zero
        fct = im**(n - nt) * sqrt(real((2 * n + 1) * (2 * nt + 1),O))
        fct = fct * invwno
        do ml = 1, 2
          m = - 3 + 2 * ml
          k = IdentifyIndex (m, n, Nrank)                 
          do mtl = 1, 2
            mt = - 3 + 2 * mtl
            kt = IdentifyIndex (mt, nt, Nrank)
            MTh  = M2inf(i,k)                                  
            MPh  = M3inf(i,k)
            MThC = conjg(M2inf(i,kt))                              
            MPhC = conjg(M3inf(i,kt))
            Ptt  = Ptt + m * mt * MTh * MThC * fct
            Ptp  = Ptp + im * m * MTh * MPhC * fct
            Ppt  = Ppt - im * mt * MPh * MThC * fct
            Ppp  = Ppp + MPh * MPhC * fct
          end do
        end do
        S(1,i) = S(1,i) + A(n,nt) * Ptt + C(n,nt) * Ptp +                           &
                          conjg(C(nt,n)) * Ppt + B(n,nt) * Ppp
        S(2,i) = S(2,i) + A(n,nt) * Ptp + C(n,nt) * Ptt +                           &
                          conjg(C(nt,n)) * Ppp + B(n,nt) * Ppt
        S(3,i) = S(3,i) + A(n,nt) * Ppp + C(n,nt) * Ppt +                           &
                          conjg(C(nt,n)) * Ptp + B(n,nt) * Ptt                            
      end do 
    end do
  end do
  deallocate (M2inf, M3inf, A, B, C)  
end subroutine SijSCpqSPH    
! **********************************************************************************
subroutine AVSCTCOEF (wavenumber, ind_ref, amin, amax, TypeDist, Npar, par, Cnst,   &
           Nrank, Nint, NP, A, B, C)
  use parameters
  implicit none
  integer      :: TypeDist, Npar, Nrank, Nint, NP
  real(O)      :: wavenumber, amin, amax, Cnst, par(Npar)
  complex(O)   :: ind_ref, A(NP,NP), B(NP,NP), C(NP,NP)
!
  integer      :: n, nt, p
  real(O)      :: r, f, PDF 
  real(O),allocatable    :: x(:), w(:)
  complex(O),allocatable :: cv(:)

  allocate (cv(2*Nrank))
  allocate (x(Nint), w(Nint))
  call Simpson (amin, amax, Nint, x, w)
  do n = 1, Nrank
    do nt = 1, Nrank
      A(n,nt) = zero
      B(n,nt) = zero
      C(n,nt) = zero
    end do
  end do
  do p = 1, Nint    
    r = x(p)
    f = w(p) * Cnst * PDF(TypeDist, Npar, par, r)
    call coefficients_fg_m (wavenumber, r, ind_ref, 1, Nrank, Nrank, cv)
    do n = 1, Nrank
      do nt = 1, Nrank      
        A(n,nt) = A(n,nt) + cv(n) * conjg(cv(nt)) * f
        B(n,nt) = B(n,nt) + cv(Nrank+n) * conjg(cv(Nrank+nt)) * f
        C(n,nt) = C(n,nt) + cv(n) * conjg(cv(Nrank+nt)) * f
      end do
    end do
  end do
  deallocate (cv, x, w)
end subroutine AVSCTCOEF
! **********************************************************************************
subroutine AvCscatCextSPH (wavenumber, snorm, ind_ref, amin, amax, TypeDist, Npar,  &
           par, Cnst, Nrank, Nint, Cext, Cscat, Qext, Qscat, AsymPar)
!-----------------------------------------------------------------------------------     
! The routine computes the average cross sections, efficiencies and asymmetry      !
! parameter.                                                                       !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer      :: TypeDist, Npar, Nrank, Nint
  real(O)      :: wavenumber, amin, amax, Cnst, par(Npar), Cext, Cscat, Qext,       &
                  Qscat, AsymPar, snorm
  complex(O)   :: ind_ref  
!      
  integer       :: n
  real(O)       :: sum, fct, wavenumber2, invwno
  complex(O)    :: sum1, sum2
  real(O),allocatable    :: ccav(:)               
  complex(O),allocatable :: cav(:), cav1(:), cav2(:)
!
  wavenumber2 = wavenumber * wavenumber
  invwno      = 1._O / wavenumber2
  allocate (cav(2*Nrank), ccav(2*Nrank), cav1(Nrank-1), cav2(Nrank))
  call AVSCTCOEF1 (wavenumber, ind_ref, amin, amax, TypeDist, Npar, par, Cnst,      &
       Nrank, Nint, cav, ccav, cav1, cav2)
  sum1 = zero
  sum  = 0._O
  sum2 = zero
  do n = 1, Nrank
    fct  = real(2 * n + 1,O)
    sum1 = sum1 + fct * (cav(n)  + cav(Nrank+n))
    sum  = sum  + fct * (ccav(n) + ccav(Nrank+n))
    fct  = fct / real(n * (n + 1),O)
    sum2 = sum2 + fct * cav2(n)
  end do        
  Cext  = - 2._O * Pi * real(sum1,O) * invwno
  Cscat =   2._O * Pi * sum * invwno
  sum1 = zero
  do n = 1, Nrank - 1
    fct  = real(n * (n + 2),O) / real(n + 1,O)
    sum1 = sum1 + fct * cav1(n)
  end do   
  AsymPar = 4._O * Pi * real(sum1 + sum2,O) * invwno / Cscat
  Qscat   = Cscat * wavenumber2 / snorm
  Qext    = Cext * wavenumber2 / snorm
  deallocate (cav, ccav, cav1, cav2)
end subroutine AvCscatCextSPH
! **********************************************************************************
subroutine AVSCTCOEF1 (wavenumber, ind_ref, amin, amax, TypeDist, Npar, par, Cnst,  &
           Nrank, Nint, cav, ccav, cav1, cav2)
  use parameters
  implicit none
  integer      :: TypeDist, Npar, Nrank, Nint
  real(O)      :: wavenumber, amin, amax, Cnst, par(Npar), ccav(2*Nrank) 
  complex(O)   :: ind_ref, cav(2*Nrank), cav1(Nrank-1), cav2(Nrank)
!
  integer      :: n, p
  real(O)      :: r, f, PDF 
  real(O),allocatable    :: x(:), w(:)
  complex(O),allocatable :: c(:)
!
  allocate (c(2*Nrank))
  allocate (x(Nint), w(Nint))
  call Simpson (amin, amax, Nint, x, w)
  do n = 1, 2*Nrank
    cav(n)  = zero
    ccav(n) = 0._O
  end do
  do n = 1, Nrank - 1
    cav1(n) = zero
  end do
  do n = 1, Nrank
    cav2(n) = zero
  end do
  do p = 1, Nint
    r = x(p)
    f = w(p) * Cnst * PDF(TypeDist, Npar, par, r)
    call coefficients_fg_m (wavenumber, r, ind_ref, 1, Nrank, Nrank, c)
    do n = 1, 2*Nrank
      cav (n) = cav(n)  + c(n) * f
      ccav(n) = ccav(n) + abs(c(n)) * abs(c(n)) * f      
    end do
    do n = 1, Nrank - 1
      cav1(n) = cav1(n) + (c(n) * conjg(c(n+1)) +                                   &
      c(Nrank+n) * conjg(c(Nrank+n+1))) * f
    end do
    do n = 1, Nrank
      cav2(n) = cav2(n) + c(n) * conjg(c(Nrank+n)) * f 
    end do
  end do
  deallocate (c, x, w)
end subroutine AVSCTCOEF1
! **********************************************************************************
subroutine GeometryPars (amin, amax, TypeDist, Npar, par, Nint, Cnst, AvrgArea,     &
           EffRadius, AvrgRadius, AvrgVolume, VWAvrgRadius, EffVar)
  use parameters
  implicit none
  integer      :: TypeDist, Npar, Nint
  real(O)      :: amin, amax, par(Npar), Cnst, AvrgArea, EffRadius, AvrgRadius,     &
                  AvrgVolume, VWAvrgRadius, EffVar
!
  integer      :: p
  real(O)      :: r, r2, r3, r4, dr2, sum, sum1, sum2, sum3, sum4, PDF 
  real(O),allocatable    :: x(:), w(:)
!
  allocate (x(Nint), w(Nint))
  call Simpson (amin, amax, Nint, x, w)  
  sum  = 0._O
  sum1 = 0._O
  sum2 = 0._O
  sum3 = 0._O
  sum4 = 0._O
  do p = 1, Nint
    r    = x(p) 
    r2   = r  * r
    r3   = r2 * r
    r4   = r2 * r2
    sum  = sum  + w(p) * PDF(TypeDist, Npar, par, r)
    sum1 = sum1 + w(p) * r  * PDF(TypeDist, Npar, par, r)
    sum2 = sum2 + w(p) * r2 * PDF(TypeDist, Npar, par, r)
    sum3 = sum3 + w(p) * r3 * PDF(TypeDist, Npar, par, r)  
    sum4 = sum4 + w(p) * r4 * PDF(TypeDist, Npar, par, r)         
  end do
  Cnst = 1._O / sum
  AvrgArea     = Pi * sum2 * Cnst
  EffRadius    = Pi * sum3 * Cnst / AvrgArea
  AvrgRadius   = sum1 * Cnst
  AvrgVolume   = 4._O * Pi * sum3 * Cnst / 3._O
  VWAvrgRadius = 4._O * Pi * sum4 * Cnst / (3._O * AvrgVolume)
  sum  = 0._O
  do p = 1, Nint
    r   = x(p)  
    r2  = r * r             
    dr2 = (r - EffRadius) * (r - EffRadius)
    sum = sum + w(p) * dr2 * r2 * PDF(TypeDist, Npar, par, r)
  end do
  EffVar = Pi * sum * Cnst / (AvrgArea * EffRadius * EffRadius)  
  deallocate (x, w)
end subroutine GeometryPars
!***********************************************************************************
function PDF (TypeDist, Npar, par, a) result(func)
  use parameters
  implicit none
  integer,intent(in) :: TypeDist, Npar
  real(O),intent(in) :: par(Npar), a
  real(O) :: func
!
  real(O) :: alpha, gamma, ac, ag, sg, ag1, sg1, lnf, lnf1
!
  select case (TypeDist)
  case (1)
!   --- modified gamma distribution ---
    alpha = par(1)
    gamma = par(2)  
    ac    = par(3)
    lnf   = alpha * log(a) - alpha * (a / ac)**gamma / gamma    
    func  = exp(lnf)
  case (2)
!   --- log normal distribution ---
    ag    = par(1)
    sg    = par(2)
    lnf   = - log(a) - 0.5_O * (log(a/ag) / log(sg))**2 
    func  = exp(lnf)
  case (3)
!   --- gamma distribution ---
    alpha = par(1)
    gamma = par(2)
    lnf   = (1._O - 3._O * gamma) * log(a) / gamma - a / alpha / gamma
    func  = exp(lnf)
  case (4)
!   --- power law distribution ---    
    alpha =   par(1)
    lnf   = - alpha * log(a)
    func  =   exp(lnf)
  case (5)
!   --- modified bimodal log normal distribution ---
    ag    = par(1)   
    sg    = par(2)
    ag1   = par(3)
    sg1   = par(4)
    gamma = par(5)
    lnf   = - 4._O * log(a) - 0.5_O * (log(a/ag)  / log(sg))**2 
    lnf1  = - 4._O * log(a) - 0.5_O * (log(a/ag1) / log(sg1))**2 
    func  = exp(lnf) + gamma * exp(lnf1)     
  end select
end function PDF
