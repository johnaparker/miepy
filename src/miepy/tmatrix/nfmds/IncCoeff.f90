! **********************************************************************************
! *            ROUTINES FOR COMPUTING THE INCIDENT FIELD COEFFICIENTS              *
! *    ------------------------------------------------------------------------    *
! *    Partial list of subroutines:                                                *
! *      PWcoefficients_ab,                    PWcoefficients_CircPolab(UNUSED),   *
! *      PWcoefficients_Polab,                 PWcoefficients_ab_m,                *
! *      PWcoefficients_CircPolab_m(UNUSED),   PWcoefficients_Polab_m,             *
! *      PWcoefficientsPARTSUB,                PWcoefficientsPARTSUBrefl,          *
! *      PWcoefficientsPARTSUBtrans,           Kcor,                               *
! *      Fresnel_aer_sub,                      Fresnel_sub_aer,                    *
! *      GBcoefficients_ab,                    GBcoefficients_ab_m,                *
! *      product_rotcoefficients,              GBcoefficients_ab_local             *
! **********************************************************************************
subroutine PWcoefficients_ab (thetaGI, phiGI, alpha, beta, gamma, alphap, Mrank,   &
           Nrank, Nmax, c) 
!-----------------------------------------------------------------------------------
! The routine computes the expansion coefficients of a linearly polarized plane    !
! wave for the azimuthal modes m = 0,1,...,Mrank.                                  ! 
!----------------------------------------------------------------------------------- 
  use parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax
  real(O)    :: thetaGI, phiGI, alpha, beta, gamma, alphap      
  complex(O) :: c(2*Nmax)
!                 
  integer    :: k, m, n, N0, ml, l      
  real(O)    :: nm, thetaLI, phiLI, e0eT, e0eP, mlr, arg
  complex(O) :: fact, factc, factp, factt
  real(O),allocatable :: Pnm(:), dPnm(:), pinm(:), taunm(:)  
!
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))
  call T_spherical_global_local (thetaGI, phiGI, alpha, beta, gamma, thetaLI, phiLI)
  call parameters_coefficients_ab (thetaGI, phiGI, thetaLI, phiLI, alpha, beta,     &
       gamma, alphap, e0eT, e0eP)
  do m = 0, Mrank        
    call leg_normalized (thetaLI, m, Nrank, Pnm, dPnm, pinm, taunm)      
    if (m == 0) then
      do k = 1, Nrank
        n     = k               
        nm    = real(2 * n * (n + 1),O)         
        nm    = 1._O / sqrt(nm)            
        fact  = 4._O * im**n * nm          
        factt = fact * taunm(n)     
        c(k)      = - factt * e0eP
        c(Nmax+k) = - im * factt * e0eT               
      end do
    else        
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      ml = m        
      do l = 1, 2
        mlr  = real(ml,O)
        arg  = mlr * phiLI
        fact = exp(- im * arg)                  
        do k = 1, Nrank - m + 1
          n     = m + k - 1
          nm    = real(2 * n * (n + 1),O)               
          nm    = 1._O / sqrt(nm)             
          factc = 4._O * im**n * fact * nm
          factp = factc * im * mlr * pinm(n)
          factt = factc * taunm(n)          
          c(N0+k)      = - factp * e0eT - factt * e0eP
          c(Nmax+N0+k) = - im * (factt * e0eT - factp * e0eP) 
        end do
        N0 =   N0 + Nrank - m + 1
        ml = - ml
      end do                  
    end if
  end do
  deallocate (Pnm, dPnm, pinm, taunm)
end subroutine PWcoefficients_ab
! **********************************************************************************
subroutine PWcoefficients_CircPolab (thetaGI, phiGI, alpha, beta, gamma, right,     &
           Mrank, Nrank, Nmax, c) 
!-----------------------------------------------------------------------------------
! The routine computes the expansion coefficients of a circularly polarized        !
! plane wave for the azimuthal modes m = 0,1,...,Mrank.                            ! 
!----------------------------------------------------------------------------------- 
  use parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax
  real(O)    :: thetaGI, phiGI, alpha, beta, gamma
  complex(O) :: c(2*Nmax)
  logical    :: right
!                 
  integer    :: k, m, n, N0, ml, l      
  real(O)    :: norm, nm, thetaLI, phiLI, mlr, arg
  complex(O) :: epol_theta, epol_phi, e0eT, e0eP, fact, factc, factp, factt
  real(O),allocatable :: Pnm(:), dPnm(:), pinm(:), taunm(:)  
!
  norm = sqrt(2._O) / 2._O 
  if (right) then
    epol_theta = norm * one
    epol_phi   = norm * im
  else
    epol_theta =   norm * one
    epol_phi   = - norm * im
  end if
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))
  call T_spherical_global_local (thetaGI, phiGI, alpha, beta, gamma, thetaLI, phiLI)
  call parameters_coefficients_Polab (thetaGI, phiGI, thetaLI, phiLI, alpha, beta,  &
       gamma, epol_theta, epol_phi, e0eT, e0eP) 
  do m = 0, Mrank        
    call leg_normalized (thetaLI, m, Nrank, Pnm, dPnm, pinm, taunm)      
    if (m == 0) then
      do k = 1, Nrank
        n     = k               
        nm    = real(2 * n * (n + 1),O)         
        nm    = 1._O / sqrt(nm)            
        fact  = 4._O * im**n * nm          
        factt = fact * taunm(n)     
        c(k)      = - factt * e0eP
        c(Nmax+k) = - im * factt * e0eT               
      end do
    else        
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      ml = m        
      do l = 1, 2
        mlr  = real(ml,O)
        arg  = mlr * phiLI
        fact = exp(- im * arg)                  
        do k = 1, Nrank - m + 1
          n     = m + k - 1
          nm    = real(2 * n * (n + 1),O)               
          nm    = 1._O / sqrt(nm)             
          factc = 4._O * im**n * fact * nm
          factp = factc * im * mlr * pinm(n)
          factt = factc * taunm(n)          
          c(N0+k)      = - factp * e0eT - factt * e0eP
          c(Nmax+N0+k) = - im * (factt * e0eT - factp * e0eP) 
        end do
        N0 =   N0 + Nrank - m + 1
        ml = - ml
      end do                  
    end if
  end do
  deallocate (Pnm, dPnm, pinm, taunm)
end subroutine PWcoefficients_CircPolab
! **********************************************************************************
subroutine PWcoefficients_Polab (thetaGI, phiGI, alpha, beta, gamma, epol_theta,    &
           epol_phi, Mrank, Nrank, Nmax, c) 
!-----------------------------------------------------------------------------------
! The routine computes the expansion coefficients of an arbitrarily polarized      !
! plane wave for the azimuthal modes m = 0,1,...,Mrank.                            ! 
!----------------------------------------------------------------------------------- 
  use parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax
  real(O)    :: thetaGI, phiGI, alpha, beta, gamma
  complex(O) :: epol_theta, epol_phi, c(2*Nmax)
!                 
  integer    :: k, m, n, N0, ml, l      
  real(O)    :: nm, thetaLI, phiLI, mlr, arg
  complex(O) :: e0eT, e0eP, fact, factc, factp, factt
  real(O),allocatable :: Pnm(:), dPnm(:), pinm(:), taunm(:)  
!  
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))
  call T_spherical_global_local (thetaGI, phiGI, alpha, beta, gamma, thetaLI, phiLI)
  call parameters_coefficients_Polab (thetaGI, phiGI, thetaLI, phiLI, alpha, beta,  &
       gamma, epol_theta, epol_phi, e0eT, e0eP) 
  do m = 0, Mrank        
    call leg_normalized (thetaLI, m, Nrank, Pnm, dPnm, pinm, taunm)      
    if (m == 0) then
      do k = 1, Nrank
        n     = k               
        nm    = real(2 * n * (n + 1),O)         
        nm    = 1._O / sqrt(nm)            
        fact  = 4._O * im**n * nm          
        factt = fact * taunm(n)     
        c(k)      = - factt * e0eP
        c(Nmax+k) = - im * factt * e0eT               
      end do
    else        
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      ml = m        
      do l = 1, 2
        mlr  = real(ml,O)
        arg  = mlr * phiLI
        fact = exp(- im * arg)                  
        do k = 1, Nrank - m + 1
          n     = m + k - 1
          nm    = real(2 * n * (n + 1),O)               
          nm    = 1._O / sqrt(nm)             
          factc = 4._O * im**n * fact * nm
          factp = factc * im * mlr * pinm(n)
          factt = factc * taunm(n)          
          c(N0+k)      = - factp * e0eT - factt * e0eP
          c(Nmax+N0+k) = - im * (factt * e0eT - factp * e0eP) 
        end do
        N0 =   N0 + Nrank - m + 1
        ml = - ml
      end do                  
    end if
  end do
  deallocate (Pnm, dPnm, pinm, taunm)
end subroutine PWcoefficients_Polab
! **********************************************************************************
subroutine PWcoefficients_ab_m (thetaGI, phiGI, alpha, beta, gamma, alphap, m,      &
           Nrank, Nmax, c)
!-----------------------------------------------------------------------------------
! The routine computes the expansion coefficients of a linearly polarized plane    !
! wave for the azimuthal mode m.                                                   !
!----------------------------------------------------------------------------------- 
  use parameters
  implicit none
  integer    :: m, Nrank, Nmax
  real(O)    :: thetaGI, phiGI, alpha, beta, gamma, alphap      
  complex(O) :: c(2*Nmax)
!      
  integer    :: k, n     
  real(O)    :: nm, thetaLI, phiLI, e0eT, e0eP, mr, arg
  complex(O) :: fact, factc, factp, factt
  real(O),allocatable :: Pnm(:), dPnm(:), pinm(:), taunm(:)  
! 
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))
  call T_spherical_global_local (thetaGI, phiGI, alpha, beta, gamma, thetaLI, phiLI)
  call parameters_coefficients_ab (thetaGI, phiGI, thetaLI, phiLI, alpha, beta,     &
       gamma, alphap, e0eT, e0eP)
  call leg_normalized (thetaLI, abs(m), Nrank, Pnm, dPnm, pinm, taunm)
  mr   = real(m,O)
  arg  = mr * phiLI
  fact = exp(- im * arg)
  do k = 1, Nmax
    if (m == 0) then
      n = k
    else
      n = abs(m) + k - 1
    end if
    nm    = real(2 * n * (n + 1),O)             
    nm    = 1._O / sqrt(nm)     
    factc = 4._O * im**n * fact * nm        
    factp = factc * im * mr * pinm(n)
    factt = factc * taunm(n)        
    c(k)      = - factp * e0eT - factt * e0eP
    c(Nmax+k) = - im * (factt * e0eT - factp * e0eP)
  end do                    
  deallocate (Pnm, dPnm, pinm, taunm)
end subroutine PWcoefficients_ab_m 
! **********************************************************************************
subroutine PWcoefficients_CircPolab_m (thetaGI, phiGI, alpha, beta, gamma, right,   &
           m, Nrank, Nmax, c)
!-----------------------------------------------------------------------------------
! The routine computes the expansion coefficients of a circularly polarized        !
! plane wave for the azimuthal mode m.                                             !
!----------------------------------------------------------------------------------- 
  use parameters
  implicit none
  integer    :: m, Nrank, Nmax
  real(O)    :: thetaGI, phiGI, alpha, beta, gamma     
  complex(O) :: c(2*Nmax)
  logical    :: right  
!      
  integer    :: k, n     
  real(O)    :: norm, nm, thetaLI, phiLI, mr, arg
  complex(O) :: epol_theta, epol_phi, e0eT, e0eP, fact, factc, factp, factt
  real(O),allocatable :: Pnm(:), dPnm(:), pinm(:), taunm(:)  
! 
  norm = sqrt(2._O) / 2._O 
  if (right) then
    epol_theta = norm * one
    epol_phi   = norm * im
  else
    epol_theta =   norm * one
    epol_phi   = - norm * im
  end if
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))
  call T_spherical_global_local (thetaGI, phiGI, alpha, beta, gamma, thetaLI, phiLI)
  call parameters_coefficients_Polab (thetaGI, phiGI, thetaLI, phiLI, alpha, beta,  &
       gamma, epol_theta, epol_phi, e0eT, e0eP)        
  call leg_normalized (thetaLI, abs(m), Nrank, Pnm, dPnm, pinm, taunm)
  mr   = real(m,O)
  arg  = mr * phiLI
  fact = exp(- im * arg)
  do k = 1, Nmax
    if (m == 0) then
      n = k
    else
      n = abs(m) + k - 1
    end if
    nm    = real(2 * n * (n + 1),O)             
    nm    = 1._O / sqrt(nm)     
    factc = 4._O * im**n * fact * nm        
    factp = factc * im * mr * pinm(n)
    factt = factc * taunm(n)        
    c(k)      = - factp * e0eT - factt * e0eP
    c(Nmax+k) = - im * (factt * e0eT - factp * e0eP)
  end do                    
  deallocate (Pnm, dPnm, pinm, taunm)
end subroutine PWcoefficients_CircPolab_m 
! **********************************************************************************
subroutine PWcoefficients_Polab_m (thetaGI, phiGI, alpha, beta, gamma, epol_theta,  &
           epol_phi, m, Nrank, Nmax, c)
!-----------------------------------------------------------------------------------
! The routine computes the expansion coefficients of an arbitrarily polarized      !
! plane wave for the azimuthal mode m.                                             !
!----------------------------------------------------------------------------------- 
  use parameters
  implicit none
  integer    :: m, Nrank, Nmax
  real(O)    :: thetaGI, phiGI, alpha, beta, gamma
  complex(O) :: epol_theta, epol_phi, c(2*Nmax)
!      
  integer    :: k, n     
  real(O)    :: nm, thetaLI, phiLI, mr, arg
  complex(O) :: e0eT, e0eP, fact, factc, factp, factt
  real(O),allocatable :: Pnm(:), dPnm(:), pinm(:), taunm(:)  
!   
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))
  call T_spherical_global_local (thetaGI, phiGI, alpha, beta, gamma, thetaLI, phiLI)
  call parameters_coefficients_Polab (thetaGI, phiGI, thetaLI, phiLI, alpha, beta,  &
       gamma, epol_theta, epol_phi, e0eT, e0eP)        
  call leg_normalized (thetaLI, abs(m), Nrank, Pnm, dPnm, pinm, taunm)
  mr   = real(m,O)
  arg  = mr * phiLI
  fact = exp(- im * arg)
  do k = 1, Nmax
    if (m == 0) then
      n = k
    else
      n = abs(m) + k - 1
    end if
    nm    = real(2 * n * (n + 1),O)             
    nm    = 1._O / sqrt(nm)     
    factc = 4._O * im**n * fact * nm        
    factp = factc * im * mr * pinm(n)
    factt = factc * taunm(n)        
    c(k)      = - factp * e0eT - factt * e0eP
    c(Nmax+k) = - im * (factt * e0eT - factp * e0eP)
  end do                    
  deallocate (Pnm, dPnm, pinm, taunm)
end subroutine PWcoefficients_Polab_m 
! **********************************************************************************
subroutine PWcoefficientsPARTSUB (beta0, alphap, m, Nrank, Nmax, c)
!------------------------------------------------------------------------------------
! The routine computes the expansion coefficients of a plane wave (traveling in     !
! the Oxz plane) for the azimuthal mode m.                                          !
!------------------------------------------------------------------------------------ 
  use parameters
  implicit none
  integer    :: m, Nrank, Nmax
  real(O)    :: beta0, alphap
  complex(O) :: c(2*Nmax)    
!
  integer    :: k, n      
  real(O)    :: nm, mr, ca, sa
  complex(O) :: factc, factp, factt
  real(O),allocatable :: Pnm(:), dPnm(:), pinm(:), taunm(:)
!
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))
  call leg_normalized (beta0, abs(m), Nrank, Pnm, dPnm, pinm, taunm)
  mr = real(m,O)
  ca = cos(alphap)
  sa = sin(alphap)
  do k = 1, Nmax
    if (m == 0) then
      n = k
    else
      n = abs(m) + k - 1
    end if
    nm    = real(2 * n * (n + 1),O)         
    nm    = 1._O / sqrt(nm)     
    factc = 4._O * im**n * nm * cos(mr * Pi)        
    factp = factc * im * mr * pinm(n)
    factt = factc * taunm(n)
    c(k)      = - factp * ca - factt * sa
    c(Nmax+k) = - im * (factt * ca - factp * sa)
  end do
  deallocate (Pnm, dPnm, pinm, taunm)     
end subroutine PWcoefficientsPARTSUB
! **********************************************************************************
subroutine PWcoefficientsPARTSUBrefl (beta0, alphap, z0, wavenumber, ind_ref, m,    &
           Nrank, Nmax, c)
!-----------------------------------------------------------------------------------
! The routine computes the expansion coefficients of the reflected plane wave      !
! (traveling in the Oxz plane) for the azimuthal mode m.                           !
!----------------------------------------------------------------------------------- 
  use parameters
  implicit none
  integer    :: m, Nrank, Nmax
  real(O)    :: beta0, alphap, z0, wavenumber      
  complex(O) :: c(2*Nmax), ind_ref
!
  integer    :: k, n      
  real(O)    :: nm, betaM0, mr
  complex(O) :: Rpar, Rperp, Tpar, Tperp, et, ep, cosb, sinb, cosb0, factc,         &
                factp, factt, phase
  real(O),allocatable :: Pnm(:), dPnm(:), pinm(:), taunm(:)
!
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))
  cosb0 = cmplx(cos(beta0),0._O)
  call Fresnel_aer_sub (cosb0, ind_ref, Rpar, Rperp, Tpar, Tperp, cosb, sinb)
  phase  = exp(2._O * im * wavenumber * cosb0 * z0)
  et     = cos(alphap) * Rpar * phase
  ep     = sin(alphap) * Rperp * phase
  betaM0 = Pi - beta0
  call leg_normalized (betaM0, abs(m), Nrank, Pnm, dPnm, pinm, taunm)
  mr = real(m,O)  
  do k = 1, Nmax
    if (m == 0) then
      n = k
    else
      n = abs(m) + k - 1
    end if
    nm    = real(2 * n * (n + 1),O)         
    nm    = 1._O / sqrt(nm)     
    factc = 4._O * im**n * nm * cos(mr * Pi)        
    factp = factc * im * mr * pinm(n)
    factt = factc * taunm(n)
    c(k)      = - factp * et - factt * ep
    c(Nmax+k) = - im * (factt * et - factp * ep)
  end do      
  deallocate (Pnm, dPnm, pinm, taunm)      
end subroutine PWcoefficientsPARTSUBrefl
! **********************************************************************************
subroutine PWcoefficientsPARTSUBtrans (beta, alphap, z0, wavenumber, ind_ref, m,    &
           Nrank, Nmax, c)
!----------------------------------------------------------------------------------
! The routine computes the expansion coefficients of the transmitted plane wave   !
! (traveling in the Oxz plane) for the azimuthal mode m.                          !
! --------------------------------------------------------------------------------- 
  use parameters
  implicit none
  integer    :: m, Nrank, Nmax
  real(O)    :: beta, alphap, z0, wavenumber
  complex(O) :: c(2*Nmax), ind_ref     
!
  integer    :: k, n
  real(O)    :: nm, mr
  complex(O) :: Rpar, Rperp, Tpar, Tperp, et, ep, cosb, sinb0, cosb0, factc, factp, &
                factt, phase
  complex(O),allocatable :: Pnm(:), dPnm(:), pinm(:), taunm(:)
!
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))
  cosb = cmplx(cos(beta),0._O)   
  call Fresnel_sub_aer (cosb, ind_ref, Rpar, Rperp, Tpar, Tperp, cosb0, sinb0)  
  phase = exp(im * wavenumber * (cosb0 - ind_ref * cosb) * z0)       
  et    = cos(alphap) * Tpar * phase
  ep    = sin(alphap) * Tperp * phase
  cosb0 = - cosb0 
  call Leg_normalized_complex (sinb0, cosb0, abs(m), Nrank, Pnm, dPnm, pinm, taunm)
  mr = real(m,O) 
  do k = 1, Nmax
    if (m == 0) then
      n = k
    else
      n = abs(m) + k - 1
    end if
    nm    = real(2 * n * (n + 1),O)         
    nm    = 1._O / sqrt(nm)     
    factc = 4._O * im**n * nm
    factp = factc * im * mr * pinm(n)
    factt = factc * taunm(n)
    c(k)      = - factp * et - factt * ep
    c(Nmax+k) = - im * (factt * et - factp * ep)
  end do   
  deallocate (Pnm, dPnm, pinm, taunm)    
end subroutine PWcoefficientsPARTSUBtrans
! **********************************************************************************
subroutine Fresnel_aer_sub (cosb0, m, Rpar, Rperp, Tpar, Tperp, cosb, sinb)
  use parameters
  use derived_parameters
  implicit none
  complex(O) :: cosb0, sinb0, m, Rpar, Rperp, Tpar, Tperp, cosb, sinb        
! 
  sinb0 = sqrt(1._O - cosb0**2)
  sinb  = sinb0 / m
  cosb  = sqrt(1._O - sinb**2)
  if (aimag(m * cosb) < 0._O) cosb = - cosb 
  if (abs(m * cosb0 + cosb) < MachEps) then
    print "(/,2x,'Warning in subroutine Fresnel_aer_sub in module IncCoeff:')"
    print "(  2x,'the denominator in the expression of the Fresnel reflection')"
    print "(  2x,'coefficients is smaller than the machine precision;')"    
  end if            
  if (abs(m * cosb0 - cosb) < MachEps) then      ! redundant
    Rpar = zero
    Tpar = 1._O / m
  else
    Rpar = (m * cosb0 - cosb) / (m * cosb0 + cosb)
    Tpar = 2._O * cosb0 / (m * cosb0 + cosb)
  end if
  if (abs(cosb0 - m * cosb) < MachEps) then      ! redundant
    Rperp = zero
    Tperp = one
  else
    Rperp = (cosb0 - m * cosb) / (cosb0 + m * cosb)
    Tperp = 2._O * cosb0 / (cosb0 + m * cosb)
  end if      
end subroutine Fresnel_aer_sub
! **********************************************************************************
subroutine Fresnel_sub_aer (cosb, m, Rpar, Rperp, Tpar, Tperp, cosb0, sinb0)
  use parameters
  use derived_parameters
  implicit none
  complex(O) :: cosb, m, Rpar, Rperp, Tpar, Tperp, sinb0, cosb0, sinb
!
  sinb  = sqrt(1._O - cosb**2)
  sinb0 = m * sinb
  cosb0 = sqrt(1._O - sinb0**2)       
  if (aimag(cosb0) < 0._O) cosb0 = - cosb0  
  if (abs(m * cosb0 + cosb) < MachEps) then
    print "(/,2x,'Warning in subroutine Fresnel_sub_aer in module IncCoeff:')"
    print "(  2x,'the denominator in the expression of the Fresnel transmission')"
    print "(  2x,'coefficients is smaller than the machine precision;')"    
  end if         
  if (abs(cosb - m * cosb0) < MachEps) then      ! redundant
    Rpar = zero
    Tpar = m   
  else
    Rpar = (cosb - m * cosb0) / (m * cosb0 + cosb)
    Tpar = 2._O * m * cosb / (m * cosb0 + cosb)
  end if
  if (abs(m * cosb - cosb0) < MachEps) then      ! redundant
    Rperp = zero
    Tperp = one 
  else
    Rperp = (m * cosb - cosb0) / (cosb0 + m * cosb)
    Tperp = 2._O * m * cosb / (cosb0 + m * cosb) 
  end if      
end subroutine Fresnel_sub_aer                  
! **********************************************************************************
subroutine GBcoefficients_ab (wavenumber, x0, y0, z0, w0, thetaGI, phiGI, alpha,    &
           beta, gamma, alphap, Mrank, Nrank, Nmax, c)
!------------------------------------------------------------------------------------
! The routine computes the expansion coefficients of a Gaussian beam for the        !
! azimuthal modes m = 0,1,...,Mrank.  Note that x0, y0, z0 are the coordinates of   !
! the beam focal point with respect to a global coordinate system (having the same  !
! origin as the particle coordinate system).                                        ! 
!------------------------------------------------------------------------------------ 
  use parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax
  real(O)    :: wavenumber, x0, y0, z0, w0, thetaGI, phiGI, alpha, beta, gamma,     &
                alphap          
  complex(O) :: c(2*Nmax)
!      
  integer    :: m, Nmaxl, k, N0, l, ml
  complex(O),allocatable :: cl(:)
!
  do m = 0, Mrank                       
    if (m == 0) then
      Nmaxl = Nrank
      allocate (cl(2*Nmaxl))
      call GBcoefficients_ab_m (wavenumber, x0, y0, z0, w0, thetaGI, phiGI, alpha,  &
           beta, gamma, alphap, m, Nmaxl, cl)
      do k = 1, Nmaxl        
        c(k)      = cl(k)
        c(Nmax+k) = cl(k+Nmaxl)
      end do
      deallocate (cl)
    else        
      N0    = Nrank + (m - 1) * (2 * Nrank - m + 2)
      Nmaxl = Nrank - m + 1
      allocate (cl(2*Nmaxl))           
      ml = m
      do l = 1, 2
        call GBcoefficients_ab_m (wavenumber, x0, y0, z0, w0, thetaGI, phiGI,       &
             alpha, beta, gamma, alphap, ml, Nmaxl, cl)
        do k = 1, Nmaxl                           
          c(N0+k)      = cl(k)
          c(Nmax+N0+k) = cl(k+Nmaxl)                
        end do
        N0 =   N0 + Nrank - m + 1
        ml = - ml
      end do 
      deallocate (cl)                 
    end if
  end do     
end subroutine GBcoefficients_ab     
! **********************************************************************************
subroutine GBcoefficients_ab_m (wavenumber, x0, y0, z0, w0, thetaGI, phiGI, alpha,  &
           beta, gamma, alphap, m, Nmax, c)
!------------------------------------------------------------------------------------
! The routine computes the expansion coefficients of a Gaussian beam for the        !
! azimuthal mode m. Note that x0, y0, z0 are the coordinates of the beam focal point!
! with respect to a global coordinate system (having the same origin as the         !
! particle coordinate system).                                                      !
!------------------------------------------------------------------------------------ 
  use parameters
  implicit none
  integer    :: m, Nmax
  real(O)    :: wavenumber, x0, y0, z0, w0, thetaGI, phiGI, alpha, beta, gamma,     &
                alphap          
  complex(O) :: c(2*Nmax)
!      
  integer    :: k, n, m1, p
  real(O)    :: xpG, ypG, zpG, xpL, ypL, zpL
  complex(O) :: suma, sumb
  complex(O),allocatable :: a(:), b(:), d(:)
!
  xpG = - x0
  ypG = - y0
  zpG = - z0
  call T_cartesian_global_local (xpG, ypG, zpG, phiGI, thetaGI, alphap, xpL, ypL, zpL)            
  do k = 1, Nmax
    if (m == 0) then
      n = k
    else
      n = abs(m) + k - 1
    end if
    allocate (a(2*n+1), b(2*n+1), d(2*n+1)) 
    call GBcoefficients_ab_local (wavenumber, xpL, ypL, zpL, w0, n, a, b)   
    call product_rotcoefficients (thetaGI, phiGI, alpha, beta, gamma, alphap, m, n, d) 
    suma = zero
    sumb = zero
    do m1 = 0, n
      if (m1 == 0) then        
        suma = suma + a(1) * d(1)
        sumb = sumb + b(1) * d(1)
      else
        p    = 2 * m1
        suma = suma + a(p) * d(p) + a(p+1) * d(p+1)
        sumb = sumb + b(p) * d(p) + b(p+1) * d(p+1)
      end if
    end do
    c(k)     = suma
    c(k+Nmax)= sumb        
    deallocate (a, b, d)
  end do
end subroutine GBcoefficients_ab_m      
! **********************************************************************************
subroutine product_rotcoefficients (thetaGI, phiGI, alpha, beta, gamma, alphap,     &
           m, n, d)
  use parameters
  implicit none
  integer     :: m, n
  real(O)     :: thetaGI, phiGI, alpha, beta, gamma, alphap
  complex(O)  :: d(2*n+1)   
!      
  integer     :: m1, m2, k       
  complex(O)  :: sum, Dnm1m2, Dnm2m, coef_rotation
!
  do m1 = 0, n
    if (m1 == 0) then     
      sum = zero
      do m2 = - n, n
        Dnm1m2 = coef_rotation (- alphap, - thetaGI, - phiGI, n, m1, m2)
        Dnm2m  = coef_rotation (  alpha, beta, gamma, n, m2, m)        
        sum    = sum + Dnm1m2 * Dnm2m
      end do
      d(1) = sum
    else
      k = 2 * m1 
!     --- m1 > 0 ---
      sum = zero         
      do m2 =  - n, n
        Dnm1m2 = coef_rotation (- alphap, - thetaGI, - phiGI, n, m1, m2)
        Dnm2m  = coef_rotation ( alpha, beta, gamma, n, m2, m)
        sum    = sum + Dnm1m2 * Dnm2m
      end do
      d(k) = sum
!     --- m1 < 0 ---
      sum = zero          
      do m2 = - n, n
        Dnm1m2 = coef_rotation (- alphap, - thetaGI, - phiGI, n, - m1, m2)
        Dnm2m  = coef_rotation ( alpha, beta, gamma, n, m2, m)
        sum    = sum + Dnm1m2 * Dnm2m
      end do
      d(k+1) = sum
    end if
  end do      
end subroutine product_rotcoefficients                    
! **********************************************************************************
subroutine GBcoefficients_ab_local (wavenumber, x0, y0, z0, w0, n, a, b)
  use parameters
  use derived_parameters
  implicit none
  integer     :: n
  real(O)     :: wavenumber, x0, y0, z0, w0      
  complex(O)  :: a(2*n+1), b(2*n+1)
!      
  integer     :: m, k
  real(O)     :: l, ro0, phi0, ron, argz, argpm, argpp, sig, w02, dz
  complex(O)  :: Q, PSI, arg, fact, Kcor, factz, factm, factp, imQ, dzc
  complex(O),allocatable :: Jbes(:)
!
  if (abs(w0) < MachEps) then
    print "(/,2x,'Warning in subroutine GBcoefficients_ab_local in module IncCoeff:')"
    print "(  2x,'the waist radius of the Gaussian beam is smaller than the machine')"
    print "(  2x,'precision;')"    
  end if
  l   = wavenumber * w0 * w0  
  dz  = 2._O * z0 / l
  dzc = cmplx(dz,0.0,O)
  Q   = 1._O / (im - dzc)
  ro0 = sqrt(x0 * x0 + y0 * y0)     
  if (ro0 < MachEps) then
    phi0 = 0._O                  
  else    
    phi0 = atan2(y0,x0)
    if (phi0 < 0._O) phi0 = 2._O * Pi + phi0     
  end if 
  ron = (real(n,O) + 0.5_O) / wavenumber
  w02 = w0 * w0
  imQ = im * Q
  PSI = imQ * exp(-imQ * (ro0 * ro0 + ron * ron) / w02)
  arg = 2._O * Q * ro0 * ron / w02
  allocate (Jbes(0:n+1))  
  call bes_J (arg, n + 1, Jbes)     
  argz  = wavenumber * z0
  factz = exp(im * argz)   
  do m = 0, n
    if (m == 0) then              
      fact =   2._O * factz * Jbes(1) * Kcor(m,n) * PSI       
      a(1) = - cos(phi0) * fact       
      b(1) =   im * sin(phi0) * fact
    else
      k     = 2 * m   
      argpm = real(m - 1,O) * phi0
      argpp = real(m + 1,O) * phi0
      fact  = factz * Kcor(m,n) * PSI
!     --- m > 0 ---      
      factm = exp(im * argpm) * Jbes(m-1)
      factp = exp(im * argpp) * Jbes(m+1)             
      a(k)  = fact * (factm - factp)          
      b(k)  = fact * (factm + factp) 
!     --- m < 0 ---
      sig    = (- 1._O)**(m - 1)
      factm  = sig * exp(- im * argpm) * Jbes(m-1)
      factp  = sig * exp(- im * argpp) * Jbes(m+1)           
      a(k+1) = fact * (factp - factm)       
      b(k+1) = fact * (factp + factm)  
    end if
  end do
  deallocate (Jbes)   
end subroutine GBcoefficients_ab_local    
! **********************************************************************************
function Kcor (m, n) result (K)
  use parameters
  implicit none
  integer      :: m, n
  complex(O)   :: K
!      
  integer      :: p
  real(O)      :: fact, sum, lnx  
!
  if (m == 0) then
    fact = 2._O * sqrt(real(n * (n + 1),O) / real(2 * n + 1,O))      
    K    = im**n * fact
  else
    sum = 0._O
    do p = 1, 2 * abs(m)
      sum = sum + log(real(n - abs(m) + p,O))
    end do
    lnx = real(abs(m) - 1,O) * log(2._O / real(2 * n + 1,O))
    lnx = lnx + 0.5_O * (log(real(2 * n + 1,O)) - log(real(n * (n + 1),O)) + sum)     
    K   = im**(n + abs(m)) * (-1._O)**abs(m) * exp(lnx)
  end if     
end function Kcor

