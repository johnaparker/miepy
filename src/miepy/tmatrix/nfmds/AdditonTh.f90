! **********************************************************************************
! *             ROUTINES FOR COMPUTING THE TRANSLATION ADDITION                    *
! *        COEFFICIENTS, ROTATION FUNCTIONS AND COUPLING COEFFICIENTS              *
! *     ----------------------------------------------------------------------     *
! *     Partial list of subroutines:                                               *
! *       Jacobi_n0,                  dmm1n,                   dmm1nVect,          *
! *       coef_rotation,              coef_rotationVect,       coef_rotationMatr,  *
! *       MatTransCmnn1,              MatTransAB_mnn1,         MatTransAB_mn_m1n1, *
! *       matrix_inverse,             MatRot_R_mn_m1n1_single, MatRot_R_mn_m1n1,   *
! *       MatTransRot_TR_mn_m1n1,     MatRotTrans_RT_mn_m1n1,                      *
! *       MatRotTransRot_RTR_mn_m1n1, identity_transformation, lnfactorial,        *  
! *       gumnkl,                     CouplingCoef                                 *
! ********************************************************************************** 
function Jacobi_n0 (beta, m, m1) result (Pmm1n0)
!-----------------------------------------------------------------------------------
! The routine computes the Jacobi polynom for m >= 0, m1 >= 0, and n0 = max(m,m1). !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: m, m1
  real(O)    :: beta, x
  complex(O) :: Pmm1n0
!
  integer    :: p, maxm, minm
  real(O)    :: prod, fact  
!
  x = cos(beta)  
  if (m == m1) then  
    fact   = ((1._O + x) / 2._O)**m
    Pmm1n0 = cmplx(fact,0.0,O)
  else
    maxm = max(m,m1)
    minm = min(m,m1)  
    prod = 1._O
    do p = 1, maxm - minm
      fact = sqrt(real(m + m1 + p,O) / real(p,O))
      prod = prod * fact
    end do  
    fact   = sqrt((1._O - x) / 2._O)
    Pmm1n0 = (-im)**(maxm - minm) * prod * fact**(maxm - minm)
    fact   = sqrt((1._O + x) / 2._O)
    Pmm1n0 = Pmm1n0 * fact**(maxm + minm)
  end if
!..................................................................................!
!       Take into account that: cos(beta/2) = -((1 + cos(beta))/2) for beta > 180  !
!       and sin(beta/2) = -((1 - cos(beta))/2) for beta < 0.                       !
!..................................................................................!
  if (beta > Pi .or. beta < 0) then
    Pmm1n0 = (-1._O)**(m + m1) * Pmm1n0
  end if 
end function Jacobi_n0
!***********************************************************************************
function dmm1n (beta, n, m, m1) result (d)
!-----------------------------------------------------------------------------------
!The routine computes the dmm1n coeffcients for m >= 0, m1 >= 0, and n >= max(m,m1)! 
!----------------------------------------------------------------------------------- 
  use parameters
  implicit none
  integer     :: n, m, m1
  real(O)     :: beta
  complex(O)  :: d
!
  integer     :: n0, l
  real(O)     :: x, fact1, fact2
  complex(O)  :: Pmm1n, Pmm1_m, Pmm1, Pmm1_p, Jacobi_n0
!
  x  = cos(beta)
  n0 = max(m,m1)
  d  = zero
  if ((m <= n) .and. (m1 <= n)) then
    if (n0 == 0) then
      if (n == 0) then
        Pmm1n = (1._O,0._O)
      else
        Pmm1_m = zero
        Pmm1   = (1._O,0._O)
        do l = 0, n - 1
          fact1  = real(2 * l + 1,O) * x / real(l + 1,O)
          fact2  = real(l,O) / real(l + 1,O)
          Pmm1_p = fact1 * Pmm1 - fact2 * Pmm1_m
          Pmm1_m = Pmm1
          Pmm1   = Pmm1_p
        end do  
        Pmm1n = Pmm1
      end if
    else
      if (n == n0) then
        Pmm1n = Jacobi_n0 (beta, m, m1)    
      else 
        Pmm1_m = zero
        Pmm1   = Jacobi_n0 (beta, m, m1)
        do l = n0, n - 1                        
          fact1  = real(l * (l + 1),O) * x - real(m * m1,O)
          fact1  = fact1 / sqrt(real((l + 1)**2 - m**2,O))
          fact1  = fact1 / sqrt(real((l + 1)**2 - m1**2,O))
          fact1  = fact1 * real(2 * l + 1,O) / real(l,O)
          fact2  = sqrt(real(l**2 - m**2,O)) * sqrt(real(l**2 - m1**2,O))
          fact2  = fact2 / sqrt(real((l + 1)**2 - m**2,O))
          fact2  = fact2 / sqrt(real((l + 1)**2 - m1**2,O))
          fact2  = fact2 * real(l + 1,O) / real(l,O)       
          Pmm1_p = fact1 * Pmm1 - fact2 * Pmm1_m
          Pmm1_m = Pmm1
          Pmm1   = Pmm1_p
        end do  
        Pmm1n = Pmm1
      end if
    end if
    d = im**(m1 - m) * Pmm1n
  end if  
end function dmm1n
!***********************************************************************************
subroutine dmm1nVect (beta, Nrank, m, m1, dmm1)
!-----------------------------------------------------------------------------------
! The routine computes the dmm1 vector coeffcients for m >= 0, m1 >= 0, and        !
! n = 0,...,Nrank.                                                                 !
!-----------------------------------------------------------------------------------  
  use parameters
  implicit none
  integer     :: Nrank, m, m1
  real(O)     :: beta
  complex(O)  :: dmm1(0:Nrank)
!
  integer     :: n0, n
  real(O)     :: x, fact1, fact2
  complex(O)  :: Jacobi_n0
!
  x  = cos(beta)
  n0 = max(m,m1)
  if (m <= Nrank .and. m1 <= Nrank) then
    if (n0 == 0) then     
      dmm1(0) = one      
      dmm1(1) = x * dmm1(0)
      do n = 1, Nrank - 1
        fact1 = real(2 * n + 1,O) * x / real(n + 1,O)
        fact2 = real(n,O) / real(n + 1,O)
        dmm1(n+1) = fact1 * dmm1(n) - fact2 * dmm1(n-1)
      end do
    else
      do n = 0, n0 - 1
        dmm1(n) = zero
      end do
      dmm1(n0) = im**(m1 - m) * Jacobi_n0 (beta, m, m1)      
      do n = n0, Nrank - 1                      
        fact1 = real(n * (n + 1),O) * x - real(m * m1,O)
        fact1 = fact1 / sqrt(real((n + 1)**2 - m**2,O))
        fact1 = fact1 / sqrt(real((n + 1)**2 - m1**2,O))
        fact1 = fact1 * real(2 * n + 1,O) / real(n,O)
        fact2 = sqrt(real(n**2 - m**2,O)) * sqrt(real(n**2 - m1**2,O))
        fact2 = fact2 / sqrt(real((n + 1)**2 - m**2,O))
        fact2 = fact2 / sqrt(real((n + 1)**2 - m1**2,O))
        fact2 = fact2 * real(n + 1,O) / real(n,O)
        dmm1(n+1) = fact1 * dmm1(n) - fact2 * dmm1(n-1)
      end do
    end if  
  else
    do n = 0, Nrank
      dmm1(n) = zero
    end do
  end if            
end subroutine dmm1nVect
!***********************************************************************************         
function coef_rotation (alfa, beta, gama, n, m, m1) result (Dnmm1)
!-----------------------------------------------------------------------------------
! The routine computes the rotation coefficients Dnmm1(alfa,beta,gama).            !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: n, m, m1
  real(O)    :: alfa, beta, gama
  complex(O) :: Dnmm1
!
  real(O)    :: arga, argg       
  complex(O) :: d, dmm1n
!
  if (m >= 0 .and. m1 >= 0) then
    d = dmm1n (beta, n, m, m1)  
  else if (m < 0 .and. m1 >= 0) then
    d = (-1._O)**n * dmm1n (beta + Pi, n, - m, m1)
  else if (m >= 0 .and. m1 < 0) then
    d = (-1._O)**n * dmm1n (beta + Pi, n, m, - m1)    
  else if (m < 0 .and. m1 < 0) then
    d = dmm1n (beta, n, - m, - m1)
  end if
  arga  = real(m,O) * alfa
  argg  = real(m1,O) * gama
  Dnmm1 = (-1._O)**(m + m1) * d * exp(im * arga) * exp(im * argg)
end function coef_rotation
!***********************************************************************************         
subroutine coef_rotationVect (alfa, beta, gama, Nrank, m, m1, Dmm1)
!------------------------------------------------------------------------------------
! The routine computes the vector rotation coefficients Dmm1(alfa,beta,gama)        !
! for n = 0,...,Nrank.                                                              !
!------------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: Nrank, m, m1
  real(O)    :: alfa, beta, gama
  complex(O) :: Dmm1(0:Nrank)
!
  integer    :: n
  real(O)    :: arga, argg
  complex(O) :: fact     
  complex(O),allocatable :: d(:)
!
  allocate (d(0:Nrank))
  if (m >= 0 .and. m1 >= 0) then
    call dmm1nVect (beta, Nrank, m, m1, d)
  else if (m < 0 .and. m1 >= 0) then  
    call dmm1nVect (beta + Pi, Nrank, - m, m1, d)
    do n = 0, Nrank
      d(n) = (-1._O)**n * d(n)
    end do
  else if (m >= 0 .and. m1 < 0) then
    call dmm1nVect (beta + Pi, Nrank, m, - m1, d)
    do n = 0, Nrank
      d(n) = (-1._O)**n * d(n)
    end do
  else if (m < 0 .and. m1 < 0) then
    call dmm1nVect (beta, Nrank, - m, - m1, d)
  end if
  arga  = real(m,O) * alfa
  argg  = real(m1,O) * gama
  fact  = (-1._O)**(m + m1) * exp(im * arga) * exp(im * argg)
  do n = 0, Nrank
    Dmm1(n) = fact * d(n)
  end do
  deallocate (d)
end subroutine coef_rotationVect
!***********************************************************************************         
subroutine coef_rotationMatr (alfa, beta, gama, Nrank, Mrank, Mrank1, D)
!-----------------------------------------------------------------------------------
! The routine computes the matrix rotation coefficients D(alfa,beta,gama)          !
! for m = -Mrank,Mrank, m1 = -Mrank1,Mrank1, n = 0,...,Nrank.                      !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: Nrank, Mrank, Mrank1
  real(O)    :: alfa, beta, gama
  complex(O) :: D(-Mrank:Mrank,-Mrank1:Mrank1,0:Nrank)
!
  integer    :: m, m1, n
  complex(O),allocatable :: Dmm1(:)
!
  allocate (Dmm1(0:Nrank))
  do m = - Mrank, Mrank
    do m1 = - Mrank1, Mrank1
      call coef_rotationVect (alfa, beta, gama, Nrank, m, m1, Dmm1)
      do n = 0, Nrank
        D(m,m1,n) = Dmm1(n)
      end do
    end do
  end do 
  deallocate (Dmm1)
end subroutine coef_rotationMatr
!***********************************************************************************
subroutine MatTransCmnn1 (index, wavenumber, z0, m, NN, c, NNP)
!-----------------------------------------------------------------------------------
! The routine computes the axial-translation matrix Cm(z0) for spherical wave      !
! functions and the azimuthal mode m. For index = 1, the translation coefficients  !
! involve the Bessel functions, while for index = 3, the translation coefficients  !
! involve the Hankel functions of the first order. Note that if index = 3, the     !
! routine must be called with positive values of the axial distance z0.            !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer     :: m, NN, NNP, index
  real(O)     :: z0, wavenumber
  complex(O)  :: c(0:NNP,0:2*NNP+1)
!                 
  integer     :: n, n1
  real(O)     :: f1, f2, f3, fm, fn ,fnm
  complex(O)  :: x 
  complex(O), allocatable :: j(:),jd(:)
!
  allocate (j(0:2*NN+1), jd(0:2*NN+1))
  if (m == 0) then 
    do n = 0, NN
      do n1 = 0, 2*NN + 1
        c(n,n1) = zero
      end do  
    end do
  end if        
  x = cmplx (wavenumber * z0,0.0,O)
  if (m == 0) then
    if (index == 1) then
      call besel_j (x, 2*NN + 1, j, jd)
    else if (index == 3) then
      call besel_h (x, 2*NN + 1, j, jd)      
    end if
    do n1 = 0, 2*NN + 1
      f1 = sqrt(real(2 * n1 + 1,O))       
      c(0,n1) = (-1._O)**n1 * f1 * j(n1)
    end do  
  else
    fm = sqrt(real(2 * m + 1,O) / 2._O / real(m,O))
    do n1 = m, 2*NN + 1 - m       
      f1 = fm * sqrt(real(n1 + m - 1,O) / real(2 * n1 - 1,O))
      f1 = f1 * sqrt(real(n1 + m,O) / real(2 * n1 + 1,O))     
      f2 = fm * sqrt(real(n1 - m + 1,O) / real(2 * n1 + 1,O))                
      f2 = f2 * sqrt(real(n1 - m + 2,O) / real(2 * n1 + 3,O))
      c(m,n1) = f1 * c(m-1,n1-1) + f2 * c(m-1,n1+1)
    end do  
  end if
  fm = sqrt(real(2 * m + 3,O))
  do n = m + 1, NN
    fn  = sqrt(real(2 * n - 1,O) * real(2 * n + 1,O))
    fnm = sqrt(real(n - m,O) * real(n + m,O))       
    if (n == m + 1) then          
      f3 = fn / fm
      f3 = f3 / fnm   
      c(n,m) = - f3 * c(n-1,m+1)
      do n1 = m + 1, 2*NN + 1 - n
        f2 = fn /  sqrt(real(2 * n1 - 1,O)  * real(2 * n1 + 1,O))
        f2 = f2 * (sqrt(real(n1 - m,O) * real(n1 + m,O)) / fnm)         
        f3 = fn /  sqrt(real(2 * n1 + 1,O) * real(2 * n1 + 3,O))
        f3 = f3 * (sqrt(real(n1 - m + 1,O) * real(n1 + m + 1,O)) / fnm)         
        c(n,n1) = f2 * c(n-1,n1-1) - f3 * c(n-1,n1+1)
      end do  
    else          
      f1 = sqrt(real(2 * n + 1,O) / real(2 * n - 3,O))
      f1 = f1 * (sqrt(real(n - m - 1,O) * real(n + m - 1,O)) / fnm)   
      f3 = fn / fm
      f3 = f3 / fnm   
      c(n,m) = f1 * c(n-2,m) - f3 * c(n-1,m+1)
      do n1 = m + 1, 2*NN + 1 - n           
        f2 = fn /  sqrt(real(2 * n1 - 1,O)  * real(2 * n1 + 1,O))
        f2 = f2 * (sqrt(real(n1 - m,O) * real(n1 + m,O)) / fnm) 
        f3 = fn /  sqrt(real(2 * n1 + 1,O) * real(2 * n1 + 3,O))
        f3 = f3 * (sqrt(real(n1 - m + 1,O) * real(n1 + m + 1,O)) / fnm)         
        c(n,n1) = f1 * c(n-2,n1) + f2 * c(n-1,n1-1) - f3 * c(n-1,n1+1)
      end do  
    end if
  end do   
  deallocate (j, jd)    
end subroutine MatTransCmnn1
!***********************************************************************************
subroutine MatTransAB_mnn1 (wavenumber, direct, z0, m, c, NNP, Nmax, Nmax1,         &
           AA, NP, MP)
!-----------------------------------------------------------------------------------
! The routine computes the axial-translation matrix Tm(z0) for vector spherical    !
! wave functions and the azimuthal mode m. For direct = t, the direct              !
! transformation matrix is computed, while for direct = f, the inverse             !
! transformation matrix is computed. The structure of the matrix is:               !
!                                 |  A    -B |                                     !                                             
!                            T =  |          |.                                    !                          
!                                 | -B     A |                                     !                            
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: m, NNP, Nmax, Nmax1, NP, MP
  real(O)    :: z0, wavenumber
  complex(O) :: AA(2*NP,2*MP), c(0:NNP,0:2*NNP+1)               
  logical    :: direct
!
  integer    :: n, n1, k, k1
  real(O)    :: xr, f1, f2, f3, f, s
  complex(O) :: A, B, x
!
  do k = 1, 2*Nmax
    do k1 = 1, 2*Nmax1
      AA(k,k1) = zero
    end do      
  end do
  xr = wavenumber * z0
  x  = cmplx(xr,0.0,O)  
  do k1 = 1, Nmax1
    if (m == 0) then
      n1 = k1
    else
      n1 = m + k1 - 1
    end if            
    f1 = sqrt(real(n1 - m + 1,O) / real(2 * n1 + 1,O))
    f1 = f1 * sqrt(real(n1 + m + 1,O) / real(2 * n1 + 3,O))
    f1 = f1 * xr / real(n1 + 1,O)
    if (m <= n1 - 1) then
      f2 = sqrt(real(n1 - m,O) / real(2 * n1 + 1,O))
      f2 = f2 * sqrt(real(n1 + m,O) / real(2 * n1 - 1,O))
      f2 = f2 * xr / real(n1,O)    
    end if
    f3  = real(m,O) / real(n1 * (n1 + 1),O)
    do k = 1, Nmax
      if (m == 0) then
        n = k
      else
        n = m + k - 1
      end if 
      if (m <= n1 - 1) then
        A = c(n,n1) + f1 * c(n,n1+1) + f2 * c(n,n1-1)
      else
        A = c(n,n1) + f1 * c(n,n1+1)
      end if              
      B = im * x * f3 * c(n,n1)   
      f = sqrt(real(n1 * (n1 + 1),O) / real(n * (n + 1),O))
      A = A * f
      B = B * f
      if (direct) then       
        AA(k,k1)       =   A
        AA(k,k1+Nmax1) = - B 
      else 
        s = (-1._O)**(n + n1) 
        AA(k,k1)       = s * A
        AA(k,k1+Nmax1) = s * B
      end if
    end do  
  end do      
  do k = 1, Nmax
    do k1 = 1, Nmax1
      AA(k+Nmax,k1)       = AA(k,k1+Nmax1)
      AA(k+Nmax,k1+Nmax1) = AA(k,k1)
    end do
  end do                
end subroutine MatTransAB_mnn1
!***********************************************************************************
subroutine MatTransAB_mn_m1n1 (index, wavenumber, x0, y0, z0, Mrank, Nrank, Nmax,   &
           Mrank1, Nrank1, Nmax1, AB, NP, MP)
!-----------------------------------------------------------------------------------
! The routine computes the axial-translation matrix T(x0,y0,z0) for vector         !
! spherical wave functions. For index = 1, the translation coefficients involve the!
! Bessel functions, while for index = 3, the translation coefficients involve the  !
! Hankel functions of the first order. The structure of the matrix is              !
!                                |  A     B |                                      ! 
!                            T = |          |.                                     !                           
!                                |  B     A |                                      ! 
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax, Mrank1, Nrank1, Nmax1, NP, MP, index
  real(O)    :: x0, y0, z0, wavenumber
  complex(O) :: AB(2*NP,2*MP)   
!
  integer    :: n, n1, m, m1, k, k1, m2, l, ml, l1, m1l, N_1, N_0, NN, MM
  real(O)    :: r, alfa, beta, z, xr, f1, f2, f3, f4, f, fn, fm, fnm
  complex(O) :: A, B, Dn1m2m1, Dn1_m2m1, AA, BB, Dnmm2, Dnm_m2, x, f3i
  complex(O), allocatable :: j(:), jd(:), c(:,:), DP(:,:,:), DM(:,:,:)
!
  NN = max(Nrank,Nrank1) + 5
  MM = max(Nrank,Nrank1)    
  do k = 1, 2*Nmax
    do k1 = 1, 2*Nmax1
      AB(k,k1) = zero
    end do  
  end do     
  allocate (c(0:NN,0:2*NN+1), j(0:2*NN+1), jd(0:2*NN+1))  
  do n = 0, NN
    do n1 = 0, 2*NN + 1
      c(n,n1) = zero
    end do  
  end do
  call T_cartesian_spherical (x0, y0, z0, r, beta, alfa)  
  xr = wavenumber * r
  x  = cmplx(xr,0.0,O)
  z  = 0._O  
  allocate (DP(-MM:MM,-MM:MM,0:NN), DM(-MM:MM,-MM:MM,0:NN))
  call coef_rotationMatr (alfa, beta, z, NN, MM, MM, DP)
  call coef_rotationMatr (z, - beta, - alfa, NN, MM, MM, DM)
  do m2 = 0, MM
    if (m2 == 0) then
      if (index == 1) then
        call besel_j (x, 2*NN+1, j, jd)         
      else if (index == 3) then
        call besel_h(x, 2*NN+1, j, jd)      
      end if                                     
      do n1 = 0, 2*NN + 1 
        f1 = sqrt(real(2 * n1 + 1,O))
        c(0,n1) = (-1._O)**n1 * f1 * j(n1)
      end do  
    else
      fm = sqrt(real(2 * m2 + 1,O) / 2._O / real(m2,O))
      do n1 = m2, 2*NN + 1 - m2     
        f1 = fm * sqrt(real(n1 + m2 - 1,O) / real(2 * n1 - 1,O))
        f1 = f1 * sqrt(real(n1 + m2,O) / real(2 * n1 + 1,O))            
        f2 = fm * sqrt(real(n1 - m2 + 1,O) / real(2 * n1 + 1,O))
        f2 = f2 * sqrt(real(n1 - m2 + 2,O) / real(2 * n1 + 3,O))
        c(m2,n1) = f1 * c(m2-1,n1-1) + f2 * c(m2-1,n1+1)
      end do  
    end if
    fm = sqrt(real(2 * m2 + 3,O))      
    do n = m2 + 1, NN
      fn  = sqrt(real(2 * n - 1,O) * real(2 * n + 1,O))
      fnm = sqrt(real(n - m2,O) * real(n + m2,O))
      if (n == m2 + 1) then         
        f3 = fn / fm
        f3 = f3 / fnm                       
        c(n,m2) = - f3 * c(n-1,m2+1)
        do n1 = m2 + 1, 2*NN + 1 - n
          f2 = fn /  sqrt(real(2 * n1 - 1,O) * real(2 * n1 + 1,O))
          f2 = f2 * (sqrt(real(n1 - m2,O) * real(n1 + m2,O)) / fnm)               
          f3 = fn /  sqrt(real(2 * n1 + 1,O) * real(2 * n1 + 3,O))
          f3 = f3 * (sqrt(real(n1 - m2 + 1,O) * real(n1 + m2 + 1,O)) / fnm)               
          c(n,n1) = f2 * c(n-1,n1-1) - f3 * c(n-1,n1+1)
        end do  
      else              
        f1 = sqrt(real(2 * n + 1,O) / real(2 * n - 3,O))
        f1 = f1 * (sqrt(real(n - m2 - 1,O) * real(n + m2 - 1,O)) / fnm)                                 
        f3 = fn / fm
        f3 = f3 / fnm           
        c(n,m2) = f1 * c(n-2,m2) - f3 * c(n-1,m2+1)
        do n1 = m2 + 1, 2*NN + 1 - n
          f2 = fn /  sqrt(real(2 * n1 - 1,O) * real(2 * n1 + 1,O))
          f2 = f2 * (sqrt(real(n1 - m2,O) * real(n1 + m2,O)) / fnm)               
          f3 = fn /  sqrt(real(2 * n1 + 1,O) * real(2 * n1 + 3,O))
          f3 = f3 * (sqrt(real(n1 - m2 + 1,O) * real(n1 + m2 + 1,O)) / fnm)                               
          c(n,n1) = f1 * c(n-2,n1) + f2 * c(n-1,n1-1) - f3 * c(n-1,n1+1)
        end do  
      end if
    end do      
    do m1 = 0, Mrank1     
      if (m1 == 0) then   
        do k1 = 1, Nrank1
          n1 = k1
          if (m2 <= n1) then            
            Dn1m2m1  = DM( m2, m1, n1)
            Dn1_m2m1 = DM(-m2, m1, n1)
            f1 = sqrt(real(n1 - m2 + 1,O) / real(2 * n1 + 1,O))
            f1 = f1 * sqrt(real(n1 + m2 + 1,O) / real(2 * n1 + 3,O))
            f1 = f1 * xr / real(n1 + 1,O)           
            if (m2 <= n1 - 1) then                              
              f2 = sqrt(real(n1 - m2,O) / real(2 * n1 + 1,O))
              f2 = f2 * sqrt(real(n1 + m2,O) / real(2 * n1 - 1,O))
              f2 = f2 * xr / real(n1,O)
            end if
            f3  = real(m2,O) / real(n1 * (n1 + 1),O)
            f3i = im * x * f3
            f4  = sqrt(real(n1 * (n1 + 1),O))
            do m = 0, Mrank                             
              if (m == 0) then                      
                do k = 1, Nrank
                  n = k                 
                  Dnmm2  = DP(m,  m2, n)
                  Dnm_m2 = DP(m,- m2, n)
                  if (m2 <= n1 - 1) then
                    A = c(n,n1) + f1 * c(n,n1+1) + f2 * c(n,n1-1)
                  else                                                    
                    A = c(n,n1) + f1 * c(n,n1+1)
                  end if                                  
                  B = f3i * c(n,n1)                               
                  f = f4 / sqrt(real(n * (n + 1),O))
                  A = A * f
                  B = B * f
                  if (m2 == 0) then
                    AA = Dnmm2 * A * Dn1m2m1
                    BB = Dnmm2 * B * Dn1m2m1
                  else
                    AA = Dnmm2 * A * Dn1m2m1 + Dnm_m2 * A * Dn1_m2m1
                    BB = Dnmm2 * B * Dn1m2m1 - Dnm_m2 * B * Dn1_m2m1
                  end if
                  AB(k,k1)       = AB(k,k1) + AA
                  AB(k,k1+Nmax1) = AB(k,k1+Nmax1) + BB    
                end do  
              else
                N_0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
                ml  = m
                do l = 1, 2                          
                  do k = 1, Nrank - m + 1
                    n = m + k - 1                    
                    Dnmm2  = DP(- ml,  m2, n)
                    Dnm_m2 = DP(- ml,- m2, n)
                    if (m2 <= n1 - 1) then
                      A = c(n,n1) + f1 * c(n,n1+1) + f2 * c(n,n1-1)
                    else                                                          
                      A = c(n,n1) + f1 * c(n,n1+1)
                    end if                                
                    B = f3i * c(n,n1)                             
                    f = f4 / sqrt(real(n * (n + 1),O))
                    A = A * f
                    B = B * f
                    if (m2 == 0) then
                      AA = Dnmm2 * A * Dn1m2m1
                      BB = Dnmm2 * B * Dn1m2m1
                    else
                      AA = Dnmm2 * A * Dn1m2m1 + Dnm_m2 * A * Dn1_m2m1
                      BB = Dnmm2 * B * Dn1m2m1 - Dnm_m2 * B * Dn1_m2m1
                    end if
                    AB(k+N_0,k1)       = AB(k+N_0,k1) + AA
                    AB(k+N_0,k1+Nmax1) = AB(k+N_0,k1+Nmax1) + BB
                  end do  
                  N_0 =   N_0 + Nrank - m + 1
                  ml  = - m
                end do      
              end if
            end do                            
          end if
        end do  
      else          
        N_1 = Nrank1 + (m1 - 1) * (2 * Nrank1 - m1 + 2)
        m1l = m1
        do l1 = 1, 2             
          do k1 = 1, Nrank1 - m1 + 1
            n1 = m1 + k1 - 1
            if (m2 <= n1) then             
              Dn1m2m1  = DM(  m2, - m1l, n1)
              Dn1_m2m1 = DM(- m2, - m1l, n1)          
              f1 = sqrt(real(n1 - m2 + 1,O) / real(2 * n1 + 1,O))
              f1 = f1 * sqrt(real(n1 + m2 + 1,O) / real(2 * n1 + 3,O))
              f1 = f1 * xr / real(n1 + 1,O)         
              if (m2 <= n1 - 1) then                            
                f2 = sqrt(real(n1 - m2,O) / real(2 * n1 + 1,O))
                f2 = f2 * sqrt(real(n1 + m2,O) / real(2 * n1 - 1,O))
                f2 = f2 * xr / real(n1,O)
              end if
              f3  = real(m2,O) / real(n1 * (n1 + 1),O)
              f3i = im * x * f3
              f4  = sqrt(real(n1 * (n1 + 1),O))
              do m = 0, Mrank
                if (m == 0) then                                  
                  do k = 1, Nrank
                    n = k                    
                    Dnmm2  = DP(m,  m2, n)
                    Dnm_m2 = DP(m,- m2, n)
                    if (m2 <= n1 - 1) then
                      A = c(n,n1) + f1 * c(n,n1+1) + f2 * c(n,n1-1)
                    else                                                          
                      A = c(n,n1) + f1 * c(n,n1+1)
                    end if                                
                    B = f3i * c(n,n1)                             
                    f = f4 / sqrt(real(n * (n + 1),O))
                    A = A * f
                    B = B * f
                    if (m2 == 0) then
                      AA = Dnmm2 * A * Dn1m2m1
                      BB = Dnmm2 * B * Dn1m2m1
                    else
                      AA = Dnmm2 * A * Dn1m2m1 + Dnm_m2 * A * Dn1_m2m1
                      BB = Dnmm2 * B * Dn1m2m1 - Dnm_m2 * B * Dn1_m2m1
                    end if
                    AB(k,k1+N_1)       = AB(k,k1+N_1) + AA
                    AB(k,k1+N_1+Nmax1) = AB(k,k1+N_1+Nmax1) + BB    
                  end do  
                else
                  N_0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
                  ml  = m
                  do l = 1, 2                                   
                    do k = 1, Nrank - m + 1
                      n = m + k - 1                      
                      Dnmm2  = DP(- ml,  m2, n)
                      Dnm_m2 = DP(- ml,- m2, n)
                      if (m2 <= n1 - 1) then
                        A = c(n,n1) + f1 * c(n,n1+1) + f2 * c(n,n1-1)
                      else                                                        
                        A = c(n,n1) + f1 * c(n,n1+1)
                      end if                              
                      B = f3i * c(n,n1)                           
                      f = f4 / sqrt(real(n * (n + 1),O))
                      A = A * f
                      B = B * f
                      if (m2 == 0) then
                        AA = Dnmm2 * A * Dn1m2m1
                        BB = Dnmm2 * B * Dn1m2m1
                      else
                        AA = Dnmm2 * A * Dn1m2m1 + Dnm_m2 * A * Dn1_m2m1
                        BB = Dnmm2 * B * Dn1m2m1 - Dnm_m2 * B * Dn1_m2m1
                      end if
                      AB(k+N_0,k1+N_1)       = AB(k+N_0,k1+N_1) + AA
                      AB(k+N_0,k1+N_1+Nmax1) = AB(k+N_0,k1+N_1+Nmax1) + BB
                    end do  
                    N_0 =   N_0 + Nrank - m + 1
                    ml  = - m
                  end do      
                end if
              end do                   
            end if
          end do  
          N_1 =   N_1 + Nrank1 - m1 + 1
          m1l = - m1
        end do          
      end if
    end do  
  end do  
  do k = 1, Nmax
    do k1 = 1, Nmax1
      AB(k+Nmax,k1)       = AB(k,k1+Nmax1)
      AB(k+Nmax,k1+Nmax1) = AB(k,k1)
     end do  
  end do
  deallocate (c, j, jd)  
  deallocate (DM, DP)   
end subroutine MatTransAB_mn_m1n1
!***********************************************************************************
subroutine matrix_inverse (Mrank, Nrank, Nmax, Mrank1, Nrank1, Nmax1, AB, NP, MP)
!-----------------------------------------------------------------------------------
! The routine computes the inverse T(-x0,-y0,-z0) of the general translation       !
! matrix T(x0,y0,z0) for vector spherical wave functions.                          !
!-----------------------------------------------------------------------------------       
  use parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax, Mrank1, Nrank1, Nmax1, NP, MP
  complex(O) :: AB(2*NP,2*MP)
!
  integer    :: n, n1, m, m1, k, k1, l, l1, N_1, N_0
  real(O)    :: s
!
  do m = 0, Mrank
    if (m == 0) then
      do k = 1, Nrank
        n = k    
        do m1 = 0, Mrank1
          if (m1 == 0) then
            do k1 = 1, Nrank1
              n1 = k1 
              s  = (- 1._O)**(n + n1)       
              AB(k,k1)       =   s * AB(k,k1)
              AB(k,k1+Nmax1) = - s * AB(k,k1+Nmax1)
            end do  
          else        
            N_1 = Nrank1 + (m1 - 1) * (2 * Nrank1 - m1 + 2)
            do l1 = 1, 2
              do k1 = 1, Nrank1 - m1 + 1
                n1 = m1 + k1 - 1  
                s  = ( -1._O)**(n + n1)      
                AB(k,k1+N_1)       =   s * AB(k,k1+N_1)
                AB(k,k1+N_1+Nmax1) = - s * AB(k,k1+N_1+Nmax1)
              end do  
              N_1 = N_1 + Nrank1 - m1 + 1
            end do          
          end if
        end do  
      end do  
    else
      N_0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      do l = 1, 2
        do k = 1, Nrank - m + 1
          n = m + k - 1    
          do m1 = 0, Mrank1
            if (m1 == 0) then
              do k1 = 1, Nrank1
                n1 = k1
                s  = ( -1._O)**(n + n1)
                AB(k+N_0,k1)       =   s * AB(k+N_0,k1)
                AB(k+N_0,k1+Nmax1) = - s * AB(k+N_0,k1+Nmax1)
              end do  
            else        
              N_1 = Nrank1 + (m1 - 1) * (2 * Nrank1 - m1 + 2)
              do l1 = 1, 2
                do k1 = 1, Nrank1 - m1 + 1
                  n1 = m1 + k1 - 1      
                  s  = ( -1._O)**(n + n1) 
                  AB(k+N_0,k1+N_1)       =   s * AB(k+N_0,k1+N_1)
                  AB(k+N_0,k1+N_1+Nmax1) = - s * AB(k+N_0,k1+N_1+Nmax1)
                end do  
                N_1 = N_1 + Nrank1 - m1 + 1
              end do           
            end if
          end do  
        end do  
        N_0 = N_0 + Nrank - m + 1
      end do
    end if
  end do  
  do k = 1, Nmax
    do k1 = 1, Nmax1
      AB(k+Nmax,k1)       = AB(k,k1+Nmax1)
      AB(k+Nmax,k1+Nmax1) = AB(k,k1)
     end do      
  end do  
end subroutine matrix_inverse
!***********************************************************************************
subroutine MatRot_R_mn_m1n1_single (alfa, beta, gama, Mrank, Nrank, Nmax, Mrank1,   &
           Nrank1, Nmax1, R, NP, MP, matrix)
!-----------------------------------------------------------------------------------
! The routine computes the individual rotation matrix R(alfa,beta,gama) for        !
! vector spherical wave functions.                                                 !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax, Mrank1, Nrank1, Nmax1, NP, MP  
  real(O)    :: alfa, beta, gama
  complex(O) :: R(NP,MP)
  logical    :: matrix
!    
  integer    :: n, m, m1, k, k1, l, ml, l1, m1l, N_1, N_0
  complex(O) :: coef_rotation
  complex(O),allocatable :: D(:,:,:)
!
  if (matrix) then
    allocate (D(-Mrank:Mrank,-Mrank1:Mrank1,0:Nrank1))
    call coef_rotationMatr (alfa, beta, gama, Nrank1, Mrank, Mrank1, D)
  end if
  do k = 1, Nmax
    do k1 = 1, Nmax1
      R(k,k1) = zero
    end do  
  end do    
  do m = 0, Mrank
    if (m == 0) then
      do k = 1, Nrank
        n = k    
        do m1 = 0, Mrank1
          if (m1 == 0) then                 
            k1 = n
            if (k1 <= Nrank1) then
              if (matrix) then                        
                R(k,k1) = D(m, m1, n)
              else
                R(k,k1) = coef_rotation (alfa, beta, gama, n, m, m1)
              end if     
            end if
          else        
            N_1 = Nrank1 + (m1 - 1) * (2 * Nrank1 - m1 + 2)
            m1l = m1
            do l1 = 1, 2                          
              k1 = n - m1 + 1
              if (1 <= k1 .and. k1 <= Nrank1 - m1 + 1) then
                if (matrix) then                        
                  R(k,k1+N_1) = D(m, - m1l, n)
                else
                  R(k,k1+N_1) = coef_rotation (alfa, beta, gama, n, m, - m1l)
                end if
              end if
              N_1 =   N_1 + Nrank1 - m1 + 1
              m1l = - m1
            end do          
          end if
        end do  
      end do  
    else
      N_0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      ml  = m
      do l = 1, 2
        do k = 1, Nrank - m + 1
          n = m + k - 1    
          do m1 = 0, Mrank1
            if (m1 == 0) then                     
              k1 = n
              if (k1 <= Nrank1) then
                if (matrix) then                        
                  R(k+N_0,k1) = D(- ml, m1, n)
                else
                  R(k+N_0,k1) = coef_rotation (alfa, beta, gama, n, - ml, m1)
                end if
              end if           
            else        
              N_1 = Nrank1 + (m1 - 1) * (2 * Nrank1 - m1 + 2)
              m1l = m1
              do l1 = 1, 2                
                k1 = n - m1 + 1
                if (1 <= k1 .and. k1 <= Nrank1 - m1 + 1) then
                  if (matrix) then                               
                    R(k+N_0,k1+N_1) = D(- ml, - m1l, n)
                  else
                    R(k+N_0,k1+N_1) = coef_rotation (alfa, beta, gama, n, - ml, - m1l)
                  end if
                end if                
                N_1 =   N_1 + Nrank1 - m1 + 1
                m1l = - m1                
              end do         
            end if
          end do  
        end do  
        N_0 =   N_0 + Nrank - m + 1
        ml  = - m
      end do     
    end if
  end do
  if (matrix) deallocate (D)         
end subroutine MatRot_R_mn_m1n1_single
!***********************************************************************************
subroutine MatRot_R_mn_m1n1 (alfa, beta, gama, Mrank, Nrank, Nmax, Mrank1, Nrank1,  &
           Nmax1, R, NP, MP, matrix)
!-----------------------------------------------------------------------------------
! The routine computes the rotation matrix R(alfa,beta,gama) for vector spherical  !
! wave functions. The structure of the matrix is:                                  !
!                                     |  Ri   0  |                                 !                                              
!                                 R = |          |,                                !                           
!                                     |  0    Ri |                                 ! 
! where Ri is the individual rotation matrix.                                      !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax, Mrank1, Nrank1, Nmax1, NP, MP  
  real(O)    :: alfa, beta, gama
  complex(O) :: R(2*NP,2*MP)
  logical    :: matrix
!      
  integer    :: n, m, m1, k, k1, l, ml, l1, m1l, N_1, N_0
  complex(O) :: coef_rotation
  complex(O),allocatable :: D(:,:,:)
!
  if (matrix) then
    allocate (D(-Mrank:Mrank,-Mrank1:Mrank1,0:Nrank1))
    call coef_rotationMatr (alfa, beta, gama, Nrank1, Mrank, Mrank1, D)
  end if
  do k = 1, 2*Nmax
    do k1 = 1, 2*Nmax1
      R(k,k1) = zero
    end do  
  end do  
  do m = 0, Mrank
    if (m == 0) then
      do k = 1, Nrank
        n = k    
        do m1 = 0, Mrank1
          if (m1 == 0) then                 
            k1 = n
            if (k1 <= Nrank1) then  
              if (matrix) then               
                R(k,k1) = D(m, m1, n)
              else
                R(k,k1) = coef_rotation (alfa, beta, gama, n, m, m1)
              end if                                                        
            end if             
          else        
            N_1 = Nrank1 + (m1 - 1) * (2 * Nrank1 - m1 + 2)
            m1l = m1
            do l1 = 1, 2                          
              k1 = n - m1 + 1
              if (1 <= k1 .and. k1 <= Nrank1 - m1 + 1) then                    
                if (matrix) then
                  R(k,k1+N_1) = D(m, - m1l, n)
                else
                  R(k,k1+N_1) = coef_rotation (alfa, beta, gama, n, m, - m1l)
                end if
              end if                
              N_1 =   N_1 + Nrank1 - m1 + 1
              m1l = - m1
            end do 
          end if 
        end do  
      end do  
    else
      N_0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      ml  = m
      do l = 1, 2
        do k = 1, Nrank - m + 1
          n = m + k - 1    
          do m1 = 0, Mrank1
            if (m1 == 0) then                     
              k1 = n
              if (k1 <= Nrank1) then
                if (matrix) then                                                  
                  R(k+N_0,k1) = D(- ml, m1, n)
                else
                  R(k+N_0,k1) = coef_rotation (alfa, beta, gama, n, - ml, m1)
                end if
              end if              
            else        
              N_1 = Nrank1 + (m1 - 1) * (2 * Nrank1 - m1 + 2)
              m1l = m1
              do l1 = 1, 2                          
                k1 = n - m1 + 1
                if (1 <=k1 .and. k1 <= Nrank1 - m1 + 1) then
                  if (matrix) then                                                                
                    R(k+N_0,k1+N_1) = D(- ml, - m1l, n)
                  else
                    R(k+N_0,k1+N_1) = coef_rotation (alfa, beta, gama, n, - ml, - m1l)
                  end if
                end if                   
                N_1 =   N_1 + Nrank1 - m1 + 1
                m1l = - m1
              end do
            end if
          end do  
        end do  
        N_0 =   N_0 + Nrank - m + 1
        ml  = - m
      end do
    end if
  end do  
  do k = 1, Nmax
    do k1 = 1, Nmax1
      R(k+Nmax,k1+Nmax1) = R(k,k1)
    end do   
  end do  
  if (matrix) deallocate (D)
end subroutine MatRot_R_mn_m1n1
!***********************************************************************************
subroutine MatTransRot_TR_mn_m1n1 (index, wavenumber, x0, y0, z0, alfa, beta, gama,  &
           Mrank, Nrank, Nmax, Mrank1, Nrank1, Nmax1, TR, NP, MP)            
!-----------------------------------------------------------------------------------
! The routine computes the transformation matrix                                   !
!                    TR = T(x0,y0,z0) * R(alfa,beta,gama),                         !
! where T is a translation matrix and R is a rotation matrix.                      !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer     :: index, Mrank, Nrank, Nmax, Mrank1, Nrank1, Nmax1, NP, MP
  complex(O)  :: TR(2*NP,2*MP)
  real(O)     :: wavenumber, x0, y0, z0, alfa, beta, gama
!
  integer     :: i, j, k, Nmaxim
  complex(O)  :: sum, sum1
  complex(O), allocatable :: R(:,:)
!
  Nmaxim = max(Nmax,Nmax1)
  allocate (R(Nmaxim,Nmaxim))
  if (Nmax1 >= Nmax) then
    call MatTransAB_mn_m1n1 (index, wavenumber, x0, y0, z0, Mrank, Nrank,           &
         Nmax, Mrank1, Nrank1, Nmax1, TR, NP, MP)  
    call MatRot_R_mn_m1n1_single (alfa, beta, gama,Mrank1, Nrank1, Nmax1,           &
         Mrank1, Nrank1, Nmax1, R, Nmaxim, Nmaxim, .true.)   
    do i = 1, Nmax
      do j = 1, Nmax1
        sum  = zero
        sum1 = zero
        do k = 1, Nmax1
          sum  = sum  + TR(i,k+Nmax1) * R(k,j)
          sum1 = sum1 + TR(i+Nmax,k+Nmax1) * R(k,j)
        end do
        TR(i,j)      = sum1  
        TR(i+Nmax,j) = sum          
      end do  
    end do
  else
    call MatTransAB_mn_m1n1 (index, wavenumber, x0, y0, z0, Mrank, Nrank,           &
         Nmax, Mrank, Nrank, Nmax, TR, NP, MP)  
    call MatRot_R_mn_m1n1_single (alfa, beta, gama, Mrank, Nrank, Nmax,             &
         Mrank1, Nrank1, Nmax1, R, Nmaxim, Nmaxim, .true.)    
    do i = 1, Nmax
      do j = 1, Nmax1
        sum  = zero
        sum1 = zero
        do k = 1, Nmax
          sum  = sum  + TR(i,k+Nmax) * R(k,j)
          sum1 = sum1 + TR(i+Nmax,k+Nmax) * R(k,j)
        end do  
        TR(i,j)      = sum1
        TR(i+Nmax,j) = sum          
      end do  
    end do
  end if
  do i = 1, Nmax
    do j = 1, Nmax1
      TR(i,j+Nmax1)      = TR(i+Nmax,j)
      TR(i+Nmax,j+Nmax1) = TR(i,j)
    end do  
  end do
  deallocate (R)  
end subroutine MatTransRot_TR_mn_m1n1
!***********************************************************************************                               
subroutine MatRotTrans_RT_mn_m1n1 (index, wavenumber, x0, y0, z0, alfa, beta, gama, &
           Mrank, Nrank, Nmax, Mrank1, Nrank1, Nmax1, RT, NP, MP)           
!-----------------------------------------------------------------------------------
! The routine computes the transformation matrix                                   !
!           RT = R(-gama,-beta,-alfa) * T(-x0,-y0,-z0),                            !
! where T is a translation matrix and R is a rotation matrix.                      !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer     :: index, Mrank, Nrank, Nmax, Mrank1, Nrank1, Nmax1, NP, MP
  real(O)     :: wavenumber, x0, y0, z0, alfa, beta, gama
  complex(O)  :: RT(2*NP,2*MP)
!   
  integer         :: i, j, k, Nmaxim
  complex(O)  :: sum, sum1
  complex(O), allocatable :: R(:,:)
!
  Nmaxim = max(Nmax,Nmax1)
  allocate (R(Nmaxim,Nmaxim))
  if (Nmax >= Nmax1) then
    call MatRot_R_mn_m1n1_single (- gama, - beta, - alfa, Mrank, Nrank, Nmax,       &
         Mrank, Nrank, Nmax, R, Nmaxim, Nmaxim, .true.)
    call MatTransAB_mn_m1n1 (index, wavenumber, - x0, - y0, - z0, Mrank, Nrank,     &
         Nmax, Mrank1, Nrank1, Nmax1, RT, NP, MP)         
    do i = 1, Nmax
      do j = 1, Nmax1
        sum  = zero
        sum1 = zero
        do k = 1, Nmax
          sum  = sum  + R(i,k) * RT(k,j+Nmax1)
          sum1 = sum1 + R(i,k) * RT(k+Nmax,j+Nmax1)
        end do  
        RT(i,j)      = sum1
        RT(i+Nmax,j) = sum          
      end do  
    end do
  else
    call MatRot_R_mn_m1n1_single (- gama, - beta, - alfa, Mrank, Nrank, Nmax,       &
         Mrank1, Nrank1, Nmax1, R, Nmaxim, Nmaxim, .true.)
    call MatTransAB_mn_m1n1 (index, wavenumber, - x0, - y0, - z0, Mrank1, Nrank1,   &
         Nmax1, Mrank1, Nrank1, Nmax1, RT, NP, MP)              
    do i = 1, Nmax
      do j = 1, Nmax1
        sum  = zero
        sum1 = zero
        do k = 1, Nmax1
          sum  = sum  + R(i,k) * RT(k,j+Nmax1)
          sum1 = sum1 + R(i,k) * RT(k+Nmax1,j+Nmax1)
        end do  
        RT(i,j)      = sum1
        RT(i+Nmax,j) = sum          
      end do  
    end do
  end if         
  do i = 1, Nmax
    do j = 1, Nmax1
      RT(i,j+Nmax1)      = RT(i+Nmax,j)
      RT(i+Nmax,j+Nmax1) = RT(i,j)
    end do  
  end do
  deallocate (R)  
end subroutine MatRotTrans_RT_mn_m1n1
!***********************************************************************************
subroutine MatRotTransRot_RTR_mn_m1n1 (index, wavenumber, x, y, z, alfa, beta, gama, &
           Mrank, Nrank, Nmax, x1, y1, z1, alfa1, beta1, gama1, Mrank1, Nrank1,      &
           Nmax1, RTR, NP, MP)                         
!-----------------------------------------------------------------------------------               
! The routine computes the transformation matrix                                   ! 
!             RTR = R(-gama,-beta,-alfa) * T(r1-r) *  R1(alfa1,beta1,gama1),       !
! where T is a translation matrix and R is a rotation matrix.                      !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer     :: index, Mrank, Nrank, Nmax, Mrank1, Nrank1, Nmax1, NP, MP
  real(O)     :: wavenumber, x, y, z, alfa, beta, gama, x1, y1, z1, alfa1,          &
                 beta1, gama1
  complex(O)  :: RTR(2*NP,2*MP)
!         
  integer         :: i, j, k, l, Nmaxim
  real(O)     :: dx, dy, dz
  complex(O)  :: sum, sum1
  complex(O), allocatable :: R(:,:),R1(:,:)     
!
  Nmaxim = max(Nmax,Nmax1)
  allocate (R(Nmaxim,Nmaxim), R1(Nmaxim,Nmaxim))
  if (Nmax >= Nmax1) then
    call MatRot_R_mn_m1n1_single (- gama, - beta, - alfa, Mrank, Nrank, Nmax,       &
         Mrank, Nrank, Nmax, R, Nmaxim, Nmaxim, .true.)
    dx = x1 - x
    dy = y1 - y
    dz = z1 - z
    call MatTransAB_mn_m1n1 (index, wavenumber, dx, dy, dz, Mrank, Nrank,           &
         Nmax, Mrank, Nrank, Nmax, RTR, NP, MP)
    call MatRot_R_mn_m1n1_single (alfa1, beta1, gama1, Mrank, Nrank, Nmax,          &
         Mrank1, Nrank1, Nmax1, R1, Nmaxim, Nmaxim, .true.)           
    do i = 1, Nmax
      do j = 1, Nmax1
        sum  = zero
        sum1 = zero
        do l = 1, Nmax
          do k = 1, Nmax
            sum  = sum  + R(i,k) * RTR(k,l+Nmax) * R1(l,j)
            sum1 = sum1 + R(i,k) * RTR(k+Nmax,l+Nmax) * R1(l,j)
          end do  
        end do
        RTR(i,j)      = sum1
        RTR(i+Nmax,j) = sum         
      end do  
    end do 
  else
    call MatRot_R_mn_m1n1_single (- gama, - beta, - alfa, Mrank, Nrank, Nmax,       &
         Mrank1, Nrank1, Nmax1, R, Nmaxim, Nmaxim, .true.)
    dx = x1 - x
    dy = y1 - y
    dz = z1 - z
    call MatTransAB_mn_m1n1 (index, wavenumber, dx, dy, dz, Mrank1, Nrank1,         &
         Nmax1, Mrank1, Nrank1, Nmax1, RTR, NP, MP)
    call MatRot_R_mn_m1n1_single (alfa1, beta1, gama1, Mrank1, Nrank1, Nmax1,       &
         Mrank1, Nrank1, Nmax1, R1, Nmaxim, Nmaxim, .true.)          
    do i = 1, Nmax
      do j = 1, Nmax1
        sum  = zero
        sum1 = zero
        do l = 1, Nmax1
          do k = 1, Nmax1
            sum  = sum  + R(i,k) * RTR(k,l+Nmax1) * R1(l,j)
            sum1 = sum1 + R(i,k) * RTR(k+Nmax1,l+Nmax1) * R1(l,j)
          end do  
        end do
        RTR(i,j)      = sum1
        RTR(i+Nmax,j) = sum         
      end do  
    end do 
  end if
  do i = 1, Nmax
    do j = 1, Nmax1
      RTR(i,j+Nmax1)      = RTR(i+Nmax,j)
      RTR(i+Nmax,j+Nmax1) = RTR(i,j)
    end do  
  end do
  deallocate (R, R1)  
end subroutine MatRotTransRot_RTR_mn_m1n1
!***********************************************************************************
subroutine identity_transformation (Mrank, Nrank, Nmax, Mrank1, Nrank1, Nmax1,      &
           T, NP, MP)
  use parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax, Mrank1, Nrank1, Nmax1, NP, MP   
  complex(O) :: T(2*NP,2*MP)
!
  integer    :: i, j, m, m1, n, n1, k, k1, N_0, N_1, l
!
  do i = 1, 2*Nmax
    do j = 1, 2*Nmax1
      T(i,j) = zero
    end do
  end do
  do m = 0, Mrank
    if (m == 0) then
      m1 = m
      do k = 1, Nrank
        n = k                           
        do k1 = 1, Nrank1
          n1 = k1
          if (n == n1) then                
            T(k,k1)            = one
            T(k+Nmax,k1+Nmax1) = one
          end if 
        end do          
      end do                                
    else
      N_0 = Nrank + (m - 1) * (2 * Nrank - m + 2)        
      m1  = m
      N_1 = Nrank1 + (m1 - 1) * (2 * Nrank1 - m1 + 2)
      if (m1 <= Mrank1) then        
        do l = 1, 2     
          do k = 1, Nrank - m + 1
            n = m + k - 1                                 
            do k1 = 1, Nrank1 - m1 + 1
              n1 = m1 + k1 - 1
              if (n == n1) then
                T(k+N_0,k1+N_1)            = one
                T(k+N_0+Nmax,k1+N_1+Nmax1) = one
              end if
            end do                
          end do          
          N_1 = N_1 + Nrank1 - m1 + 1
          N_0 = N_0 + Nrank - m + 1      
        end do
      end if
    end if
  end do  
end subroutine identity_transformation    
!***********************************************************************************
function lnfactorial (n) result (lnfact)
  use parameters
  implicit none
  integer :: n, k
  real(O) :: lnfact, sum, kr
!
  if (n < 0) then
    print "(/,2x,'Error in subroutine lnfactorial in module AdditonTh:')" 
    print "(  2x,'the argument of the factorial function is negative; ')"              
    stop    
  else if (n > 1) then
    sum = 0._O
    do k = 1,n
      kr  = real(k,O) 
      sum = sum + log(kr)
    end do
    lnfact = sum
  else
    lnfact = 0._O
  end if 
end function lnfactorial
!***********************************************************************************
function gumnkl (u, m, n, k, l) result (g)
  use parameters
  implicit none
  integer :: u, m, n, k, l
  real(O) :: g
!  
  integer :: int
  real(O) :: lnfactorial, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, f
!
  int = u + m + k
  a1  = lnfactorial (int)
  int = u - m - k
  a2  = lnfactorial (int)
  int = n + l - u
  a3  = lnfactorial (int)
  int = u + n - l
  a4  = lnfactorial (int)
  int = u + l - n
  a5  = lnfactorial (int)
  int = n - m
  a6  = lnfactorial (int)
  int = n + m
  a7  = lnfactorial (int)
  int = l - k
  a8  = lnfactorial (int)
  int = l + k
  a9  = lnfactorial (int)
  int = n + l + u + 1
  a10 = lnfactorial (int)
  g   = 0.5_O * (a1 + a2 + a3 + a4 + a5 - a6 - a7 - a8 - a9 - a10)
  f   = sqrt(real(2 * u + 1,O))
  g   = f * exp(g)
end function gumnkl
!***********************************************************************************
subroutine CouplingCoef (m, n, k, l, Cmnkl)
!-----------------------------------------------------------------------------------
! The routine computes the coupling coefficients Cmnkl.                            !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer :: m, n, k, l
  real(O) :: Cmnkl(0:n+l)
!
  integer :: u, Linf
  real(O) :: bu, cu, Su, Sup1, Sum1, gumnkl
!
  Linf = max(abs(m + k), abs(n - l))
  if (Linf > 0) then
    do u = 0, Linf - 1
      Cmnkl(u) = 0._O
    end do
  end if    
  Cmnkl(n+l) = gumnkl (n + l, m, n, k, l)   
  if (n + l > Linf) then  
    u  = n + l
    bu = real((m - k) * u * (u + 1) - (m + k) * (n * (n + 1) - l * (l + 1)),O)
    bu = bu / real(n + l - u + 1,O) / real(n + l + u + 1,O)
    bu = bu * real(2 * u + 1,O) / real(u + 1,O)
    Cmnkl(n+l-1) = bu * gumnkl (n + l - 1, m, n, k, l) 
    if (n + l - 1 > Linf) then
      Sup1 = 1._O
      Su   = bu
      do u = n + l - 1, Linf + 1, -1    
        bu   = real((m - k) * u * (u + 1) - (m + k) * (n * (n + 1) - l * (l + 1)),O)
        bu   = bu / real(n + l - u + 1,O) / real(n + l + u + 1,O)
        bu   = bu * real(2 * u + 1,O) / real(u + 1,O)  
        cu   = real((u + m + k + 1) * (u - m - k + 1),O)
        cu   = cu * real(u,O) / real(u + 1,O)
        cu   = cu * real(u + n - l + 1,O) / real(n + l - u + 1,O)
        cu   = - cu * real(u + l - n + 1,O) / real(n + l + u + 1,O)
        Sum1 = bu * Su + cu * Sup1
        Cmnkl(u-1) = gumnkl(u - 1, m, n, k, l) * Sum1
        Sup1 = Su
        Su   = Sum1
      end do
    end if
  end if  
end subroutine CouplingCoef
