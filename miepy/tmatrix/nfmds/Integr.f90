! **********************************************************************************
! *                      ROUTINES FOR NUMERICAL INTEGRATION                        *
! *      ----------------------------------------------------------------------    *
! *      Partial list of subroutines:                                              *
! *        Gauss_Legendre,       Laguerre,         Simpson,        Trapez,         *
! *        Gauss_Legendre1,      Laguerre1,        Lgroot,         Lgrecr,         *
! *        gaussq,               class,            gausq2,         readinputIntegr *
! **********************************************************************************
subroutine Gauss_Legendre (x1, x2, Np, wp, xp)
!------------------------------------------------------------------------------------ 
! The routine computes the weights and nodes for the Gauss-Legendre quadrature      !
! method.                                                                           !
!                                                                                   !
! Input parameters:                                                                 !
! - x1, x2 (real variables) - lower and upper integration limits.                   !
! - Np (integer) - number of quadrature points.                                     !
!                                                                                   !
! Output parameters:                                                                !
! - wp (real array) - weights.                                                      !
! - xp (real array) - quadrature points (nodes).                                    !
!                                                                                   !
! The following parameters (specified in the group statement "Integration"          !
! from the input file "../INPUTFILES/Input.dat") control the integration            !
! process:                                                                          !
!                                                                                   !
! - TypeIntegr (character array) - specifies the source of the numerical            !
!   integration routines. The permissive values are:                                !
!   - 'MET1' - routines from Numerical Recipes,                                     !
!   - 'MET2' - modified routines from Slatec library.                               !
!   The recommended value of the parameter TypeIntegr is 'MET1'.                    !
!                                                                                   !
! - epsGauss (real) - tolerance for computing the roots of the Legendre polynomials.!
!   The accuracy of the T-matrix calculation strongly depends on epsGauss. The      !
!   value of this parameter should be set AS SMALL AS POSSIBLE. For double          !
!   precision arithmetic, epsGauss = 1.e-15,...,1.e-12.                             !
!------------------------------------------------------------------------------------
  use parameters
  implicit none
  integer       :: Np, i
  real(O)       :: x1, x2, wp(Np), xp(Np), epsGauss, epsLaguerre, xs, xd
  character(20) :: TypeIntegr
  real(O),allocatable :: work(:)
!   
  call readinputIntegr ( TypeIntegr, epsGauss, epsLaguerre ) 
  call check_Integration (TypeIntegr) 
  if (TypeIntegr(1:4) == 'MET1') then
    call Gauss_Legendre1 (x1, x2, Np, wp, xp, epsGauss)
  else if (TypeIntegr(1:4) == 'MET2') then   
    allocate (work(Np))
    call gaussq (1, Np, work, xp, wp)   
    xd = (x2 - x1) / 2._O
    xs = (x2 + x1) / 2._O
    do i = 1, Np
      wp(i) = xd * wp(i)
      xp(i) = xd * xp(i) + xs
    end do
    deallocate (work)
  end if
end subroutine Gauss_Legendre
! **********************************************************************************
subroutine Laguerre (n, x, a)
!------------------------------------------------------------------------------------ 
! The routine computes the weights and nodes for the Laguerre quadrature method.    !
!                                                                                   !
! Input parameters:                                                                 !
! - n (integer) - number of quadrature points.                                      !
!                                                                                   !
! Output parameters:                                                                !
! - x (real array) - quadrature points (nodes).                                     !
! - a (real array) - weights.                                                       !
!                                                                                   !
! The following parameters (specified in the group statement "Integration"          !
! from the input file "../INPUTFILES/Input.dat") control the integration            !
! process:                                                                          !
!                                                                                   !
! - TypeIntegr (character array) - specifies the source of the numerical            !
!   integration routines. The permissive values are:                                !
!   - 'MET1' - routines from Numerical Recipes,                                     !
!   - 'MET2' - modified routines from Slatec library.                               !
!   The recommended value of the parameter TypeIntegr is 'MET1'.                    !
!                                                                                   !
! - epsLaguerre (real) - tolerance for computing the roots of the Laguerre          !
!   polynoms. The default value is 1.e-10.                                          !
!------------------------------------------------------------------------------------
  use parameters
  implicit none      
  integer :: n
  real(O) :: x(n), a(n), epsGauss, epsLaguerre
  character(20) :: TypeIntegr
  real(O),allocatable :: work(:)    
!         
  call readinputIntegr ( TypeIntegr, epsGauss, epsLaguerre )  
  call check_Integration (TypeIntegr)  
  close (unit = iInput)
  if (TypeIntegr(1:4) == 'MET1') then
    call  Laguerre1 (n, x, a, epsLaguerre)
  else if (TypeIntegr(1:4) == 'MET2') then   
    allocate (work(n))
    call gaussq (2, n, work, x, a)      
    deallocate (work)
  end if
end subroutine Laguerre
! **********************************************************************************
subroutine readinputIntegr ( TypeIntegr, epsGauss, epsLaguerre )
  use parameters
  implicit none
  integer       :: ios
  real(O)       :: epsGauss, epsLaguerre
  character(20) :: TypeIntegr
  character(80) :: string
  logical       :: XFindPar
!
  open (unit = iInput, file = FileInput, status = "old", position = "rewind")
  TypeIntegr  = 'MET1'  
  epsGauss    = 1.e-10_O  
  epsLaguerre = 1.e-10_O
  string      = 'Integration'
  if (XFindPar (iInput, string)) then
    read (iInput, *, iostat = ios) TypeIntegr
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable TypeIntegr;')"
      stop
    end if
    read (iInput, *, iostat = ios) epsGauss
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsGauss;')"
      stop
    end if
    read (iInput, *, iostat = ios) epsLaguerre
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsLaguerre;')"
      stop
    end if
  else
    print "(/,2x,'Group name Integration not found;')"
    stop  
  end if  
  call check_Integration (TypeIntegr)  
  close (unit = iInput) 
end subroutine readinputIntegr  
!***********************************************************************************
subroutine Simpson (a, b, N, x, w)
!------------------------------------------------------------------------------------ 
! The routine computes the weights and nodes for the Simpson quadrature method.     !
!                                                                                   !
! Input parameters:                                                                 !
! - a, b (real variables) - lower and upper integration limits.                     !
! - N (integer) - number of quadrature points.                                      !
!                                                                                   !
! Output parameters:                                                                !
! - w (real array) - weights.                                                       !
! - x (real array) - quadrature points (nodes).                                     !
!                                                                                   !
! Note: The number of integration points must be an odd number.                     !
!------------------------------------------------------------------------------------
  use parameters
  implicit none
  integer  :: N, i
  real(O)  :: a, b, h, x(N), w(N), N1r, i1r
!         
  N1r = real(N - 1,O)
  h   = (b - a) / N1r
  do i = 1, N
    i1r  = real(i - 1,O)
    x(i) = a + i1r * h
  end do  
  if (mod(N - 1,2) == 0) then
    w(1) = h / 3._O
    w(N) = h / 3._O
    do i = 2, N - 1, 2
      w(i) = 4._O * h / 3._O
    end do
    do i = 3, N - 2, 2
      w(i) = 2._O * h / 3._O
    end do
  else
    print "(/,2x,'Warning in subroutine Simpson in module Integr:')"
    print "(  2x,'the number of integration points is an even number and therefore')"
    print "(  2x,'the last interval integration is performed with the trapez rule;')"
    w(1)   = h / 3._O
    w(N-1) = 5._O * h / 6._O
    w(N)   = h / 2._O
    do i = 2, N - 2, 2
      w(i) = 4._O * h / 3._O
    end do
    do i = 3, N - 3, 2
      w(i) = 2._O * h / 3._O
    end do       
  end if
end subroutine Simpson
!***********************************************************************************
subroutine Trapez (a, b, N, x, w)
!------------------------------------------------------------------------------------ 
! The routine computes the weights and nodes for the trapez integration method.     !
!                                                                                   !
! Input parameters:                                                                 !
! - a, b (real variables) - lower and upper integration limits.                     !
! - N (integer) - number of quadrature points.                                      !
!                                                                                   !
! Output parameters:                                                                !
! - w (real array) - weights.                                                       !
! - x (real array) - quadrature points (nodes).                                     !
!------------------------------------------------------------------------------------
  use parameters
  implicit none
  integer  :: N, i
  real(O)  :: a, b, h, x(N), w(N), N1r, i1r
!         
  N1r = real(N - 1,O)
  h   = (b - a) / N1r
  do i = 1, N
    i1r  = real(i - 1,O)
    x(i) = a + i1r * h
  end do  
  w(1) = 0.5_O * h
  w(N) = 0.5_O * h
  do i = 2, N - 1
    w(i) = h
  end do 
end subroutine Trapez
! **********************************************************************************
! *                  ROUTINES FOR COMPUTING THE NODES AND WEIGHTS                  *
! *              FOR GAUSSIAN-TYPE QUADRATURE RULES - NUMERICAL RECIPES            *
! **********************************************************************************
subroutine Gauss_Legendre1 (x1, x2, Np, wp, xp, epsGauss)
!-----------------------------------------------------------------------------------
! The significance of the input and ouput parameters is as in the routine          !
! Gauss_Legendre. The following parmeters control the computation:                 !
! - TolRootPol (real) - if epsGauss < TolRootPol, then epsGauss is setted to       !
!   TolRootPol. Default value: TolRootPol = MachEps.                               !
! - NIterPol (integer) - maximum number of iteration for computing the roots       !
!   with the desired accuracy. Default value: NIterPol = 100000.                   !
! These parameters are specified in the routine "DrvParameters" from the file      !
! '../TMATROUTINES/MachParam.f90'.                                                 !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer :: Np, i, j, M, iter
  real(O) :: x1, x2, wp(Np), xp(Np), xm, xl, z, p1, p2, p3, pp, z1,                 &
             jr, j1r, j2r, Npr, ir, epsGauss
  logical :: more
!  
  if (epsGauss < TolRootPol) then
    epsGauss = TolRootPol 
    print "(/,2x,'Warning in subroutine Gauss_Legendre in module Integr:')"
    print "(  2x,'the tolerance epsGauss is too low and epsGauss has been setted')"    
    print "(  2x,'to TolRootPol, where TolRootPol =',1pe13.4,';')", TolRootPol    
  end if
  Npr = real(Np,O)
  M   = int(0.5_O * (Npr + 1._O))
  xm  = 0.5_O * (x1 + x2)
  xl  = 0.5_O * (x2 - x1)
  do i = 1, M
    ir   = real(i,O)
    z    = cos(Pi * (ir - 0.25_O) / (Npr + 0.5_O))
    more =.true.
    iter = 0
    do while (more)
      iter = iter + 1     
      p1   = 1._O
      p2   = 0._O
      do j = 1,Np
        p3  = p2
        p2  = p1
        jr  = real(j,O)
        j1r = real(j - 1,O)
        j2r = real(2 * j - 1,O)
        p1  = (j2r * z * p2 - j1r * p3) / jr
      end do
      pp = Npr * (z * p1 - p2) / (z * z - 1._O)
      z1 = z
      z  = z1 - p1 / pp      
      if (abs(z - z1) <= epsGauss .or. iter == NIterPol) more = .false.
    end do 
    if (abs(z - z1) > epsGauss ) then        
      print "(/,2x,'Error in subroutine Gauss_Legendre in module Integr:')"
      print "(  2x,'the  root was not determined with the prescribed accuracy because ')"
      print "(  2x,'the tolerance epsGauss or the iteration number NIterPol specified ')"
      print "(  2x,'in the subroutine MachParam are too low;')"
      stop 
    end if  
    xp(i)      = xm - xl * z
    xp(Np+1-i) = xm + xl * z
    wp(i)      = 2._O * xl / ((1._O - z * z) * pp * pp)
    wp(Np+1-i) = wp(i)
  end do      
end subroutine Gauss_Legendre1
! **********************************************************************************
subroutine Laguerre1 (n, x, a, epsLaguerre)
!-----------------------------------------------------------------------------------
! The significance of the input and ouput parameters is as in the routine          !
! Laquerre. The parameter NIterPol  specifies the maximum number of iteration      !
! for computing the roots with the desired accuracy. NIterPol is specified in      !
! the routine "DrvParameters" from the file '../TMATROUTINES/MachParam.f90' and    !
! the default value is NIterPol = 100000.                                          !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer :: n
  real(O) :: x(n), a(n), epsLaguerre
!      
  integer :: i, exp, j
  real(O) :: csx, csa, r1, dpn, pn1, xt, nr, i2r, tmp
!    
  nr  = real(n,O) 
  csx = 0._O
  csa = 0._O
  do i = 1, n
    if (i == 1) xt = 3._O / (1._O + 2.4_O * nr)        
    if (i == 2) xt = xt + 15._O / (1._O + 2.5_O * nr)        
    if (i > 2) then
      i2r = real(i - 2,O)
      r1  = (1._O + 2.55_O * i2r) / (1.9_O * i2r)
      xt  = xt + r1 * (xt - x(i-2))
    end if
    call Lgroot (xt, n, dpn, pn1, epsLaguerre, exp)
    x(i) = xt
    tmp  = 1._O / dpn / pn1 / nr
    if (exp > 0) then
      do j = 1, 2 * exp
        tmp = tmp * 0.1_O
      end do
    end if
    a(i) = tmp 
    csx  = csx + xt
    csa  = csa + a(i)
  end do      
end subroutine Laguerre1
! **********************************************************************************
subroutine Lgroot (x, n, dpn, pn1, epsLaguerre, exp)
  use parameters
  use derived_parameters
  implicit none  
  integer  :: n, exp
  real(O)  :: x, dpn, pn1, epsLaguerre
!      
  integer  :: iter
  real(O)  :: p, d, dp
  logical  :: more
!
  iter  = 0
  more  = .true.
  do while (more)
    iter = iter + 1
    call Lgrecr (p, dp, pn1, x, n, exp)
    d = p / dp
    x = x - d  
    if (abs(d / x) <= epsLaguerre .or. iter == NIterPol) more = .false.                         
  end do
  if (abs(d / x) > epsLaguerre ) then
    print "(/,2x,'Error in subroutine Laguerre in module Integr:')"
    print "(  2x,'the root was not determined with the prescribed accuracy because the')"
    print "(  2x,'tolerance  epsLaguerre  or the iteration number  NIterPol  specified')"
    print "(  2x,'in the subroutine MachParam are too low;')"
    stop
  end if
  dpn = dp                                 
end subroutine Lgroot
! ******** *************************************************************************
subroutine Lgrecr (pn, dpn, pn1, x, n, exp)
  use parameters
  implicit none
  integer  :: n, exp
  real(O)  :: pn, dpn, pn1, x
!      
  integer  :: j, dexp
  real(O)  :: p1, p, dp1, dp, q, dq, jr, j1r, j2r, big, low
!
  p1   = 1._O
  p    = x - 1._O
  dp1  = 0._O
  dp   = 1._O
  exp  = 0
  big  = 1.e+3_O
  low  = 1.e-3_O 
  dexp = int(log10(big))
  do j = 2, n
    jr  = real(j,O)
    j1r = real(j - 1,O)
    j2r = real(2*j - 1,O)
    q   = (x - j2r) *  p / jr - j1r * p1 / jr
    dq  = (x - j2r) * dp / jr + p / jr - j1r * dp1 / jr
    if (q > big .or. dq > big) then
      q   =  q * low
      dq  = dq * low
      p   =  p * low
      dp  = dp * low
      exp = exp + dexp
    end if
    p1  = p
    p   = q
    dp1 = dp
    dp  = dq
  end do
  pn  = p
  dpn = dp
  pn1 = p1        
end subroutine Lgrecr
! **********************************************************************************
! *         ROUTINES FOR COMPUTING THE NODES AND WEIGHTS FOR GAUSSIAN-TYPE         *
! *           QUADRATURE RULES - SIMPLIFIED VERSION FROM SLATEC LIBRARY            *
! **********************************************************************************
subroutine gaussq (kindt, n, b, t, w)                  
!-----------------------------------------------------------------------------------
! Input parameters:                                                                !
! - kindt (integer) - specifies the type of quadrature rule, i.e.,                 !
!   kindt = 1 for Legendre quadrature with w(x) = 1 on (-1, 1), and                !
!   kindt = 2 for the generalized Laguerre quadrature with w(x) = exp(-x) on       !
!   (0, +infinity).                                                                !
! - n (integer) - number of points used for the quadrature rule.                   ! 
! - b (real array) - scratch array of length n.                                    ! 
!                                                                                  !
! Output parameters:                                                               !
! - t (real array) - quadrature nodes.                                             !
! - w (real array) - weights.                                                      !
!-----------------------------------------------------------------------------------
  use parameters 
  implicit none
  integer  :: kindt, n
  real(O)  :: b(n), t(n), w(n)
!  
  integer  :: i 
  real(O)  :: muzero      
!                 
  call class (kindt, n, b, t, muzero)   
  w(1) = 1._O
  do i = 2, n
    w(i) = 0._O
  end do
  call gausq2 (kindt, n, t, b, w)
  do i = 1, n
    w(i) = muzero * w(i) * w(i)
  end do                              
end subroutine gaussq
!***********************************************************************************
subroutine class (kindt, n, b, a, muzero)
  use parameters
  implicit none
  integer  :: kindt, n
  real(O)  :: a(n), b(n), muzero
!      
  integer  :: nm1, i     
  real(O)  :: abi
!
  nm1 = n - 1
  if (kindt == 1) then
    muzero = 2._O
    do i = 1, nm1
      a(i) = 0._O
      abi  = i
      b(i) = abi / sqrt(4._O * abi * abi - 1._O)
    end do
    a(n) = 0._O  
  else if (kindt == 2) then      
    muzero = 1._O 
    do i = 1, nm1
      a(i) = 2._O * i - 1._O 
      b(i) = real(i,O)       
    end do
    a(n) = 2._O * n - 1 
  end if             
end subroutine class
!***********************************************************************************
subroutine gausq2 (kindt, n, d, e, z)
  use parameters
  use derived_parameters
  implicit none
  integer :: kindt, n
  real(O) :: d(n), e(n), z(n)
!
  integer :: i, j, k, l, m, ii, mml
  real(O) :: b, c, f, g, p, r, s
!    
  e(n) = 0._O
  do l = 1, n
    j = 0
    do while (j < NIterPol)
      j = j + 1    
      m = l - 1
      do while (m < n)
        m = m + 1           
        if (m < n .and. abs(e(m)) <= MachEps * (abs(d(m)) + abs(d(m+1)))) exit
      end do
      p = d(l)
      if (m == l)  exit                     
      g   = (d(l+1) - p) / (2._O * e(l))
      r   = sqrt(g * g + 1._O)
      g   = d(m) - p + e(l) / (g + sign(r, g))
      s   = 1._O
      c   = 1._O
      p   = 0._O
      mml = m - l                
      do ii = 1, mml
        i = m - ii
        f = s * e(i)
        b = c * e(i)
        if (abs(f) < abs(g)) then
          s = f / g
          r = sqrt(s * s + 1._O)
          e(i+1) = g * r
          c = 1._O / r
          s = s * c         
        else
          c = g / f
          r = sqrt(c * c + 1._O)
          e(i+1) = f * r
          s = 1._O / r
          c = c * s
        end if              
        g = d(i+1) - p
        r = (d(i) - g) * s + 2._O * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
        f = z(i+1)
        z(i+1) = s * z(i) + c * f
        z(i)   = c * z(i) - s * f
      end do                                                                     
      d(l) = d(l) - p
      e(l) = g
      e(m) = 0._O
    end do
    if (j == NIterPol) then
      if (kindt == 1) then
        print "(/,2x,'Error in subroutine gausq2 in module Integr:')"
        print "(  2x,'the Legendre quadratures were not determined with the prescribed ')"
        print "(  2x,'accuracy, because the iteration number NIterPol specified in the')"
        print "(  2x,'subroutine MachParam is too low;')"
      else if (kindt == 2) then
        print "(/,2x,'Error in subroutine gausq2 in module Integr:')"
        print "(  2x,'the Laguerre quadratures were not determined with the prescribed ')"
        print "(  2x,'accuracy, because the iteration number NIterPol specified in the')"
        print "(  2x,'subroutine MachParam is too low;')"
      end if              
      stop
    end if  
  end do  
  do ii = 2, n
    i = ii - 1
    k = i
    p = d(i)
    do j = ii, n
      if (d(j) < p) then
        k = j
        p = d(j)
      end if  
    end do
    if (k /= i) then
      d(k) = d(i)
      d(i) = p
      p    = z(i)
      z(i) = z(k)
      z(k) = p
    end if  
  end do     
end subroutine gausq2
