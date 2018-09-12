! **********************************************************************************
! *                  ROUTINES FOR COMPUTING MACHINE CONSTANTS                      *
! *    -----------------------------------------------------------------------     *
! *    List of subroutines:                                                        *
! *      DrvParameters,       machr.                                               *  
! **********************************************************************************
subroutine DrvParameters
!-----------------------------------------------------------------------------------
! The routine computes the following machine constants:                            !
! - NBaseDig (integer),        NIterBes (integer),         NIterPol (integer),     !
! - LargestPosNumber (real),   SmallestPosNumber (real),   MachEps (real),         !               
! - ZeroCoord (real),          TolJ0Val (real),            TolRootPol (real),      !
! - InitBesVal (real),         FactNBes (real),            LargestBesVal (real),   !
! - UpperBoundSeq (real),      LowerBoundSeq (real),       ZeroSinXX (real),       ! 
! - MaxArgBes (real),          ZeroLUVal (real),           LargestSplineVal (real) !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none  
  integer  :: ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp
  real(O)  :: eps, epsneg, xmin, xmax, x, step, fct
!  
  call machr ( ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp,          &
       eps, epsneg, xmin, xmax )    
  LargestPosNumber  = xmax
  SmallestPosNumber = xmin
  MachEps  = eps
  NBaseDig = it
  NIterBes = 5000 
  NIterPol = 100000    
!  
  ZeroCoord  = eps**0.666_O
  if (O == kind(1.e0)) then
    fct = 10._O
  else
    fct = 1._O
  end if
  TolJ0Val   = fct * eps**0.666_O
  TolRootPol = eps
  InitBesVal = 1.e-35_O 
  FactNBes   = 400._O
  LargestBesVal = xmax**0.333_O 
  MaxArgBes     = 1.e+4_O  
  UpperBoundSeq = 1.e+10_O
  LowerBoundSeq = 1.e-10_O  
  x    = 1._O
  step = 0.5_O
  fct  = sin(x) / x
  do while (fct < 1._O)
    x   = step * x
    fct = sin(x) / x
  end do
  ZeroSinXX = x 
  ZeroLUVal = 1.e-20_O
  x = sqrt (xmax)
  LargestSplineVal = sqrt(x)                  
end subroutine DrvParameters
! **********************************************************************************
subroutine machr ( ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp,      &
           eps, epsneg, xmin, xmax)
!----------------------------------------------------------------------------------- 
! The routine computes the parameters of the floating-point arithmetic system.     !
!                                                                                  !
! Output parameters:                                                               !
! - ibeta (integer) - the radix for the floating-point representation.             !
! - it (integer) - the number of base ibeta digits in the floating-point           !
!   significand.                                                                   !
! - irnd (integer) - specifies if floating-point addition chops or rounds          !
!   appear.                                                                        ! 
! - ngrd (integer) - the number of guard digits for multiplication with            !
!   truncating arithmetic.                                                         !
! - machep (integer) - the largest negative integer such that                      !
!   1.0 + real(ibeta) ** machep /= 1.0.                                            !
! - negeps (integer) - the largest negative integer such that                      ! 
!   1.0 - real(ibeta) ** negeps /= 1.0.                                            !
! - iexp (integer) - the number of bits reserved for the representation of         !
!   the exponent of a floating-point number.                                       !
! - minexp (integer) - the largest in magnitude negative integer such that         !
!   real(ibeta) ** minexp is positive and normalized.                              !
! - maxexp (integer) - the smallest positive power of beta that overflows.         !
! - eps (real) - the smallest positive floating-point number such that             !
!   1.0 + eps /= 0.0.                                                              !
! - epsneg (real) - the smallest positive floating-point number such that          !
!   1.0 - epsneg /= 0.0.                                                           !
! - xmin (real) - the smallest non-vanishing normalized floating-point power.      !
!   of the radix: xmin = real(ibeta) ** minexp.                                    !
! - xmax (real) - the largest finite floating-point number.                        !
!   xmax = (1.0 - epsneg) * real(ibeta) ** maxexp.                                 ! 
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer  :: ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp
  real(O)  :: eps, epsneg, xmin, xmax
!  
  integer  :: i, itemp, iz, j, k, mx, nxres
  real(O)  :: a, b, beta, betah, betain, oner, t, temp, temp1, tempa, twor,          &
              y, z, zeror 
!
  oner  = 1._O
  twor  = oner + oner
  zeror = oner - oner
  a = oner
  do 
    a = a + a
    temp  = a + oner
    temp1 = temp - a
    if (temp1 - oner /= zeror ) exit    
  end do
  b = oner
  do 
    b = b + b
    temp  = a + b
    itemp = int (temp - a)
    if (itemp /= 0) exit    
  end do  
  ibeta = itemp
  beta  = real(ibeta,O)
  it = 0
  b  = oner
  do
    it = it + 1
    b  = b * beta
    temp  = b + oner
    temp1 = temp - b
    if (temp1 - oner /= zeror) exit
  end do                 
  irnd  = 0
  betah = beta / twor
  temp  = a + betah
  if (temp - a /= zeror) irnd = 1
  tempa = a + beta
  temp  = tempa + betah
  if (irnd == 0 .and. temp - tempa /= zeror) irnd = 2
  negep  = it +3
  betain = oner / beta
  a = oner
  do i = 1, negep
    a = a * betain
  end do
  b = a
  do 
    temp = oner - a
    if (temp - oner /= zeror) exit
    a = a * beta
    negep = negep - 1
  end do  
  negep  = - negep
  epsneg = a 
  if (ibeta /= 2 .and. irnd /= 0) then
    a = (a * (oner + a)) / twor
    temp = oner - a 
    if (temp - oner /= zeror) epsneg = a
  end if
  machep = -it - 3
  a = b
  do
    temp = oner + a
    if (temp - oner /= zeror) exit
    a = a * beta
    machep = machep + 1
  end do
  eps  = a
  temp = tempa + beta * (oner + eps)
  if (ibeta /= 2 .and. irnd /= 0) then
    a = (a * (oner + a)) / twor
    temp = oner + a
    if (temp - oner /= zeror) eps = a
  end if
  ngrd = 0
  temp = oner + eps
  if (irnd == 0 .and. temp * oner - oner /= zeror) ngrd = 1
  i = 0
  k = 1
  z = betain
  t = oner + eps
  nxres = 0
  do 
    y = z
    z = y * y
    a = z * oner
    temp = z * t
    if (a + a == zeror .or. abs(z) >= y) exit
    temp1 = temp * betain
    if (temp1 * beta == z) exit
    i = i + 1
    k = k + k
  end do
  if (ibeta /= 10) then
    iexp = i + 1
    mx   = k + k
  else
    iexp = 2
    iz   = ibeta
    do 
      if (k < iz) exit
      iz   = iz * ibeta
      iexp = iexp + 1        
    end do
    mx = iz + iz - 1
  end if
  do
    xmin = y
    y = y * betain    
    a = y * oner
    temp = y * t
    if ( a + a == zeror .or. abs(y) >= xmin) exit
    k = k + 1
    temp1 = temp * betain
    if (temp1 * beta == y) then
      nxres = 3
      xmin  = y
      exit
    end if
  end do
  minexp = - k  
  if (mx <= k + k - 3 .and. ibeta /= 10) then
    mx   = mx + mx
    iexp = iexp + 1    
  end if
  maxexp = mx + minexp
  irnd   = irnd + nxres
  if (irnd == 2 .or. irnd == 5) maxexp = maxexp - 2
  if (irnd == 3 .or. irnd == 4) maxexp = maxexp - it
  i = maxexp + minexp
  if (ibeta ==2 .and. i == 0) maxexp = maxexp - 1
  if (i > 20) maxexp = maxexp - 1
  if (a /= y) maxexp = maxexp - 2
  xmax = oner - epsneg
  if (xmax * oner /= xmax) xmax = oner - beta * epsneg
  xmax = xmax / ( beta * beta * beta * xmin)
  i = maxexp + minexp + 3
  do j = 1, i
    if (ibeta == 2) then
      xmax = xmax + xmax
    else
      xmax = xmax * beta
    end if
  end do
end subroutine machr
