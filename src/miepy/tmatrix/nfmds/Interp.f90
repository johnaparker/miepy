! **********************************************************************************
! *                           INTERPOLATION ROUTINES                               *
! *      --------------------------------------------------------------------      *
! *      Partial list of subroutines:                                              *
! *        Interp,  linterp,  spline,  splint,  dpchim,  dpchfe,  dchfev,  dpchst, *
! *        readinputInterp                                                         *
! ********************************************************************************** 
subroutine Interp (n, x, y, xpoint, ypoint)
!------------------------------------------------------------------------------------ 
! Given the pairs (x_1,y_1), (x_2,y_2),..., (x_n,y_n), the routine computes the     !
! value ypoint at xpoint using linear, spline or cubic Hermite interpolation.       !
!                                                                                   !
! Input parameters:                                                                 !
! - n (integer)    - number of data points.                                         !
! - x (real array) - set of discrete points (abscissas).                            !
! - y (real array) - function values at discrete points.                            ! 
! - xpoint (real)  - evaluation point.                                              !
!                                                                                   !
! Output parameters:                                                                !
! - ypoint (real) - value of the function at xpoint.                                !
!                                                                                   !
! The parameter TypeInterp (specified in the group statement "Interpolation"        !
! from the input file "../INPUTFILES/Input.dat") provides the type of               !
! interpolation method.                                                             ! 
!------------------------------------------------------------------------------------
  use parameters
  implicit none
  integer   :: n
  real(O)   :: x(n), y(n), xpoint, ypoint
!
  integer       :: i
  real(O)       :: yd1, ydn
  character(20) :: TypeInterp
  logical       :: ascending
  real(O),allocatable :: x1(:), y1(:), yd(:)             
!   
  call readinputInterp ( TypeInterp )    
  allocate (x1(n), y1(n))            
  ascending = .true.
  do i = 1, n - 1
    if (x(i) > x(i+1)) then
      ascending = .false.
      exit
    end if
  end do
  do i = 1, n
    if (ascending) then
      x1(i) = x(i)
      y1(i) = y(i)
    else
      x1(i) = x(n-i+1)
      y1(i) = y(n-i+1)
    end if
  end do   
  if (xpoint < x1(1) .and. xpoint > x1(n)) then
    print "(/,2x,'Error in subroutine Interp in file Interp.f90:')"
    print "(  2x, a)",                                                               &
   'the interpolation point does not belong to the set of discrete points;'    
    stop    
  end if                                             
  if (TypeInterp(1:6) == 'LINEAR') then
    call linterp (x1, y1, n, xpoint, ypoint)  
  else if (TypeInterp(1:6) == 'SPLINE') then    
    allocate (yd(n))
    yd1 = 0._O
    ydn = 0._O     
    call spline (x1, y1, n, yd1, ydn, yd)  
    call splint (x1, y1, yd, n, xpoint, ypoint)
    deallocate (yd)
  else if (TypeInterp(1:7) == 'HERMITE') then   
    allocate ( yd(n))                  
    call dpchim (n, x1, y1, yd)        
    call dpchfe (n, x1, y1, yd, xpoint, ypoint)
    deallocate  (yd)
  end if
  deallocate  (x1, y1)  
end subroutine Interp
! ************************************************************************************
subroutine readinputInterp ( TypeInterp ) 
  use parameters
  implicit none
  integer       :: ios
  character(20) :: TypeInterp
  character(80) :: string
  logical       :: XFindPar
!
  open (unit = iInput, file = FileInput, status = "old", position = "rewind")  
  TypeInterp = 'LINEAR'  
  string     = 'Interpolation'
  if (XFindPar (iInput, string)) then
    read (iInput, *, iostat = ios) TypeInterp
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable TypeInterp;')"
      stop
    end if
  else
    print "(/,2x,'Group name Interpolation not found;')"
    stop  
  end if  
  call check_Interpolation (TypeInterp) 
  close (unit = iInput)   
end subroutine readinputInterp   
! ************************************************************************************
! *                           LINEAR INTERPOLATION ROUTINE                           *
! ************************************************************************************
subroutine linterp (xa, ya, n, x, y)
  use parameters
  use derived_parameters
  implicit none
  integer  :: n
  real(O)  :: xa(n), ya(n), x, y
!  
  integer  :: i,ix
  real(O)  :: xs, xd, csi, dx
  logical  :: found  
! 
  found = .false.   
  do i = 1, n - 1
    xs = xa(i)
    xd = xa(i+1)
    dx = xd - xs
    if (abs(dx) < MachEps) then
      print "(/,2x,'Error in subroutine linterp in file Interp.f90:')"
      print "(  2x,'coincident knots;')"
      stop      
    end if      
    if (x >= xs .and. x < xd) then
      ix    = i
      csi   = (x - xs) / dx
      found = .true.
      exit
    end if                   
  end do
  if (found) then
    y = (1._O - csi) * ya(ix) + csi * ya(ix+1)
  else
    y = ya(n)
  end if    
end subroutine linterp  
! ************************************************************************************
! *         SPLINE INTERPOLATION ROUTINE - F90-VERSION FROM NUMERICAL RECIPES        *
! ************************************************************************************
subroutine spline (x, y, n, yp1, ypn, y2)
!-----------------------------------------------------------------------------------
! The routine computes the derivatives needed to determine a spline interpolant    !
! to given data.                                                                   !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none  
  integer  :: n
  real(O)  :: x(n), y(n), yp1, ypn, y2(n)
!    
  integer  :: i, k
  real(O)  :: sig, p, qn, un, dy, dx, dyp, dym, dxp, dxm
  real(O),allocatable :: u(:)
!  
  allocate (u(n))
  if (yp1 > LargestSplineVal) then
    y2(1) = 0._O
    u (1) = 0._O
  else
    y2(1) = - 0.5_O
    dy    = y(2) - y(1)
    dx    = x(2) - x(1)
    if (abs(dx) < MachEps) then
      print "(/,2x,'Error in subroutine spline in file Interp.f90:')"
      print "(  2x,'coincident knots;')"
      stop      
    end if 
    u (1) = (3._O / dx ) * (dy / dx - yp1)    
  end if
  do i = 2, n - 1
    sig   = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
    p     = sig * y2(i-1) + 2._O
    y2(i) = (sig - 1._O) / p
    dyp   = y(i+1) - y(i)
    dym   = y(i) - y(i-1)
    dxp   = x(i+1) - x(i)
    dxm   = x(i) - x(i-1)
    dx    = x(i+1) - x(i-1)
    if (abs(dxp) < MachEps ) then
      print "(/,2x,'Error in subroutine spline in file Interp.f90:')"
      print "(  2x,'coincident knots;')"
      stop      
    end if 
    u(i)  = (6._O * (dyp / dxp - dym / dxm) / dx - sig * u(i-1)) / p    
  end do
  if (ypn > LargestSplineVal) then
    qn = 0._O
    un = 0._O
  else
    qn = 0.5_O
    dy = y(n) - y(n-1)
    dx = x(n) - x(n-1)
    un = (3._O / dx) * (ypn - dy / dx)    
  end if
  y2(n) = (un - qn * u(n-1)) / (qn * y2(n-1) + 1._O)
  do k = n - 1, 1, - 1
    y2(k) = y2(k) * y2(k+1) + u(k)
  end do
  deallocate (u)
end subroutine spline
! *********************************************************************************** 
subroutine splint (xa, ya, y2a, n, x, y)   
!------------------------------------------------------------------------------------
! The routine evaluates a spline interpolant at a point.                            !
!------------------------------------------------------------------------------------    
  use parameters
  implicit none
  integer  :: n
  real(O)  :: xa(n), ya(n), y2a(n), x, y
!  
  integer  :: klo, khi, k
  real(O)  :: h, a, b 
  logical  :: more
!  
  klo = 1
  khi = n
  more = .true.
  do while (more)
    k = (khi + klo) / 2
    if (xa(k) > x) then
      khi = k
    else
      klo = k
    end if
    if (khi - klo <= 1) more = .false.
  end do
  h = xa(khi) - xa(klo)
  a = (xa(khi) - x) / h
  b = (x - xa(klo)) / h
  y = a * ya(klo) + b * ya(khi) + ((a**3 - a) * y2a(klo) +                          &
     (b**3 - b) * y2a(khi)) * (h**2) / 6._O  
end subroutine splint    
! ************************************************************************************
! *           HERMITE INTERPOLATION ROUTINE - F90-VERSION FROM SLATEC LIBRARY        *
! ************************************************************************************            
subroutine dpchim (n, x, f, d)
!-----------------------------------------------------------------------------------
! The routine computes the derivatives needed to determine a monotone piecewise    !
! cubic Hermite interpolant to given data.                                         !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer  :: n
  real(O)  :: x(n), f(n), d(n)
!
  integer  :: i, nless1
  real(O)  :: del1, del2, dmax, dmin, drat1, drat2, h1, h2, hsum, hsumt3,           &
              w1, w2, dpchst                
!              
  if (n < 2) then
    print "(/,2x,'Error in subroutine dpchfe in file Interp.f90:')"
    print "(  2x,'number of data points less than two;')"
    stop      
  end if             
  nless1 = n - 1
  h1     = x(2) - x(1)
  del1   = (f(2) - f(1)) / h1  
  if (nless1 > 1) then
    h2    =   x(3) - x(2)
    del2  =   (f(3) - f(2)) / h2
    hsum  =   h1 + h2
    w1    =   (h1 + hsum) / hsum
    w2    = - h1 / hsum
    d(1)  =   w1 * del1 + w2 * del2
    if (dpchst (d(1), del1) <= 0._O)  then
      d(1) = 0._O
    else if (dpchst (del1, del2) < 0._O)  then
      dmax = 3._O * del1
      if (abs(d(1)) > abs(dmax))  d(1) = dmax
    end if
    do i = 2, nless1          
      if (i /= 2) then
        h1   = h2
        h2   = x(i+1) - x(i)
        hsum = h1 + h2
        del1 = del2
        del2 = (f(i+1) - f(i)) / h2
      end if          
      d(i) = 0._O
      if (dpchst (del1, del2) > 0._O) then
        hsumt3 = hsum + hsum + hsum
        w1     = (hsum + h1) / hsumt3
        w2     = (hsum + h2) / hsumt3
        dmax   = max(abs(del1), abs(del2) )
        dmin   = min(abs(del1), abs(del2) )
        drat1  = del1 / dmax
        drat2  = del2 / dmax
        d(i)   = dmin / (w1 * drat1 + w2 * drat2)
      end if      
    end do   
    w1   = - h2 / hsum
    w2   =  (h2 + hsum) / hsum
    d(n) = w1 * del1 + w2 * del2
    if (dpchst (d(n), del2) <= 0._O)  then
      d(n) = 0._O
    else if (dpchst (del1, del2) < 0._O)  then
      dmax = 3._O * del2
      if (abs(d(n)) > abs(dmax)) d(n) = dmax
    end if
  else
    d(1) = del1
    d(n) = del1
  end if       
end subroutine dpchim
!************************************************************************************
subroutine dpchfe (n, x, f, d, xe, fe)
!------------------------------------------------------------------------------------
! The routine evaluates a cubic polynomial given in Hermite form at a point.        !
!------------------------------------------------------------------------------------ 
  use parameters  
  implicit none
  integer  :: n
  real(O)  :: x(n), f(n), d(n), xe, fe
!
  integer  :: i, ix
  real(O)  :: xs, xd, x1, x2, f1, f2, d1, d2
  logical  :: found
!  
  found = .false.
  do i = 1, n - 1
    xs = x(i)
    xd = x(i+1)            
    if (xe >= xs .and. xe < xd) then
      ix    = i
      found = .true.      
      exit             
    end if                   
  end do
  if (found) then  
    x1 = x(ix)
    x2 = x(ix+1)
    f1 = f(ix)
    f2 = f(ix+1)
    d1 = d(ix)
    d2 = d(ix+1)          
    call dchfev (x1, x2, f1, f2, d1, d2, xe, fe)
  else
    fe = f(n)
  end if      
end subroutine dpchfe
! ***********************************************************************************
subroutine dchfev (x1, x2, f1, f2, d1, d2, xe, fe) 
  use parameters
  use derived_parameters  
  implicit none  
  real(O)  :: x1, x2, f1, f2, d1, d2, xe, fe
!
  real(O)  :: c2, c3, del1, del2, delta, h, x     
!    
  h =   x2 - x1
  if (abs(h) < MachEps) then
    print "(/,2x,'Error in subroutine dchfev in file Interp.f90:')"
    print "(  2x,'coincident knots;')"
    stop      
  end if          
  delta   =   (f2 - f1) / h
  del1    =   (d1 - delta) / h
  del2    =   (d2 - delta) / h                                          
  c2      = - (del1 + del1 + del2)
  c3      =   (del1 + del2) / h
  x       =   xe - x1
  fe      =   f1 + x * (d1 + x * (c2 + x * c3))        
end subroutine dchfev
! ***********************************************************************************
function dpchst (arg1, arg2) result (func)
  use parameters
  implicit none
  real(O),intent(in) :: arg1, arg2
  real(O)            :: func
!  
  func = 0._O
  if ((arg1 /= 0._O) .and. (arg2 /= 0._O)) then
    func = sign(1._O,arg1) * sign(1._O,arg2)
  end if 
end function dpchst
