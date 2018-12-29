! **********************************************************************************
! *       ROUTINES FOR COMPUTING THE SPHERICAL BESSEL AND HANKEL FUNCTIONS,        *
! *    THE CYLINDRICAL BESSEL FUNCTIONS AND THE ASSOCIATED LEGENDRE FUNCTIONS      *
! *    -------------------------------------------------------------------------   *
! *    Partial list of subroutines:                                                *
! *      besel_j,                 besel_h,              besel_h_complete,          *
! *      J01,                     bes_J,                Pmm,                       *
! *      pimm,                    Leg,                  P_mm_real,                 *
! *      pi_mm_real,              tau_mm_real,          Leg_normalized,            *
! *      P_mm,                    pi_mm,                tau_mm,                    *
! *      Leg_normalized_complex                                                    * 
! *    Partial list of alternative subroutines:                                    * 
! *      besel_jA,                besel_hA,             besel_h_completeA,         *
! *      besel_jB,                besel_jC,             bes_JA,                    *
! *      MSTA1,                   MSTA2,                ENVJ                       *
! **********************************************************************************
subroutine besel_j (z, N, j, jd)
!-----------------------------------------------------------------------------------
! The routine computes the spherical Bessel functions j(z) and the derivatives     !
! jd(z) = d(z*j(z))/dz. The order of the Bessel functions varies between 0 and N,  !
! and the backward recurrence relation                                             !
!                     j(n-1,z) + j(n+1,z) = [(2n+1) / z] * j(n,z).                 !
! is used for computation. Here j(n,z) denotes the Bessel function of order n and  !
! argument z. The following parameters control the computation:                    !                                
! - ZeroSinXX (real) - for |z| < ZeroSinXX, we set j(0,z) = sin(z)/z = 1 and       !
!   j(n,z) = 0, for n = 1,2,...,N.                                                 !
! - FactNBes (real) - determines the order of the Bessel functions Nstart,         !
!   at which the backward recurrence is started, i.e., Nstart = FactNBes * N.      !
!   Default value: FactNBes = 400._O.                                              !
! - InitBesVal (real) - initial value of the Bessel function at the beginning of   !
!   the backward recurrence. Default value: InitBesVal = 1.e-35_O.                 !
! - MaxArgBes (real) - for |z| > MaxArgBes, the asymptotic expressions of          !
!   spherical Bessel functions for large arguments are used.                       !
! - UpperBoundSeq (real) - upper bound of the sequence of Bessel functions         !
!   (if a Bessel function value exceeds this bound the entire sequence is scaled). !
!   Default value: UpperBoundSeq = 1.e+10_O.                                       !
! - LowerBoundSeq (real) - lower bound of the sequence of Bessel functions.        !
!   Default value: LowerBoundSeq = 1.e-10_O.                                       !
! These parameters are specified in the routine "DrvParameters" from the file      !
! '../TMATROUTINES/MachParam.f90'.                                                 !
!                                                                                  !
! Input parameters:                                                                ! 
! - z (complex) - argument of the spherical Bessel functions.                      !
! - N (integer) - maximum order of the spherical Bessel functions.                 !
!                                                                                  !
! Output parameters:                                                               !
! - j  (complex array) - spherical Bessel functions.                               !
! - jd (complex array) - derivatives of the spherical Bessel functions.            !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: N
  complex(O) :: z, j(0:N), jd(0:N)
!
  complex(O) :: j0, j1, f2, f1, f, kc, scale, argc
  real(O)    :: a0, arg
  integer    :: M, Ma, Mn, k, l
!  
  a0 = abs(z)
  do k = 0, N 
    j(k)  = zero
    jd(k) = zero
  end do
  if (a0 < ZeroSinXX) then
    j(0) = one
  else if (a0 >= ZeroSinXX .and. a0 < MaxArgBes) then
    j(0) = sin(z) / z
    j(1) = j(0) / z - cos(z) / z
    if (N >= 2) then
      j0 = j(0)
      j1 = j(1)
      Ma = nint(a0 + 4._O * a0**0.33_O + 2._O + sqrt(101._O + a0)) + 20  
      Mn = N + nint(sqrt(FactNBes * N))
      M  = max(Ma,Mn)  
      f2 = zero
      f1 = cmplx(InitBesVal,0.0,O)
      do k = M, 0, -1
        kc = cmplx(2 * k + 3,0.0,O)
        f  = kc * f1 / z - f2   
        if (k <= N) j(k) = f
        f2 = f1
        f1 = f
        if (abs(f1) > UpperBoundSeq) then
          f2 = f2 * LowerBoundSeq
          f1 = f1 * LowerBoundSeq
          do l = k, N
            j(l) = j(l) * LowerBoundSeq
          end do
        else if (abs(f1) < LowerBoundSeq) then
          f2 = f2 * UpperBoundSeq
          f1 = f1 * UpperBoundSeq
          do l = k, N
            j(l) = j(l) * UpperBoundSeq
          end do
        end if
      end do        
      if (abs(f1) > abs(f2)) then
        scale = j0 / f1
      else
        scale = j1 / f2
      end if                             
      do k = 0, N
        j(k) = scale * j(k)
      end do
    end if 
    do k = 1, N
      kc    = cmplx(k,0.0,O)
      jd(k) = z * j(k-1) - kc * j(k)
    end do
  else
    do k = 0, N
      arg  = 0.5_O * real(k + 1,O) * pi
      argc = cmplx(arg,0.0,O) 
      j(k) = cos (z - argc) / z
    end do 
    do k = 1, N
      kc    = cmplx(k,0.0,O)
      jd(k) = z * j(k-1) - kc * j(k)
    end do
  end if 
end subroutine besel_j
!**********************************************************************************
subroutine besel_h (z, N, h, hd)
!-----------------------------------------------------------------------------------
! The routine computes the spherical Hankel functions of the first kind            !
! h(z) = j(z) + i*y(z) and the derivatives hd(z) = d(z*h(z))/dz, where j(z) and    !
! y(z) are the spherical Bessel and Neumann functions, respectively. The order of  !
! the functions varies between 0 and N, and the control parameters ZeroSinXX,      !
! FactNBes, InitBesVal, UpperBoundSeq and LowerBoundSeq are specified in the       !
! routine "DrvParameters" from the file '../TMATROUTINES/MachParam.f90'.           !
!                                                                                  !
! Input parameters:                                                                ! 
! - z (complex) - argument of the spherical Hankel functions.                      !
! - N (integer) - maximum order of the spherical Hankel functions.                 !
!                                                                                  !
! Output parameters:                                                               !
! - h  (complex array) - spherical Hankel functions.                               !
! - hd (complex array) - derivatives of the spherical Hankel functions.            !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: N
  complex(O) :: z, h(0:N), hd(0:N)
!
  integer    :: Ma, Mn, M, k, l
  real(O)    :: a0, arg
  complex(O) :: j0, j1, f, f1, f2, kc, z2, z3, scale, argc
  complex(O),allocatable :: j(:), y(:)
!                   
  z2 = z * z
  z3 = z2 * z
  a0 = abs(z)
  allocate (j(0:N), y(0:N))
  do k = 0, N 
    j(k)  = zero
    y(k)  = zero
    h(k)  = zero
    hd(k) = zero
  end do
  if (a0 < ZeroSinXX) then
    j(0) = one
  else if (a0 >= ZeroSinXX .and. a0 < MaxArgBes) then
    j(0) = sin(z) / z
    j(1) = j(0) / z - cos(z) / z
    if (N >= 2) then
      j0 = j(0)
      j1 = j(1)
      Ma = nint(a0 + 4._O * a0**0.33_O + 2._O + sqrt(101._O + a0)) + 20  
      Mn = N + nint(sqrt(FactNBes * N))
      M  = max(Ma,Mn)  
      f2 = zero
      f1 = cmplx(InitBesVal,0.0,O)
      do k = M, 0, -1
        kc = cmplx(2 * k + 3,0.0,O)
        f  = kc * f1 / z - f2   
        if (k <= N) j(k) = f
        f2 = f1
        f1 = f
        if (abs(f1) > UpperBoundSeq) then
          f2 = f2 * LowerBoundSeq
          f1 = f1 * LowerBoundSeq
          do l = k, N
            j(l) = j(l) * LowerBoundSeq
          end do
        else if (abs(f1) < LowerBoundSeq) then
          f2 = f2 * UpperBoundSeq
          f1 = f1 * UpperBoundSeq
          do l = k, N
            j(l) = j(l) * UpperBoundSeq
          end do
        end if
      end do        
      if (abs(f1) > abs(f2)) then
        scale = j0 / f1
      else
        scale = j1 / f2
      end if                             
      do k = 0, N
        j(k) = scale * j(k)
      end do      
    end if               
    y(0) = - cos(z) / z
    y(1) = (y(0) - sin(z)) / z
    if (N >= 2) then   
      do k = 2, N
        if (abs(j(k-1)) > abs(j(k-2))) then               
          y(k) = (j(k) * y(k-1) - one / z2) / j(k-1)
        else 
          kc   = cmplx(2 * k - 1,0.0,O)
          y(k) = (j(k) * y(k-2) - kc / z3) / j(k-2)
        end if 
      end do 
    end if    
    do k = 0, N
      h(k) = j(k) + im * y(k)
    end do
    do k = 1, N
      kc    = cmplx(k,0.0,O)
      hd(k) = z * h(k-1) - kc * h(k)
    end do
  else
    do k = 0, N 
      arg  = 0.5_O * real(k + 1,O) * pi
      argc = cmplx(arg,0.0,O)    
      h(k) = exp(im * (z - argc)) / z      
    end do 
    do k = 1, N
      kc    = cmplx(k,0.0,O)
      hd(k) = z * h(k-1) - kc * h(k)
    end do
  end if
  deallocate (j, y)
end subroutine besel_h
!**********************************************************************************
subroutine besel_h_complete (z, N, j, y, h, hd)
!-----------------------------------------------------------------------------------
! The routine computes the spherical Bessel and Neumann functions j(z) and y(z),   !
! the spherical Hankel functions of the first kind h(z) = j(z) + i*y(z) and the    !
! derivatives hd(z) = d(z*h(z))/dz. The order of the functions varies between 0    !
! and N, and the control parameters ZeroSinXX, FactNBes, InitBesVal, UpperBoundSeq !
! and LowerBoundSeq are specified in the routine "DrvParameters" from the file     !
! '../TMATROUTINES/MachParam.f90'.                                                 !
!                                                                                  !
! Input parameters:                                                                ! 
! - z (complex) - argument of the spherical Bessel, Neumann and Hankel functions.  !
! - N (integer) - maximum order of the spherical functions.                        !
!                                                                                  !
! Output parameters:                                                               !
! - j  (complex array) - spherical Bessel functions.                               !
! - y  (complex array) - spherical Neumann functions.                              !
! - h  (complex array) - spherical Hankel functions.                               !
! - hd (complex array) - derivatives of the spherical Hankel functions.            !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: N
  complex(O) :: z, j(0:N), y(0:N), h(0:N), hd(0:N)
!
  integer    :: Ma, Mn, M, k, l
  real(O)    :: a0, arg
  complex(O) :: j0, j1, f, f1, f2, kc, z2, z3, scale, argc
!                   
  z2 = z * z
  z3 = z2 * z
  a0 = abs(z)
  do k = 0, N 
    j(k)  = zero
    y(k)  = zero
    h(k)  = zero
    hd(k) = zero
  end do
  if (a0 < ZeroSinXX) then
    j(0) = one
  else if (a0 >= ZeroSinXX .and. a0 < MaxArgBes) then    
    j(0) = sin(z) / z
    j(1) = j(0) / z - cos(z) / z
    if (N >= 2) then
      j0 = j(0)
      j1 = j(1)
      Ma = nint(a0 + 4._O * a0**0.33_O + 2._O + sqrt(101._O + a0)) + 20  
      Mn = N + nint(sqrt(FactNBes * N))
      M  = max(Ma,Mn)  
      f2 = zero
      f1 = cmplx(InitBesVal,0.0,O)
      do k = M, 0, -1
        kc = cmplx(2 * k + 3,0.0,O)
        f  = kc * f1 / z - f2   
        if (k <= N) j(k) = f
        f2 = f1
        f1 = f
        if (abs(f1) > UpperBoundSeq) then
          f2 = f2 * LowerBoundSeq
          f1 = f1 * LowerBoundSeq
          do l = k, N
            j(l) = j(l) * LowerBoundSeq
          end do
        else if (abs(f1) < LowerBoundSeq) then
          f2 = f2 * UpperBoundSeq
          f1 = f1 * UpperBoundSeq
          do l = k, N
            j(l) = j(l) * UpperBoundSeq
          end do
        end if
      end do        
      if (abs(f1) > abs(f2)) then
        scale = j0 / f1
      else
        scale = j1 / f2
      end if                             
      do k = 0, N
        j(k) = scale * j(k)
      end do      
    end if               
    y(0) = - cos(z) / z
    y(1) = (y(0) - sin(z)) / z
    if (N >= 2) then   
      do k = 2, N
        if (abs(j(k-1)) > abs(j(k-2))) then               
          y(k) = (j(k) * y(k-1) - one / z2) / j(k-1)
        else 
          kc   = cmplx(2 * k - 1,0.0,O)
          y(k) = (j(k) * y(k-2) - kc / z3) / j(k-2)
        end if 
      end do 
    end if    
    do k = 0, N
      h(k) = j(k) + im * y(k)
    end do
    do k = 1, N
      kc    = cmplx(k,0.0,O)
      hd(k) = z * h(k-1) - kc * h(k)
    end do
  else    
    do k = 0, N 
      arg  = 0.5_O * real(k + 1,O) * pi
      argc = cmplx(arg,0.0,O)     
      j(k) = cos (z - argc) / z
      y(k) = sin (z - argc) / z
      h(k) = j(k) + im * y(k) 
    end do   
    do k = 1, N
      kc    = cmplx(k,0.0,O)
      hd(k) = z * h(k-1) - kc * h(k)
    end do
  end if
end subroutine besel_h_complete
!***********************************************************************************
subroutine J01 (z, J0, J1)
  use parameters
  use derived_parameters
  implicit none
  complex(O) :: z, J0, J1
!
  complex(O) :: z2, z1, r, t1, p0, q0, u, t2, p1, q1, kc
  real(O)    :: a0, rp2, A(12), B(12), A1(12), B1(12)
  integer    :: k, k0 
!
  rp2 = 2._O / pi
  a0  = abs(z)
  z2  = z * z
  z1  = z
  if (a0 == 0._O) then
    J0 = one
    J1 = zero          
  else 
    if (real(z) < 0._O) z1 = - z
    if (a0 <= 12._O) then
      J0 = one
      r  = one
      do k = 1, 40
        kc = cmplx(k * k,0.0,O)             
        r  = - 0.25_O * r * z2 / kc
        J0 = J0 + r
        if (abs(r) < MachEps * abs(J0)) exit
      end do 
      J1 = one
      r  = one
      do k = 1, 40
        kc = cmplx(k * (k + 1),0.0,O)
        r  = - 0.25_O * r * z2 / kc
        J1 = J1 + r
        if (abs(r) < MachEps * abs(J1)) exit
      end do 
      J1 = 0.5_O * z1 * J1      
    else 
      A  = (/-.703125e-01_O          , .1121520996093750e+00_O,                     & 
             -.5725014209747314e+00_O, .6074042001273483e+01_O,                     &
             -.1100171402692467e+03_O, .3038090510922384e+04_O,                     &
             -.1188384262567832e+06_O, .6252951493434797e+07_O,                     &
             -.4259392165047669e+09_O, .3646840080706556e+11_O,                     &
             -.3833534661393944e+13_O, .4854014686852901e+15_O/)
      B  = (/ .732421875e-01_O       ,-.2271080017089844e+00_O,                     &
              .1727727502584457e+01_O,-.2438052969955606e+02_O,                     &
              .5513358961220206e+03_O,-.1825775547429318e+05_O,                     &
              .8328593040162893e+06_O,-.5006958953198893e+08_O,                     &
              .3836255180230433e+10_O,-.3649010818849833e+12_O,                     &
              .4218971570284096e+14_O,-.5827244631566907e+16_O/)
      A1 = (/ .1171875e+00_O         ,-.1441955566406250e+00_O,                     &
              .6765925884246826e+00_O,-.6883914268109947e+01_O,                     &
              .1215978918765359e+03_O,-.3302272294480852e+04_O,                     &
              .1276412726461746e+06_O,-.6656367718817688e+07_O,                     &
              .4502786003050393e+09_O,-.3833857520742790e+11_O,                     &
              .4011838599133198e+13_O,-.5060568503314727e+15_O/)
      B1 = (/-.1025390625e+00_O      , .2775764465332031e+00_O,                     &
             -.1993531733751297e+01_O, .2724882731126854e+02_O,                     &
             -.6038440767050702e+03_O, .1971837591223663e+05_O,                     &
             -.8902978767070678e+06_O, .5310411010968522e+08_O,                     &
             -.4043620325107754e+10_O, .3827011346598605e+12_O,                     &
             -.4406481417852278e+14_o, .6065091351222699e+16_O/)
      k0 = 12
      if (a0 >= 35._O) k0 = 10
      if (a0 >= 50._O) k0 = 8
      t1 = z1 - 0.25_O * pi
      p0 = one
      do k = 1, k0
        p0 = p0 + A(k) * z1**(- 2 * k)
      end do
      q0 = - 0.125_O / z1
      do k = 1, k0
        q0 = q0 + B(k) * z1**(- 2 * k - 1)
      end do
      u  = sqrt(rp2 / z1)
      J0 = u * (p0 * cos(t1) - q0 * sin(t1))      
      t2 = z1 - 0.75_O * pi
      p1 = one
      do k = 1, k0
        p1 = p1 + A1(k) * z1**(- 2 * k)
      end do
      q1 = 0.375_O / z1
      do k = 1, k0
        q1 = q1 + B1(k) * z1**(- 2 * k - 1)
      end do
      J1 = u * (p1 * cos(t2) - q1 * sin(t2))      
    end if 
    if (real(z) < 0._O) then
      J1 = - J1
    end if
  end if     
end subroutine J01
!***********************************************************************************
subroutine bes_J (z, N, J)
!-----------------------------------------------------------------------------------
! The routine computes the cylindrical Bessel function J(z). The order of the      !
! Bessel function varies between 0 and N, and the control parameters ZeroSinXX,    !
! FactNBes, InitBesVal, UpperBoundSeq and LowerBoundSeq are specified in the       !
! routine "DrvParameters" from the file '../TMATROUTINES/MachParam.f90'.           !
!                                                                                  !
! Input parameters:                                                                ! 
! - z (complex) - argument of the cylindrical Bessel functions.                    !
! - N (integer) - maximum order of cylindrical Bessel functions.                   !
!                                                                                  !
! Output parameters:                                                               !
! - J (complex array) - cylindrical Bessel functions.                              !  
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: N
  complex(O) :: z, J(0:N)
!
  complex(O) :: J0, J1, S, S0, f2, f1, f, kc
  real(O)    :: y0, a0
  integer    :: k, M, Ma, Mn, l
!   
  y0 = abs(aimag(z))
  a0 = abs(z)
  do k = 0, N
    J(k) = zero
  end do
  if (a0 < ZeroSinXX) then      
    J(0) = one
  else 
    call J01 (z, J0, J1)      
    J(0) = J0      
    J(1) = J1
    if (N > 1) then   
      if (a0 <= 300._O .or. N > int(0.25_O * a0)) then
        Ma = nint(a0 + 4._O * a0**0.33_O + 2._O + sqrt(101._O + a0)) + 10 
        Mn = N + nint(sqrt(FactNBes * N))
        M  = max(Ma,Mn)                                   
        S  = zero     
        f2 = zero
        f1 = cmplx(InitBesVal,0.0,O) 
        do k = M, 0, -1
          kc = cmplx(k + 1,0.0,O)
          f  = 2._O * kc / z * f1 - f2
          if (k <= N) J(k) = f
          if (k == 2 * int(k / 2) .and. k /= 0) then
            if (y0 <= 1._O) then
              S = S + 2._O * f
            else 
              S = S + (- 1._O)**(k / 2) * 2._O * f
            end if 
          end if    
          f2 = f1
          f1 = f
          if (abs(f1) > UpperBoundSeq) then
            f2 = f2 * LowerBoundSeq
            f1 = f1 * LowerBoundSeq
            S  = S * LowerBoundSeq
            do l = k, N
              J(l) = J(l) * LowerBoundSeq
            end do
          else if (abs(f1) < LowerBoundSeq) then
            f2 = f2 * UpperBoundSeq
            f1 = f1 * UpperBoundSeq
            S  = S * UpperBoundSeq
            do l = k, N
              J(l) = J(l) * UpperBoundSeq
            end do
          end if
        end do 
        if (y0 <= 1._O) then
          S0 = S + f1
        else 
          S0 = (S + f1) / cos(z)
        end if 
        do k = 0, N
          J(k) = J(k) / S0
        end do 
      else      
        do k = 2, N
          kc   = cmplx(k - 1,0.0,O)
          J(k) = 2._O * kc / z * J(k-1) - J(k-2)  
        end do
      end if 
    end if
  endif   
end subroutine bes_J
!***********************************************************************************
subroutine Pmm (theta, m, Pm)
!-----------------------------------------------------------------------------------
! The routine computes the Legendre function Pmm(theta) for m > 0.                 !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none  
  integer :: k, m
  real(O) :: theta, Pm, tamp, f, sth
!
  sth  = sin(theta)
  tamp = 1._O
  do k = 1, m
    f    = real(2 * k - 1,O)
    tamp = f * sth * tamp  
  end do
  Pm = tamp 
end subroutine Pmm
!***********************************************************************************
subroutine pimm (theta, m, pim)
!-----------------------------------------------------------------------------------
! The routine computes the auxiliary function pimm,                                !
!                    pimm(theta) = Pmm[cos(theta)] / sin(theta)                    !
! for m > 1.                                                                       !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer :: k, m
  real(O) :: theta, pim, tamp, f, sth  
!
  sth  = sin(theta)
  tamp = 1._O
  do k = 1, m - 1
    f    = real(2 * k - 1,O)
    tamp = f * sth * tamp  
  end do
  f   = real(2 * m - 1,O)
  pim = f * tamp
end subroutine pimm
!***********************************************************************************
subroutine Leg (theta, m, Nmax, Pnm, dPnm, pinm, taunm)
!-----------------------------------------------------------------------------------
! The routine computes the associated Legendre functions Pnm[cos(theta)], the      !
! derivatives dPn0[cos(theta)] = dPn0[cos(theta)]/d[cos(theta)], and the auxiliary !
! functions pinm(theta) = Pnm[cos(theta)]/sin(theta), and taunm(theta) =           !
! dPnm[cos(theta)]/d(theta) of degree m. The order of the associated Legendre      !
! functions varies between 0 and Nmax, and note that for m > 0, the first m values !
! of Pnm, pinm and taunm (corresponding to the orders 0, 1, ..., m - 1) are zero.  !
! The derivatives dPn0 are computed for the azimuthal mode m = 0 (for m >0, the    !
! Nmax values of dPn0 are zero).                                                   ! 
!                                                                                  !
! Input parameteres:                                                               !
! - theta (real) - polar angle (argument of the associated Legendre functions).    !
! - m (integer) - degree of the associated Legendre functions.                     !
! - Nmax (integer) - maximum order of the associated Legendre functions.           !
!                                                                                  ! 
! Output parameters:                                                               !
! - Pnm (real array) - associated Legendre functions.                              ! 
! - dPnm (real array) - derivatives of the associated Legendre functions           !
!   for m = 0.                                                                     !
! - pinm and taunm (real arrays) - auxiliary functions.                            !  
!----------------------------------------------------------------------------------- 
  use parameters
  implicit none
  real(O) :: theta
  integer :: m, Nmax
  real(O) :: Pnm(0:Nmax), dPnm(0:Nmax), pinm(0:Nmax), taunm(0:Nmax)
!
  real(O) :: Pm, pim, f1, f2, cth, sth
  integer :: n
!
  cth = cos(theta)
  sth = sin(theta)
  if (m == 0) then
    Pnm(0) = 1._O
    Pnm(1) = cth
    do n = 1, Nmax - 1
      f1 = real(2 * n + 1,O) / real(n + 1,O)
      f2 = real(n,O) / real(n + 1,O)
      Pnm(n+1) = f1 * cth * Pnm(n) - f2 * Pnm(n-1)
    end do
  else
    do n = 0, m - 1
      Pnm(n) = 0._O
    end do
    call Pmm (theta, m, Pm)
    Pnm(m) = Pm
    do n = m, Nmax - 1
      f1 = real(2 * n + 1,O) / real(n + 1 - m,O)
      f2 = real(n + m,O) / real(n + 1 - m,O)
      Pnm(n+1) = f1 * cth * Pnm(n) - f2 * Pnm(n-1)
    end do
  end if
  if (m == 0) then
    dPnm(0) = 0._O
    do n = 1, Nmax
      f1 = real(n,O)
      dPnm(n) = f1 * Pnm(n-1) + cth * dPnm(n-1)
    end do
  else
    do n = 0, Nmax
      dPnm(n) = 0._O
    end do
  end if
  if (m == 0) then
    do n = 0, Nmax
      pinm(n)  =   0._O
      taunm(n) = - sth * dPnm(n)
    end do
  else
    if (m == 1) then
      pinm(0)  = 0._O     
      taunm(0) = 0._O
      pinm(1)  = 1._O
    else
      do n = 0, m - 1
        pinm(n)  = 0._O
        taunm(n) = 0._O
      end do
      call pimm (theta, m, pim)
      pinm(m) = pim
    end if
    do n = m, Nmax - 1
      f1 = real(2 * n + 1,O) / real(n + 1 - m,O)
      f2 = real(n + m,O) / real(n + 1 - m,O)
      pinm(n+1) = f1 * cth * pinm(n) - f2 * pinm(n-1)
    end do
    do n = m, Nmax
      f1 = real(n,O)
      f2 = real(n + m,O)
      taunm(n) = f1 * cth * pinm(n) - f2 * pinm(n-1)
    end do
  end if   
end subroutine Leg
!***********************************************************************************
subroutine P_mm_real (theta, m, Pm)
!-----------------------------------------------------------------------------------
! The routine computes the normalized Legendre function Pmm(theta) of real         !
! argument (real angle).                                                           !  
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer :: k, m
  real(O) :: theta, Pm, tamp, f, cth, sth
!
  cth = cos(theta)
  sth = sin(theta) 
  if (m == 0) then
    Pm = sqrt(3._O / 2._O) * cth
  else
    tamp = 1._O
    do k = 1, m
      f    = sqrt(real(m + k,O) / 4._O / real(k,O))
      tamp = f * sth * tamp        
    end do
    f  = sqrt(real(2 * m + 1,O) / 2._O)
    Pm = f * tamp
  end if
end subroutine P_mm_real
!***********************************************************************************
subroutine pi_mm_real (theta, m, pim)
!-----------------------------------------------------------------------------------
! The routine computes the normalized auxiliary function pimm,                     !
!                    pimm(theta) = Pmm[cos(theta)] / sin(theta)                    !
! of real argument (real angle).                                                   !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer :: k, m
  real(O) :: theta, pim, tamp, f, sth
!
  sth = sin(theta)
  if (m == 0) then
    pim = 0._O
  else
    tamp = 1._O
    do k = 1, m - 1
      f    = sqrt(real(m + k,O) / 4._O / real(k,O))
      tamp = f * sth * tamp        
    end do   
    f   = sqrt(real(2 * m + 1,O)) / 2._O
    pim = f * tamp
  end if
end subroutine pi_mm_real
!***********************************************************************************
subroutine tau_mm_real (theta, m, taum)
!-----------------------------------------------------------------------------------
! The routine computes the normalized auxiliary function taumm,                    !
!                  taumm(theta) = dPmm[cos(theta)] / d(theta)                      !
! of real argument (real angle).                                                   !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer :: k, m
  real(O) :: theta, taum, tamp, f, cth, sth  
!
  cth = cos(theta)
  sth = sin(theta)
  if (m == 0) then
    taum = - sqrt(3._O / 2._O) * sth
  else
    tamp = 1._O
    do k = 1, m - 1
      f    = sqrt(real(m + k,O) / 4._O / real(k,O)) 
      tamp = f * sth * tamp       
    end do
    f    = real(m,O) * sqrt(real(2 * m + 1,O)) / 2._O
    taum = f * cth * tamp
  end if
end subroutine tau_mm_real
!***********************************************************************************
subroutine Leg_normalized (theta, m, Nmax, Pnm, dPnm, pinm, taunm)
!------------------------------------------------------------------------------------
! The routine computes the normalized associated Legendre functions                 !
! Pnm[cos(theta)], the derivatives dPn0[cos(theta)] = dPn0[cos(theta)]/d[cos(theta)]!
! and the auxiliary functions pinm(theta) = Pnm[cos(theta)]/sin(theta), and         !
! taunm(theta) = dPnm[cos(theta)]/d(theta) of degree m. The order of the normalized !
! associated Legendre functions varies between 0 and Nmax, and note that for        !
! m > 0, the first m values of Pnm, pinm and taunm (corresponding to the orders     !
! 0, 1, ..., m - 1)are zero. The derivatives dPn0 are computed for the azimuthal    !
! mode m = 0 (for m >0, the Nmax values of dPn0 are zero).                          ! 
!                                                                                   !
! Input parameteres:                                                                !
! - theta (real) - polar angle (argument of the normalized associated Legendre      !
!   functions).                                                                     !
! - m (integer) - degree of the normalized associated Legendre functions.           !
! - Nmax (integer) - maximum order of the normalized associated Legendre            !
!   functions.                                                                      !
!                                                                                   ! 
! Output parameters:                                                                !
! - Pnm (real array) - normalized associated Legendre functions.                    ! 
! - dPnm (real array) - derivatives of the normalized associated Legendre           !
!   functions for m = 0.                                                            !
! - pinm and taunm (real arrays) - normalized auxiliary functions.                  !  
!------------------------------------------------------------------------------------ 
  use parameters
  implicit none
  real(O) :: theta
  integer :: m, Nmax
  real(O) :: Pnm(0:Nmax), dPnm(0:Nmax), pinm(0:Nmax), taunm(0:Nmax)
!
  real(O) :: Pm, pim, f1, f2, cth, sth
  integer :: n
!
  cth = cos(theta)
  sth = sin(theta)
  if (m == 0) then
    Pnm(0) = sqrt(2._O) / 2._O
    Pnm(1) = sqrt(3._O / 2._O) * cos(theta)
    do n = 1, Nmax - 1
      f1 = sqrt(real(2 * n + 1,O) / real(n + 1,O))
      f1 = f1 * sqrt(real(2 * n + 3,O) / real(n + 1,O))
      f2 = sqrt(real(2 * n + 3,O) / real(2 * n - 1,O))
      f2 = f2 * real(n,O) / real(n + 1,O)
      Pnm(n+1) = f1 * cth * Pnm(n) - f2 * Pnm(n-1)
    end do
  else
    do n = 0, m - 1
      Pnm(n) = 0._O
    end do
    call P_mm_real (theta, m, Pm)
    Pnm(m) = Pm
    do n = m, Nmax - 1
      f1 = sqrt(real(2 * n + 1,O) / real(n + 1 - m,O))
      f1 = f1 * sqrt(real(2 * n + 3,O) / real(n + 1 + m,O))
      f2 = sqrt(real(2 * n + 3,O) / real(2 * n - 1,O))
      f2 = f2 * sqrt(real(n - m,O) / real(n - m + 1,O))
      f2 = f2 * sqrt(real(n + m,O) / real(n + m + 1,O))
      Pnm(n+1) = f1 * cth * Pnm(n) - f2 * Pnm(n-1)
    end do
  end if
  if (m == 0) then
    dPnm(0) = 0._O
    do n = 1, Nmax
      f2 = sqrt(real(2 * n + 1,O) / real(2 * n - 1,O))
      f1 = real(n,O) * f2
      dPnm(n) = f1 * Pnm(n-1) + f2 * cth * dPnm(n-1)
    end do
  else
    do n = 0, Nmax
      dPnm(n) = 0._O
    end do
  end if
  if (m == 0) then
    do n = 0, Nmax
      pinm(n)  =   0._O
      taunm(n) = - sth * dPnm(n)
    end do
  else
    if (m == 1) then
      pinm(0)  = 0._O     
      taunm(0) = 0._O
      pinm(1)  = sqrt(3._O) / 2._O
    else
      do n = 0, m - 1
        pinm(n)  = 0._O
        taunm(n) = 0._O
      end do
      call pi_mm_real (theta, m, pim)
      pinm(m) = pim
    end if
    do n = m, Nmax - 1
      f1 = sqrt(real(2 * n + 1,O) / real(n + 1 - m,O))
      f1 = f1 * sqrt(real(2 * n + 3,O) / real(n + 1 + m,O))
      f2 = sqrt(real(2 * n + 3,O) / real(2 * n - 1,O))
      f2 = f2 * sqrt(real(n - m,O) / real(n - m + 1,O))
      f2 = f2 * sqrt(real(n + m,O) / real(n + m + 1,O))
      pinm(n+1) = f1 * cth * pinm(n) - f2 * pinm(n-1)
    end do
    do n = m, Nmax
      f1 = real(n,O)
      f2 = sqrt(real(2 * n + 1,O) / real(2 * n - 1,O))
      f2 = f2 * sqrt(real(n - m,O) / real(n + m,O))
      f2 = f2 * real(n + m,O)
      taunm(n) = f1 * cth * pinm(n) - f2 * pinm(n-1)
    end do
  end if   
end subroutine Leg_normalized
!***********************************************************************************
subroutine P_mm (sint, cost, m, Pm)
!-----------------------------------------------------------------------------------
! The routine computes the normalized Legendre function Pmm(theta) of complex      !
! argument (complex angle).                                                        !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer     :: k, m
  real(O)     :: f
  complex(O)  :: sint, cost, Pm, tamp  
!
  if (m == 0) then
    Pm = sqrt(3._O/2._O) * cost
  else
    tamp = one
    do k = 1, m
      f    = sqrt(real(m + k,O) / 4._O / real(k,O)) 
      tamp = f * sint * tamp  
    end do
    f  = sqrt(real(2 * m + 1,O) / 2._O)
    Pm = f * tamp
  end if 
end subroutine P_mm
!***********************************************************************************
subroutine pi_mm (sint, m, pim)
!-----------------------------------------------------------------------------------
! The routine computes the normalized auxiliary function                           !
!                     pimm(theta) = Pmm[cos(theta)] / sin(theta)                   !
! of complex argument (complex angle).                                             !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer     :: k, m
  real(O)     :: f
  complex(O)  :: sint, pim, tamp  
!
  if (m == 0) then
    pim = zero
  else
    tamp = one
    do k = 1, m - 1
      f    = sqrt(real(m + k,O) / 4._O / real(k,O))
      tamp = sint * f * tamp  
    end do    
    f   = sqrt(real(2 * m + 1,O)) / 2._O
    pim = f * tamp
  end if 
end subroutine pi_mm
!***********************************************************************************
subroutine tau_mm (sint, cost, m, taum)
!-----------------------------------------------------------------------------------
! The routine computes the normalized auxiliary function                           !
!                      taumm(theta) = dPmm[cos(theta)] / d(theta)                  !
! of complex argument (complex angle).                                             !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer     :: k, m
  real(O)     :: f
  complex(O)  :: sint, cost, taum, tamp  
!
  if (m == 0) then
    taum = - sqrt(3._O / 2._O) * sint
  else
    tamp = one
    do k = 1, m - 1
      f    = sqrt(real(m + k,O) / 4._O / real(k,O))
      tamp = f * sint * tamp   
    end do    
    f    = real(m,O) * sqrt(real(2 * m + 1,O)) / 2._O
    taum = f * cost * tamp
  end if 
end subroutine tau_mm
!***********************************************************************************
subroutine Leg_normalized_complex (sint, cost, m, Nmax, Pnm, dPnm, pinm, taunm)
!------------------------------------------------------------------------------------
! The routine computes the normalized associated Legendre functions                 !
! Pnm[cos(theta)], the derivatives dPn0[cos(theta)] = dPn0[cos(theta)]/d[cos(theta)]!
! and the auxiliary functions pinm(theta) = Pnm[cos(theta)]/sin(theta), and         !
! taunm(theta) = dPnm[cos(theta)]/d(theta) of degree m and complex argument.        !
! The order of the normalized  associated Legendre functions varies between 0 and   !
! Nmax, and note that for m > 0, the first m values of Pnm, pinm and taunm          !
! (corresponding to the orders 0, 1, ..., m - 1) are zero. The derivatives dPn0     !
! are computed for the azimuthal mode m = 0 (for m >0, the Nmax values of dPn0      !
! are zero).                                                                        ! 
!                                                                                   !
! Input parameteres:                                                                !
! - sint (complex) - sint = sin(theta), where theta is the complex polar angle.     !
! - cost (complex) - cost = cos(theta), where theta is the complex polar angle.     !
! - m (integer) - degree of the normalized associated Legendre functions.           !
! - Nmax (integer) - maximum order of the normalized associated Legendre            !
!   functions.                                                                      !
!                                                                                   ! 
! Output parameters:                                                                !
! - Pnm (complex array) - normalized associated Legendre functions.                 ! 
! - dPnm (complex array) - derivatives of the normalized associated Legendre        !
!   functions for m = 0.                                                            !
! - pinm and taunm (complex arrays) - normalized auxiliary functions.               !  
!------------------------------------------------------------------------------------ 
  use parameters
  implicit none
  integer     :: m, Nmax
  complex(O)  :: sint, cost
  complex(O)  :: Pnm(0:Nmax), dPnm(0:Nmax), pinm(0:Nmax), taunm(0:Nmax)
!
  complex(O)  :: Pm, pim
  real(O)     :: f1, f2
  integer     :: n
!
  if (m == 0) then
    Pnm(0) = cmplx(sqrt(2._O) / 2._O,0.0,O)
    Pnm(1) = sqrt(3._O / 2._O) * cost
    do n = 1, Nmax - 1
      f1 = sqrt(real(2 * n + 1,O) / real(n + 1,O))
      f1 = f1 * sqrt(real(2 * n + 3,O) / real(n + 1,O))
      f2 = sqrt(real(2 * n + 3,O) / real(2 * n - 1,O))
      f2 = f2 * real(n,O) / real(n + 1,O)
      Pnm(n+1) = f1 * cost * Pnm(n) - f2 * Pnm(n-1)
    end do
  else
    do n = 0, m - 1
      Pnm(n) = zero
    end do
    call P_mm (sint, cost, m, Pm)
    Pnm(m) = Pm
    do n = m, Nmax - 1
      f1 = sqrt(real(2 * n + 1,O) / real(n + 1 - m,O))
      f1 = f1 * sqrt(real(2 * n + 3,O) / real(n + 1 + m,O))
      f2 = sqrt(real(2 * n + 3,O) / real(2 * n - 1,O))
      f2 = f2 * sqrt(real(n - m,O) / real(n - m + 1,O))
      f2 = f2 * sqrt(real(n + m,O) / real(n + m + 1,O))
      Pnm(n+1) = f1 * cost * Pnm(n) - f2 * Pnm(n-1)
    end do
  end if
  if (m == 0) then
    dPnm(0) = zero
    do n = 1, Nmax
      f2 = sqrt(real(2 * n + 1,O) / real(2 * n - 1,O))
      f1 = real(n,O) * f2
      dPnm(n) = f1 * Pnm(n-1) + f2 * cost * dPnm(n-1)
    end do
  else
    do n = 0, Nmax
      dPnm(n) = zero
    end do
  end if
  if (m == 0) then
    do n = 0, Nmax
      pinm(n)  =   zero
      taunm(n) = - sint * dPnm(n)
    end do
  else
    if (m == 1) then
      pinm(0)  = zero
      taunm(0) = zero
      pinm(1)  = cmplx(sqrt(3._O) / 2._O,0.0,O)
    else
      do n = 0, m - 1
        pinm(n)  = zero
        taunm(n) = zero
      end do
      call pi_mm (sint, m, pim)     
      pinm(m) = pim
    end if
    do n = m, Nmax - 1
      f1 = sqrt(real(2 * n + 1,O) / real(n + 1 - m,O))
      f1 = f1 * sqrt(real(2 * n + 3,O) / real(n + 1 + m,O))
      f2 = sqrt(real(2 * n + 3,O) / real(2 * n - 1,O))
      f2 = f2 * sqrt(real(n - m,O) / real(n - m + 1,O))
      f2 = f2 * sqrt(real(n + m,O) / real(n + m + 1,O))
      pinm(n+1) = f1 * cost * pinm(n) - f2 * pinm(n-1)
    end do
    do n = m, Nmax
      f1 = real(n,O)
      f2 = sqrt(real(2 * n + 1,O) / real(2 * n - 1,O))
      f2 = f2 * sqrt(real(n - m,O) / real(n + m,O))
      f2 = f2 * real(n + m,O)
      taunm(n) = f1 * cost * pinm(n) - f2 * pinm(n-1)
    end do    
  end if   
end subroutine Leg_normalized_complex
! **********************************************************************************
! *       ALTERNATIVE ROUTINES FOR COMPUTING THE SPHERICAL BESSEL AND HANKEL       *
! *               FUNCTIONS AND THE CYLINDRICAL BESSEL FUNCTIONS                   *
! **********************************************************************************
subroutine besel_jA (z, N, j, jd)
!-----------------------------------------------------------------------------------
! The routine computes the spherical Bessel functions j(z) and the derivatives     !
! jd(z) = d(z*j(z))/dz. The order of the Bessel functions varies between 0 and N,  !
! and the backward recurrence relation                                             !
!          csi(n,z) = 1.0 / ((2n+1) / z - csi(n+1,z)), n = Nstart,...,1,           !
! and the forward recurrence relations                                             !
!          j(n,z)   = csi(n,z) * j(n-1,z)                                          !
!          jd(n,z)  = ( z - n * csi(n,z)) * j(n-1,z),    n = 1,...,N,              !
! with j(0,z) = sin(z) / z, are used for computation. Here Nstart is the downward  !
! starting index and j(n,z) denotes the Bessel function of order n and argument z. !
! Note that:                                                                       !
!          csi(n,z) = j(n,z) / j(n-1,z).                                           ! 
! The following parmeters control the computation:                                 !                                
! - ZeroSinXX (real) - for |z| < ZeroSinXX, we set j(0,z) = sin(z)/z = 1 and       !
!   j(n,z) = 0, for n = 1,2,...,N.                                                 !
! - NIterBes (integer) - maximum number of backward recurrences. Default value     !
!   NIterBes = 5000.                                                               !
! - TolJ0Val (real) - if abs(j(0,z) - sin(z)/z) > TolJ0Val, the starting index     !
!   Nstart (of the backward recurrence) is increased by 20 and the recurrence is   !
!   restarted. Default value: TolJ0Val = MachEps**0.666_O                          !
! These parameters are specified in the routine "DrvParameters" from the file      !
! '../TMATROUTINES/MachParam.f90'.                                                 !
!                                                                                  !
! Input parameters:                                                                ! 
! - z (complex) - argument of the spherical Bessel functions.                      !
! - N (integer) - maximum order of the spherical Bessel functions.                 !
!                                                                                  !
! Output parameters:                                                               !
! - j  (complex array) - spherical Bessel functions.                               !
! - jd (complex array) - derivatives of the spherical Bessel functions.            !
!                                                                                  !
! Note: This routine does not use the asymptotic expressions of spherical Bessel   !
! functions for large arguments.                                                   !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: N
  complex(O) :: z, j(0:N), jd(0:N)
!  
  integer    :: M, k
  real(O)    :: a0, deps
  complex(O) :: csin1, csin, fct, j0, kc
  logical    :: more
  complex(O),allocatable :: csi(:)
!                   
  a0 = abs(z)   
  allocate (csi(N))
  do k = 0, N 
    j(k)  = zero
    jd(k) = zero
  end do
  if (a0 < ZeroSinXX) then
    j(0) = one
  else
    j(0)  = sin(z) / z
    j(1)  = j(0) / z - cos(z) / z
    jd(1) = z * j(0) - j(1)
    if (N >= 2) then
      M = N + nint(a0 + 4._O * a0**0.33_O + 2._O + sqrt(101._O + a0)) + 10 
      more = .true.       
      do while (more .and. M < NIterBes)
        csin1 = zero
        do k = M, 1, -1        
          fct  = cmplx(2 * k + 1,0.0,O) / z
          csin = one / (fct - csin1)
          if (k <= N) csi(k) = csin
          csin1 = csin
        end do
        j(0) = cos(z) / (one - z * csi(1))
        j0   = sin(z) / z 
        deps = abs (j0 - j(0))                          
        if (deps < TolJ0Val) then
          more = .false.
        else
          M = M + 20
        end if
      end do
      if (more) then
        print "(/,2x,'Error in subroutine besel_j:')"
        print "(2x,'the  csi  sequence can not be computed with the desired accuracy.')"
        print "(2x,'The tolerance TolJ0Val or the iteration number NIterBes specified')"
        print "(2x,'in the subroutine MachParam are too low.                         ')"
        stop 
      end if              
      j(1)  = csi(1) * j(0)
      jd(1) = (z - csi(1)) * j(0)
      do k = 2, N
        kc    = cmplx(k,0.0,O)
        j(k)  = csi(k) * j(k-1)
        jd(k) = (z - kc * csi(k)) * j(k-1)
      end do      
    end if     
  end if 
  deallocate (csi)
end subroutine besel_jA
!**********************************************************************************
subroutine besel_hA (z, N, h, hd)
!-----------------------------------------------------------------------------------
! The routine computes the spherical Hankel functions of the first kind            !
! h(z) = j(z) + i*y(z) and the derivatives hd(z) = d(z*h(z))/dz, where j(z) and    !
! y(z) are the spherical Bessel and Neumann functions, respectively. The order of  !
! the functions varies between 0 and N, and the control parameters ZeroSinXX,      !
! NIterBes and TolJ0Val are specified in the routine "DrvParameters" from the file !
! '../TMATROUTINES/MachParam.f90'.                                                 !
!                                                                                  !
! Input parameters:                                                                ! 
! - z (complex) - argument of the spherical Hankel functions.                      !
! - N (integer) - maximum order of the spherical Hankel functions.                 !
!                                                                                  !
! Output parameters:                                                               !
! - h  (complex array) - spherical Hankel functions.                               !
! - hd (complex array) - derivatives of the spherical Hankel functions.            !
!                                                                                  !
! Note: This routine does not use the asymptotic expressions of spherical Bessel   !
! functions for large arguments.                                                   !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: N
  complex(O) :: z, h(0:N), hd(0:N)
!
  integer    :: M, k
  real(O)    :: a0, deps
  complex(O) :: csin1, csin, fct, j0, kc, z2, z3
  logical    :: more
  complex(O),allocatable :: csi(:), j(:), y(:)
!                   
  z2 = z * z
  z3 = z2 * z
  a0 = abs(z)
  allocate (csi(N), j(0:N), y(0:N))
  do k = 0, N 
    j(k)  = zero
    y(k)  = zero
    h(k)  = zero
    hd(k) = zero
  end do
  if (a0 < ZeroSinXX) then
    j(0) = one
  else
    j(0) = sin(z) / z
    j(1) = j(0) / z - cos(z) / z
    if (N >= 2) then
      M = N + nint(a0 + 4._O * a0**0.33_O + 2._O + sqrt(101._O + a0)) + 10
      more = .true.       
      do while (more .and. M < NIterBes)
        csin1 = zero
        do k = M, 1, -1        
          fct  = cmplx(2 * k + 1,0.0,O) / z
          csin = one / (fct - csin1)
          if (k <= N) csi(k) = csin
          csin1 = csin
        end do
        j(0) = cos(z) / (one - z * csi(1))
        j0   = sin(z) / z 
        deps = abs (j0 - j(0))                          
        if (deps < TolJ0Val) then
          more = .false.
        else
          M = M + 20
        end if
      end do  
      if (more) then
        print "(/,2x,'Error in subroutine besel_h:')"
        print "(2x,'the  csi  sequence can not be computed with the desired accuracy.')"
        print "(2x,'The tolerance TolJ0Val or the iteration number NIterBes specified')"
        print "(2x,'in the subroutine MachParam are too low.                         ')"                
        stop 
      end if              
      j(1) = csi(1) * j(0)    
      do k = 2, N       
        j(k) = csi(k) * j(k-1)              
      end do      
    end if               
    y(0) = - cos(z) / z
    y(1) = (y(0) - sin(z)) / z
    if (N >= 2) then   
      do k = 2, N
        if (abs(j(k-1)) > abs(j(k-2))) then               
          y(k) = (j(k) * y(k-1) - one / z2) / j(k-1)
        else 
          kc   = cmplx(2 * k - 1,0.0,O)
          y(k) = (j(k) * y(k-2) - kc / z3) / j(k-2)
        end if 
      end do 
    end if    
    do k = 0, N
      h(k) = j(k) + im * y(k)
    end do
    do k = 1, N
      kc    = cmplx(k,0.0,O)
      hd(k) = z * h(k-1) - kc * h(k)
    end do
  end if
  deallocate (csi, j, y)
end subroutine besel_hA
!**********************************************************************************
subroutine besel_h_completeA (z, N, j, y, h, hd)
!-----------------------------------------------------------------------------------
! The routine computes the spherical Bessel and Neumann functions j(z) and y(z),   !
! the spherical Hankel functions of the first kind h(z) = j(z) + i*y(z) and the    !
! derivatives hd(z) = d(z*h(z))/dz. The order of the functions varies between 0    !
! and N, and the control parameters ZeroSinXX, NIterBes and TolJ0Val are specified !
! in the routine "DrvParameters" from the file '../TMATROUTINES/MachParam.f90'.    !
!                                                                                  !
! Input parameters:                                                                ! 
! - z (complex) - argument of the spherical Bessel, Neumann and Hankel functions.  !
! - N (integer) - maximum order of the spherical functions.                        !
!                                                                                  !
! Output parameters:                                                               !
! - j  (complex array) - spherical Bessel functions.                               !
! - y  (complex array) - spherical Neumann functions.                              !
! - h  (complex array) - spherical Hankel functions.                               !
! - hd (complex array) - derivatives of the spherical Hankel functions.            !
!                                                                                  !
! Note: This routine does not use the asymptotic expressions of spherical Bessel   !
! functions for large arguments.                                                   !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: N
  complex(O) :: z, j(0:N), y(0:N), h(0:N), hd(0:N)
!
  integer    :: M, k
  real(O)    :: a0, deps
  complex(O) :: csin1, csin, fct, j0, kc, z2, z3
  logical    :: more
  complex(O),allocatable :: csi(:)
!                   
  z2 = z * z
  z3 = z2 * z
  a0 = abs(z)  
  allocate (csi(N))
  do k = 0, N 
    j(k)  = zero
    y(k)  = zero
    h(k)  = zero
    hd(k) = zero
  end do
  if (a0 < ZeroSinXX) then
    j(0) = one
  else
    j(0) = sin(z) / z
    j(1) = j(0) / z - cos(z) / z
    if (N >= 2) then
      M = N + nint(a0 + 4._O * a0**0.33_O + 2._O + sqrt(101._O + a0)) + 10
      more = .true.       
      do while (more .and. M < NIterBes)
        csin1 = zero
        do k = M, 1, -1        
          fct  = cmplx(2 * k + 1,0.0,O) / z
          csin = one / (fct - csin1)
          if (k <= N) csi(k) = csin
          csin1 = csin
        end do
        j(0) = cos(z) / (one - z * csi(1))
        j0   = sin(z) / z 
        deps = abs (j0 - j(0))                          
        if (deps < TolJ0Val) then
          more = .false.
        else
          M = M + 20
        end if
      end do
      if (more) then
        print "(/,2x,'Error in subroutine besel_h_complete:')"
        print "(2x,'the  csi  sequence can not be computed with the desired accuracy.')"
        print "(2x,'The tolerance TolJ0Val or the iteration number NIterBes specified')"
        print "(2x,'in the subroutine MachParam are too low.                         ')"                
        stop 
      end if              
      j(1) = csi(1) * j(0)    
      do k = 2, N       
        j(k) = csi(k) * j(k-1)              
      end do      
    end if               
    y(0) = - cos(z) / z
    y(1) = (y(0) - sin(z)) / z
    if (N >= 2) then   
      do k = 2, N
        if (abs(j(k-1)) > abs(j(k-2))) then               
          y(k) = (j(k) * y(k-1) - one / z2) / j(k-1)
        else 
          kc   = cmplx(2 * k - 1,0.0,O)
          y(k) = (j(k) * y(k-2) - kc / z3) / j(k-2)
        end if 
      end do 
    end if    
    do k = 0, N
      h(k) = j(k) + im * y(k)
    end do
    do k = 1, N
      kc    = cmplx(k,0.0,O)
      hd(k) = z * h(k-1) - kc * h(k)
    end do
  end if
  deallocate (csi)
end subroutine besel_h_completeA
! **********************************************************************************
subroutine besel_jB (z, N, j, y)
!-----------------------------------------------------------------------------------
! The routine computes the spherical Bessel functions j(z) and the spherical       !
! Neumann functions y(z). The order of the functions varies between 0 and N, and   !
! a backward recurrence relation is used to compute the Bessel functions, while    !
! an upward recurrence relation is employed to compute the spherical Neumann       !
! functions. The control parameters ZeroSinXX, FactNBes, InitBesVal, UpperBoundSeq !
! and LowerBoundSeq are specified in the routine "DrvParameters" from the file     !
! '../TMATROUTINES/MachParam.f90'.                                                 !
!                                                                                  !
! Input parameters:                                                                ! 
! - z (complex) - argument of the spherical Bessel and Neumann functions.          !
! - N (integer) - maximum order of the spherical functions.                        !
!                                                                                  !
! Output parameters:                                                               !
! - j (complex array) - spherical Bessel functions.                                !
! - y (complex array) - spherical Neumann functions.                               !
!                                                                                  !
! Note: This routine does not use the asymptotic expressions of spherical Bessel   !
! functions for large arguments.                                                   !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: N
  complex(O) :: z, j(0:N), y(0:N)
!
  complex(O) :: j0, j1, j2, f2, f1, f, kc, scale, z2, z3
  real(O)    :: a0, ajscale
  integer    :: M, Ma, Mn, k, l
!
  z2 = z * z
  z3 = z2 * z
  a0 = abs(z) 
  do k = 0, N 
    j(k) = zero
    y(k) = zero
  end do
  if (a0 < ZeroSinXX) then
    j(0) = one
  else
    j(0) = sin(z) / z
    j(1) = j(0) / z - cos(z) / z
    if (N >= 2) then
      j0 = j(0)
      j1 = j(1)
      j2 = 3._O * j1 / z - j0
      Ma = nint(a0 + 4._O * a0**0.33_O + 2._O + sqrt(101._O + a0)) + 10 
      Mn = N + nint(sqrt(FactNBes * N))
      M  = max(Ma,Mn)  
      f2 = zero
      f1 = cmplx(InitBesVal,0.0,O)
      do k = M - 1, N, -1
        kc = cmplx(2 * k + 1,0.0,O)
        f  = kc * f1 / z - f2               
        f2 = f1
        f1 = f
        if (abs(f1) > UpperBoundSeq) then
          f2 = f2 * LowerBoundSeq
          f1 = f1 * LowerBoundSeq             
        else if (abs(f1) < LowerBoundSeq) then
          f2 = f2 * UpperBoundSeq
          f1 = f1 * UpperBoundSeq            
        end if        
      end do
      j(N)   = f2
      j(N-1) = f1
      do k = N - 1, 1, -1
        kc = cmplx(2 * k + 1,0.0,O)
        j(k-1) = kc * j(k)/z - j(k+1)
        if (abs(j(k-1)) > UpperBoundSeq) then           
          do l = k - 1, N
            j(l) = j(l) * LowerBoundSeq
          end do
        else if (abs(j(k-1)) < LowerBoundSeq) then          
          do l = k - 1, N
            j(l) = j(l) * UpperBoundSeq
          end do 
        end if
      end do
      ajscale = max(abs(j0), abs(j1), abs(j2))
      if (ajscale == abs(j0)) then
        scale = j0 / j(0)
      else if (ajscale == abs(j1)) then
        scale = j1 / j(1)
      else
        scale = j2 / j(2)           
      end if
      do k = 0, N
        j(k) = scale * j(k)
      end do
    end if 
    y(0) = - cos(z) / z
    y(1) = (y(0) - sin(z)) / z
    if (N >= 2) then   
      do k = 2, N
        if (abs(j(k-1)) > abs(j(k-2))) then
          y(k) = (j(k) * y(k-1) - 1._O / z2) / j(k-1)
        else 
          kc   = cmplx(2 * k - 1,0.0,O)
          y(k) = (j(k) * y(k-2) - kc / z3) / j(k-2)
        end if 
      end do 
    end if
  end if 
end subroutine besel_jB
!***********************************************************************************
subroutine besel_jC (z, N, j, y)
!-----------------------------------------------------------------------------------
! The routine computes the spherical Bessel functions j(z) and the spherical       !
! Neumann functions y(z). The order of the functions varies between 0 and N, and   !
! a backward recurrence relation is used to compute the Bessel functions, while    !
! an upward recurrence relation is employed to compute the spherical Neumann       !
! functions. The parameters which control the computation are ZeroSinXX,           !
! InitBesVal, UpperBoundSeq, LowerBoundSeq and LargestBesVal. LargestBesVal        !
! is real and denotes the largest (accepted) values of the spherical Neumann       !
! functions (default value: LargestBesVal = LargestPosNumber**0.333). These        !
! parameters are specified in the routine "DrvParameters" from the file            !
! '../TMATROUTINES/MachParam.f90'.                                                 !
!                                                                                  !
! Input parameters:                                                                ! 
! - z (complex) - argument of the spherical Bessel functions.                      !
! - N (integer) - maximum order of the spherical Bessel functions.                 !
!                                                                                  !
! Output parameters:                                                               !
! - j (complex array) - spherical Bessel functions.                                !
! - y (complex array) - spherical Neumann functions.                               !
!                                                                                  !
! Note: This routine does not use the asymptotic expressions of spherical Bessel   !
! functions for large arguments.                                                   !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: N
  complex(O) :: z, j(0:N), y(0:N)
!
  complex(O) :: j0, j1, j2, f2, f1, f, y0, y1, g0, g1, yk, h2, h1, h0, yl2, yl1,    &
                ylk, P11, P12, P21, P22, scale, z2, kc
  real(O)    :: a0, yak, ya1, wa, ya0, ajscale
  integer    :: NM, M, k, l, LB, LB0, MSTA1, MSTA2
!
  z2 = z * z
  a0 = abs(z) 
  NM = N         
  do k = 0, N 
    j(k) = zero
    y(k) = zero
  end do
  if (a0 < ZeroSinXX) then
    j(0) = one
  else
    j(0) = sin(z) / z
    j(1) = (j(0) - cos(z)) / z
    if (NM >= 2) then
       j0 = j(0)
       j1 = j(1)
       j2 = 3._O * j1 / z - j0      
       M  = MSTA1 (a0, 200)
       if (M < N) then
         NM = M       
       else 
         M = MSTA2 (A0, N, 15)
       end if                                   
       f2 = zero
       f1 = cmplx(InitBesVal,0.0,O)
       do k = M - 1, NM, -1
         kc = cmplx(2 * k + 1,0.0,O)
         f  = kc * f1 / z - f2   
         f2 = f1
         f1 = f
         if (abs(f1) > UpperBoundSeq) then
           f2 = f2 * LowerBoundSeq
           f1 = f1 * LowerBoundSeq            
         else if (abs(f1) < LowerBoundSeq) then
           f2 = f2 * UpperBoundSeq
           f1 = f1 * UpperBoundSeq           
         end if         
       end do
       j(NM)   = f2
       j(NM-1) = f1
       do k = NM - 1, 1, -1
         kc     = cmplx(2 * k + 1,0.0,O) 
         j(k-1) = kc * j(k) / z - j(k+1)
         if (abs(j(k-1)) > UpperBoundSeq) then      
           do l = k - 1, NM
             j(l) = j(l) * LowerBoundSeq
           end do
         else if (abs(j(k-1)) < LowerBoundSeq) then         
           do l = k - 1, NM
             j(l) = j(l) * UpperBoundSeq
           end do 
         end if 
       end do
       ajscale = max(abs(j0), abs(j1), abs(j2))
       if (ajscale == abs(j0)) then
         scale = j0 / j(0)
       else if (ajscale == abs(j1)) then
         scale = j1 / j(1)
       else
         scale = j2 / j(2)
       end if          
       do k = 0, NM
         j(k) = scale * j(k)
       end do
    end if 
    y(0) = - cos(z) / z
    y(1) = (y(0) - sin(z)) / z    
    if (NM >= 2) then
       y0  = y(0)
       y1  = y(1)   
       ya0 = abs(y0)
       LB  = 0
       g0  = y0
       g1  = y1
       do k = 2, NM
         kc = cmplx(2 * k - 1,0.0,O)
         yk = kc / z * g1 - g0
         if (abs (yk) > LargestBesVal) exit  
         yak = abs(yk)
         ya1 = abs(g0)
         if ((yak < ya0) .and. (yak < ya1)) LB = k
         y(k) = yk
         g0   = g1
         g1   = yk
       end do       
       do while ((LB > 4) .and. (aimag(z) /= 0._O) .and. (LB /= LB0))     
         h2  = one
         h1  = zero
         LB0 = LB
         do k = LB, 1, -1
           kc = cmplx(2 * k + 1,0.0,O)
           h0 = kc / z * h1 - h2
           h2 = h1
           h1 = h0
         end do
         P12 = h0
         P22 = h2
         h2  = zero
         h1  = one
         do k = LB, 1, -1
           kc = cmplx(2 * k + 1,0.0,O)
           h0 = kc / z * h1 - h2
           h2 = h1
           h1 = h0
         end do
         P11 = h0
         P21 = h2
         if (LB == NM) j(LB+1) = (2._O * LB + 1._O) / z * j(LB) - j(LB-1)
         if (abs(j(0)) > abs(j(1))) then                   
           y(LB+1) = (j(LB+1) * y0 - P11 / z2) / j(0)
           y(LB)   = (j(LB) * y0   + P12 / z2) / j(0)
         else 
           y(LB+1) = (j(LB+1) * y1 - P21 / z2) / j(1)
           y(LB)   = (j(LB) * y1   + P22 / z2) / j(1)
         end if 
         yl2 = y(LB+1)
         yl1 = y(LB)
         do k = LB - 1, 0, -1
           kc = cmplx(2 * k + 3,0.0,O)
           ylk = kc / z * yl1 - yl2
           y(k)= ylk 
           yl2 = yl1
           yl1 = ylk 
         end do
         yl1 = y(LB)
         yl2 = y(LB+1)
         do k = LB + 1, NM - 1
           kc     = cmplx(2 * k + 1,0.0,O)
           ylk    = kc /z * yl2 - yl1 
           y(k+1) = ylk 
           yl1    = yl2
           yl2    = ylk
         end do
         do k = 2, NM
           wa = abs(y(k))
           IF (wa < abs(y(k-1))) LB = k
         end do 
       end do     
    end if
  end if 
end subroutine besel_jC
! **********************************************************************************
subroutine bes_JA (z, N, J)
!-----------------------------------------------------------------------------------
! The routine computes the cylindrical Bessel function J(z). The order of the      !
! Bessel functions varies between 0 and N, and the control parameters ZeroSinXX,   !
! InitBesVal, UpperBoundSeq and LowerBoundSeq are specified in the routine         !
! "DrvParameters" from the file '../TMATROUTINES/MachParam.f90'.                   !
!                                                                                  !
! Input parameters:                                                                ! 
! - z (complex) - argument of the cylindrical Bessel functions.                    !
! - N (integer) - maximum order of cylindrical Bessel functions.                   !
!                                                                                  !
! Output parameters:                                                               !
! - J (complex array) - cylindrical Bessel functions.                              !  
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none 
  integer    :: N
  complex(O) :: z, J(0:N)
!
  complex(O) :: J0, J1, f2, f1, f, scale, kc
  real(O)    :: a0
  integer    :: NM, M, k, l, MSTA1, MSTA2 
!
  a0 = abs(z)
  NM = N  
  do k = 0, N
    J(k) = zero
  end do
  if (a0 < ZeroSinXX) then      
    J(0) = one
  else 
    call J01 (z, J0, J1)
    J(0) = J0
    J(1) = J1   
    if (N > 1) then
      if (N < int(0.25_O * a0)) then        
        do k = 2, N
          kc   = cmplx(2 * (k - 1),0.0,O)
          J(k) = kc * J(k-1) / z - J(k-2)    
        end do 
      else 
        M = MSTA1 (a0, 200)
        if (M < N) then
          NM = M              
        else 
          M = MSTA2 (A0, N, 15)
        end if  
        f2 = zero
        f1 = cmplx(InitBesVal,0.0,O) 
        do k = M, 0, -1
          kc = cmplx(2 * (k + 1),0.0,O)           
          f  = kc * f1 / z - f2           
          if (k <= NM) J(k) = f
          f2 = f1
          f1 = f
          if (abs(f1) > UpperBoundSeq) then
            f2 = f2 * LowerBoundSeq
            f1 = f1 * LowerBoundSeq
            do l = k, NM
              J(l) = J(l) * LowerBoundSeq
            end do
          else if (abs(f1) < LowerBoundSeq) then
            f2 = f2 * UpperBoundSeq
            f1 = f1 * UpperBoundSeq
            do l = k, NM
              J(l) = J(l) * UpperBoundSeq
            end do       
          end if 
        end do
        if (abs(f1) > abs(f2)) then
          scale = J0 / f1
        else 
          scale = J1 / f2
        end if 
        do k = 0, NM
          J(k) = scale * J(k)
        end do 
      end if
    end if
  end if    
end subroutine bes_JA
!***********************************************************************************
function MSTA1 (x, MP) result (MSTA)
  use parameters
  implicit none 
  real(O), intent(in) :: x
  integer,intent(in)  :: MP  
  integer             :: MSTA, N0, N1, it, NN 
  real(O)             :: a0, f0, f1, f, ENVJ    
!
  a0 = abs(x)
  if (a0 < 1._O) a0 = 1._O
  N0 = int(1.1_O * a0) + 1
  f0 = ENVJ(N0,a0) - MP
  N1 = N0 + 5
  f1 = ENVJ(N1,a0) - MP
  do it = 1, 50   
    NN = N1 - int((N1 - N0) / (1._O - f0 / f1))   
    f  = ENVJ(NN,a0) - MP
    if (abs(NN - N1) < 1) exit
    N0 = N1
    f0 = f1
    N1 = NN
    f1 = f
  end do 
  MSTA = NN  
end function MSTA1
!***********************************************************************************
function MSTA2 (x, N, MP) result (MSTA)
  use parameters
  implicit none 
  real(O), intent(in) :: x
  integer,intent(in)  :: N, MP 
  integer             :: MSTA, N0, N1, it, NN
  real(O)             :: a0, hmp, ejn, obj, f0, f1, f, ENVJ  
!
  a0 = abs(X)
  if (a0 < 1._O) a0 = 1._O
  hmp = 0.5_O * MP
  ejn = ENVJ(N,a0)
  if (ejn <= hmp) then
    obj = 1._O * MP
    N0  = int(1.1_O * a0)
  else
    obj = hmp + ejn
    N0  = N
  end if
  f0 = ENVJ(N0,a0) - obj
  N1 = N0 + 5
  f1 = ENVJ(N1,A0) - obj
  do it = 1, 50
    NN = N1 - int((N1 - N0) / (1._O - f0 / f1))
    f  = ENVJ(NN,a0) - obj
    if (abs(NN - N1) < 1) exit
    N0 = N1
    f0 = f1
    N1 = NN
    f1 = f
  end do
  MSTA = NN + 10
end function MSTA2
!***********************************************************************************
function ENVJ (N, x) result (ENV)  
  use parameters
  implicit none 
  real(O), intent(in) :: x
  integer,intent(in)  :: N 
  real(O)             :: ENV
!
  ENV = 0.5_O * log10(6.28_O * N) - N * log10(1.36_O * x / N) 
end function ENVJ



