! **********************************************************************************
! *    ROUTINES FOR COMPUTING THE SCATTERED FIELD COEFFICIENTS FOR SPHERES AND     *
! *    THE Q MATRIX FOR UNIAXIAL ANISOTROPIC PARTICLES AND PARTICLES ON OR NEAR    *
! *                              A PLANE SURFACE                                   * 
! *    -----------------------------------------------------------------------     *
! *    Partial list of subroutines:                                                *
! *      coefficients_fg,                  coefficients_fg_m,                      *
! *      coefficients_fg_m_coated,         coefficients_fg_m_rec,                  *
! *      vector_Q_sphere,                  vector_Q_sphere1,                       *
! *      vector_Q_sphere1_m,               An,      Bn,                            *
! *      matrix_S1,                        matrix_Q_anis,                          *
! *      ConvergenceMatrixElemAnis,        matrixPARTSUB                           *
! *************************************
subroutine coefficients_fg (wavenumber, r, ind_ref, Mrank, Nrank, Nmax, c)
!-----------------------------------------------------------------------------------
! The routine computes the scattered field coefficients for a spherical particle   !
! and the azimuthal modes m = 0,1,..., Mrank.                                      !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none      
  integer    :: Mrank, Nrank, Nmax
  real(O)    :: r, wavenumber
  complex(O) :: ind_ref, c(2*Nmax)
!
  integer    :: k, n, m, l, N0
  real(O)    :: nr
  complex(O) :: x, f, f1, f2
  complex(O),allocatable :: j(:), y(:), h(:), hd(:), A(:)
!
  allocate (j(0:Nrank), y(0:Nrank), h(0:Nrank), hd(0:Nrank), A(Nrank))
  x = cmplx(wavenumber * r,0.0,O)
  if (abs(x) < MachEps) x = cmplx(MachEps,MachEps,O)
  call besel_h_complete (x, Nrank, j, y, h, hd)
  call An (x, ind_ref, Nrank, A)
  do m = 0, Mrank         
    if (m == 0) then
      do k = 1, Nrank
        n  = k
        nr = real(n,O)
        f  = ind_ref * A(n) + nr / x                                                                                             
        f1 = f * j(n) - j(n-1)
        f2 = f * h(n) - h(n-1)
        c(k) = - f1 / f2
        f  = A(n) / ind_ref + nr / x
        f1 = f * j(n) - j(n-1)
        f2 = f * h(n) - h(n-1)
        c(Nmax+k) = - f1 / f2 
      end do
    else         
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      do l = 1, 2
        do k = 1, Nrank - m + 1
          n  = m + k - 1
          nr = real(n,O)
          f  = ind_ref * A(n) + nr / x                                                                                           
          f1 = f * j(n) - j(n-1)
          f2 = f * h(n) - h(n-1)
          c(N0+k) = - f1 / f2 
          f  = A(n) / ind_ref + nr / x
          f1 = f * j(n) - j(n-1)
          f2 = f * h(n) - h(n-1)
          c(Nmax+N0+k) = - f1 / f2 
        end do
        N0 = N0 + Nrank - m + 1
      end do       
    end if
  end do      
  deallocate (j, y, h, hd, A)      
end subroutine coefficients_fg
!***********************************************************************************
subroutine coefficients_fg_m (wavenumber, r, ind_ref, m, Nrank, Nmax, c)
!-----------------------------------------------------------------------------------
! The routine computes the scattered field coefficients for a spherical particle   !
! and the azimuthal modes m.                                                       !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none  
  integer    :: m, Nrank, Nmax
  real(O)    :: r, wavenumber
  complex(O) :: ind_ref, c(2*Nmax)
!    
  integer    :: k, n
  real(O)    :: nr
  complex(O) :: x, f, f1, f2
  complex(O),allocatable :: j(:), y(:), h(:), hd(:), A(:)
!
  allocate (j(0:Nrank), y(0:Nrank), h(0:Nrank), hd(0:Nrank), A(Nrank))
  x = cmplx(wavenumber * r,0.0,O)
  if (abs(x) < MachEps) x = cmplx(MachEps,MachEps,O)
  call besel_h_complete (x, Nrank, j, y, h, hd)
  call An (x, ind_ref, Nrank, A)
  if (m == 0) then
    do k = 1, Nmax
      n  = k
      nr = real(n,O)
      f  = ind_ref * A(n) + nr / x                                                                                   
      f1 = f * j(n) - j(n-1)
      f2 = f * h(n) - h(n-1)
      c(k) = - f1 / f2
      f  = A(n) / ind_ref + nr / x 
      f1 = f * j(n) - j(n-1)
      f2 = f * h(n) - h(n-1)
      c(Nmax+k) = - f1 / f2 
    end do
  else
    do k = 1, Nmax
      n  = m + k - 1
      nr = real(n,O)
      f  = ind_ref * A(n) + nr / x                                                                                           
      f1 = f * j(n) - j(n-1)
      f2 = f * h(n) - h(n-1)
      c(k) = - f1 / f2 
      f  = A(n) / ind_ref + nr / x 
      f1 = f * j(n) - j(n-1)
      f2 = f * h(n) - h(n-1)
      c(Nmax+k) = - f1 / f2             
    end do                                                 
  end if        
  deallocate (j, y, h, hd, A)      
end subroutine coefficients_fg_m
!***********************************************************************************
subroutine coefficients_fg_m_coated (wavenumber, Npart, r, ind_ref, m, Nrank,       &
           Nmax, c)
!-----------------------------------------------------------------------------------
! The routine computes the scattered field coefficients for a coated spherical     !
! particle and the azimuthal modes m.                                              !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none  
  integer    :: m, Npart, Nrank, Nmax
  real(O)    :: wavenumber, r(Npart) 
  complex(O) :: ind_ref(Npart), c(2*Nmax)   
!
  integer    :: p
  complex(O) :: x, ind_refC 
!
  if (Npart == 1) then
    x = cmplx(wavenumber * r(1),0.0,O)
    ind_refC = ind_ref(1)
  else
    x = ind_ref(Npart-1) * wavenumber * r(Npart)
    ind_refC = ind_ref(Npart) / ind_ref(Npart-1)
  end if
  if (abs(x) < MachEps) x = cmplx(MachEps,MachEps,O)
  call coefficients_fg_m_rec (.true., x, ind_refC, m, Nrank, Nmax, c)
  if (Npart > 1) then
    do p = Npart - 1, 1, -1
      if (p > 1) then
        x = ind_ref(p-1) * wavenumber * r(p)
        ind_refC = ind_ref(p) / ind_ref(p-1)
      else
        x = cmplx(wavenumber * r(1),0.0,O)
        ind_refC = ind_ref(1)
      end if
      if (abs(x) < MachEps) x = cmplx(MachEps,MachEps,O)
      call coefficients_fg_m_rec (.false., x, ind_refC, m, Nrank, Nmax, c)
    end do
  end if
end subroutine coefficients_fg_m_coated    
!***********************************************************************************
subroutine coefficients_fg_m_rec (firstcall, x, ind_ref, m, Nrank, Nmax, c)
!-----------------------------------------------------------------------------------
! The routine computes the scattered field coefficients for a coated sphere and    !
! the azimuthal mode m using a reccurence relation.                                !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none  
  integer    :: m, Nrank, Nmax
  complex(O) :: x, ind_ref, c(2*Nmax)
  logical    :: firstcall
!    
  integer    :: k, n
  real(O)    :: nr
  complex(O) :: xm, f, f1, f2, b1, b2, pn, sn, qn, rn
  complex(O),allocatable :: j(:), y(:), h(:), hd(:), A(:), B(:), jm(:), ym(:),      &
                            hm(:), hdm(:)
!
  allocate (j(0:Nrank), y(0:Nrank), h(0:Nrank), hd(0:Nrank), A(Nrank))   
  call besel_h_complete (x, Nrank, j, y, h, hd)
  call An (x, ind_ref, Nrank, A)
  if (.not. firstcall) then
    xm = x * ind_ref
    allocate (jm(0:Nrank), ym(0:Nrank), hm(0:Nrank), hdm(0:Nrank), B(Nrank)) 
    call besel_h_complete (xm, Nrank, jm, ym, hm, hdm)  
    call Bn (x, ind_ref, Nrank, B)
  end if
  if (m == 0) then
    do k = 1, Nmax
      n  = k
      nr = real(n,O)
      if (firstcall) then
        pn = one
        sn = one
        qn = one
        rn = one
      else
        b1 = hm(n) / jm(n)
        b2 = B(n) / A(n)
        pn = one + c(k) * b1
        sn = one + c(Nmax+k) * b1
        qn = one + c(Nmax+k) * b1 * b2
        rn = one + c(k) * b1 * b2
      end if
      f  = rn * ind_ref * A(n) + pn * nr / x
      f1 = f * j(n) - pn * j(n-1)
      f2 = f * h(n) - pn * h(n-1)
      c(k) = - f1 / f2
      f  = qn * A(n) / ind_ref + sn * nr / x 
      f1 = f * j(n) - sn * j(n-1)
      f2 = f * h(n) - sn * h(n-1)
      c(Nmax+k) = - f1 / f2 
    end do
  else
    do k = 1, Nmax
      n  = m + k - 1
      nr = real(n,O)
      if (firstcall) then
        pn = one
        sn = one
        qn = one
        rn = one
      else
        b1 = hm(n) / jm(n)
        b2 = B(n) / A(n)
        pn = one + c(k) * b1
        sn = one + c(Nmax+k) * b1
        qn = one + c(Nmax+k) * b1 * b2
        rn = one + c(k) * b1 * b2
      end if
      f  = rn * ind_ref * A(n) + pn * nr / x
      f1 = f * j(n) - pn * j(n-1)
      f2 = f * h(n) - pn * h(n-1)
      c(k) = - f1 / f2
      f  = qn * A(n) / ind_ref + sn * nr / x 
      f1 = f * j(n) - sn * j(n-1)
      f2 = f * h(n) - sn * h(n-1)
      c(Nmax+k) = - f1 / f2             
    end do                                                 
  end if        
  deallocate (j, y, h, hd, A)
  if (.not.firstcall) deallocate (jm, ym, hm, hdm, B)      
end subroutine coefficients_fg_m_rec
!***********************************************************************************
subroutine vector_Q_sphere (index1, index2, wavenumber, r, ind_ref, Mrank, Nrank,   &
           Nmax, c)     
!-----------------------------------------------------------------------------------
! The routine computes the Q vector for a spherical particle and the azimuthal     !
! modes m = 0,1,...,Mrank.                                                         !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none  
  integer    :: Mrank, Nrank, Nmax, index1, index2
  real(O)    :: r, wavenumber
  complex(O) :: ind_ref, c(2*Nmax)
!
  integer    :: k, n, m, l, N0
  complex(O) :: x, mx, f, f1, f2
  complex(O),allocatable :: j(:), jd(:), h(:), hd(:)
!
  allocate (j(0:Nrank), jd(0:Nrank), h(0:Nrank), hd(0:Nrank)) 
  x = cmplx(wavenumber * r,0.0,O)
  if (abs(x) < MachEps) x = cmplx(MachEps,MachEps,O)
  mx = x * ind_ref
  if (index1 == 3 .and. index2 == 1) then
    call besel_h (x, Nrank, h, hd)
    call besel_j (mx, Nrank, j, jd)
  else if (index1 == 3 .and. index2 == 3) then
    call besel_h (x, Nrank, h, hd)
    call besel_h (mx, Nrank, j, jd)
  else if (index1 == 1 .and. index2 == 1) then
    call besel_j (x, Nrank, h, hd)
    call besel_j (mx, Nrank, j, jd)
  else if (index1 == 1 .and. index2 == 3) then
    call besel_j (x, Nrank, h, hd)
    call besel_h (mx, Nrank, j, jd)   
  end if
  f = im * x
  do m = 0, Mrank         
    if (m == 0) then
      do k = 1, Nrank
        n  = k  
        f1 = f * j(n) * hd(n)
        f2 = f * h(n) * jd(n)                                                                                    
        c(k)      = f1 - f2 
        c(Nmax+k) = (ind_ref**2 * f1 - f2) / ind_ref
      end do
    else         
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      do l = 1, 2
        do k = 1, Nrank - m + 1
          n = m + k - 1
          f1 = f * j(n) * hd(n)
          f2 = f * h(n) * jd(n)
          c(N0+k)      = f1 - f2 
          c(Nmax+N0+k) = (ind_ref**2 * f1 - f2) / ind_ref               
        end do                                             
        N0 = N0 + Nrank - m + 1
      end do      
    end if
  end do
  deallocate (j, jd, h, hd)     
end subroutine vector_Q_sphere
!***********************************************************************************
subroutine vector_Q_sphere1 (index1, index2, wavenumber, r, ind_ref, Mrank, Nrank,  &
           Nmax, c)
!-----------------------------------------------------------------------------------
! The routine computes the Q vector for a spherical particle and the azimuthal     !
! modes m = 0,1,..., Mrank. In contrast to the subroutine "vector_Q_sphere" the    !
! derivative of the Bessel and Hankel functions are expressed in terms of the      !
! An and Bn coefficients.                                                          !
!-----------------------------------------------------------------------------------
  use parameters        
  use derived_parameters
  implicit none      
  integer    :: Mrank, Nrank, Nmax, index1, index2
  real(O)    :: r, wavenumber
  complex(O) :: ind_ref, c(2*Nmax)
!
  integer    :: k, n, m, l, N0
  real(O)    :: nr
  complex(O) :: x, mx, f, f1, f2, fn
  complex(O),allocatable :: j(:), jd(:), h(:), hd(:), A(:)
!
  allocate (j(0:Nrank), jd(0:Nrank), h(0:Nrank), hd(0:Nrank), A(Nrank)) 
  x = cmplx(wavenumber * r,0._O)
  if (abs(x) < MachEps) x = cmplx(MachEps,MachEps,O)
  mx = x * ind_ref
  if (index2 == 1) then
    call An (x, ind_ref, Nrank, A)
  else
    call Bn (x, ind_ref, Nrank, A)
  end if
  if (index1 == 3 .and. index2 == 1) then
    call besel_h (x, Nrank, h, hd)
    call besel_j (mx, Nrank, j, jd)
  else if(index1 == 3 .and. index2 == 3) then
    call besel_h (x, Nrank, h, hd)
    call besel_h (mx, Nrank, j, jd)
  else if(index1 == 1 .and. index2 == 1) then
    call besel_j (x, Nrank, h, hd)
    call besel_j (mx, Nrank, j, jd)
  else if (index1 == 1 .and. index2 == 3) then
    call besel_j (x, Nrank, h, hd)
    call besel_h (mx, Nrank, j, jd)   
  end if
  f = - im * x * x
  do m = 0, Mrank         
    if (m == 0) then
      do k = 1, Nrank
        n  = k                          
        nr = real(n,O)
        f1 = ind_ref * A(n) + nr / x
        f2 = A(n) / ind_ref + nr / x                                                             
        fn = f * j(n)
        c(k) = fn * (f1*h(n) - h(n-1))
        c(Nmax+k) = ind_ref * fn * (f2 * h(n) - h(n-1))
      end do
    else
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      do l = 1, 2
        do k = 1, Nrank - m + 1
          n = m + k - 1
          nr = real(n,O)
          f1 = ind_ref * A(n) + nr / x
          f2 = A(n) / ind_ref + nr / x                                                           
          fn = f * j(n)
          c(N0+k) = fn * (f1 * h(n) - h(n-1))
          c(Nmax+N0+k) = ind_ref * fn * (f2 * h(n) - h(n-1))            
        end do                                             
        N0 = N0 + Nrank - m + 1
      end do      
    end if
  end do
  deallocate (j, jd, h, hd, A)      
end subroutine vector_Q_sphere1
!***********************************************************************************
subroutine vector_Q_sphere1_m (index1, index2, wavenumber, r,ind_ref, m, Nrank,     &
           Nmax, c)
!-----------------------------------------------------------------------------------
! The routine computes the Q vector for a spherical particle and the azimuthal     !
! mode m.                                                                          !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none  
  integer    :: m, Nrank, Nmax, index1, index2
  real(O)    :: wavenumber, r
  complex(O) :: ind_ref, c(2*Nmax)
!
  integer    :: k, n
  real(O)    :: nr
  complex(O) :: x, mx, f, f1, f2, fn
  complex(O),allocatable :: j(:), jd(:), h(:), hd(:), A(:)
!
  allocate (j(0:Nrank), jd(0:Nrank), h(0:Nrank), hd(0:Nrank), A(Nrank)) 
  x = cmplx(wavenumber * r,0.0,O)
  if (abs(x) < MachEps) x = cmplx(MachEps,MachEps,O)
  mx = x * ind_ref
  if (index2 == 1) then
    call An (x, ind_ref, Nrank, A)
  else
    call Bn (x, ind_ref, Nrank, A)
  end if
  if (index1 == 3 .and. index2 == 1) then
    call besel_h (x, Nrank, h, hd)
    call besel_j (mx, Nrank, j, jd)
  else if (index1 == 3 .and. index2 == 3) then
    call besel_h (x, Nrank, h, hd)
    call besel_h (mx, Nrank, j, jd)
  else if (index1 == 1 .and. index2 == 1) then
    call besel_j (x, Nrank, h, hd)
    call besel_j (mx, Nrank, j, jd)
  else if (index1 == 1 .and. index2 == 3) then
    call besel_j (x, Nrank, h, hd)
    call besel_h (mx, Nrank, j, jd)   
  end if 
  f = - im * x * x     
  if (m == 0) then
    do k = 1, Nmax
      n = k
      nr = real(n,O)
      f1 = ind_ref * A(n) + nr / x
      f2 = A(n) / ind_ref + nr / x                                                           
      fn = f * j(n)                                                                                          
      c(k)      = fn * (f1 * h(n) - h(n-1))
      c(Nmax+k) = ind_ref * fn * (f2 * h(n) - h(n-1))
    end do
  else
    do k = 1, Nmax
      n = m + k - 1
      nr = real(n,O)
      f1 = ind_ref * A(n) + nr / x
      f2 = A(n) / ind_ref + nr / x                                                           
      fn = f * j(n)
      c(k)      = fn * (f1 * h(n) - h(n-1))
      c(Nmax+k) = ind_ref * fn * (f2 * h(n) - h(n-1))           
    end do                                                 
  end if          
  deallocate (j, jd, h, hd, A)      
end subroutine vector_Q_sphere1_m       
!***********************************************************************************
subroutine An (x, m, Nrank, A)
!-----------------------------------------------------------------------------------
! The routine computes                                                             ! 
!                  A(n) = d[ln(x*jn(x))] = [x*jn(x)]'/[x*jn(x)],                   !
! where jn are the Bessel functions.                                               !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none 
  integer    :: Nrank
  complex(O) :: x, m, A(Nrank)
!       
  real(O)    :: DNc, DNmax, n1r
  integer    :: Nmax, n
  complex(O) :: mx, fp, f, tamp
!
  DNc   = abs(x + 4._O * x**0.33_O + 2._O)
  DNmax = max(DNc,abs(m * x)) + 100._O
  Nmax  = int(DNmax)
  mx = m * x
  fp = zero
  do n = Nmax, 1, -1
    n1r  = real(n + 1,O)
    tamp = n1r / mx
    f = tamp - 1._O / (fp + tamp)
    if (n <= Nrank) A(n) = f
    fp = f            
  end do        
end subroutine An
!***********************************************************************************
subroutine Bn (x, m, Nrank, B)
!-----------------------------------------------------------------------------------
! The routine computes                                                             !
!                    B(n) = d[ln(x*hn(x))] = [x*hn(x)]'/[x*hn(x)],                 !
! where hn are the Hankel functions.                                               !  
!-----------------------------------------------------------------------------------
  use parameters
  implicit none 
  integer    :: Nrank
  complex(O) :: x, m, B(Nrank)
!       
  integer    :: n
  real(O)    :: n1r
  complex(O) :: mx, smx, cmx, j0, j1, y0, y1, h0, h1, tamp
!
  mx  = m * x
  smx = sin(mx)
  cmx = cos(mx)
  j0  = smx / mx
  j1  = j0 / mx - cmx / mx
  y0  = - cmx / mx
  y1  = y0 / mx - smx / mx
  h0  = j0 + im * y0
  h1  = j1 + im * y1
  B(1)= h0 / h1 - 1._O / mx
  do n = 1, Nrank - 1
    n1r  = real(n + 1,O)          
    tamp = n1r / mx 
    B(n+1) = - tamp + 1._O / (tamp - B(n))         
  end do        
end subroutine Bn
!***********************************************************************************
subroutine matrix_S1 (Nmax, Nmax1, S2, n2p, m2p, S1, n1p, m1p)
  use parameters
  implicit none
  integer    :: Nmax, Nmax1, n1p, m1p, n2p, m2p, k, k1  
  complex(O) :: S1(2*n1p,2*m1p), S2(2*n2p,2*m2p)
!
  do k = 1, Nmax
    do k1 = 1, Nmax1
      S1(k,k1) = S2(k1,k)
      S1(k,k1+Nmax1) = - S2(k1+Nmax1,k)
      S1(k+Nmax,k1)  =   S1(k,k1+Nmax1)
      S1(k+Nmax,k1+Nmax1) = S1(k,k1)
    end do     
  end do      
end subroutine matrix_S1
!***********************************************************************************
subroutine matrix_Q_anis (FileGeom, TypeGeom, index1, index2, k, ind_ref, ind_refZ, &
           alphaP, betaP, gamaP, Nsurf, surf, rpX, npX, areaX, Nface, Mrank, Nrank, &
           Nmax, Nbeta, NintAL, Nparam, Nintparam, paramG1, paramG2, weightsG,      &
           A, nap, map)
!-----------------------------------------------------------------------------------
! The routine compute the Q matrix of an uniaxial anisotropic particle.            ! 
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer     :: TypeGeom, index1, index2, Mrank, Nrank, Nmax, Nsurf, Nface,         &
                 NintAL, Nparam, Nintparam(Nparam), Nbeta, nap, map
  real(O)     :: k, alphaP, betaP, gamaP, surf(Nsurf), paramG1(Nparam,NintAL*NintAL),&
                 paramG2(Nparam,NintAL*NintAL), weightsG(Nparam,NintAL*NintAL),      &
                 rpX(3,NfacePD), npX(3,NfacePD), areaX(NfacePD)
  complex(O)  :: ind_ref, ind_refZ, A(2*nap,2*map)
  logical     :: FileGeom
!
  integer     :: i, j, pint, iparam, Nintl
  real(O)     :: r, theta, phi, dA, param1, param2, pondere, sign, n(3), x, y, z,    &
                 xl, yl, zl, thetal, phil, rl, nr, nt, np, nx, ny, nz, nxl, nyl, nzl 
  complex(O)  :: ki, zkl, ml(3), nl(3), xec(3), yec(3), v1, v2, fact, f, xhc(3),     &
                 yhc(3), mixt_product
  complex(O),allocatable :: Xe(:,:), Ye(:,:), m3(:,:), n3(:,:), Xh(:,:), Yh(:,:)
!
  allocate (Xe(3,Nmax), Ye(3,Nmax), Xh(3,Nmax), Yh(3,Nmax), m3(3,Nmax), n3(3,Nmax)) 
  do i = 1, 2*Nmax
    do j = 1, 2*Nmax
      A(i,j) = zero
    end do
  end do      
  ki   = k * ind_ref
  sign = 1._O
  if (index1 == 3 .and. index2 == 1) sign = - 1._O
  f = sign * im * k * k / Pi
  if (.not. FileGeom) then
    do iparam = 1, Nparam
      Nintl = Nintparam(iparam)
      do pint = 1, Nintl          
        param1  = paramG1(iparam,pint)
        param2  = paramG2(iparam,pint)
        pondere = weightsG(iparam,pint)
        call elem_geom3D (TypeGeom, Nsurf, surf, param1, param2, iparam, r, theta,   &
             phi, dA, n, .false., 1)                                                      
        call MN_anisotrop (ki, ind_ref, ind_refZ, r, theta, phi, alphaP, betaP,      &
             gamaP, Mrank, Nrank, Nmax, Nbeta, Xe, Ye, Xh, Yh)
!       --- compute the localized SWVF in the principal coordinate system ----
        x = r * sin(theta) * cos(phi)
        y = r * sin(theta) * sin(phi)
        z = r * cos(theta) 
        call T_cartesian_global_local (x, y, z, alphaP, betaP, gamaP, xl, yl, zl)  
        call T_cartesian_spherical (xl, yl, zl, rl, thetal, phil) 
        zkl = cmplx(k * rl,0.0,O) 
        if (index1 == 3 .and. index2 == 1) then       
          call MN_complete (3, zkl, thetal, phil, Mrank, Nrank, Nmax,.false.,.false.,&
               m3, n3)
        else if (index1 == 1 .and. index2 == 1) then            
          call MN_complete (1, zkl, thetal, phil, Mrank, Nrank, Nmax,.false.,.false.,&
               m3, n3)
        end if
!       --- transform the unit normal vector in the principal coordinate system ---
        nr = n(1)
        nt = n(2)
        np = n(3)
        nx = sin(theta) * cos(phi)  * nr + cos(theta) * cos(phi)  * nt - sin(phi) * np 
        ny = sin(theta) * sin(phi)  * nr + cos(theta) * sin(phi)  * nt + cos(phi) * np
        nz =             cos(theta) * nr -             sin(theta) * nt 
        call T_cartesian_global_local (nx, ny, nz, alphaP, betaP, gamaP, nxl, nyl, nzl)
        n(1) =   sin(thetal) * cos(phil) * nxl + sin(thetal) * sin(phil) * nyl +     &
                 cos(thetal) * nzl
        n(2) =   cos(thetal) * cos(phil) * nxl + cos(thetal) * sin(phil) * nyl -     &
                 sin(thetal) * nzl
        n(3) =               - sin(phil) * nxl +               cos(phil) * nyl      
        fact = f * dA * pondere
        do i = 1, Nmax
          ml(1) = m3(1,i)
          ml(2) = m3(2,i)
          ml(3) = m3(3,i)
          nl(1) = n3(1,i)
          nl(2) = n3(2,i)
          nl(3) = n3(3,i)
          do j = 1, Nmax
            xec(1) = Xe(1,j)
            xec(2) = Xe(2,j)
            xec(3) = Xe(3,j)
            yec(1) = Ye(1,j)
            yec(2) = Ye(2,j)
            yec(3) = Ye(3,j)
            xhc(1) = Xh(1,j)
            xhc(2) = Xh(2,j)
            xhc(3) = Xh(3,j)
            yhc(1) = Yh(1,j)
            yhc(2) = Yh(2,j)
            yhc(3) = Yh(3,j)
            v1 = mixt_product (n, xec, nl)
            v2 = mixt_product (n, xhc, ml)                                                                             
            A(i,j) = A(i,j) + (v1 + ind_ref * v2) * fact
            v1 = mixt_product (n, yec, nl)
            v2 = mixt_product (n, yhc, ml)
            A(i,j+Nmax) = A(i,j+Nmax) + (v1 + ind_ref * v2) * fact
            v1 = mixt_product (n, xec, ml)
            v2 = mixt_product (n, xhc, nl)
            A(i+Nmax,j) = A(i+Nmax,j) + (v1 + ind_ref * v2) * fact
            v1 = mixt_product (n, yec, ml)
            v2 = mixt_product (n, yhc, nl)
            A(i+Nmax,j+Nmax) = A(i+Nmax,j+Nmax) + (v1 + ind_ref * v2) * fact          
          end do
        end do    
      end do
    end do
  else
    do pint = 1, Nface
      x = rpX(1,pint)
      y = rpX(2,pint)
      z = rpX(3,pint)
      call T_cartesian_spherical (x, y, z, r, theta, phi)
      dA = areaX(pint)
      nr =   sin(theta) * cos(phi) * npX(1,pint) +                                   &
             sin(theta) * sin(phi) * npX(2,pint) + cos(theta) * npX(3,pint) 
      nt =   cos(theta) * cos(phi) * npX(1,pint) +                                   &
             cos(theta) * sin(phi) * npX(2,pint) - sin(theta) * npX(3,pint)
      np = - sin(phi) * npX(1,pint) + cos(phi) * npX(2,pint)
      call MN_anisotrop (ki, ind_ref, ind_refZ, r, theta, phi, alphaP, betaP, gamaP, &
           Mrank, Nrank, Nmax, Nbeta, Xe, Ye, Xh, Yh)
!     --- compute the localized SWVF in the principal coordinate system ----       
      call T_cartesian_global_local (x, y, z, alphaP, betaP, gamaP, xl, yl, zl)  
      call T_cartesian_spherical (xl, yl, zl, rl, thetal, phil) 
      zkl = cmplx(k * rl,0.0,O) 
      if (index1 == 3 .and. index2 == 1) then       
        call MN_complete (3, zkl, thetal, phil, Mrank, Nrank, Nmax,.false.,.false.,  &
             m3, n3)
      else if (index1 == 1 .and. index2 == 1) then      
        call MN_complete (1, zkl, thetal, phil, Mrank, Nrank, Nmax,.false.,.false.,  &
             m3, n3)
      end if
!     --- transform the unit normal vector in the principal coordinate system ---     
      nx = sin(theta) * cos(phi)  * nr + cos(theta) * cos(phi)  * nt - sin(phi) * np 
      ny = sin(theta) * sin(phi)  * nr + cos(theta) * sin(phi)  * nt + cos(phi) * np
      nz =             cos(theta) * nr -             sin(theta) * nt 
      call T_cartesian_global_local (nx, ny, nz, alphaP, betaP, gamaP, nxl, nyl, nzl)
      n(1) =   sin(thetal) * cos(phil) * nxl + sin(thetal) * sin(phil) * nyl +       &
               cos(thetal) * nzl
      n(2) =   cos(thetal) * cos(phil) * nxl + cos(thetal) * sin(phil) * nyl -       &
               sin(thetal) * nzl
      n(3) =               - sin(phil) * nxl +               cos(phil) * nyl      
      fact = f * dA
      do i = 1, Nmax
        ml(1) = m3(1,i)
        ml(2) = m3(2,i)
        ml(3) = m3(3,i)
        nl(1) = n3(1,i)
        nl(2) = n3(2,i)
        nl(3) = n3(3,i)
        do j = 1, Nmax
          xec(1) = Xe(1,j)
          xec(2) = Xe(2,j)
          xec(3) = Xe(3,j)
          yec(1) = Ye(1,j)
          yec(2) = Ye(2,j)
          yec(3) = Ye(3,j)
          xhc(1) = Xh(1,j)
          xhc(2) = Xh(2,j)
          xhc(3) = Xh(3,j)
          yhc(1) = Yh(1,j)
          yhc(2) = Yh(2,j)
          yhc(3) = Yh(3,j)
          v1 = mixt_product (n, xec, nl)
          v2 = mixt_product (n, xhc, ml)                                                                               
          A(i,j) = A(i,j) + (v1 + ind_ref * v2) * fact
          v1 = mixt_product (n, yec, nl)
          v2 = mixt_product (n, yhc, ml)
          A(i,j+Nmax) = A(i,j+Nmax) + (v1 + ind_ref * v2) * fact
          v1 = mixt_product (n, xec, ml)
          v2 = mixt_product (n, xhc, nl)
          A(i+Nmax,j) = A(i+Nmax,j) + (v1 + ind_ref * v2) * fact
          v1 = mixt_product (n, yec, ml)
          v2 = mixt_product (n, yhc, nl)
          A(i+Nmax,j+Nmax) = A(i+Nmax,j+Nmax) + (v1 + ind_ref * v2) * fact          
        end do
      end do          
    end do
  end if                            
  deallocate (Xe, Ye, Xh, Yh, m3, n3)
end subroutine matrix_Q_anis 
!***********************************************************************************
subroutine ConvergenceMatrixElemAnis (TypeGeom, index1, index2, k, ind_ref,         &
           ind_refZ, alphaP, betaP, gamaP, Nsurf, surf, mline, mcol, Nrank,         &
           NintMax, Nint1, Nint2, Nbeta, Nparam)
  use parameters
  use derived_parameters
  implicit none
  integer     :: TypeGeom, index1, index2, mline, mcol, Nrank, NintMax, Nint1,      &
                 Nint2, Nsurf, Nparam, Nbeta
  real(O)     :: k, alphaP, betaP, gamaP, surf(Nsurf)
  complex(O)  :: ind_ref, ind_refZ
!
  integer     :: NmaxL, NmaxC, ml_minus, i, j, pint, iparam, Nint1l, Nint2l,        &
                 Nintl, NintAL, iint, imn, jmn
  real(O)     :: r, theta, phi, dA, param1, param2, pondere, sign, n(3), x, y, z,   &
                 xl, yl, zl, thetal, phil, rl, nr, nt, np, nx, ny, nz, nxl, nyl,    &
                 nzl, MaxA(2), MinA(2)
  complex(O)  :: ki, zkl, ml(3), nl(3), xec(3), yec(3), v1, v2, fact, f, xhc(3),    &
                 yhc(3), mixt_product, A11(2,2), A12(2,2), A21(2,2), A22(2,2),      &
                 expfctL
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG1(:,:), paramG2(:,:), weightsG(:,:) 
  complex(O),allocatable :: Xe(:,:), Ye(:,:), m3(:,:), n3(:,:), Xh(:,:), Yh(:,:)
!
  NmaxL = Nrank - mline + 1  ! mline >= 0
  NmaxC = Nrank - mcol  + 1  ! mcol  >= 0
  ml_minus = - mline
  allocate (Xe(3,NmaxC), Ye(3,NmaxC), Xh(3,NmaxC), Yh(3,NmaxC),                     &
            m3(3,NmaxL), n3(3,NmaxL))      
  ki   = k * ind_ref
  sign = 1._O
  if (index1 == 3 .and. index2 == 1) sign = - 1._O
  f = sign * im * k * k / Pi
  if (mline == 1) then 
    print "(2x,'Integral Test over Matrix Elements of Blocks (1,Mrank):')"
  else 
    print "(2x,'Integral Test over Matrix Elements of Blocks (Mrank,Mrank):')"
  end if
  print "(2x,'Nint1   Nint2          max-min (A11,A22)             max-min (A12,A21) ')"    
  do iint = 1, NintMax
    Nint1l = Nint1 + iint - 1
    Nint2l = Nint2 + iint - 1
    NintAL = max(Nint1l,Nint2l)                      
    allocate (paramG1(Nparam,NintAL*NintAL), paramG2(Nparam,NintAL*NintAL),         &
              weightsG(Nparam,NintAL*NintAL))
    allocate (Nintparam(Nparam))
    call interpolation_list3D (TypeGeom, Nsurf, surf, Nint1l, Nint2l, NintAL,       &
         Nparam, Nintparam, paramG1, paramG2, weightsG, .false., 1) 
    do i = 1, 2
      do j = 1, 2
        A11(i,j) = zero
        A12(i,j) = zero
        A21(i,j) = zero
        A22(i,j) = zero
      end do
    end do   
    do iparam = 1, Nparam
      Nintl = Nintparam(iparam)
      do pint = 1, Nintl          
        param1  = paramG1(iparam,pint)
        param2  = paramG2(iparam,pint)
        pondere = weightsG(iparam,pint)
        call elem_geom3D (TypeGeom, Nsurf, surf, param1, param2, iparam, r, theta,   &
             phi, dA, n, .false., 1)                                                      
        call MN_m_anisotrop (ki, ind_ref, ind_refZ, r, theta, phi, alphaP, betaP,    &
             gamaP, mcol, Nrank, NmaxC, Nbeta, Xe, Ye, Xh, Yh)        
!       --- compute the localized SWVF in the principal coordinate system ----
        x = r * sin(theta) * cos(phi)
        y = r * sin(theta) * sin(phi)
        z = r * cos(theta) 
        call T_cartesian_global_local (x, y, z, alphaP, betaP, gamaP, xl, yl, zl)  
        call T_cartesian_spherical (xl, yl, zl, rl, thetal, phil) 
        zkl = cmplx(k * rl,0.0,O) 
        if (index1 == 3 .and. index2 == 1) then  
          call MN (3, zkl, thetal, ml_minus, Nrank, NmaxL, m3, n3)                                 
        else if (index1 == 1 .and. index2 == 1) then            
          call MN (1, zkl, thetal, ml_minus, Nrank, NmaxL, m3, n3)          
        end if
!       --- transform the unit normal vector in the principal coordinate system ---
        nr = n(1)
        nt = n(2)
        np = n(3)
        nx = sin(theta) * cos(phi)  * nr + cos(theta) * cos(phi)  * nt - sin(phi) * np 
        ny = sin(theta) * sin(phi)  * nr + cos(theta) * sin(phi)  * nt + cos(phi) * np
        nz =             cos(theta) * nr -             sin(theta) * nt 
        call T_cartesian_global_local (nx, ny, nz, alphaP, betaP, gamaP, nxl, nyl, nzl)
        n(1) = sin(thetal) * cos(phil) * nxl + sin(thetal) * sin(phil) * nyl + cos(thetal) * nzl
        n(2) = cos(thetal) * cos(phil) * nxl + cos(thetal) * sin(phil) * nyl - sin(thetal) * nzl
        n(3) =             - sin(phil) * nxl +               cos(phil) * nyl      
        fact = f * dA * pondere
        expfctL =  exp(im * ml_minus * phil)            
        do i = 1, 2
          if (i == 1) then
            imn = 1
          else
            imn = NmaxL
          end if            
          ml(1) = m3(1,imn) * expfctL
          ml(2) = m3(2,imn) * expfctL
          ml(3) = m3(3,imn) * expfctL
          nl(1) = n3(1,imn) * expfctL
          nl(2) = n3(2,imn) * expfctL
          nl(3) = n3(3,imn) * expfctL
          do j = 1, 2
            if (j == 1) then
              jmn = 1
            else
              jmn = NmaxC
            end if           
            xec(1) = Xe(1,jmn)
            xec(2) = Xe(2,jmn)
            xec(3) = Xe(3,jmn)
            yec(1) = Ye(1,jmn)
            yec(2) = Ye(2,jmn)
            yec(3) = Ye(3,jmn)
            xhc(1) = Xh(1,jmn)
            xhc(2) = Xh(2,jmn)
            xhc(3) = Xh(3,jmn)
            yhc(1) = Yh(1,jmn)
            yhc(2) = Yh(2,jmn)
            yhc(3) = Yh(3,jmn)
            v1 = mixt_product (n, xec, nl)
            v2 = mixt_product (n, xhc, ml)                                                                             
            A11(i,j) = A11(i,j) + (v1 + ind_ref * v2) * fact
            v1 = mixt_product (n, yec, nl)
            v2 = mixt_product (n, yhc, ml)
            A12(i,j) = A12(i,j) + (v1 + ind_ref * v2) * fact
            v1 = mixt_product (n, xec, ml)
            v2 = mixt_product (n, xhc, nl)
            A21(i,j) = A21(i,j) + (v1 + ind_ref * v2) * fact
            v1 = mixt_product (n, yec, ml)
            v2 = mixt_product (n, yhc, nl)
            A22(i,j) = A22(i,j) + (v1 + ind_ref * v2) * fact          
          end do
        end do    
      end do
    end do 
    MaxA(1) = max(abs(A11(1,1)), abs(A11(1,2)), abs(A11(2,1)), abs(A11(2,2)),       &
                  abs(A22(1,1)), abs(A22(1,2)), abs(A22(2,1)), abs(A22(2,2)))
    MinA(1) = min(abs(A11(1,1)), abs(A11(1,2)), abs(A11(2,1)), abs(A11(2,2)),       &
                  abs(A22(1,1)), abs(A22(1,2)), abs(A22(2,1)), abs(A22(2,2)))
        
    MaxA(2) = max(abs(A12(1,1)), abs(A12(1,2)), abs(A12(2,1)), abs(A12(2,2)),       &
                  abs(A21(1,1)), abs(A21(1,2)), abs(A21(2,1)), abs(A21(2,2)))
    MinA(2) = min(abs(A12(1,1)), abs(A12(1,2)), abs(A12(2,1)), abs(A12(2,2)),       &
                  abs(A21(1,1)), abs(A21(1,2)), abs(A21(2,1)), abs(A21(2,2)))
    print "(3x,i3,5x,i3,4x,1pe14.6,1x,1pe14.6,1x,1pe14.6,1x,1pe14.6)", Nint1l,      &
            Nint2l, MaxA(1), MinA(1), MaxA(2), MinA(2)
    deallocate (paramG1, paramG2, weightsG, Nintparam)
  end do
  deallocate (Xe, Ye, Xh, Yh, m3, n3)
end subroutine  ConvergenceMatrixElemAnis       
!***********************************************************************************
subroutine matrixPARTSUB (k, ind_ref, z0, m, Nrank, Nmax, Nx, x, pond, A, na, ma)
!-----------------------------------------------------------------------------------
! The routine computes the reflection matrix for the scattering of a particle      !
! near or on a plane surface.                                                      ! 
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: m, Nrank, Nmax, Nx, na, ma
  real(O)    :: k, z0, x(Nx), pond(Nx)
  complex(O) :: ind_ref, A(2*na,2*ma)
!
  integer    :: i, j, pint, n, n1
  real(O)    :: D, nm, n1m, mr
  complex(O) :: fact, t, Rpar, Rperp, Tpar, Tperp, q, cosb, sinb, cosbM0,           &
                sinbM0, sinb0, cosb0, f, fpp, ftp, fpt, ftt
  complex(O),allocatable :: Pnm(:), dPnm(:), pinm(:), taunm(:),                     &
                            PnmM(:), dPnmM(:), pinmM(:), taunmM(:)      
!
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))
  allocate (PnmM(0:Nrank), dPnmM(0:Nrank), pinmM(0:Nrank), taunmM(0:Nrank))
  do i = 1, 2*Nmax
    do j = 1, 2*Nmax
      A(i,j) = zero
    end do
  end do                           
  q = cmplx(k * z0,0.0,O)
  q = 2._O * im * q   
  if (abs(q) < MachEps) q = cmplx(MachEps,MachEps,O)
  mr = real(m,O)
  do pint = 1, Nx
    t     = one - x(pint) / q
    fact  = 4._O * exp(q) * pond(pint) / q    
    cosb0 = t
    sinb0 = sqrt(one - t * t)
    call Leg_normalized_complex (sinb0, cosb0, abs(m), Nrank, Pnm, dPnm, pinm, taunm)                     
    call Fresnel_aer_sub (cosb0, ind_ref, Rpar, Rperp, Tpar, Tperp, cosb, sinb)     
    cosbM0 = - t
    sinbM0 = sqrt(one - t * t)
    call Leg_normalized_complex (sinbM0, cosbM0, abs(m), Nrank, PnmM, dPnmM,        &
         pinmM, taunmM)                                                 
    do i = 1, Nmax
      if (m == 0) then
        n1 = i
      else
        n1 = abs(m) + i - 1
      end if
      n1m = real(2 * n1 * (n1 + 1),O)
      n1m = 1._O / sqrt(n1m)
      do j = 1, Nmax
        if (m == 0) then
          n = j
        else
          n = abs(m) + j - 1
        end if
        nm  = real(2 * n * (n + 1),O)
        nm  = 1._O / sqrt(nm)
        D   = nm * n1m
        f   = fact * D * im**(n1 - n)
        fpp = f * pinm(n) * pinmM(n1)
        ftt = f * taunm(n) * taunmM(n1)
        fpt = f * pinm(n) * taunmM(n1)
        ftp = f * taunm(n) * pinmM(n1)        
        A(i,j) = A(i,j) + (mr**2 * fpp * Rpar + ftt * Rperp)        
        A(i+Nmax,j) = A(i+Nmax,j) + mr * (fpt * Rpar + ftp * Rperp)
        A(i,j+Nmax) = A(i,j+Nmax) + mr *(ftp * Rpar + fpt * Rperp)
        A(i+Nmax,j+Nmax) = A(i+Nmax,j+Nmax) + (mr**2 * fpp * Rperp + ftt * Rpar)
      end do
    end do
  end do           
  deallocate (Pnm, dPnm, pinm, taunm, PnmM, dPnmM, pinmM, taunmM)        
end subroutine matrixPARTSUB 

! ***********************************************************************************
subroutine matrix_Q_Biaxial2 (k, mx, my, mz, alfaP, betaP, gamaP, Mrank, Nrank,  &
           Nmax, Nbeta, Nalpha, Q31, Q11, nap, map, mesh)
!-----------------------------------------------------------------------------------
!            The routine compute the Q matrix of a biaxial particle.               ! 
!            A - Q31, B - Q11                                                      !
!-----------------------------------------------------------------------------------
  use parameters
  use surface
  implicit none
  integer     :: Mrank, Nrank, Nmax, Nbeta, Nalpha, nap, map
  real(O)     :: k, alfaP, betaP, gamaP
  complex(O)  :: mx, my, mz
  complex(O), optional  :: Q31(2*nap,2*map), Q11(2*nap,2*map)
  type(t_mesh) :: mesh 
!
  integer     :: pint
  real(O)     :: r, teta, phi, dA, n(3), x, y, z, xl, yl, zl, tetal, phil, rl, nr, nt, np, nx, ny, nz, nxl, nyl, nzl 
  complex(O)  :: zkl, fact, f, mixt_product, ind_ref
  complex(O),allocatable,save :: Xe(:,:), Ye(:,:), m1(:,:), n1(:,:), m3(:,:), n3(:,:), Xh(:,:), Yh(:,:)
  complex(O),allocatable,save :: Xe1(:,:), Ye1(:,:), Xh1(:,:), Yh1(:,:)
  complex(O),allocatable,save :: Q31_temp(:,:), Q11_temp(:,:)
  logical     :: is_calc_q31, is_calc_q11
  type(t_mesh) :: mesh_qsvwf
!
  !$OMP THREADPRIVATE(Xe, Ye, Xh, Yh, m1, n1, m3, n3, Q31_temp, Q11_temp)
  !$OMP THREADPRIVATE(Xe1, Ye1, Xh1, Yh1)

  is_calc_q31 = present(Q31)
  is_calc_q11 = present(Q11)
  
  !$OMP PARALLEL NUM_THREADS (OMP_thread_cnt)
  allocate (Xe(3,Nmax), Ye(3,Nmax), Xh(3,Nmax), Yh(3,Nmax), m1(3,Nmax), n1(3,Nmax), m3(3,Nmax), n3(3,Nmax)) 
!  allocate (Xe1(3,Nmax), Ye1(3,Nmax), Xh1(3,Nmax), Yh1(3,Nmax))
  !$OMP END PARALLEL

  if (is_calc_q31) then
    Q31(1:2*nap,1:2*map) = zero
    !$OMP PARALLEL NUM_THREADS (OMP_thread_cnt)
    allocate (Q31_temp(2*nap,2*map))
    Q31_temp(1:2*nap,1:2*map) = zero
    !$OMP END PARALLEL
  end if

  if (is_calc_q11) then
    Q11(1:2*nap,1:2*map) = zero
    !$OMP PARALLEL NUM_THREADS (OMP_thread_cnt)
    allocate (Q11_temp(2*nap,2*map))
    Q11_temp(1:2*nap,1:2*map) = zero
    !$OMP END PARALLEL
  end if
  
  ind_ref = sqrt( 0.5_O * ( mx * mx + my * my ) )                      
  f = im * k * k / Pi

  call mesh_clear(mesh_qsvwf)

  !$OMP PARALLEL DO NUM_THREADS (OMP_thread_cnt)
  do pint = 1, mesh%items_count
	call mesh_item_get(mesh, pint, teta, phi, r, n(1), n(2), n(3), dA)	
    call  MN_biaxial (k, mx, my, mz, r, teta, phi, Mrank, Nrank, Nmax, Nbeta, Nalpha, Xe, Ye, Xh, Yh, mesh_qsvwf)
!   --- compute the localized SWVF in the principal coordinate system ----
    x = r * sin(teta) * cos(phi)
    y = r * sin(teta) * sin(phi)
    z = r * cos(teta) 
    call T_cartesian_global_local (x, y, z, alfaP, betaP, gamaP, xl, yl, zl)  
    call T_cartesian_spherical (xl, yl, zl, rl, tetal, phil) 
    zkl = cmplx(k * rl,0.0,O) 
    call MN_complete (3, zkl, tetal, phil, Mrank, Nrank, Nmax,.false.,.false., m3, n3)
    call MN_complete (1, zkl, tetal, phil, Mrank, Nrank, Nmax,.false.,.false., m1, n1)
!   --- transform the unit normal vector in the principal coordinate system ---
    nr = n(1)
    nt = n(2)
    np = n(3)
    nx = sin(teta) * cos(phi)  * nr + cos(teta) * cos(phi)  * nt - sin(phi) * np 
    ny = sin(teta) * sin(phi)  * nr + cos(teta) * sin(phi)  * nt + cos(phi) * np
    nz =             cos(teta) * nr -             sin(teta) * nt 
    call T_cartesian_global_local (nx, ny, nz, alfaP, betaP, gamaP, nxl, nyl, nzl)
    n(1) =   sin(tetal) * cos(phil) * nxl + sin(tetal) * sin(phil) * nyl +      &
                 cos(tetal) * nzl
    n(2) =   cos(tetal) * cos(phil) * nxl + cos(tetal) * sin(phil) * nyl -      &
                 sin(tetal) * nzl
    n(3) =              - sin(phil) * nxl +              cos(phil) * nyl      
    fact = f * dA
    if (is_calc_q31) then
        call mixt_matrix_biax(Q31_temp, Xe, Ye, Xh, Yh, m3, n3, Nmax, nap, map, n, -fact, ind_ref)
    end if
    if (is_calc_q11) then
        call mixt_matrix_biax(Q11_temp, Xe, Ye, Xh, Yh, m1, n1, Nmax, nap, map, n,  fact, ind_ref)
    end if
  end do
  !$OMP END PARALLEL DO  

  !$OMP PARALLEL NUM_THREADS(OMP_thread_cnt)
  deallocate (Xe, Ye, Xh, Yh, m3, n3, m1, n1) 
  if (is_calc_q31) then
    !$OMP CRITICAL
    Q31(1:2*nap,1:2*map) = Q31(1:2*nap,1:2*map) + Q31_temp(1:2*nap,1:2*map)
    !$OMP END CRITICAL 
    deallocate (Q31_temp)
  end if
  if (is_calc_q11) then
    !$OMP CRITICAL
    Q11(1:2*nap,1:2*map) = Q11(1:2*nap,1:2*map) + Q11_temp(1:2*nap,1:2*map)
    !$OMP END CRITICAL 
    deallocate (Q11_temp)
  end if
  !$OMP END PARALLEL
end subroutine

! ***********************************************************************************
subroutine mixt_matrix_biax(A, Xe, Ye, Xh, Yh, m3, n3, Nmax, nap, map, n, fact, ind_ref)
use parameters
implicit none

integer, intent(in)   :: Nmax, nap, map
real(O),intent(in)    :: n(3)
complex(O),intent(in) :: fact, ind_ref
complex(O),intent(in) :: Xe(3,Nmax), Ye(3,Nmax), Xh(3,Nmax), Yh(3,Nmax), m3(3,Nmax), n3(3,Nmax)
complex(O),intent(inout) :: A(2*nap,2*map)
!
integer    :: i, j
complex(O) :: xec(3), yec(3), xhc(3), yhc(3), ml(3), nl(3), v1, v2, mixt_product

do i = 1, Nmax
    ml(1:3) = m3(1:3,i)
    nl(1:3) = n3(1:3,i)
    do j = 1, Nmax
        xec(1:3) = Xe(1:3,j)
        yec(1:3) = Ye(1:3,j)
        xhc(1:3) = Xh(1:3,j)
        yhc(1:3) = Yh(1:3,j)
        v1 = mixt_product (n, xec, nl)
        v2 = mixt_product (n, xhc, ml)
        A(i,j) = A(i,j) + (v1 + ind_ref * v2) * fact
        v1 = mixt_product (n, yec, nl)
        v2 = mixt_product (n, yhc, ml)
        A(i,j+Nmax) = A(i,j+Nmax) + (v1 + ind_ref * v2) * fact
        v1 = mixt_product (n, xec, ml)
        v2 = mixt_product (n, xhc, nl)
        A(i+Nmax,j) = A(i+Nmax,j) + (v1 + ind_ref * v2) * fact
        v1 = mixt_product (n, yec, ml)
        v2 = mixt_product (n, yhc, nl)
        A(i+Nmax,j+Nmax) = A(i+Nmax,j+Nmax) + (v1 + ind_ref * v2) * fact
    end do
end do
end subroutine

!************************************************************************************
subroutine MN_biaxial_old (k0, mx, my, mz, rBIG, tetaBIG, phiBIG, alfaP, betaP, gamaP,   &
           Mrank, Nrank, Nmax, Nbeta, Nalpha, Xe, Ye, Xh, Yh, mesh)
  use parameters
  use derived_parameters
  use surface
  implicit none
  integer    :: Mrank, Nrank, Nmax, Nbeta, Nalpha       
  real(O)    :: k0, rBIG, tetaBIG, phiBIG, alfaP, betaP, gamaP
  complex(O) :: mx, my, mz, Xe(3,Nmax), Ye(3,Nmax), Xh(3,Nmax), Yh(3,Nmax)
  type(t_mesh) :: mesh
!
  integer    :: i, k, pbeta, m, n, N0, ml, l, palpha, pint
  real(O)    :: Pi2, ct, st, cp, sp, stcp, stsp, xBIG, yBIG, zBIG, x, y, z, r, teta, &
                phi, cth, sth, cph, sph, sthcph, sthsph, cthcph, cthsph, abeta,      &
                bbeta, aalpha, balpha, beta, cb, sb, cb2, sb2, alpha, ca, sa, ca2,   &
                sa2, S(3,3), fact, nm, mlr, dalpha				    
  complex(O) :: mx2, my2, mz2, mav, imx, imy, imz, lkb, lka, lbb, lba, lab, laa,     &
                delta, l1, l2, sl1, sl2, dl, f, f2, if2, kn, k1x, k1y, k1z, k1r,     & 
                exp1, k2x, k2y, k2z, k2r, exp2, factnm, fif2, f2if2, a1,             &
                a2, a3, b1, b2, b3, Xe1, Xe2, Xe3, Ye1, Ye2, Ye3, Xh1, Xh2, Xh3,     &
                Yh1, Yh2, Yh3 						   
  real(O),allocatable :: Pnm(:), dPnm(:), pinm(:), taunm(:), xbeta(:), wbeta(:),     &
                         xalpha(:), walpha(:)
!
  Pi2  = 2._O * Pi  
! --- coordinates of the generic point in the principal coordinate system ---  
  ct   = cos(tetaBIG)
  st   = sin(tetaBIG)
  cp   = cos(phiBIG)
  sp   = sin(phiBIG)
  stcp = st * cp 
  stsp = st * sp 
  xBIG = rBIG * stcp
  yBIG = rBIG * stsp
  zBIG = rBIG * ct
  call T_cartesian_global_local (xBIG, yBIG, zBIG, alfaP, betaP, gamaP, x, y, z)
! --- polar angles of the generic point in the principal coordinate system ---
  call T_cartesian_spherical (x, y, z, r, teta, phi)  
  cth = cos(teta)
  sth = sin(teta)
  cph = cos(phi)
  sph = sin(phi)
  sthcph = sth * cph 
  sthsph = sth * sph 
  cthcph = cth * cph 
  cthsph = cth * sph  
! --- refractive indices --- 
  mx2 = mx * mx
  my2 = my * my
  mz2 = mz * mz
  mav = 0.5_O * ( mx2 + my2 )
  imx = 1._O / mx2
  imy = 1._O / my2  
  imz = 1._O / mz2 
! --- initialization ---
  do i = 1, 3
    do k = 1, Nmax
      Xe(i,k) = zero
      Ye(i,k) = zero
      Xh(i,k) = zero
      Yh(i,k) = zero
    end do
  end do  
! --- quadrature ---  
!$OMP CRITICAL
  if (mesh%items_count == 0) then
    allocate (xbeta(Nbeta), wbeta(Nbeta), xalpha(Nalpha), walpha(Nalpha))   
    abeta = 0._O
    bbeta = Pi
    call Gauss_Legendre (abeta, bbeta, Nbeta, wbeta, xbeta)                   
    aalpha = 0._O
    balpha = Pi2
    call Gauss_Legendre (aalpha, balpha, Nalpha, walpha, xalpha)
    !
    call mesh_clear(mesh)
    call mesh_expand(mesh, Nalpha*Nbeta)
    do palpha=1,Nalpha
        do pbeta=1,Nbeta
            beta  = xbeta(pbeta)           
            alpha = xalpha(palpha)
            fact  = wbeta(pbeta) * walpha(palpha) * sin(beta) / 4._O / Pi

            mesh%items_count = mesh%items_count + 1
            call mesh_item_set_default(mesh, mesh%items_count, beta, alpha, fact)
        end do
    end do
    !
    deallocate (xalpha, walpha, xbeta, wbeta)
  end if
!$OMP END CRITICAL 

! --- integration loop --- 
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank)) 
  do pint=1, mesh%items_count
    beta  = mesh%items(pint)%theta
    alpha = mesh%items(pint)%phi
    fact  = mesh%items(pint)%weight

    cb   = cos( beta )
    sb   = sin( beta )  
    cb2  = cb * cb
    sb2  = sb * sb 
    ca    = cos( alpha )
    sa    = sin( alpha )  
    ca2   = ca * ca
    sa2   = sa * sa 

!     --- biaxial quantities ---
      lkb = mav * ( imx * ca2 + imy * sa2 - imz ) * sb *cb           
      lka = mav * ( imy - imx ) * sa * ca * sb
      lbb = mav * ( ( imx * ca2 + imy * sa2 ) *cb2 + imz * sb2 )
      lba = mav * ( imy - imx ) * sa * ca * cb
      lab = lba
      laa = mav * ( imx * sa2 + imy * ca2 )   
!
      delta = sqrt( ( lbb - laa ) * ( lbb - laa ) + 4._O * lba * lba )
      l1    = 0.5_O * ( lbb + laa + delta )
      l2    = 0.5_O * ( lbb + laa - delta )
      sl1   = sqrt( l1 )
      sl2   = sqrt( l2 )               
!
      dl =  0.5_O * ( lbb - laa - delta )
      if ( abs(dl) < MachEps ) then
        if2   = zero
	    fif2  = zero
        f2if2 = one 	
      else
        f     = - lba / dl
        f2    =  f * f 	
	    if2   =  one  / ( one + f2 )
	    fif2  =  f    / ( one + f2 )
        f2if2 =  f2   / ( one + f2 )	
      end if                   
!
      kn   = sqrt( mav / l1 )
      k1x  = kn * sb * ca
      k1y  = kn * sb * sa
      k1z  = kn * cb
      k1r  = k1x * x + k1y * y + k1z * z
      exp1 = exp( im * k0 * k1r )      
!
      kn   = sqrt( mav / l2 )
      k2x  = kn * sb * ca
      k2y  = kn * sb * sa
      k2z  = kn * cb
      k2r  = k2x * x + k2y * y + k2z * z
      exp2 = exp( im * k0 * k2r )
      
!     --- rotation matrix ---
      call m_unitvct_spherical (beta, alpha, S)      
!     --- integration weight including the factor 4*Pi  ---
!      fact = wbeta(pbeta) * walpha(palpha) * sin(beta) / 4._O / Pi
!     --- compute vectors --- 
      do m = 0, Mrank         
        if (m == 0) then   
          call Leg_normalized (beta, m, Nrank, Pnm, dPnm, pinm, taunm)  
          do k = 1, Nrank
            n      = k
            nm     = real(2 * n * (n + 1),O)
            nm     = 1._O / sqrt(nm) 
	        factnm = nm * (-im)**(n + 1)
!           --- compute Xe ---	  
            a1 =( ( lkb * fif2  + lka * if2   ) * im * taunm(n)  * exp1 -            &	         
                  ( lkb * fif2  - lka * f2if2 ) * im * taunm(n)  * exp2 ) * factnm		      	             	                            
            a2 =(  l1 * fif2 * im * taunm(n)  * exp1 -                               &
	               l2 * fif2 * im * taunm(n)  * exp2 ) * factnm		     
            a3 =( l1 *  if2  * im * taunm(n)  * exp1 +                               &
	              l2 * f2if2 * im * taunm(n)  * exp2 ) * factnm
!            
            b1 = S(1,1) * a1 + S(2,1) * a2 + S(3,1) * a3
            b2 = S(1,2) * a1 + S(2,2) * a2 + S(3,2) * a3	    	    
            b3 = S(1,3) * a1 + S(2,3) * a2 + S(3,3) * a3
!         
            Xe1 = - b1
            Xe2 = - b2
            Xe3 = - b3 
            Xe(1,k) = Xe(1,k) + fact * Xe1
            Xe(2,k) = Xe(2,k) + fact * Xe2
            Xe(3,k) = Xe(3,k) + fact * Xe3  
!           --- compute Ye --- 	  
            a1 =( ( lkb * f2if2 + lka * fif2  ) * taunm(n) * exp1 +                  &
	              ( lkb * if2   - lka * fif2  ) * taunm(n) * exp2 ) * factnm		      	             	                            
            a2 =( l1 * f2if2 * taunm(n) * exp1 +                                     &
	              l2 * if2   * taunm(n) * exp2 ) * factnm		     
            a3 =( l1 * fif2 * taunm(n) * exp1 -                                      &
	              l2 * fif2 * taunm(n) * exp2 ) * factnm
!            
            b1 = S(1,1) * a1 + S(2,1) * a2 + S(3,1) * a3
            b2 = S(1,2) * a1 + S(2,2) * a2 + S(3,2) * a3	    	    
            b3 = S(1,3) * a1 + S(2,3) * a2 + S(3,3) * a3
!         
            Ye1 = - b1
            Ye2 = - b2
            Ye3 = - b3 
            Ye(1,k) = Ye(1,k) + fact * Ye1
            Ye(2,k) = Ye(2,k) + fact * Ye2
            Ye(3,k) = Ye(3,k) + fact * Ye3 
!           --- compute Xh ---
            a1 = zero
            a2 = ( sl1 * if2   * im * taunm(n) * exp1 +                               &
	               sl2 * f2if2 * im * taunm(n) * exp2 ) * factnm
            a3 = (-sl1 * fif2 * im * taunm(n)  * exp1 +                               &
	               sl2 * fif2 * im * taunm(n)  * exp2 ) * factnm
!            
            b1 = S(1,1) * a1 + S(2,1) * a2 + S(3,1) * a3
            b2 = S(1,2) * a1 + S(2,2) * a2 + S(3,2) * a3	    	    
            b3 = S(1,3) * a1 + S(2,3) * a2 + S(3,3) * a3
!         
            Xh1 = im * b1
            Xh2 = im * b2
            Xh3 = im * b3 
            Xh(1,k) = Xh(1,k) + fact * Xh1
            Xh(2,k) = Xh(2,k) + fact * Xh2
            Xh(3,k) = Xh(3,k) + fact * Xh3 	    	     	    	    	    
!           --- compute Yh ---
            a1 = zero
            a2 = ( sl1 * fif2 * taunm(n) * exp1 -                                     &
	               sl2 * fif2 * taunm(n) * exp2 ) * factnm		    
            a3 = (-sl1 * f2if2 * taunm(n) * exp1 -                                    &
	               sl2 * if2   * taunm(n) * exp2 ) * factnm
!            
            b1 = S(1,1) * a1 + S(2,1) * a2 + S(3,1) * a3
            b2 = S(1,2) * a1 + S(2,2) * a2 + S(3,2) * a3	    	    
            b3 = S(1,3) * a1 + S(2,3) * a2 + S(3,3) * a3
!         
            Yh1 = im * b1
            Yh2 = im * b2
            Yh3 = im * b3 
            Yh(1,k) = Yh(1,k) + fact * Yh1
            Yh(2,k) = Yh(2,k) + fact * Yh2
            Yh(3,k) = Yh(3,k) + fact * Yh3              	                	    	    	    	    	    	    	    	    	    	    	    
          end do
        else
          N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)         
          call Leg_normalized (beta, m, Nrank, Pnm, dPnm, pinm, taunm)
          ml = m
          do l = 1, 2
            mlr = real(ml,O)
            do k = 1, Nrank - m + 1
              n      = m + k - 1            
              nm     = real(2 * n * (n + 1),O)
              nm     = 1._O / sqrt(nm)
              factnm = nm * (-im)**(n + 1) * exp( im * mlr * alpha )
!             --- compute Xe ---	    	    	    	     	    	      
              a1 =(( ( lkb * f2if2 + lka * fif2  ) * mlr * pinm(n) +                    &
	                 ( lkb * fif2  + lka * if2   ) * im * taunm(n) ) * exp1 -           &
	               (-( lkb * if2   - lka * fif2  ) * mlr * pinm(n) +                    &
		             ( lkb * fif2  - lka * f2if2 ) * im * taunm(n) ) * exp2 ) * factnm		      	             	                            
              a2 =((  l1 * f2if2 * mlr * pinm(n) + l1 * fif2 * im * taunm(n) ) * exp1 - &
	               (- l2 * if2   * mlr * pinm(n) + l2 * fif2 * im * taunm(n) ) * exp2 ) * factnm		     
              a3 =((  l1 * fif2 * mlr * pinm(n) + l1 *  if2  * im * taunm(n) ) * exp1 + &
	               (- l2 * fif2 * mlr * pinm(n) + l2 * f2if2 * im * taunm(n) ) * exp2 ) * factnm	      
!            
              b1 = S(1,1) * a1 + S(2,1) * a2 + S(3,1) * a3
              b2 = S(1,2) * a1 + S(2,2) * a2 + S(3,2) * a3	    	    
              b3 = S(1,3) * a1 + S(2,3) * a2 + S(3,3) * a3              
!         
              Xe1 = - b1
              Xe2 = - b2
              Xe3 = - b3 
              Xe(1,k+N0) = Xe(1,k+N0) + fact * Xe1
              Xe(2,k+N0) = Xe(2,k+N0) + fact * Xe2
              Xe(3,k+N0) = Xe(3,k+N0) + fact * Xe3  
!             --- compute Ye --- 	    	    	    	     	    
              a1 =(( ( lkb * f2if2 + lka * fif2  ) * taunm(n) +                         &
	                 ( lkb * fif2  + lka * if2   ) * im * mlr * pinm(n) ) * exp1 -      &
	               (-( lkb * if2   - lka * fif2  ) * taunm(n) +                         &
		             ( lkb * fif2  - lka * f2if2 ) * im * mlr * pinm(n) ) * exp2 ) * factnm		      	             	                            
              a2 =((  l1 * f2if2 * taunm(n) + l1 * fif2 * im * mlr * pinm(n) ) * exp1 - &
	               (- l2 * if2   * taunm(n) + l2 * fif2 * im * mlr * pinm(n) ) * exp2 ) * factnm		     
              a3 =((  l1 * fif2 * taunm(n) + l1 *  if2  * im * mlr * pinm(n) ) * exp1 + &
	               (- l2 * fif2 * taunm(n) + l2 * f2if2 * im * mlr * pinm(n) ) * exp2 ) * factnm
!            
              b1 = S(1,1) * a1 + S(2,1) * a2 + S(3,1) * a3
              b2 = S(1,2) * a1 + S(2,2) * a2 + S(3,2) * a3	    	    
              b3 = S(1,3) * a1 + S(2,3) * a2 + S(3,3) * a3
!         
              Ye1 = - b1
              Ye2 = - b2
              Ye3 = - b3 
              Ye(1,k+N0) = Ye(1,k+N0) + fact * Ye1
              Ye(2,k+N0) = Ye(2,k+N0) + fact * Ye2
              Ye(3,k+N0) = Ye(3,k+N0) + fact * Ye3 
!             --- compute Xh ---	    	    	    	     	    
              a1 = zero
              a2 = ((  sl1 * fif2 * mlr * pinm(n) + sl1 * if2   * im * taunm(n) ) * exp1 + &
	                (- sl2 * fif2 * mlr * pinm(n) + sl2 * f2if2 * im * taunm(n) ) * exp2 ) * factnm
              a3 =(-(  sl1 * f2if2 * mlr * pinm(n) + sl1 * fif2 * im * taunm(n) ) * exp1 + &
	                (- sl2 * if2   * mlr * pinm(n) + sl2 * fif2 * im * taunm(n) ) * exp2 ) * factnm
!            
              b1 = S(1,1) * a1 + S(2,1) * a2 + S(3,1) * a3
              b2 = S(1,2) * a1 + S(2,2) * a2 + S(3,2) * a3	    	    
              b3 = S(1,3) * a1 + S(2,3) * a2 + S(3,3) * a3
!         
              Xh1 = im * b1
              Xh2 = im * b2
              Xh3 = im * b3 
              Xh(1,k+N0) = Xh(1,k+N0) + fact * Xh1
              Xh(2,k+N0) = Xh(2,k+N0) + fact * Xh2
              Xh(3,k+N0) = Xh(3,k+N0) + fact * Xh3 	    	     	    	    	    
!             --- compute Yh ---
              a1 = zero
              a2 = ((  sl1 * fif2 * taunm(n) + sl1 * if2   * im * mlr * pinm(n) ) * exp1 + &
	                (- sl2 * fif2 * taunm(n) + sl2 * f2if2 * im * mlr * pinm(n) ) * exp2 ) * factnm		    
              a3 =(-(  sl1 * f2if2 * taunm(n) + sl1 * fif2 * im * mlr * pinm(n) ) * exp1 + &
	                (- sl2 * if2   * taunm(n) + sl2 * fif2 * im * mlr * pinm(n) ) * exp2 ) * factnm
!            
              b1 = S(1,1) * a1 + S(2,1) * a2 + S(3,1) * a3
              b2 = S(1,2) * a1 + S(2,2) * a2 + S(3,2) * a3	    	    
              b3 = S(1,3) * a1 + S(2,3) * a2 + S(3,3) * a3
!         
              Yh1 = im * b1
              Yh2 = im * b2
              Yh3 = im * b3 
              Yh(1,k+N0) = Yh(1,k+N0) + fact * Yh1
              Yh(2,k+N0) = Yh(2,k+N0) + fact * Yh2
              Yh(3,k+N0) = Yh(3,k+N0) + fact * Yh3     	      	       	              	    	  	  			  	             
            end do 
            N0 =   N0 + Nrank - m + 1                   
            ml = - m           
          end do 
        end if 
      end do ! m
  end do ! pint

! --- vectors in polar coordinates of the principal coordinate system ---              
  do k = 1, Nmax
    Xe1 = Xe(1,k) 
    Xe2 = Xe(2,k) 
    Xe3 = Xe(3,k)  
    Xe(1,k) =   Xe1 * sthcph + Xe2 * sthsph + Xe3 * cth
    Xe(2,k) =   Xe1 * cthcph + Xe2 * cthsph - Xe3 * sth
    Xe(3,k) = - Xe1 * sph    + Xe2 * cph
!
    Ye1 = Ye(1,k) 
    Ye2 = Ye(2,k) 
    Ye3 = Ye(3,k)  
    Ye(1,k) =   Ye1 * sthcph + Ye2 * sthsph + Ye3 * cth
    Ye(2,k) =   Ye1 * cthcph + Ye2 * cthsph - Ye3 * sth
    Ye(3,k) = - Ye1 * sph    + Ye2 * cph
!
    Xh1 = Xh(1,k) 
    Xh2 = Xh(2,k) 
    Xh3 = Xh(3,k)  
    Xh(1,k) =   Xh1 * sthcph + Xh2 * sthsph + Xh3 * cth
    Xh(2,k) =   Xh1 * cthcph + Xh2 * cthsph - Xh3 * sth
    Xh(3,k) = - Xh1 * sph    + Xh2 * cph
!
    Yh1 = Yh(1,k) 
    Yh2 = Yh(2,k) 
    Yh3 = Yh(3,k)  
    Yh(1,k) =   Yh1 * sthcph + Yh2 * sthsph + Yh3 * cth
    Yh(2,k) =   Yh1 * cthcph + Yh2 * cthsph - Yh3 * sth
    Yh(3,k) = - Yh1 * sph    + Yh2 * cph            
  end do        
  deallocate (Pnm, dPnm, pinm, taunm)      
end subroutine MN_biaxial_old

!************************************************************************************
subroutine MN_biaxial (k0, mx, my, mz, rBIG, tetaBIG, phiBIG, &
           Mrank, Nrank, Nmax, Nbeta, Nalpha, Xe, Ye, Xh, Yh, mesh)
  use parameters
  use derived_parameters
  use surface
  implicit none
  integer    :: Mrank, Nrank, Nmax, Nbeta, Nalpha       
  real(O)    :: k0, rBIG, tetaBIG, phiBIG
  complex(O) :: mx, my, mz, Xe(3,Nmax), Ye(3,Nmax), Xh(3,Nmax), Yh(3,Nmax)
  type(t_mesh) :: mesh
!
  integer    :: i, k, pbeta, m, n, N0, ml, l, palpha, pint
  real(O)    :: ct, st, cp, sp, stcp, stsp, xBIG, yBIG, zBIG, x, y, z, r, teta, &
                phi, cth, sth, cph, sph, sthcph, sthsph, cthcph, cthsph, abeta,      &
                bbeta, aalpha, balpha, beta, cb, sb, cb2, sb2, alpha, ca, sa, ca2,   &
                sa2, S(3,3), fact, nm, mlr, dalpha				    
  complex(O) :: mx2, my2, mz2, mav, imx, imy, imz, &
                delta, l1, l2, sl1, sl2, dl, f, f2, if2, kn, k1x, k1y, k1z, k1r,     & 
                exp1, k2x, k2y, k2z, k2r, exp2, factnm, fif2, f2if2, a1,             &
                a2, a3, b1, b2, b3 						   
  real(O),allocatable :: Pnm(:), dPnm(:), pinm(:), taunm(:), xbeta(:), wbeta(:),     &
                         xalpha(:), walpha(:)
!
  complex(O) :: tensor_e(1:3,1:3), tensor_m, sq_tensor_m
  complex(O), dimension(1:3,1:3) :: tensor_l, tensor_e1, t_b
  complex(O) :: k1, k2, det
  real(O)    :: rot_op(1:3,1:3)
  complex(O), dimension(1:3) :: w_e1, w_e2, w_h1, w_h2
  integer, allocatable :: rel_ind(:,:)
  complex(O), allocatable :: dn(:)
  complex(O) :: a_x1, a_y1, a_x2, a_y2, sq_mav
  logical :: is_limit, is_first

! --- init variables
tensor_e(:,:) = 0._O
tensor_e(1,1) = mx*mx ! k
tensor_e(2,2) = my*my ! beta
tensor_e(3,3) = mz*mz ! alpha
tensor_m      = 1._O
sq_tensor_m   = sqrt(tensor_m)
mav    = (tensor_e(1,1)+tensor_e(2,2))/2._O
sq_mav = sqrt( mav )


! relation between index (M,N) and linear index
allocate(rel_ind(-Mrank:Mrank, Nrank))
rel_ind = 0
do m = 0, Mrank
    if (m == 0) then
        do k = 1, Nrank
            rel_ind(m,k) = k
	    end do
    else
        N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
        do k = 1, Nrank - m + 1
            n      = m + k - 1   
            rel_ind(m, n) = k+N0
        end do
        N0 =   N0 + Nrank - m + 1                   
        do k = 1, Nrank - m + 1
            n      = m + k - 1   
            rel_ind(-m, n) = k+N0
        end do
    end if
end do

! coefficions dn
allocate(dn(1:Nrank))
do n = 1, Nrank
    nm    = real(2 * n * (n + 1), O)
    dn(n) = (-im)**(n + 1) / sqrt(nm)
end do

! x, y, z
ct = cos(tetaBIG)
st = sin(tetaBIG)
cp = cos(phiBIG)
sp = sin(phiBIG)
x = rBIG * st * cp
y = rBIG * st * sp
z = rBIG * ct

Xe(1:3,1:Nmax) = zero
Ye(1:3,1:Nmax) = zero
Xh(1:3,1:Nmax) = zero
Yh(1:3,1:Nmax) = zero
! --- init variables

! --- tensor_e^{-1}
tensor_e1 = tensor_e
call identity_matrix(3, t_b, 3, 3)
call LU_parallel(tensor_e1, 3, 3, t_b, 3, 3, 3)
tensor_e1 = t_b
! --- tensor_e^{-1}

! --- generate quadrature
!$OMP CRITICAL
is_first = mesh%items_count == 0
if (mesh%items_count == 0) then
    allocate (xbeta(Nbeta), wbeta(Nbeta), xalpha(Nalpha), walpha(Nalpha))   
    abeta = 0._O
    bbeta = Pi
    call Gauss_Legendre (abeta, bbeta, Nbeta, wbeta, xbeta)                   
    aalpha = 0._O
    balpha = 2._O * Pi
    call Gauss_Legendre (aalpha, balpha, Nalpha, walpha, xalpha)
    !
    call mesh_clear(mesh)
    call mesh_expand(mesh, Nbeta*Nalpha)
    do palpha=1,Nalpha
        do pbeta=1,Nbeta
            beta  = xbeta(pbeta)           
            alpha = xalpha(palpha)
            fact  = wbeta(pbeta) * walpha(palpha) * sin(beta) / 4._O / Pi
            mesh%items_count = mesh%items_count + 1
            call mesh_item_set_default(mesh, mesh%items_count, beta, alpha, fact)
        end do
    end do
    !
    deallocate (xalpha, walpha, xbeta, wbeta)
end if
!$OMP END CRITICAL
! --- generate quadrature

! --- integration loop
allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank)) 
do pint=1, mesh%items_count
    beta  = mesh%items(pint)%theta
    alpha = mesh%items(pint)%phi
    fact  = mesh%items(pint)%weight
    
    ! rotation operatior R(alpha, beta)
    cb = cos(beta)
    sb = sin(beta)  
    ca = cos(alpha)
    sa = sin(alpha)
    rot_op(1,1:3) = (/ ca * sb,     ca * cb,    -sa /)
    rot_op(2,1:3) = (/ sa * sb,     sa * cb,     ca /)
    rot_op(3,1:3) = (/      cb,        - sb,   0._O /)
    
    ! tensor_lambda = R(alpha, beta)^(-1) * tensor_e^(-1) * R(alpha, beta)
    ! find solution [ tensor_e * R(alpha, beta) ] * tensor_lambda = [ R(alpha, beta) ]
    tensor_l = MATMUL(tensor_e1, rot_op)
    t_b = rot_op
    call transpose(t_b, 3, 3, 3)
    tensor_l = MATMUL(t_b, tensor_l) * mav
    
    ! find eigennumbers - solution det[(tensor_lambda - tensor_m * (ko/k)^2)_{1:2,1:2} ] = 0
    det = (tensor_l(2,2) - tensor_l(3,3))*(tensor_l(2,2) - tensor_l(3,3)) + 4._O * tensor_l(3,2) * tensor_l(2,3)
    det = sqrt( det )
    l1 = sqrt( 0.5_O * (tensor_l(3,3) + tensor_l(2,2) + det) )
    l2 = sqrt( 0.5_O * (tensor_l(3,3) + tensor_l(2,2) - det) )
    k1 = sq_tensor_m * sq_mav / l1 
    k2 = sq_tensor_m * sq_mav / l2
    
    ! solution eigenvectors - k1: v1=(f, 1) k2: v2=(-1, f)
    dl = 0.5_O * (tensor_l(2,2) - tensor_l(3,3) - det)
    if (abs(dl) < MachEps) then
        is_limit = .true.
        f   = 0._O
        if2 = 1._O
    else
        is_limit = .false.
        f   = - tensor_l(2,3) / dl
        if2 = 1._O / (1._O + f*f)
     end if
    
    ! exponents
    dl = cmplx(0._O, k0 * ( sb * ca * x + sb * sa * y + cb * z )) ! im * k0 * ( sb * ca * x + sb * sa * y + cb * z )
    exp1 = exp( k1 * dl)
    exp2 = exp( k2 * dl)

    ! basis vectors for E
    if (is_limit) then
        w_e1(1:3) = (/(0._O, 0._O), (-1._O, 0._O),  (0._O, 0._O)/) * mav
        w_e2(1:3) = (/(0._O, 0._O),  (0._O, 0._O), (-1._O, 0._O)/) * mav
    else
        w_e1(1:3) = (/(0._O, 0._O),  -f,  (-1._O, 0._O)/) * mav
        w_e2(1:3) = (/(0._O, 0._O), (1._O, 0._O),    -f/) * mav
    end if
    t_b = MATMUL(tensor_e1, rot_op)
    w_e1 = MATMUL(t_b, w_e1)
    w_e2 = MATMUL(t_b, w_e2)

    ! basis vectors for H
    if (is_limit) then
        w_h1(1:3) = (/(0._O, 0._O), (0._O, 0._O), (-1._O, 0._O)/) * l1 * im * sq_tensor_m
        w_h2(1:3) = (/(0._O, 0._O), (1._O, 0._O),  (0._O, 0._O)/) * l2 * im * sq_tensor_m
    else
        w_h1(1:3) = (/(0._O, 0._O), (1._O, 0._O),  -f/) * l1 * im * sq_tensor_m
        w_h2(1:3) = (/(0._O, 0._O),   f, (1._O, 0._O)/) * l2 * im * sq_tensor_m
    end if
    w_h1 = MATMUL(rot_op, w_h1)
    w_h2 = MATMUL(rot_op, w_h2)

    ! compute vectors --- 
    do m = 0, Mrank
        call Leg_normalized (beta, m, Nrank, Pnm, dPnm, pinm, taunm)  
        ml = m
        if (m == 0) then
            l = 1
        else
            l = 2
        end if
        do while (l .ge. 1)
            mlr = real(ml, O)
            if (is_limit) then
                factnm = exp( cmplx(0._O, mlr * alpha) )
            else
                factnm = exp( cmplx(0._O, mlr * alpha) ) * if2
            end if
            do n = max(1, m), Nrank
                if (is_limit) then
                    a_x1 =             mlr * pinm(n)
                    a_x2 = cmplx(0._O, taunm(n))      ! im * taunm(n)
                    a_y1 =             taunm(n)
                    a_y2 = cmplx(0._O, mlr * pinm(n)) ! im * mlr * pinm(n)
                else
                    a_x1 =  f * mlr * pinm(n) +     cmplx(0._O, taunm(n)) ! im * taunm(n)
                    a_x2 =    - mlr * pinm(n) + f * cmplx(0._O, taunm(n)) ! im * taunm(n)
                    a_y1 =       f * taunm(n) +     cmplx(0._O, mlr * pinm(n)) ! im * mlr * pinm(n)
                    a_y2 =         - taunm(n) + f * cmplx(0._O, mlr * pinm(n)) ! im * mlr * pinm(n)
                end if
            
                dl = dn(n) * factnm * fact

                a_x1 = a_x1 * exp1 * dl
                a_x2 = a_x2 * exp2 * dl
                a_y1 = a_y1 * exp1 * dl
                a_y2 = a_y2 * exp2 * dl

                k = rel_ind(ml,n)
                Xe(1:3,k) = Xe(1:3,k) + (a_x1 * w_e1(1:3) + a_x2 * w_e2(1:3))
                Ye(1:3,k) = Ye(1:3,k) + (a_y1 * w_e1(1:3) + a_y2 * w_e2(1:3))
                Xh(1:3,k) = Xh(1:3,k) + (a_x1 * w_h1(1:3) + a_x2 * w_h2(1:3))
                Yh(1:3,k) = Yh(1:3,k) + (a_y1 * w_h1(1:3) + a_y2 * w_h2(1:3))
	        end do ! n
	        ml = -ml
	        l  = l - 1
	    end do
    end do ! m
end do
deallocate (Pnm, dPnm, pinm, taunm)
! --- integration loop

! --- vectors in polar coordinates of the principal coordinate system
cth = cos(tetaBIG)
sth = sin(tetaBIG)
cph = cos(phiBIG)
sph = sin(phiBIG)
rot_op(1,1:3) = (/ cph * sth,     sph * sth,    cth /)
rot_op(2,1:3) = (/ cph * cth,     sph * cth,   -sth /)
rot_op(3,1:3) = (/      -sph,           cph,   0._O /)

do k = 1, Nmax
    Xe(1:3,k) = MATMUL(rot_op, Xe(1:3,k))
    Ye(1:3,k) = MATMUL(rot_op, Ye(1:3,k))
    Xh(1:3,k) = MATMUL(rot_op, Xh(1:3,k))
    Yh(1:3,k) = MATMUL(rot_op, Yh(1:3,k))
end do

! --- deallocate variables
deallocate (dn, rel_ind)
! --- deallocate variables

end subroutine MN_biaxial
