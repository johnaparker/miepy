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
! **********************************************************************************
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
