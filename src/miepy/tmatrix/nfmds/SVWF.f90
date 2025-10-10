! **********************************************************************************
! *          ROUTINES FOR COMPUTING THE SPHERICAL VECTOR WAVE FUNCTIONS            *
! *    -----------------------------------------------------------------------     *
! *    Partial list of subroutines:                                                *
! *      MN,                  MN_infinit,      MN_DS,          MN_poles_COMP,      *
! *      MN_poles_COMP1,      MN_DS_COMP,      MN_DS_COMP1,    MN_poles_LAY,       *
! *      MN_DS_LAY,           MN_complete,     MN_infinit_complete,                *
! *      MN_infinit_reflection_complete,       MN_left_right,  MN_anisotrop,       *
! *      MN_m_anisotrop.                                                           *
! **********************************************************************************
subroutine MN (index, z, theta, m, Nrank, Nmax, MV, NV)
!----------------------------------------------------------------------------------- 
! The routine computes the localized vector spherical wave functions (excepting    !
! the factor exp(j*m*phi)) for the azimuthal mode m.                               !
!                                                                                  !
! Input parameters:                                                                !
! - index (integer) - if index = 1, the regular vector wave functions are          !
!   computed. If index = 3, the radiating vector wave functions are computed.      !
! - z (complex,) - z = k * r, where k is the wave number and r is the module of    !
!   the position vector.                                                           !
! - theta (real) - zenith angle.                                                   !
! - m (integer) - azimuthal mode.                                                  !
! - Nrank (integer) - maximum expansion order.                                     !
! - Nmax (integer) - specifies the dimension of the vectors.                       !
!   For m = 0, Nmax = Nrank, while for m > 0, Nmax = Nrank - |m| + 1.              !
!                                                                                  !
! Output parameters:                                                               !
! - MV, NV (complex arrays) - vector spherical wave functions.                     !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: index, m, Nrank, Nmax  
  real(O)    :: theta
  complex(O) :: z, MV(3,Nmax), NV(3,Nmax)
!
  integer    :: k, n
  real(O)    :: nm, mr, nr, fp, ft, fl  
  complex(O) :: zc
  real(O),allocatable    :: Pnm(:), dPnm(:), pinm(:), taunm(:)
  complex(O),allocatable :: jh(:), jhd(:)
! 
  allocate (jh(0:Nrank), jhd(0:Nrank))
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))
  zc = z
  if (abs(zc) < MachEps) zc = cmplx(MachEps,MachEps,O) 
  if (index == 1) then
    call besel_j (zc, Nrank, jh, jhd) 
  else if (index == 3) then
    call besel_h (zc, Nrank, jh, jhd)      
  endif
  call leg_normalized (theta, abs(m), Nrank, Pnm, dPnm, pinm, taunm)
  mr = real(m,O)  
  do k = 1, Nmax
    if (m == 0) then
      n = k
    else
      n = abs(m) + k - 1
    end if
    nr = real(n * (n + 1),O)
    nm = 1._O / sqrt(2._O * nr)
    fp = mr * pinm(n) * nm
    ft = taunm(n) * nm
    fl = nr * Pnm(n) * nm
    MV(1,k) =   zero
    MV(2,k) =   jh(n)  * im * fp   
    MV(3,k) = - jh(n)  * ft         
    NV(1,k) =   jh(n)  * fl / zc
    NV(2,k) =   jhd(n) * ft / zc
    NV(3,k) =   jhd(n) * im * fp / zc           
  end do  
  deallocate (jh, jhd)
  deallocate (Pnm, dPnm, pinm, taunm)      
end subroutine MN
!***********************************************************************************
subroutine MN_infinit (theta, phi, m, Nrank, Nmax, Minf, Ninf)
!-----------------------------------------------------------------------------------
! The routine computes the localized vector spherical wave functions at infinity   !
! (excepting the factor exp(j*k*R)/(k*R)) for the azimuthal mode m.                !
!                                                                                  !
! Input parameters:                                                                !
! - theta, phi (real variables) - polar angles.                                    !
! - m (integer) - azimuthal mode.                                                  !
! - Nrank (integer) - maximum expansion order.                                     !
! - Nmax (integer) - specifies the dimension of the vectors.                       !
!   For m = 0, Nmax = Nrank, while for m > 0, Nmax = Nrank - |m| + 1.              !
!                                                                                  !
! Output parameters:                                                               !
! - Minf, Ninf (complex arrays) - vector spherical wave functions at infinity.     !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: m, Nrank, Nmax
  real(O)    :: theta, phi
  complex(O) :: Minf(3,Nmax), Ninf(3,Nmax)
!
  integer    :: k, n  
  real(O)    :: arg, nm, mr
  complex(O) :: fact, factc, factt, factp
  real(O),allocatable :: Pnm(:), dPnm(:), pinm(:), taunm(:)
!
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))
  call leg_normalized (theta, abs(m), Nrank, Pnm, dPnm, pinm, taunm)
  mr   = real(m,O)
  arg  = mr * phi
  fact = exp(im * arg)  
  do k = 1, Nmax
    if (m == 0) then
      n = k
    else
      n = abs(m) + k - 1
    end if
    nm    = real(2 * n * (n + 1),O)         
    nm    = 1._O / sqrt(nm)     
    factc = fact  * (-im)**(n + 1) * nm     
    factp = factc * mr * pinm(n)
    factt = factc * taunm(n)
    Minf(1,k) =   zero
    Minf(2,k) =   im * factp
    Minf(3,k) = - factt
    Ninf(1,k) =   zero
    Ninf(2,k) =   im * factt
    Ninf(3,k) = - factp
  end do
  deallocate (Pnm, dPnm, pinm, taunm)     
end subroutine MN_infinit
!***********************************************************************************
subroutine MN_DS (index, k, r, theta, zRe, zIm, m, Nrank, MV, NV)
!-----------------------------------------------------------------------------------
! The routine computes the distributed vector spherical wave functions (excepting  !
! the factor exp(j*m*phi)) for the azimuthal mode m.                               !
!                                                                                  ! 
! Input parameters:                                                                !
! - index (integer) - if index = 1, the regular vector wave functions are          !
!   computed. If index = 3, the radiating vector wave functions are computed.      !
! - k (complex) - wave number.                                                     !
! - r (real) - module of the position vector.                                      !
! - theta (real) - zenith angle.                                                   !
! - zRe, zIm (real arrays) - coordinates of distributed sources in the complex     !
!   plane.                                                                         !                                      
! - m (integer) - azimuthal mode.                                                  !
! - Nrank (integer) - maximum expansion order or the number of distributed         !
!   sources.                                                                       !
!                                                                                  ! 
! Output parameters:                                                               !
! - MV, NV (complex arrays) - distributed vector spherical wave functions.         !
!-----------------------------------------------------------------------------------                
  use parameters
  use derived_parameters
  implicit none
  integer    :: index, m, Nrank  
  real(O)    :: r, theta, zRe(Nrank), zIm(Nrank)
  complex(O) :: k, MV(3,Nrank), NV(3,Nrank) 
!         
  integer    :: p, n, Nbes     
  real(O)    :: ro, z, nm, mr, nr, cth, sth, dz
  complex(O) :: RR, sint, cost, argJ, sinc, cosc, Pmm, pimm, taumm, fp, ft, fl,     &
                factp, factt, factl, roc, dzc
  complex(O),allocatable :: jh(:), jhd(:)
!
  if (m == 0) then
    n = 1
  else
    n = abs(m)
  end if
  Nbes = n + 1
  mr   = real(m,O) 
  nr   = real(n * (n + 1),O)
  nm   = 1._O / sqrt(2._O * nr)   
  allocate (jh(0:Nbes), jhd(0:Nbes))
  sth = sin(theta)
  cth = cos(theta)
  ro  = r * sth
  z   = r * cth
  roc = cmplx(ro,0.0,O)
  do p = 1, Nrank
    dz  = z - zRe(p)
    dzc = cmplx(dz,0.0,O)    
    RR  = sqrt(roc * roc + (dzc - im * zIm(p))**2)
    if (abs(RR) < MachEps) RR = cmplx(MachEps,MachEps,O)
    sint = ro / RR
    cost = (dzc - im * zIm(p)) / RR
    argJ = k * RR
    if (index == 1) then
      call besel_j (argJ, Nbes, jh, jhd)
    else if (index == 3) then
      call besel_h (argJ, Nbes, jh, jhd)
    end if
    call P_mm (sint, cost, abs(m), Pmm)
    call pi_mm (sint, abs(m), pimm)
    call tau_mm (sint, cost, abs(m), taumm)
    sinc = sth * cost - cth * sint
    cosc = cth * cost + sth * sint
    fp   = im * mr * pimm * nm
    ft   = taumm * nm
    fl   = nr * Pmm * nm
    factp   =   jh(n)  * fp     
    MV(1,p) =   factp  * sinc
    MV(2,p) =   factp  * cosc
    MV(3,p) = - jh(n)  * ft             
    factl   =   jh(n)  * fl               
    factt   =   jhd(n) * ft
    NV(1,p) =   (   factl * cosc + factt * sinc) / argJ
    NV(2,p) =   ( - factl * sinc + factt * cosc) / argJ
    NV(3,p) =   jhd(n) * fp / argJ    
  end do
  deallocate (jh, jhd)      
end subroutine MN_DS
!***********************************************************************************
subroutine MN_poles_COMP (index, kk, r, theta, m, Npart, zpart, Nrankp, Nmax, MV, NV)
!-----------------------------------------------------------------------------------
! The routine computes the localized vector spherical wave functions (excepting    !
! the factor exp(j*m*phi)) for a composite particle and the azimuthal mode m.      !
!                                                                                  !
! Input parameters:                                                                !
! - index (integer) - if index = 1, the regular vector wave functions are          !
!   computed. If index = 3, the radiating vector wave functions are computed.      !
! - kk (complex) - wave number.                                                    !
! - r (real) - module of the position vector.                                      !
! - theta (real) - zenith  angle.                                                  !
! - m (integer) - azimuthal mode.                                                  !
! - Npart (integer) - number of homogeneous regions.                               !
! - zpart (real array) -  axial coordinates of the region origins.                 !
! - Nrankp (integer array) - maximum expansion orders of the regions.              !
! - Nmax (integer) - specifies the dimension of the vectors.                       !
!                                                                                  !
! Output parameters:                                                               !
! - MV, NV (complex arrays) - vector spherical wave functions.                     !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: index, m, Npart, Nmax, Nrankp(Npart)  
  real(O)    :: r, theta, zpart(Npart)
  complex(O) :: kk, MV(3,Nmax), NV(3,Nmax)
!
  integer    :: k, n, p, Nstart, Nmaxp, Nbes
  real(O)    :: r1, theta1, sinc, cosc, nm, mr, nr, fp, ft, fl, cth  
  complex(O) :: zz, nvr, nvtheta, nvphi, mvtheta
  real(O),allocatable    :: Pnm(:), dPnm(:), pinm(:), taunm(:)
  complex(O),allocatable :: jh(:), jhd(:)
!
  cth  = cos(theta)
  Nbes = 1
  do p = 1, Npart
    if (Nrankp(p) > Nbes) then
      Nbes = Nrankp(p)
    end if 
  end do
  allocate (jh(0:Nbes), jhd(0:Nbes))
  allocate (Pnm(0:Nbes), dPnm(0:Nbes), pinm(0:Nbes), taunm(0:Nbes))      
  mr = real(m,O)
  Nstart = 0
  do p = 1, Npart       
    if (Nrankp(p) >= abs(m)) then
      if (m == 0) then
        Nmaxp = Nrankp(p)
      else
        Nmaxp = Nrankp(p) - abs(m) + 1
      end if
      r1 = sqrt(r**2 + zpart(p)**2 - 2._O * r * zpart(p) * cth)
      if (r1 < MachEps) r1 = MachEps
      theta1 = acos((r * cth - zpart(p)) / r1)
      zz    = kk * r1
      if (index == 1) then
        call besel_j (zz, Nrankp(p), jh, jhd)
      else if (index == 3) then
        call besel_h (zz, Nrankp(p), jh, jhd)
      end if
      call leg_normalized (theta1, abs(m), Nrankp(p), Pnm, dPnm, pinm, taunm)
      sinc = sin(theta - theta1)
      cosc = cos(theta - theta1)
      do k = 1, Nmaxp
        if (m == 0) then
          n = k
        else
          n = abs(m) + k - 1
        end if
        nr = real(n * (n + 1),O)
        nm = 1._O / sqrt(2._O * nr)
        fp = mr * pinm(n) * nm
        ft = taunm(n) * nm
        fl = nr * Pnm(n) * nm
        mvtheta = jh(n) * im * fp 
        MV(1,Nstart+k) =    mvtheta * sinc
        MV(2,Nstart+k) =    mvtheta * cosc
        MV(3,Nstart+k) = - jh(n)   * ft
        nvr     = fl * jh(n) / zz
        nvtheta = jhd(n) * ft / zz
        nvphi   = jhd(n) * im * fp / zz
        NV(1,Nstart+k) =   nvr * cosc + nvtheta * sinc
        NV(2,Nstart+k) = - nvr * sinc + nvtheta * cosc
        NV(3,Nstart+k) =   nvphi                                                   
      end do    
      Nstart = Nstart + Nmaxp          
    endif
  end do
  deallocate (jh, jhd)
  deallocate (Pnm, dPnm, pinm, taunm)      
end subroutine MN_poles_COMP
!***********************************************************************************
subroutine MN_poles_COMP1 (index, ipart, kk, r, theta, m, Npart, zpart, Nrankp,     &
           Nmax, MV, NV)
!-----------------------------------------------------------------------------------
! The routine computes the localized vector spherical wave functions (excepting    !
! the factor exp(j*m*phi)) for the homogeneous region ipart (which is a part of    !
! the composite particle) and the azimuthal mode m.                                !
!                                                                                  ! 
! Input parameters:                                                                !
! - index (integer) - if index = 1, the regular vector wave functions are          !
!   computed. If index = 3, the radiating vector wave functions are computed.      !
! - ipart (integer) - index of the homogeneous region.                             ! 
! - kk (complex) - wave number.                                                    !
! - r (real) - module of the position vector.                                      !
! - theta (real) - zenith angle.                                                   !
! - m (integer) - azimuthal mode.                                                  !
! - Npart (integer) - number of homogeneous regions.                               !
! - zpart (real array) - axial coordinates of the region origins.                  !
! - Nrankp (integer array) - maximum expansion orders of the regions.              !
! - Nmax (integer) - specifies the dimension of the vectors.                       !
!                                                                                  !
! Output parameters:                                                               !
! - MV, NV (complex arrays) - vector spherical wave functions.                     !
!----------------------------------------------------------------------------------- 
  use parameters
  use derived_parameters
  implicit none
  integer    :: index, m, Npart, Nmax, Nrankp(Npart), ipart     
  real(O)    :: r, theta, zpart(Npart)
  complex(O) :: kk, MV(3,Nmax), NV(3,Nmax)
!
  integer    :: k, n, p, Nstart, Nmaxp, Nbes
  real(O)    :: r1, theta1, sinc, cosc, nm, mr, nr, fp, ft, fl, cth  
  complex(O) :: zz, nvr, nvtheta, nvphi, mvtheta
  real(O),allocatable    :: Pnm(:), dPnm(:), pinm(:), taunm(:)
  complex(O),allocatable :: jh(:), jhd(:)
!
  cth  = cos(theta)
  Nbes = 1
  do p = 1, Npart
    if (Nrankp(p) > Nbes) then
      Nbes = Nrankp(p)
    end if 
  end do
  allocate (jh(0:Nbes), jhd(0:Nbes))
  allocate (Pnm(0:Nbes), dPnm(0:Nbes), pinm(0:Nbes), taunm(0:Nbes))
  mr = real(m,O)
  do k = 1, 3
    do p = 1, Nmax
      MV(k,p) = zero
      NV(k,p) = zero
    end do
  end do      
  Nstart = 0
  do p = 1, Npart       
    if (Nrankp(p) >= abs(m)) then
      if (m == 0) then
        Nmaxp = Nrankp(p)
      else
        Nmaxp = Nrankp(p) - abs(m) + 1
      end if
      if (p == ipart) then
        r1 = sqrt(r**2 + zpart(p)**2 - 2._O * r * zpart(p) * cth)
        if (r1 < MachEps) r1 = MachEps
        theta1 = acos((r * cth - zpart(p)) / r1)
        zz     = kk * r1
        if (index == 1) then
          call besel_j (zz, Nrankp(p), jh, jhd)
        else if (index == 3) then
          call besel_h (zz, Nrankp(p), jh, jhd)
        end if
        call leg_normalized (theta1, abs(m), Nrankp(p), Pnm, dPnm, pinm, taunm)
        sinc = sin(theta - theta1)
        cosc = cos(theta - theta1)
        do k = 1, Nmaxp
          if (m == 0) then
            n = k
          else
            n = abs(m) + k - 1
          end if
          nr = real(n * (n + 1),O)
          nm = 1._O / sqrt(2._O * nr)
          fp = mr * pinm(n) * nm
          ft = taunm(n) * nm
          fl = nr * Pnm(n) * nm
          mvtheta = jh(n) * im * fp 
          MV(1,Nstart+k) =   mvtheta * sinc
          MV(2,Nstart+k) =   mvtheta * cosc
          MV(3,Nstart+k) = - jh(n)  * ft
          nvr     = fl * jh(n) / zz
          nvtheta = jhd(n) * ft / zz
          nvphi   = jhd(n) * im * fp / zz
          NV(1,Nstart+k) =   nvr * cosc + nvtheta * sinc
          NV(2,Nstart+k) = - nvr * sinc + nvtheta * cosc
          NV(3,Nstart+k) =   nvphi
        end do
      end if     
      Nstart = Nstart + Nmaxp          
    endif
  end do      
  deallocate (jh, jhd)
  deallocate (Pnm, dPnm, pinm, taunm)      
end subroutine MN_poles_COMP1
!***********************************************************************************
subroutine MN_DS_COMP (index, k, r, theta, Npart, Nrankpmax, Nrankp, zRe, zIm, m,   &
           Nrank, MV, NV)
!-----------------------------------------------------------------------------------
! The routine computes the distributed vector spherical wave functions (excepting  !
! the factor exp(j*m*phi)) for a composite particle and the azimuthal mode m.      !
!                                                                                  !
! Input parameters:                                                                !
! - index (integer) - if index = 1, the regular vector wave functions are          !
!   computed. If index = 3, the radiating vector wave functions are computed.      !
! - k (complex) - wave number.                                                     !
! - r (real) - module of the position vector.                                      !
! - theta (real) - zenith angle.                                                   !
! - m (integer) - azimuthal mode.                                                  !
! - Npart (integer) - number of homogeneous regions.                               !
! - zRe, zIm (real array) - coordinates of the region origins in the complex       !
!   plane.                                                                         !
! - Nrankpmax (integer) - specifies the dimension of zRe and zIm.                  !
! - Nrankp (integer array) - maximum expansion orders of the homogeneous regions.  !
! - Nrank (integer) - specifies the dimension of the vectors.                      !
!                                                                                  !
! Output parameters:                                                               !
! - MV, NV (complex arrays) - vector spherical wave functions.                     !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: index, Npart, Nrankpmax, Nrankp(Npart), m, Nrank           
  real(O)    :: r, theta, zRe(Npart,Nrankpmax), zIm(Npart,Nrankpmax)
  complex(O) :: k, MV(3,Nrank), NV(3,Nrank)
!
  integer    :: p, n, Nbes, Nstart, Nrankpl, ipart      
  real(O)    :: ro, z, nm, mr, nr, sth, cth, dz
  complex(O) :: RR, sint, cost, argJ, sinc, cosc, Pmm, pimm, taumm, fp, ft, fl,     &
                factp, factt, factl, roc, dzc
  complex(O),allocatable :: jh(:), jhd(:)
!
  if (m == 0) then
    n = 1
  else
    n = abs(m)
  end if
  Nbes = n + 1
  mr   = real(m,O)
  nr   = real(n * (n + 1),O)
  nm   = 1._O / sqrt(2._O * nr)
  allocate (jh(0:Nbes), jhd(0:Nbes))
  sth = sin(theta)
  cth = cos(theta)
  ro  = r * sth
  z   = r * cth
  roc = cmplx(ro,0.0,O)
  Nstart = 0
  do ipart = 1, Npart
    Nrankpl = Nrankp(ipart)
    do p = 1, Nrankpl
      dz  = z - zRe(ipart,p)
      dzc = cmplx(dz,0.0,O)      
      RR  = sqrt(roc * roc + (dzc - im * zIm(ipart,p))**2)
      if (abs(RR) < MachEps) RR = cmplx(MachEps,MachEps,O)
      sint = ro / RR
      cost = (dzc - im * zIm(ipart,p)) / RR
      argJ = k * RR
      if (index == 1) then
        call besel_j (argJ, Nbes, jh, jhd)
      else if (index == 3) then
        call besel_h (argJ, Nbes, jh, jhd)
      end if 
      call P_mm (sint, cost, abs(m), Pmm)
      call pi_mm (sint, abs(m), pimm)
      call tau_mm (sint, cost, abs(m), taumm)
      sinc  = sth * cost - cth * sint
      cosc  = cth * cost + sth * sint
      fp    = im * mr * pimm * nm
      ft    = taumm * nm
      fl    = nr * Pmm * nm
      factp = jh(n) * fp
      MV(1,Nstart+p) =   factp * sinc
      MV(2,Nstart+p) =   factp * cosc
      MV(3,Nstart+p) = - jh(n) * ft
      factt = jhd(n) * ft
      factl = jh(n)  * fl             
      NV(1,Nstart+p) = (  factl * cosc + factt * sinc) / argJ
      NV(2,Nstart+p) = (- factl * sinc + factt * cosc) / argJ
      NV(3,Nstart+p) = jhd(n) * fp / argJ 
    end do
    Nstart = Nstart + Nrankpl
  end do
  deallocate (jh, jhd)       
end subroutine MN_DS_COMP
!***********************************************************************************
subroutine MN_DS_COMP1 (index, ipart, k, r, theta, Npart, Nrankpmax, Nrankp, zRe,   &
           zIm, m, Nrank, MV, NV)
!-----------------------------------------------------------------------------------
! The routine computes the localized vector spherical wave functions (excepting    !
! the factor exp(j*m*phi)) for the homogeneous region ipart (which is a part of    !
! a composite paticle) and the azimuthal mode m.                                   !
!                                                                                  ! 
! Input parameters:                                                                !
! - index (integer) - if index = 1, the regular vector wave functions are          !
!   computed. If index = 3, the radiating vector wave functions are computed.      !
! - ipart (integer) - index of the homogeneous region.                             !
! - k (complex) - wave number.                                                     !
! - r (real) - module of the position vector.                                      !
! - theta (real) - zenith angle.                                                   !
! - m (integer) - azimuthal mode.                                                  !
! - Npart (integer) - number of homogeneous regions.                               !
! - zRe, zIm (real array) - coordinates of the region origins in the complex       !
!   plane.                                                                         ! 
! - Nrankpmax (integer) - specifies the dimension of zRe and zIm.                  !
! - Nrankp (integer array) - maximum expansion orders of the homogeneous regions.  ! 
! - Nrank (integer) - specifies the dimension of the vectors.                      !
!                                                                                  !
! Output parameters:                                                               !
! - MV, NV (complex arrays) - vector spherical wave functions.                     !
!----------------------------------------------------------------------------------- 
  use parameters
  use derived_parameters
  implicit none
  integer    :: index, ipart, Npart, Nrankpmax, Nrankp(Npart), m, Nrank               
  real(O)    :: r, theta, zRe(Npart,Nrankpmax), zIm(Npart,Nrankpmax)
  complex(O) :: k, MV(3,Nrank), NV(3,Nrank)
!
  integer    :: l, p, n, Nbes, Nstart, Nrankpl, iipart      
  real(O)    :: ro, z, nm, mr, nr, sth, cth, dz
  complex(O) :: RR, sint, cost, argJ, sinc, cosc, Pmm, pimm, taumm, fp, ft, fl,     &
                factp, factl, factt, roc, dzc
  complex(O),allocatable :: jh(:), jhd(:)
!
  if (m == 0) then
    n = 1
  else
    n = abs(m)
  end if
  Nbes = n + 1
  mr   = real(m,O)
  nr   = real(n * (n + 1),O)
  nm   = 1._O / sqrt(2._O * nr)  
  allocate (jh(0:Nbes), jhd(0:Nbes))
  do l = 1, 3
    do p = 1, Nrank
      MV(l,p) = zero
      NV(l,p) = zero
    end do
  end do
  sth = sin(theta)
  cth = cos(theta)
  ro  = r * sth
  z   = r * cth 
  roc = cmplx(ro,0.0,O) 
  Nstart = 0
  do iipart = 1, Npart
    Nrankpl = Nrankp(iipart)
    if (iipart == ipart) then
      do p = 1, Nrankpl
        dz  = z - zRe(ipart,p)
        dzc = cmplx(dz,0.0,O)              
        RR  = sqrt(roc * roc + (dzc - im * zIm(ipart,p))**2)
        if (abs(RR) < MachEps) RR = cmplx(MachEps,MachEps,O) 
        sint = ro / RR
        cost = (dzc - im * zIm(ipart,p)) / RR
        argJ = k * RR
        if (index == 1) then
          call besel_j (argJ, Nbes, jh, jhd)
        else if (index == 3) then
          call besel_h (argJ, Nbes, jh, jhd)
        end if         
        call P_mm (sint, cost, abs(m), Pmm)
        call pi_mm (sint, abs(m), pimm)
        call tau_mm (sint, cost, abs(m), taumm)
        sinc  = sth * cost - cth * sint
        cosc  = cth * cost + sth * sint
        fp    = im * mr * pimm * nm
        ft    = taumm * nm
        fl    = nr * Pmm * nm
        factp = jh(n) * fp
        MV(1,Nstart+p) =   factp * sinc
        MV(2,Nstart+p) =   factp * cosc
        MV(3,Nstart+p) = - jh(n) * ft
        factt = jhd(n) * ft
        factl = jh(n)  * fl               
        NV(1,Nstart+p) = (  factl * cosc + factt * sinc) / argJ
        NV(2,Nstart+p) = (- factl * sinc + factt * cosc) / argJ
        NV(3,Nstart+p) = jhd(n) * fp / argJ 
      end do
    endif
    Nstart = Nstart + Nrankpl
  end do
  deallocate (jh, jhd)      
end subroutine MN_DS_COMP1
!***********************************************************************************
subroutine MN_poles_LAY (index, kk, r, theta, m, ipart, Npart, Nrankp, Nmaxpmax,    &
           Nmaxp, zpart, MV, NV)
!-----------------------------------------------------------------------------------
! The routine computes the localized vector spherical wave functions (excepting    !
! the factor exp(j*m*phi)) for the layer ipart and the azimuthal mode m.           !
!                                                                                  !
! Input Parameters:                                                                !
! - index (integer) - if index = 1, the regular vector wave functions are          !
!   computed. If index = 3, the radiating vector wave functions are computed.      !
! - kk (complex) - wave number.                                                    !
! - r (real) - module of the position vector.                                      !
! - theta (real) - zenith angle.                                                   !
! - m (integer) - azimuthal mode.                                                  !
! - ipart (integer) - index of the layer.                                          !
! - Npart (integer) - number of layers.                                            !
! - zpart (real array) - axial coordinates of the layer origins.                   !
! - Nrankp (integer array) - maximum expansion orders of the layers.               !
! - Nmaxp (integer array) - real dimensions of the vector functions corresponding  !
!   to a specific layer.                                                           !
! - Nmaxpmax (integer) - specifies the physical dimension of the vectors.          !
!                                                                                  !
! Output parameters:                                                               !
! - MV, NV (complex arrays) - vector spherical wave functions.                     !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: index, m, ipart, Npart, Nrankp(Npart), Nmaxpmax, Nmaxp(Npart)
  real(O)    :: r, theta, zpart(Npart)
  complex(O) :: kk, MV(3,Nmaxpmax), NV(3,Nmaxpmax)
!      
  integer    :: l, k, n, p, Nbes
  real(O)    :: r1, theta1, sinc, cosc, nm, mr, nr, fp, ft, fl, cth  
  complex(O) :: zz, nvr, nvtheta, nvphi, mvtheta
  real(O),allocatable    :: Pnm(:), dPnm(:), pinm(:), taunm(:)
  complex(O),allocatable :: jh(:), jhd(:)
!
  cth = cos(theta)
  do l = 1, 3
    do p = 1, Nmaxpmax
      MV(l,p) = zero
      NV(l,p) = zero
    end do
  end do
  Nbes = Nrankp(ipart)   
  allocate (jh(0:Nbes), jhd(0:Nbes))
  allocate (Pnm(0:Nbes), dPnm(0:Nbes), pinm(0:Nbes), taunm(0:Nbes))      
  mr = real(m,O)
  r1 = sqrt(r**2 + zpart(ipart)**2 - 2._O * r * zpart(ipart) * cth)
  if (r1 < MachEps) r1 = MachEps
  theta1 = acos((r * cth - zpart(ipart)) / r1)
  zz    = kk * r1
  if (index == 1) then  
    call besel_j (zz, Nrankp(ipart), jh, jhd)
  else if (index == 3) then
    call besel_h (zz, Nrankp(ipart), jh, jhd)
  end if
  call leg_normalized (theta1, abs(m), Nrankp(ipart), Pnm, dPnm, pinm, taunm)
  sinc = sin(theta - theta1)
  cosc = cos(theta - theta1)
  do k = 1, Nmaxp(ipart)
    if (m == 0) then
      n = k
    else
      n = abs(m) + k - 1
    end if
    nr = real(n * (n + 1),O)
    nm = 1._O / sqrt(2._O * nr)
    fp = mr * pinm(n) * nm
    ft = taunm(n) * nm
    fl = nr * Pnm(n) * nm
    mvtheta = jh(n) * im * fp 
    MV(1,k) =   mvtheta * sinc
    MV(2,k) =   mvtheta * cosc
    MV(3,k) = - jh(n)  * ft
    nvr     =   fl * jh(n) / zz
    nvtheta =   jhd(n) * ft / zz
    nvphi   =   jhd(n) * im * fp / zz
    NV(1,k) =   nvr * cosc + nvtheta * sinc
    NV(2,k) = - nvr * sinc + nvtheta * cosc
    NV(3,k) =   nvphi
  end do        
  deallocate (jh, jhd)
  deallocate (Pnm, dPnm, pinm, taunm)  
end subroutine MN_poles_LAY
!***********************************************************************************
subroutine MN_DS_LAY (index, k, r, theta, ipart, Npart, Nrankpmax, Nrankp,          &
           zRe, zIm, m, MV, NV)
!-----------------------------------------------------------------------------------
! The routine computes the distributed vector spherical wave functions (excepting  !
! the factor exp(j*m*phi)) for the layer ipart and the azimuthal mode m.           !
!                                                                                  !
! Input parameters:                                                                !
! - index (integer) - if index = 1, the regular vector wave functions are          !
!   computed. If index = 3, the radiating vector wave functions are computed.      !
! - k (complex) - wave number.                                                     !
! - r (real) - module of the position vector.                                      !
! - theta (real) - zenith angle.                                                   !
! - m (integer) - azimuthal mode.                                                  !
! - ipart (integer) - index of the layer.                                          !
! - Npart (integer) - number of layers.                                            !
! - zRe, zIm (real array) - coordinates of the region origins in the complex       !
!   plane.                                                                         !
! - Nrankp (integer array) - maximum expansion orders of the layers.               !
! - Nrankpmax (integer) - specifies the physical dimension of the vectors.         !
!                                                                                  !
! Output parameters:                                                               !
! - MV, NV (complex arrays) - vector spherical wave functions.                     !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: index, ipart, Npart, Nrankpmax, Nrankp(Npart), m      
  real(O)    :: r, theta
  real(O)    :: zRe(Npart,Nrankpmax), zIm(Npart,Nrankpmax)
  complex(O) :: k, MV(3,Nrankpmax), NV(3,Nrankpmax)
!
  integer    :: p, n, Nbes, Nrankpl, l      
  real(O)    :: ro, z, nm, mr, nr, sth, cth, dz
  complex(O) :: Pmm, pimm, taumm, RR, sint, cost, argJ, sinc, cosc, fp, ft, fl,     &
                factp, factl, factt, roc, dzc
  complex(O),allocatable :: jh(:), jhd(:)
!
  if (m == 0) then
    n = 1
  else
    n = abs(m)
  end if
  Nbes = n + 1
  mr   = real(m,O)
  nr   = real(n * (n + 1),O)
  nm   = 1._O / sqrt(2._O * nr)   
  allocate (jh(0:Nbes), jhd(0:Nbes))
  do l = 1, 3
    do p = 1, Nrankpmax
      MV(l,p) = zero
      NV(l,p) = zero
    end do
  end do
  sth = sin(theta)
  cth = cos(theta)
  ro  = r * sth
  z   = r * cth  
  roc = cmplx(ro,0.0,O) 
  Nrankpl = Nrankp(ipart)
  do p = 1, Nrankpl 
    dz  = z - zRe(ipart,p)
    dzc = cmplx(dz,0.0,O)   
    RR  = sqrt(roc * roc + (dzc - im * zIm(ipart,p))**2)
    if (abs(RR) < MachEps) RR = cmplx(MachEps,MachEps,O) 
    sint = ro / RR
    cost = (dzc - im * zIm(ipart,p)) / RR
    argJ = k * RR
    if (index == 1) then
      call besel_j (argJ, Nbes, jh, jhd)
    else if (index == 3) then
      call besel_h (argJ, Nbes, jh, jhd)
    end if
    call P_mm (sint, cost, abs(m), Pmm)
    call pi_mm (sint, abs(m), pimm)
    call tau_mm (sint, cost, abs(m), taumm)
    sinc  = sth * cost - cth * sint
    cosc  = cth * cost + sth * sint
    fp    = im * mr * pimm * nm
    ft    = taumm * nm
    fl    = nr * Pmm * nm
    factp = jh(n) * fp
    MV(1,p) =   factp * sinc
    MV(2,p) =   factp * cosc
    MV(3,p) = - jh(n) * ft
    factt   = jhd(n) * ft
    factl   = jh(n) * fl                
    NV(1,p) = (  factl * cosc + factt * sinc) / argJ
    NV(2,p) = (- factl * sinc + factt * cosc) / argJ
    NV(3,p) = jhd(n) * fp / argJ 
  end do   
  deallocate (jh, jhd)
end subroutine MN_DS_LAY
!***********************************************************************************
subroutine MN_complete (index, z, theta, phi, Mrank, Nrank, Nmax, plus, sym, MV, NV)
!-----------------------------------------------------------------------------------
! The routine computes the localized vector spherical wave functions for the       !
! azimuthal modes m = 0,1,...,Mrank.                                               !
!                                                                                  !  
! Input parameters:                                                                !
! - index (integer) - if index = 1, the regular vector wave functions are          !
!   computed. If index = 3, the radiating vector wave functions are computed.      !
! - z (complex) - z = k*r, where k is the wave number and r is the module of the   !
!   position vector.                                                               !
! - theta, phi (real variables) - polar angles.                                    !
! - Mrank (integer) - number of azimuthal modes.                                   !
! - Nrank (integer) - maximum expansion order.                                     !
! - Nmax (integer) - specifies the dimension of the vectors.                       !
! - plus (logical) -  if plus = true the vector functions are computed for the     !
!   azimuthal modes m = 0,+1,-1,...,+Mrank,-Mrank. If plus = false the vector      !
!   functions are computed for the azimuthal modes m = 0,-1,+1,...,-Mrank,+Mrank.  !
! - sym (logical) - if sym = true the azimuthal dependence exp(j*m*phi) is         !
!   omitted. If sym = false, this dependence is taken into account.                !
!                                                                                  ! 
! Output parameters:                                                               !
! - MV, NV (complex arrays) - vector spherical wave functions.                     !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: index, Mrank, Nrank, Nmax  
  real(O)    :: theta, phi
  complex(O) :: z, MV(3,Nmax), NV(3,Nmax)
  logical    :: plus, sym
!
  integer    :: k, m, l, ml, N0, n  
  real(O)    :: nr, nm, mlr, ft, fl, arg
  complex(O) :: zc, fact, factp, factt, factl  
  real(O),allocatable    :: Pnm(:), dPnm(:), pinm(:), taunm(:)
  complex(O),allocatable :: jh(:), jhd(:)
!
  allocate (jh(0:Nrank), jhd(0:Nrank))
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))
  zc = z
  if (abs(zc) < MachEps) zc = cmplx(MachEps,MachEps,O) 
  if (index == 1) then
    call besel_j (zc, Nrank, jh, jhd)
  else if (index == 3) then
    call besel_h (zc, Nrank, jh, jhd)
  end if                    
  do m = 0, Mrank
    call leg_normalized (theta, m, Nrank, Pnm, dPnm, pinm, taunm)    
    if (m == 0) then      
      do k = 1, Nrank
        n  = k
        nr = real(n * (n + 1),O)
        nm = 1._O / sqrt(2._O * nr)     
        ft = taunm(n) * nm
        fl = nr * Pnm(n) * nm        
        MV(1,k) =   zero
        MV(2,k) =   zero
        MV(3,k) = - jh(n)  * ft
        NV(1,k) =   jh(n)  * fl / zc
        NV(2,k) =   jhd(n) * ft / zc
        NV(3,k) =   zero
      end do
    else    
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      ml = m
      if (.not. plus) ml = - m                                 
      do l = 1, 2
        mlr = real(ml,O)
        if (sym) then
          fact = one
        else 
          arg  = mlr * phi
          fact = exp(im * arg) 
        end if 
        do k = 1, Nrank - m + 1
          n     = m + k - 1
          nr    = real(n * (n + 1),O)
          nm    = 1._O / sqrt(2._O * nr)
          factp = im * mlr * pinm(n) * nm * fact
          factt = taunm(n) * nm * fact
          factl = nr * Pnm(n) * nm * fact
          MV(1,N0+k) =   zero
          MV(2,N0+k) =   jh(n)  * factp   
          MV(3,N0+k) = - jh(n)  * factt         
          NV(1,N0+k) =   jh(n)  * factl / zc
          NV(2,N0+k) =   jhd(n) * factt / zc
          NV(3,N0+k) =   jhd(n) * factp / zc
        end do
        N0 =   N0 + Nrank - m + 1
        ml = - ml
      end do
    end if
  end do
  deallocate (jh, jhd)
  deallocate (Pnm, dPnm, pinm, taunm)   
end subroutine MN_complete
!***********************************************************************************
subroutine MN_infinit_complete (theta, phi, Mrank, Nrank, Nmax, plus, Minf, Ninf)
!-----------------------------------------------------------------------------------
! The routine computes the localized vector spherical wave functions at infinity   !
! (excepting the factor exp(j*k*R)/(k*R)) for the azimuthal modes                  !
!  m = 0,1,...,Mrank.                                                              !
!                                                                                  !
! Input parameters:                                                                !
! - theta, phi (real variables) - polar angles.                                    !
! - Mrank (integer) - number of azimuthal modes.                                   !
! - Nrank (integer) - maximum expansion order.                                     !
! - Nmax (integer) - specifies the dimension of the vectors.                       !
! - plus (logical) - if plus = true the vector functiosn are computed for the      !
!   azimuthal modes m = 0,+1,-1,...,+Mrank,-Mrank. If plus = false, the vector     !
!   functions are computed for the azimuthal modes m = 0,-1,+1,...,-Mrank,+Mrank.  !
!                                                                                  !
! Output parameters:                                                               !
! - Minf, Ninf (complex arrays) - vector spherical wave functions.                 !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none  
  integer    :: Mrank, Nrank, Nmax
  real(O)    :: theta, phi
  complex(O) :: Minf(3,Nmax), Ninf(3,Nmax)
  logical    :: plus
!
  integer    :: k, m, l, ml, N0, n  
  real(O)    :: nm, mlr, arg
  complex(O) :: fact, factc, factt, factp
  real(O),allocatable :: Pnm(:), dPnm(:), pinm(:), taunm(:)
!
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))      
  do m = 0, Mrank
    call leg_normalized (theta, m, Nrank, Pnm, dPnm, pinm, taunm)    
    if (m == 0) then      
      do k = 1, Nrank
        n     = k
        nm    = real(2 * n * (n + 1),O)         
        nm    = 1._O / sqrt(nm)             
        factc = (-im)**(n + 1) * nm
        factt = factc * taunm(n)
        Minf(1,k) =   zero
        Minf(2,k) =   zero
        Minf(3,k) = - factt
        Ninf(1,k) =   zero
        Ninf(2,k) =   im * factt
        Ninf(3,k) =   zero
      end do
    else    
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      ml = m
      if (.not. plus) ml = - m
      do l = 1, 2
        mlr  = real(ml,O)
        arg  = mlr * phi
        fact = exp(im * arg)            
        do k = 1, Nrank - m + 1
          n     = m + k - 1
          nm    = real(2 * n * (n + 1),O)               
          nm    = 1._O / sqrt(nm)             
          factc = fact * (-im)**(n + 1) * nm
          factp = factc * mlr * pinm(n)
          factt = factc * taunm(n)
          Minf(1,N0+k) =   zero
          Minf(2,N0+k) =   im * factp
          Minf(3,N0+k) = - factt
          Ninf(1,N0+k) =   zero
          Ninf(2,N0+k) =   im * factt
          Ninf(3,N0+k) = - factp
        end do
        N0 =   N0 + Nrank - m + 1
        ml = - ml
      end do      
    end if
  end do
  deallocate (Pnm, dPnm, pinm, taunm)     
end subroutine  MN_infinit_complete
!***********************************************************************************
subroutine MN_infinit_reflection_complete (wavenumber, z0, ind_ref, theta, phi,     &
           Mrank, Nrank, Nmax, Minf, Ninf)
!-----------------------------------------------------------------------------------
! The routine computes the reflected localized vector spherical wave functions     !
! at infinity (excepting the factor exp(j*k*R)/(k*R)) and the azimuthal modes      !
! m = 0,+1,-1,...,+Mrank,-Mrank.                                                   !
!                                                                                  !
! Input parameters:                                                                !
! - wavenumber (real) - wave number in the ambient medium.                         !
! - z0 (real) - axial coordinate of the substrate.                                 !
! - ind_ref (complex) - relative refractive index of the substrate.                !
! - theta, phi (real variables) - polar angles.                                    !
! - Mrank (integer) - number of azimuthal modes.                                   !
! - Nrank (integer) - maximum expansion order.                                     !
! - Nmax (integer) - specifies the dimension of the vectors.                       !
!                                                                                  !
! Output parameters:                                                               !
! - Minf, Ninf (complex arrays) - vector spherical wave functions.                 !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax, k, m, N0, ml, l, n
  real(O)    :: theta, theta0, phi, wavenumber, z0, arga, mlr, nm
  complex(O) :: ind_ref, Rpar, Rperp, Tpar, Tperp, cost, sint, cost0, fact,         &
                facta, factr, factc, factt, factp, arg, Minf(3,Nmax), Ninf(3,Nmax)
  real(O),allocatable :: Pnm(:), dPnm(:), pinm(:), taunm(:)       
!
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))
  theta0 = Pi - theta
  cost0 = cmplx(cos(theta0),0.0,O)
  call Fresnel_aer_sub (cost0, ind_ref, Rpar, Rperp, Tpar, Tperp, cost, sint)
  arg   = 2._O * cost0 * wavenumber * z0 
  factr = exp(im * arg)       
  do m = 0, Mrank                  
    call leg_normalized (theta0, m, Nrank, Pnm, dPnm, pinm, taunm)
    if (m == 0) then
      do k = 1, Nrank
        n     = k               
        nm    = real(2 * n * (n + 1),O)         
        nm    = 1._O / sqrt(nm)             
        factc = factr * (-im)**(n + 1) * nm
        factt = factc * taunm(n)                        
        Minf(1,k) =   zero
        Minf(2,k) =   zero
        Minf(3,k) = - factt * Rperp
        Ninf(1,k) =   zero
        Ninf(2,k) =   im * factt * Rpar
        Ninf(3,k) =   zero
      end do    
    else
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      ml = m
      do l = 1, 2
        mlr   = real(ml,O)
        arga  = mlr * phi
        facta = exp(im * arga)
        fact  = facta * factr
        do k = 1, Nrank - m + 1
          n     = m + k - 1
          nm    = real(2 * n * (n + 1),O)               
          nm    = 1._O / sqrt(nm)            
          factc = fact * (-im)**(n + 1) * nm    
          factp = factc * mlr * pinm(n)
          factt = factc * taunm(n)
          Minf(1,N0+k) =   zero
          Minf(2,N0+k) =   im * factp * Rpar
          Minf(3,N0+k) = - factt * Rperp
          Ninf(1,N0+k) =   zero
          Ninf(2,N0+k) =   im * factt * Rpar
          Ninf(3,N0+k) = - factp * Rperp         
        end do
        N0 =   N0 + Nrank - m + 1
        ml = - ml       
      end do
    end if                         
  end do
  deallocate (Pnm, dPnm, pinm, taunm)      
end subroutine MN_infinit_reflection_complete
!***********************************************************************************
subroutine MN_left_right (Nmax, mvl, nvl, mvr, nvr, mv, nv)
  use parameters
  implicit none
  integer    :: Nmax, i, j
  complex(O) :: mvl(3,Nmax), nvl(3,Nmax), mvr(3,Nmax), nvr(3,Nmax), mv(3,Nmax),     &
                nv(3,Nmax)
!
  do i = 1, 3
    do j = 1, Nmax
      mv(i,j) = mvl(i,j) + nvl(i,j)
      nv(i,j) = mvr(i,j) - nvr(i,j) 
    end do
  end do
end subroutine MN_left_right
!***********************************************************************************
subroutine MN_anisotrop (k1, ind_ref, ind_refZ, rBIG, thetaBIG, phiBIG, alphaP,    &
           betaP, gammaP, Mrank, Nrank, Nmax, Nbeta, Xe, Ye, Xh, Yh)
!-----------------------------------------------------------------------------------
! The routine computes the system of vector wave functions (vector quasispherical  !
! wave functions) for uniaxial anisotropic particles.                              !
!                                                                                  !
! Input parameters:                                                                !
! - k1 (complex) - wave number.                                                    ! 
! - ind_ref, ind_refZ (complex variables) - relative refractive indices of the     !
!   uniaxial anisotropic particle.                                                 !
! - rBIG, thetaBIG, phiBIG (real avriables) - spherical coordinates of the         !
!   generic point in the particle coordinate system.                               ! 
! - alphaP, betaP, gammaP (real variables) - Euler angles specifying the           !
!   orientation of the principal coordinate system with respect to the particle    !
!   coordinate system.                                                             !
! - Mrank (integer) - number of azimuthal modes.                                   !
! - Nrank (integer) - maximum expansion order.                                     !
! - Nmax (integer) - specifies the dimension of the vectors.                       !
! - Nbeta (integer) - number of integration points.                                !
!                                                                                  !
! Output parameters:                                                               !
! - Xe, Ye, Xh, Yh (complex arrays) - vector wave functions for the uniaxial       !
!   anisotropic particle.                                                          !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax, Nbeta       
  real(O)    :: rBIG, thetaBIG, phiBIG, alphaP, betaP, gammaP
  complex(O) :: k1, ind_ref, ind_refZ, Xe(3,Nmax), Ye(3,Nmax), Xh(3,Nmax), Yh(3,Nmax)
!
  integer    :: i, k, p, m, n, N0, ml, l      
  real(O)    :: xBIG, yBIG, zBIG, x, y, z, r, theta, phi, abeta, bbeta, beta, fact, &
                Pi2, ct, st, cp, sp, cth, sth, cph, sph, cb, sb, nm, mr, mrp, mrm,  &
                arg, argp, argm, s, mlr, stcp, stsp, sthcph, sthsph, cthcph, cthsph
  complex(O) :: f1, f2, x1, x2, y1, y2, fcp, fsp, f1s, ex1, ex2, f12p, f12m, o1, o2,&
                s1, c1, s2, c2, ft, fp, fJ1, fJ2, fJ1p, fJ1m, fJ2p, fJ2m, Imx1,     &
                Imx2, ImCx1, ImCx2, ImSx1, ImSx2, Xe1, Xe2, Xe3, Ye1, Ye2, Ye3,     &
                Xh1, Xh2, Xh3, Yh1, Yh2, Yh3, ind_refZc, ind_ref2, ind_refZ2,       &
                k1rsth, k1rcth, Pi2im, Piimp, Piimm
  real(O),allocatable    :: Pnm(:), dPnm(:), pinm(:), taunm(:), xbeta(:), wbeta(:) 
  complex(O),allocatable :: J1(:), J2(:)
!
  Pi2  = 2._O * Pi  
  ct   = cos(thetaBIG)
  st   = sin(thetaBIG)
  cp   = cos(phiBIG)
  sp   = sin(phiBIG)
  stcp = st * cp 
  stsp = st * sp 
  xBIG = rBIG * stcp
  yBIG = rBIG * stsp
  zBIG = rBIG * ct
  call T_cartesian_global_local (xBIG, yBIG, zBIG, alphaP, betaP, gammaP, x, y, z)  
  call T_cartesian_spherical (x, y, z, r, theta, phi)  
  cth = cos(theta)
  sth = sin(theta)
  cph = cos(phi)
  sph = sin(phi)
  sthcph = sth * cph 
  sthsph = sth * sph 
  cthcph = cth * cph 
  cthsph = cth * sph
  fcp = Pi2 * im * cph  
  fsp = Pi2 * im * sph
  k1rsth = k1 * r * sth
  k1rcth = k1 * r * cth
  do i = 1, 3
    do k = 1, Nmax
      Xe(i,k) = zero
      Ye(i,k) = zero
      Xh(i,k) = zero
      Yh(i,k) = zero
    end do
  end do  
  allocate (xbeta(Nbeta), wbeta(Nbeta))
  abeta = 0._O
  bbeta = Pi
  call Gauss_Legendre (abeta, bbeta, Nbeta, wbeta, xbeta)                        
  allocate (J1(0:Mrank+1), J2(0:Mrank+1))
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))
  ind_refZc = ind_refZ
  if (abs(ind_refZc) < MachEps) ind_refZc = cmplx(MachEps,MachEps,O)
  ind_ref2  = ind_ref * ind_ref
  ind_refZ2 = ind_refZc * ind_refZc
  do p = 1, Nbeta
    beta = xbeta(p)
    cb   = cos(beta)
    sb   = sin(beta)
    f1   = (ind_refZ2 * cb**2 + ind_ref2 * sb**2) / ind_refZ2
    f2   = (ind_refZ2 - ind_ref2) * sb * cb / ind_refZ2
    f1s  = sqrt(f1) 
    x1   = k1rsth * sb
    x2   = x1 / f1s
    y1   = k1rcth * cb 
    y2   = y1 / f1s
    ex1  = exp(im * y1)
    ex2  = exp(im * y2)
    f12p =   f1 * cb + f2 * sb
    f12m = - f1 * sb + f2 * cb      
    call bes_J (x1, Mrank + 1, J1)
    call bes_J (x2, Mrank + 1, J2)
    fact = - wbeta(p) * sb / 4._O / Pi      
    do m = 0, Mrank         
      if (m == 0) then   
        call Leg_normalized (beta, m, Nrank, Pnm, dPnm, pinm, taunm)   
        Imx1  = Pi2 * J1(0)
        Imx2  = Pi2 * J2(0)         
        ImCx1 = fcp * J1(1)
        ImCx2 = fcp * J2(1) 
        ImSx1 = fsp * J1(1)
        ImSx2 = fsp * J2(1)
        o1 = Imx1  * ex1
        o2 = Imx2  * ex2
        s1 = ImSx1 * ex1
        c1 = ImCx1 * ex1
        s2 = ImSx2 * ex2
        c2 = ImCx2 * ex2                
        do k = 1, Nrank
          n   = k
          nm  = real(2 * n * (n + 1),O)
          nm  = 1._O / sqrt(nm)                             
          ft  = taunm(n) * nm * (-im)**(n + 1) * fact            
          f1  = im * ft                     
          Xe1 = - f1 * s1
          Xe2 =   f1 * c1
          Xe3 =   zero
          Xe(1,k) = Xe(1,k) + Xe1
          Xe(2,k) = Xe(2,k) + Xe2
          Xe(3,k) = Xe(3,k) + Xe3                           
!
          f2  = f12p * ft
          Ye1 = f2   * c2
          Ye2 = f2   * s2
          Ye3 = f12m * ft * o2      
          Ye(1,k) = Ye(1,k) + Ye1
          Ye(2,k) = Ye(2,k) + Ye2
          Ye(3,k) = Ye(3,k) + Ye3     
!             
          f1  =   ft * cb                       
          Xh1 =   f1 * c1
          Xh2 =   f1 * s1
          Xh3 = - ft * sb * o1
          Xh(1,k) = Xh(1,k) + Xh1
          Xh(2,k) = Xh(2,k) + Xh2
          Xh(3,k) = Xh(3,k) + Xh3 
!                                       
          f2  =   f1s * im * ft
          Yh1 = - f2 * s2
          Yh2 =   f2 * c2
          Yh3 =   zero
          Yh(1,k) = Yh(1,k) + Yh1
          Yh(2,k) = Yh(2,k) + Yh2
          Yh(3,k) = Yh(3,k) + Yh3                 
        end do         
      else   
        N0    = Nrank + (m - 1) * (2 * Nrank - m + 2)
        mr    = real(m,O)
        mrp   = real(m + 1,O)
        mrm   = real(m - 1,O)
        arg   = mr  * phi
        argp  = mrp * phi
        argm  = mrm * phi
        Pi2im = Pi2 * im**m     
        Piimp = Pi  * im**(m + 1)
        Piimm = Pi  * im**(m - 1)           
        call Leg_normalized (beta, m, Nrank, Pnm, dPnm, pinm, taunm)
        ml = m
        s  = 1._O               
        do l = 1, 2
          mlr   = real(ml,O)
          f1    = exp(s * im * arg)               
          fJ1   = f1 * J1(m)
          fJ2   = f1 * J2(m)
          f1    = exp(s * im * argp)
          fJ1p  = f1 * J1(m+1)
          fJ2p  = f1 * J2(m+1)
          f1    = exp(s * im * argm)
          fJ1m  = f1 * J1(m-1)
          fJ2m  = f1 * J2(m-1)
          f1    = Pi2im
          Imx1  = f1 * fJ1
          Imx2  = f1 * fJ2
          f1    = Piimp * fJ1p
          f2    = Piimm * fJ1m
          ImCx1 =   f1 + f2
          ImSx1 = - im * (s * f1 - s * f2)
          f1    = Piimp * fJ2p
          f2    = Piimm * fJ2m
          ImCx2 = f1 + f2
          ImSx2 = - im * (s * f1 - s * f2)
          o1 = Imx1  * ex1
          o2 = Imx2  * ex2
          s1 = ImSx1 * ex1
          c1 = ImCx1 * ex1
          s2 = ImSx2 * ex2
          c2 = ImCx2 * ex2
          do k = 1, Nrank - m + 1
            n   = m + k - 1            
            nm  = real(2 * n * (n + 1),O)
            nm  = 1._O / sqrt(nm)            
            fp  = mlr * pinm(n) * nm * (-im)**(n + 1) * fact
            ft  = taunm(n) * nm * (-im)**(n + 1) * fact            
            f1  = im * ft
            f2  = f12p * fp    
            Xe1 = - f1 * s1 + f2 * c2
            Xe2 =   f1 * c1 + f2 * s2
            Xe3 =   f12m * fp * o2
            Xe(1,k+N0) = Xe(1,k+N0) + Xe1
            Xe(2,k+N0) = Xe(2,k+N0) + Xe2
            Xe(3,k+N0) = Xe(3,k+N0) + Xe3
!               
            f1  =   im * fp
            f2  =   f12p * ft
            Ye1 = - f1 * s1 + f2 * c2
            Ye2 =   f1 * c1 + f2 * s2
            Ye3 =   f12m * ft * o2
            Ye(1,k+N0) = Ye(1,k+N0) + Ye1
            Ye(2,k+N0) = Ye(2,k+N0) + Ye2
            Ye(3,k+N0) = Ye(3,k+N0) + Ye3             
!              
            f1  =   ft * cb
            f2  =   f1s * im * fp
            Xh1 =   f1 * c1 - f2 * s2
            Xh2 =   f1 * s1 + f2 * c2
            Xh3 = - ft * sb * o1
            Xh(1,k+N0) = Xh(1,k+N0) + Xh1
            Xh(2,k+N0) = Xh(2,k+N0) + Xh2
            Xh(3,k+N0) = Xh(3,k+N0) + Xh3
!
            f1  =   fp * cb
            f2  =   f1s * im * ft
            Yh1 =   f1 * c1 - f2 * s2
            Yh2 =   f1 * s1 + f2 * c2
            Yh3 = - fp * sb * o1     
            Yh(1,k+N0) = Yh(1,k+N0) + Yh1
            Yh(2,k+N0) = Yh(2,k+N0) + Yh2
            Yh(3,k+N0) = Yh(3,k+N0) + Yh3
          end do  
          N0 =   N0 + Nrank - m + 1                   
          ml = - m
          s  = - s
        end do        
      end if
    end do  
  end do
! --- pass to spherical coordinates ---              
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
  deallocate (J1, J2, xbeta, wbeta)
  deallocate (Pnm, dPnm, pinm, taunm)      
end subroutine MN_anisotrop
!***********************************************************************************
subroutine MN_m_anisotrop (k1, ind_ref, ind_refZ, rBIG, thetaBIG, phiBIG, alphaP,   &
           betaP, gammaP, m, Nrank, Nmax, Nbeta, Xe, Ye, Xh, Yh)

!-----------------------------------------------------------------------------------
! The routine computes the system of vector wave functions (vector quasispherical  !
! wave functions) for uniaxial anisotropic particles and the azimuthal mode m.     !
!                                                                                  !
! Input parameters:                                                                !
! - k1 (complex) - wave number.                                                    ! 
! - ind_ref, ind_refZ (complex variables) - relative refractive indices of the     !
!   uniaxial anisotropic particle.                                                 !
! - rBIG, thetaBIG, phiBIG (real avriables) - spherical coordinates of the         !
!   generic point in the particle coordinate system.                               ! 
! - alphaP, betaP, gammaP (real variables) - Euler angles specifying the           !
!   orientation of the principal coordinate system with respect to the particle    !
!   coordinate system.                                                             !
! - m (integer) - azimuthal mode. m is postive or negative.                        !
! - Nrank (integer) - maximum expansion order.                                     !
! - Nmax (integer) - specifies the dimension of the vectors.                       !
! - Nbeta (integer) - number of integration points.                                !
!                                                                                  !
! Output parameters:                                                               !
! - Xe, Ye, Xh, Yh (complex arrays) - vector wave functions for the uniaxial       !
!   anisotropic particle.                                                          !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: m, Nrank, Nmax, Nbeta   
  real(O)    :: rBIG, thetaBIG, phiBIG, alphaP, betaP, gammaP
  complex(O) :: k1, ind_ref, ind_refZ, Xe(3,Nmax), Ye(3,Nmax), Xh(3,Nmax), Yh(3,Nmax)
!
  integer    :: i, k, p, mabs, n      
  real(O)    :: xBIG, yBIG, zBIG, x, y, z, r, theta, phi, abeta, bbeta, beta, fact, &
                Pi2, ct, st, cp, sp, cth, sth, cph, sph, cb, sb, nm, mr, mrp, mrm,  &
                arg, argp, argm, s, mlr, stcp, stsp, sthcph, sthsph, cthcph, cthsph
  complex(O) :: f1, f2, x1, x2, y1, y2, fcp, fsp, f1s, ex1, ex2, f12p, f12m, o1, o2,&
                s1, c1, s2, c2, ft, fp, fJ1, fJ2, fJ1p, fJ1m, fJ2p, fJ2m, Imx1,     &
                Imx2, ImCx1, ImCx2, ImSx1, ImSx2, Xe1, Xe2, Xe3, Ye1, Ye2, Ye3,     &
                Xh1, Xh2, Xh3, Yh1, Yh2, Yh3, ind_refZc, ind_ref2, ind_refZ2,       &
                k1rsth, k1rcth, Pi2im, Piimp, Piimm
  real(O),allocatable    :: Pnm(:), dPnm(:), pinm(:), taunm(:), xbeta(:), wbeta(:) 
  complex(O),allocatable :: J1(:), J2(:)
!
  Pi2  = 2._O * Pi  
  ct   = cos(thetaBIG)
  st   = sin(thetaBIG)
  cp   = cos(phiBIG)
  sp   = sin(phiBIG)
  stcp = st * cp 
  stsp = st * sp 
  xBIG = rBIG * stcp
  yBIG = rBIG * stsp
  zBIG = rBIG * ct
  call T_cartesian_global_local (xBIG, yBIG, zBIG, alphaP, betaP, gammaP, x, y, z)  
  call T_cartesian_spherical (x, y, z, r, theta, phi)  
  cth = cos(theta)
  sth = sin(theta)
  cph = cos(phi)
  sph = sin(phi)
  sthcph = sth * cph 
  sthsph = sth * sph 
  cthcph = cth * cph 
  cthsph = cth * sph
  fcp = Pi2 * im * cph  
  fsp = Pi2 * im * sph
  k1rsth = k1 * r * sth
  k1rcth = k1 * r * cth
  do i = 1, 3
    do k = 1, Nmax
      Xe(i,k) = zero
      Ye(i,k) = zero
      Xh(i,k) = zero
      Yh(i,k) = zero
    end do
  end do  
  allocate (xbeta(Nbeta), wbeta(Nbeta))
  abeta = 0._O
  bbeta = Pi
  call Gauss_Legendre (abeta, bbeta, Nbeta, wbeta, xbeta)                        
  mabs  = abs(m)
  mr    = real(mabs,O)
  mrp   = real(mabs + 1,O)
  mrm   = real(mabs - 1,O)
  arg   = mr  * phi
  argp  = mrp * phi
  argm  = mrm * phi
  Pi2im = Pi2 * im**mabs        
  Piimp = Pi  * im**(mabs + 1)
  Piimm = Pi  * im**(mabs - 1)  
  allocate (J1(0:mabs+1), J2(0:mabs+1))
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))
  ind_refZc = ind_refZ
  if (abs(ind_refZc) < MachEps) ind_refZc = cmplx(MachEps,MachEps,O)
  ind_ref2  = ind_ref * ind_ref
  ind_refZ2 = ind_refZc * ind_refZc
  do p = 1, Nbeta
    beta = xbeta(p)
    cb   = cos(beta)
    sb   = sin(beta)
    f1   = (ind_refZ2 * cb**2 + ind_ref2 * sb**2) / ind_refZ2
    f2   = (ind_refZ2 - ind_ref2) * sb * cb / ind_refZ2
    f1s  = sqrt(f1) 
    x1   = k1rsth * sb
    x2   = x1 / f1s
    y1   = k1rcth * cb 
    y2   = y1 / f1s
    ex1  = exp(im * y1)
    ex2  = exp(im * y2)
    f12p =   f1 * cb + f2 * sb
    f12m = - f1 * sb + f2 * cb      
    call bes_J (x1, mabs + 1, J1)
    call bes_J (x2, mabs + 1, J2)
    fact = - wbeta(p) * sb / 4._O / Pi          
    if (m == 0) then     
      call Leg_normalized (beta, 0, Nrank, Pnm, dPnm, pinm, taunm)   
      Imx1  = Pi2 * J1(0)
      Imx2  = Pi2 * J2(0)           
      ImCx1 = fcp * J1(1)
      ImCx2 = fcp * J2(1) 
      ImSx1 = fsp * J1(1)
      ImSx2 = fsp * J2(1)
      o1 = Imx1  * ex1
      o2 = Imx2  * ex2
      s1 = ImSx1 * ex1
      c1 = ImCx1 * ex1
      s2 = ImSx2 * ex2
      c2 = ImCx2 * ex2              
      do k = 1, Nmax
        n   = k
        nm  = real(2 * n * (n + 1),O)
        nm  = 1._O / sqrt(nm)                               
        ft  = taunm(n) * nm * (-im)**(n + 1) * fact            
        f1  = im * ft                       
        Xe1 = - f1 * s1
        Xe2 =   f1 * c1
        Xe3 =   zero
        Xe(1,k) = Xe(1,k) + Xe1
        Xe(2,k) = Xe(2,k) + Xe2
        Xe(3,k) = Xe(3,k) + Xe3                             
!
        f2  = f12p * ft
        Ye1 = f2   * c2
        Ye2 = f2   * s2
        Ye3 = f12m * ft * o2        
        Ye(1,k) = Ye(1,k) + Ye1
        Ye(2,k) = Ye(2,k) + Ye2
        Ye(3,k) = Ye(3,k) + Ye3       
!             
        f1  =   ft * cb                 
        Xh1 =   f1 * c1
        Xh2 =   f1 * s1
        Xh3 = - ft * sb * o1
        Xh(1,k) = Xh(1,k) + Xh1
        Xh(2,k) = Xh(2,k) + Xh2
        Xh(3,k) = Xh(3,k) + Xh3 
!                                       
        f2  =   f1s * im * ft
        Yh1 = - f2 * s2
        Yh2 =   f2 * c2
        Yh3 =   zero
        Yh(1,k) = Yh(1,k) + Yh1
        Yh(2,k) = Yh(2,k) + Yh2
        Yh(3,k) = Yh(3,k) + Yh3                   
      end do         
    else                    
      call Leg_normalized (beta, mabs, Nrank, Pnm, dPnm, pinm, taunm)     
      if (m > 0) then
        s = 1._O        
      else
        s = - 1._O
      end if                      
      mlr   = real(m,O)
      f1    = exp(s * im * arg)               
      fJ1   = f1 * J1(mabs)
      fJ2   = f1 * J2(mabs)
      f1    = exp(s * im * argp)
      fJ1p  = f1 * J1(mabs+1)
      fJ2p  = f1 * J2(mabs+1)
      f1    = exp(s * im * argm)
      fJ1m  = f1 * J1(mabs-1)
      fJ2m  = f1 * J2(mabs-1)
      f1    = Pi2im
      Imx1  = f1 * fJ1
      Imx2  = f1 * fJ2
      f1    = Piimp * fJ1p
      f2    = Piimm * fJ1m
      ImCx1 =   f1 + f2
      ImSx1 = - im * (s * f1 - s * f2)
      f1    = Piimp * fJ2p
      f2    = Piimm * fJ2m
      ImCx2 = f1 + f2
      ImSx2 = - im * (s * f1 - s * f2)
      o1 = Imx1  * ex1
      o2 = Imx2  * ex2
      s1 = ImSx1 * ex1
      c1 = ImCx1 * ex1
      s2 = ImSx2 * ex2
      c2 = ImCx2 * ex2
      do k = 1, Nmax
        n   = mabs + k - 1            
        nm  = real(2 * n * (n + 1),O)
        nm  = 1._O / sqrt(nm)            
        fp  = mlr * pinm(n) * nm * (-im)**(n + 1) * fact
        ft  = taunm(n) * nm * (-im)**(n + 1) * fact            
        f1  = im * ft
        f2  = f12p * fp    
        Xe1 = - f1 * s1 + f2 * c2
        Xe2 =   f1 * c1 + f2 * s2
        Xe3 =   f12m * fp * o2
        Xe(1,k) = Xe(1,k) + Xe1
        Xe(2,k) = Xe(2,k) + Xe2
        Xe(3,k) = Xe(3,k) + Xe3
!               
        f1  =   im * fp
        f2  =   f12p * ft
        Ye1 = - f1 * s1 + f2 * c2
        Ye2 =   f1 * c1 + f2 * s2
        Ye3 =   f12m * ft * o2
        Ye(1,k) = Ye(1,k) + Ye1
        Ye(2,k) = Ye(2,k) + Ye2
        Ye(3,k) = Ye(3,k) + Ye3               
!              
        f1  =   ft * cb
        f2  =   f1s * im * fp
        Xh1 =   f1 * c1 - f2 * s2
        Xh2 =   f1 * s1 + f2 * c2
        Xh3 = - ft * sb * o1
        Xh(1,k) = Xh(1,k) + Xh1
        Xh(2,k) = Xh(2,k) + Xh2
        Xh(3,k) = Xh(3,k) + Xh3
!
        f1  =   fp * cb
        f2  =   f1s * im * ft
        Yh1 =   f1 * c1 - f2 * s2
        Yh2 =   f1 * s1 + f2 * c2
        Yh3 = - fp * sb * o1         
        Yh(1,k) = Yh(1,k) + Yh1
        Yh(2,k) = Yh(2,k) + Yh2
        Yh(3,k) = Yh(3,k) + Yh3                                                                                                     
      end do                 
    end if      
  end do
! --- pass to spherical coordinates ---              
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
  deallocate (J1, J2, xbeta, wbeta)
  deallocate (Pnm, dPnm, pinm, taunm)      
end subroutine MN_m_anisotrop


