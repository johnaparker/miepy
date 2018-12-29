! **********************************************************************************
! *         ROUTINES FOR COMPUTING THE Q MATRIX AND THE INCIDENT MATRIX FOR        * 
! *                 AXISYMMETRIC AND NONAXISYMMETRIC PARTICLES                     *
! *    ------------------------------------------------------------------------    *
! *    Partial list of subroutines:                                                *
! *      matrix_Q,                   matrix_Q_sym,          matrix_Q_m,            *
! *      incident_matrix_m,          AzimMirSym,            ConvergenceMatrixElem  *
! **********************************************************************************
subroutine matrix_Q (FileGeom, TypeGeom, index1, index2, k, ind_ref, Nsurf, surf,   &
           rp, np, area, Nface, Mrank, Nrank, Nmax, NintAL, Nparam, Nintparam,      &
           paramG1, paramG2, weightsG, mirror, Nazimutsym, perfectcond, chiral, kb, &
           A, nap, map)
!------------------------------------------------------------------------------------
! The routine computes the Q matrix of a nonaxisymmetric particle.                  !
!------------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: Nface, Mrank, Nrank, Nmax, Nsurf, NintAL, Nparam, Nintparam(Nparam),&
                TypeGeom, index1, index2, Nazimutsym, nap, map      
  real(O)    :: k, surf(Nsurf), paramG1(Nparam,NintAL*NintAL),                      &
                paramG2(Nparam,NintAL*NintAL), weightsG(Nparam,NintAL*NintAL),      &
                rp(3,NfacePD), np(3,NfacePD), area(NfacePD), kb
  complex(O) :: ind_ref, A(2*nap,2*map)
  logical    :: FileGeom, mirror, perfectcond, chiral
!
  integer    :: NmaxL, NmaxC, i, j, pint, iparam, Nintl
  real(O)    :: r, theta, phi, dA, param1, param2, pondere, nuv(3), sign, x, y, z 
  complex(O) :: ki, zc, zl, zcl, zcr, fact, f
  complex(O),allocatable :: mv1(:,:), nv1(:,:), mv3(:,:), nv3(:,:)  
!
  NmaxL = Nmax
  NmaxC = Nmax                
  allocate (mv1(3,NmaxC), nv1(3,NmaxC), mv3(3,NmaxL), nv3(3,NmaxL))
  do i = 1, 2*NmaxL
    do j = 1, 2*NmaxC
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
        call elem_geom3D (TypeGeom, Nsurf, surf, param1, param2, iparam, r, theta,  &
             phi, dA, nuv, mirror, Nazimutsym)                       
        zc = ki * r
        zl = cmplx(k * r,0.0,O)
        if (chiral) then
          zcl = zc / (1._O - kb)
          zcr = zc / (1._O + kb)
        end if      
        call mvnv (index1, index2, chiral, zl, zc, zcl, zcr, theta, phi, Mrank,     &
             Nrank, Nmax, NmaxC, NmaxL, mv3, nv3, mv1, nv1)
        fact = f * dA * pondere
        call matQ (NmaxC, NmaxL, chiral, perfectcond, ind_ref, fact, mv3, nv3,      &
             mv1, nv1, nuv, A, nap, map)      
      end do
    end do 
  else
    do pint = 1, Nface
      x = rp(1,pint)
      y = rp(2,pint)
      z = rp(3,pint)
      call T_cartesian_spherical (x, y, z, r, theta, phi)
      dA = area(pint)                        
      nuv(1) =   sin(theta) * cos(phi) * np(1,pint) +                               &
                 sin(theta) * sin(phi) * np(2,pint) + cos(theta) * np(3,pint) 
      nuv(2) =   cos(theta) * cos(phi) * np(1,pint) +                               &
                 cos(theta) * sin(phi) * np(2,pint) - sin(theta) * np(3,pint)
      nuv(3) = - sin(phi) * np(1,pint) + cos(phi) * np(2,pint)
      zc = ki * r
      zl = cmplx(k * r,0.0,O)
      if (chiral) then
        zcl = zc / (1._O - kb)
        zcr = zc / (1._O + kb)
      end if      
      call mvnv (index1, index2, chiral, zl, zc, zcl, zcr, theta, phi, Mrank,       &
           Nrank, Nmax, NmaxC, NmaxL, mv3, nv3, mv1, nv1)
      fact = f * dA
      call matQ (NmaxC, NmaxL, chiral, perfectcond, ind_ref, fact, mv3, nv3,        &
           mv1, nv1, nuv, A, nap, map)   
    end do
  end if
  call AzimMirSym (Nazimutsym, mirror, Mrank, Nrank, NmaxC, NmaxL, A, nap, map)  
  deallocate (mv1, nv1, mv3, nv3) 
end subroutine matrix_Q 
!***********************************************************************************
subroutine matrix_Q_sym (TypeGeom, index1, index2, k, ind_ref, Nsurf, surf, Mrank,  &
           Nrank, Nmax, Nint, Nparam, Nintparam, paramG, weightsG, mirror,          &
           perfectcond, A, nap, map)
!------------------------------------------------------------------------------------
! The routine computes the Q matrix of an axisymmetric particle, for the            !
! azimuthal modes m = 0, 1,..., Mrank. Only localized sources are used for          !
! calculation.                                                                      !
!------------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax, Nsurf, Nint, Nparam, Nintparam(Nparam),         &
                TypeGeom, index1, index2, nap, map
  real(O)    :: k, surf(Nsurf), paramG(Nparam,Nint), weightsG(Nparam,Nint)      
  complex(O) :: ind_ref, A(2*nap,2*map)
  logical    :: mirror, perfectcond
!      
  integer    :: i, j, pint, m, N0, N1, iparam, Nintl, ni, nj
  real(O)    :: r, theta, phi, dA, param, pondere, nuv(3), s
  complex(O) :: ki, zc, zl, mvl(3), nvl(3), mvc(3), nvc(3), v1, v2, fact, f,        &
                mixt_product
  complex(O),allocatable :: mv1(:,:), nv1(:,:), mv3(:,:), nv3(:,:)
!
  allocate (mv1(3,Nmax), nv1(3,Nmax), mv3(3,Nmax), nv3(3,Nmax))
  do i = 1, 2*Nmax
    do j = 1, 2*Nmax
      A(i,j) = zero
    end do
  end do      
  ki = k * ind_ref
  f  = im * 2._O * k * k
  do iparam = 1, Nparam
    Nintl = Nintparam(iparam)
    do pint = 1, Nintl
      param   = paramG(iparam,pint)
      pondere = weightsG(iparam,pint)    
      call elem_geomAXSYM (TypeGeom, Nsurf, surf, param, iparam, r, theta, phi,     &
           dA, nuv)        
      zc = ki * r
      zl = cmplx(k * r,0.0,O)     
      if (index1 == 3 .and. index2 == 1) then
        call MN_complete (3, zl, theta, phi, Mrank, Nrank, Nmax,.false.,.true.,     &
             mv3, nv3)          
        call MN_complete (1, zc, theta, phi, Mrank, Nrank, Nmax,.true.,.true.,      &
             mv1, nv1)
      else if (index1 == 3 .and. index2 == 3) then
        call MN_complete (3, zl, theta, phi, Mrank, Nrank, Nmax,.false.,.true.,     &
             mv3, nv3)
        call MN_complete (3, zc, theta, phi, Mrank, Nrank, Nmax,.true.,.true.,      &
             mv1, nv1)
      else if (index1 == 1 .and. index2 == 1) then
        call MN_complete (1, zl, theta, phi, Mrank, Nrank, Nmax,.false.,.true.,     &
             mv3, nv3)          
        call MN_complete (1, zc, theta, phi, Mrank, Nrank, Nmax,.true.,.true.,      &
             mv1, nv1)        
      else if (index1 == 1 .and. index2 == 3) then
        call MN_complete (1, zl, theta, phi, Mrank, Nrank, Nmax,.false.,.true.,     &
             mv3, nv3)
        call MN_complete (3, zc, theta, phi, Mrank, Nrank, Nmax,.true.,.true.,      &
             mv1, nv1)
      end if     
      fact = f * dA * pondere
      do m = 0, Mrank
        if (m == 0) then
          do i = 1, Nrank
            ni = i
            mvl(1) = mv3(1,i)
            mvl(2) = mv3(2,i)
            mvl(3) = mv3(3,i)
            nvl(1) = nv3(1,i)
            nvl(2) = nv3(2,i)
            nvl(3) = nv3(3,i)
            do j = 1, Nrank
              nj = j
              mvc(1) = mv1(1,j)
              mvc(2) = mv1(2,j)
              mvc(3) = mv1(3,j)
              nvc(1) = nv1(1,j)
              nvc(2) = nv1(2,j)
              nvc(3) = nv1(3,j)
              s = (-1._O)**(ni + nj)
              v1 = mixt_product (nuv, mvc, nvl)
              v2 = mixt_product (nuv, nvc, mvl)
              if (.not. perfectcond) then
                A(i,j) = A(i,j) + (v1 + ind_ref * v2) * fact
                if(mirror) A(i,j) = A(i,j) + s * (v1 + ind_ref * v2) * fact
              else 
                A(i,j) = A(i,j) + v2 * fact
                if (mirror) A(i,j) = A(i,j) + s * v2 * fact
              end if                                                                                                 
              v1 = mixt_product (nuv, nvc, mvl)
              v2 = mixt_product (nuv, mvc, nvl)
              if (.not. perfectcond) then
                A(i+Nmax,j+Nmax) = A(i+Nmax,j+Nmax) + (v1 + ind_ref * v2) * fact
                if (mirror) A(i+Nmax,j+Nmax) = A(i+Nmax,j+Nmax) +                   &
                                              s * (v1 + ind_ref * v2) * fact
              else 
                A(i+Nmax,j+Nmax) = A(i+Nmax,j+Nmax) + v2*fact
                if (mirror) A(i+Nmax,j+Nmax) = A(i+Nmax,j+Nmax) + s * v2 * fact
              end if                                              
            end do
          end do
        else          
          N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
          do i = 1, Nrank - m + 1
            ni = m + i - 1
            mvl(1) = mv3(1,i+N0)
            mvl(2) = mv3(2,i+N0)
            mvl(3) = mv3(3,i+N0)
            nvl(1) = nv3(1,i+N0)
            nvl(2) = nv3(2,i+N0)
            nvl(3) = nv3(3,i+N0)
            do j = 1, Nrank - m + 1
              nj = m + j - 1
              mvc(1) = mv1(1,j+N0)
              mvc(2) = mv1(2,j+N0)
              mvc(3) = mv1(3,j+N0)
              nvc(1) = nv1(1,j+N0)
              nvc(2) = nv1(2,j+N0)
              nvc(3) = nv1(3,j+N0)
              s  = (-1._O)**(ni + nj)
              v1 = mixt_product (nuv, mvc, nvl)
              v2 = mixt_product (nuv, nvc, mvl)
              if (.not. perfectcond) then
                A(i+N0,j+N0) = A(i+N0,j+N0) + (v1 + ind_ref * v2) * fact
                if(mirror) A(i+N0,j+N0) = A(i+N0,j+N0) +                            &
                                         s * (v1 + ind_ref * v2) * fact
              else 
                A(i+N0,j+N0) = A(i+N0,j+N0) + v2 * fact
                if(mirror) A(i+N0,j+N0) = A(i+N0,j+N0) + s * v2 * fact
              end if                                                                        
              v1 = mixt_product (nuv, nvc, nvl)
              v2 = mixt_product (nuv, mvc, mvl)
              if (.not. perfectcond) then
                A(i+N0,j+N0+Nmax) = A(i+N0,j+N0+Nmax) + (v1 + ind_ref * v2) * fact
                if (mirror) A(i+N0,j+N0+Nmax) = A(i+N0,j+N0+Nmax) -                 &
                                               s * (v1 + ind_ref * v2) * fact
              else 
                A(i+N0,j+N0+Nmax) = A(i+N0,j+N0+Nmax) + v2 * fact
                if (mirror) A(i+N0,j+N0+Nmax) = A(i+N0,j+N0+Nmax) - s * v2 * fact
              end if                                              
              v1 = mixt_product (nuv, mvc, mvl)
              v2 = mixt_product (nuv, nvc, nvl)
              if (.not. perfectcond) then 
                A(i+N0+Nmax,j+N0) = A(i+N0+Nmax,j+N0) + (v1 + ind_ref * v2) * fact
                if (mirror) A(i+N0+Nmax,j+N0) = A(i+N0+Nmax,j+N0) -                 &
                                               s * (v1 + ind_ref * v2) * fact
              else 
                A(i+N0+Nmax,j+N0) = A(i+N0+Nmax,j+N0) + v2 * fact
                if (mirror) A(i+N0+Nmax,j+N0) = A(i+N0+Nmax,j+N0) - s * v2 * fact
              end if                                                              
              v1 = mixt_product (nuv, nvc, mvl)
              v2 = mixt_product (nuv, mvc, nvl)
              if (.not. perfectcond) then 
                A(i+N0+Nmax,j+N0+Nmax) = A(i+N0+Nmax,j+N0+Nmax) +                   &
                                        (v1 + ind_ref * v2) * fact
                if (mirror) A(i+N0+Nmax,j+N0+Nmax) = A(i+N0+Nmax,j+N0+Nmax) +       &
                                                    s * (v1+ind_ref * v2) * fact
              else 
                A(i+N0+Nmax,j+N0+Nmax) = A(i+N0+Nmax,j+N0+Nmax) + v2 * fact
                if (mirror) A(i+N0+Nmax,j+N0+Nmax) = A(i+N0+Nmax,j+N0+Nmax) +       &
                                                    s * v2 * fact
              end if                                              
            end do
          end do
          N1 = N0 + Nrank - m + 1
          do i = 1, Nrank - m + 1
            do j = 1, Nrank - m + 1
              A(i+N1,j+N1) = A(i+N0,j+N0)
              A(i+N1,j+N1+Nmax) = - A(i+N0,j+N0+Nmax)
              A(i+N1+Nmax,j+N1) = - A(i+N0+Nmax,j+N0)
              A(i+N1+Nmax,j+N1+Nmax) = A(i+N0+Nmax,j+N0+Nmax)
            end do
          end do
        end if
      end do 
    end do
  end do
  deallocate (mv1, nv1, mv3, nv3)       
end subroutine matrix_Q_sym 
!***********************************************************************************
subroutine matrix_Q_m (FileGeom, TypeGeom, index1, index2, k, ind_ref, Nsurf, surf, &
           rp, np, area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam,  &
           paramG, weightsG, mirror, perfectcond, DS, chiral, kb, A, nap, map)
!------------------------------------------------------------------------------------
! The routine computes the Q matrix of an axisymmetric particle, and the            !
! azimuthal mode m                                                                  !
!------------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: m, Nrank, Nmax, Nsurf, Nface, Nint, Nparam, Nintparam(Nparam),      &
                TypeGeom, index1, index2, nap, map          
  real(O)    :: k, kb, surf(Nsurf), paramG(Nparam,Nint), weightsG(Nparam,Nint),     &
                zRe(Nrank), zIm(Nrank), rp(2,NfacePD), np(2,NfacePD), area(NfacePD)
  complex(O) :: ind_ref, A(2*nap,2*map)
  logical    :: FileGeom, mirror, perfectcond, DS, chiral
!      
  integer    :: NmaxL, NmaxC, i, j, pint, iparam, Nintl
  real(O)    :: r, theta, phi, dA, param, pondere, nuv(3), sign, x, z, za, thetaA
  complex(O) :: kc, ki, kil, kir, zc, zl, zcl, zcr, fact, f
  complex(O),allocatable :: mv1(:,:), nv1(:,:), mv3(:,:), nv3(:,:)
!
  NmaxL = Nmax
  NmaxC = Nmax
  if (DS) then
    NmaxC = Nrank 
    if (index1 == 3) then
      NmaxL = Nrank   
    else if (index1 == 1) then
      NmaxL = Nmax    
    end if
  end if  
  allocate (mv1(3,NmaxC), nv1(3,NmaxC), mv3(3,NmaxL), nv3(3,NmaxL))
  kc = cmplx(k,0.0,O)
  do i = 1, 2*NmaxL
    do j = 1, 2*NmaxC
      A(i,j) = zero
    end do
  end do      
  ki = k * ind_ref
  if (chiral) then
    kil = ki / (1._O - kb)
    kir = ki / (1._O + kb)
  end if
  sign = 1._O
  if (index1 == 3 .and. index2 == 1) sign = - 1._O
  f = sign * im * 2._O * k * k 
  if (.not. FileGeom) then
    do iparam = 1, Nparam
      Nintl = Nintparam(iparam)
      do pint = 1, Nintl
        param   = paramG(iparam,pint)
        pondere = weightsG(iparam,pint)   
        call elem_geomAXSYM (TypeGeom, Nsurf, surf, param, iparam, r, theta, phi,   &
             dA, nuv)  
        zc = ki * r
        zl = cmplx(k * r,0.0,O)
        if (chiral) then
          zcl = zc / (1._O - kb)
          zcr = zc / (1._O + kb)
        end if
        call mvnv_m (index1, index2, chiral, DS, zl, zc, zcl, zcr, kc, ki, kil, kir,&
             r, theta, m, Nrank, Nmax, NmaxC, NmaxL, zRe, zIm, mv3, nv3, mv1, nv1)                   
        fact = f * dA * pondere 
        call matQ_m (m, NmaxC, NmaxL, chiral, perfectcond, mirror, ind_ref, fact,   &
             mv3, nv3, mv1, nv1, nuv, A, nap, map)       
      end do
    end do
  else
    do pint = 1, Nface
      x  = rp(1,pint)     
      z  = rp(2,pint)
      r  = sqrt(x * x + z * z)      
      za = abs(z)
      if (za < MachEps) then
        theta = Pi / 2._O
      else
        thetaA = atan(x / za)   ! x > 0
        if (z >= MachEps) then
          theta = thetaA
        else
          theta = Pi - thetaA
        end if
      end if 
      dA = area(pint) 
      nuv(1) =   np(2,pint) * cos(theta) + np(1,pint) * sin(theta)
      nuv(2) = - np(2,pint) * sin(theta) + np(1,pint) * cos(theta)
      nuv(3) = 0._O
      zc = ki * r
      zl = cmplx(k * r,0.0,O)
      if (chiral) then
        zcl = zc / (1._O - kb)
        zcr = zc / (1._O + kb)
      end if
      call mvnv_m (index1, index2, chiral, DS, zl, zc, zcl, zcr, kc, ki, kil, kir,  &
           r, theta, m, Nrank, Nmax, NmaxC, NmaxL, zRe, zIm, mv3, nv3, mv1, nv1)             
      fact = f * dA  
      call matQ_m (m, NmaxC, NmaxL, chiral, perfectcond, mirror, ind_ref, fact,     &
           mv3, nv3, mv1, nv1, nuv, A, nap, map)            
    end do  
  end if            
  deallocate (mv1, nv1, mv3, nv3)       
end subroutine matrix_Q_m
!***********************************************************************************
subroutine incident_matrix_m (FileGeom, TypeGeom, k, Nsurf, surf, rp, np, area,     &
           Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,        &
           weightsG, A, nap, map)
!-----------------------------------------------------------------------------------
! The routine computes the incident matrix of an axisymmetric particle             !
! using distributed sources.                                                       !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: TypeGeom, m, Nrank, Nmax, Nint, Nsurf, Nparam, Nintparam(Nparam),   &
                nap, map, Nface
  real(O)    :: k, surf(Nsurf), zRe(Nrank), zIm(Nrank), paramG(Nparam,Nint),        &
                weightsG(Nparam,Nint), rp(2,NfacePD), np(2,NfacePD), area(NfacePD)
  complex(O) :: A(2*nap,2*map)
  logical    :: FileGeom  
!    
  integer    :: i, j, pint, iparam, Nintl
  real(O)    :: r, theta, phi, dA, param, pondere, nuv(3), x, z, za, thetaA
  complex(O) :: kc, zl, fact, f
  complex(O),allocatable :: mv1(:,:), nv1(:,:), mv3(:,:), nv3(:,:)
!
  allocate (mv1(3,Nmax), nv1(3,Nmax), mv3(3,Nrank), nv3(3,Nrank))
  kc = cmplx(k,0.0,O)
  do i = 1, 2*Nrank
    do j = 1, 2*Nmax
      A(i,j) = zero
    end do
  end do
  f = - im * 2._O * k * k
  if (.not. FileGeom) then
    do iparam = 1, Nparam
      Nintl = Nintparam(iparam)
      do pint = 1, Nintl
        param   = paramG(iparam,pint)
        pondere = weightsG(iparam,pint)
        call elem_geomAXSYM (TypeGeom, Nsurf, surf, param, iparam, r, theta, phi,   &
             dA, nuv)        
        zl = cmplx(k * r,0.0,O)
        call mvnvinc_m (zl, kc, r, theta, m, Nrank, Nmax, zRe, zIm, mv3, nv3,       &
             mv1, nv1)          
        fact = f * dA * pondere
        call matQinc_m (m, Nmax, Nrank, fact, mv3, nv3, mv1, nv1, nuv, A, nap, map)       
      end do
    end do
  else
    do pint = 1, Nface
      x  = rp(1,pint)     
      z  = rp(2,pint)
      r  = sqrt(x * x + z * z)      
      za = abs(z)
      if (za < MachEps) then
        theta = Pi / 2._O
      else
        thetaA = atan(x / za)   ! x > 0
        if (z >= MachEps) then
          theta = thetaA
        else
          theta = Pi - thetaA
        end if
      end if 
      dA = area(pint) 
      nuv(1) =   np(2,pint) * cos(theta) + np(1,pint) * sin(theta)
      nuv(2) = - np(2,pint) * sin(theta) + np(1,pint) * cos(theta)
      nuv(3) = 0._O
      zl = cmplx(k * r,0.0,O)      
      call mvnvinc_m (zl, kc, r, theta, m, Nrank, Nmax, zRe, zIm, mv3, nv3,         &
           mv1, nv1)          
      fact = f * dA 
      call matQinc_m (m, Nmax, Nrank, fact, mv3, nv3, mv1, nv1, nuv, A, nap, map) 
    end do
  end if            
  deallocate (mv1, nv1, mv3, nv3)     
end subroutine incident_matrix_m
!***********************************************************************************
subroutine AzimMirSym (Nazimutsym, mirror, Mrank, Nrank, NmaxC, NmaxL, A, nap, map)
!-----------------------------------------------------------------------------------
! The routine completes the calculation of the Q matrix of a nonaxisymmetric       !
! particle taking into account the mirror and azimuthal symmetry.                  !
! Note: The real dimensions of A are (2*NmaxL,2*NmaxC). The rows are characterized !
! by the azimuthal indices                                                         !
!                             m = 0,-1,+1,-2,+2,...,                               !
! while the columns are characterized by zhe azimuthal indices                     ! 
!                             m'= 0,+1,-1,+2,-2,...                                !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: Nazimutsym, Mrank, Nrank, NmaxC, NmaxL, nap, map  
  complex(O) :: A(2*nap,2*map)
  logical    :: mirror
!
  integer    :: m, ml, n, m1, m1l, n1, i, j, p, Nc, Nl, l, l1
  real(O)    :: dalfa, arg, sp, sm 
  complex(O) :: fact
!
  if (Nazimutsym >= 2 .or. mirror) then
    if (Nazimutsym >= 2) dalfa = 2._O * Pi / Nazimutsym
    do m = 0, Mrank
      if (m == 0) then
        ml = m
        do i = 1, Nrank
          n = i
          do m1 = 0, Mrank
            if (m1 == 0) then
              m1l = m1                    
              if (mirror) then                         
                do j = 1, Nrank
                  n1 = j
                  sp = 1._O + (-1._O)**(n + n1 + m + m1)
                  sm = 2._O - sp
                  A(i,j) = sp * A(i,j)
                  A(i,j+NmaxC) = sm * A(i,j+NmaxC)
                  A(i+NmaxL,j) = sm * A(i+NmaxL,j)
                  A(i+NmaxL,j+NmaxC) = sp * A(i+NmaxL,j+NmaxC)
                end do
              end if 
              if (Nazimutsym >= 2) then
                fact = one
                do p = 2, Nazimutsym
                  arg  = real((ml + m1l) * (p - 1),O) * dalfa
                  fact = fact + exp(im * arg)
                end do
                do j = 1, Nrank                    
                  A(i,j) = fact * A(i,j)
                  A(i,j+NmaxC) = fact * A(i,j+NmaxC)
                  A(i+NmaxL,j) = fact * A(i+NmaxL,j)
                  A(i+NmaxL,j+NmaxC) = fact * A(i+NmaxL,j+NmaxC)
                end do
              end if
            else                          
              m1l = m1                    
              Nc  = Nrank + (m1 - 1) * (2*Nrank - m1 + 2)
              do l1 = 1, 2
                if (mirror) then
                  do j = 1, Nrank - m1 + 1
                    n1 = m1 + j - 1 
                    sp = 1._O + (-1._O)**(n + n1 + m + m1)
                    sm = 2._O - sp             
                    A(i,j+Nc) = sp * A(i,j+Nc) 
                    A(i,j+Nc+NmaxC) = sm * A(i,j+Nc+NmaxC)
                    A(i+NmaxL,j+Nc) = sm * A(i+NmaxL,j+Nc)
                    A(i+NmaxL,j+Nc+NmaxC) = sp * A(i+NmaxL,j+Nc+NmaxC)
                  end do
                end if
                if (Nazimutsym >= 2) then
                  fact = one
                  do p = 2, Nazimutsym
                    arg  = real((ml + m1l) * (p - 1),O) * dalfa
                    fact = fact + exp(im * arg)
                  end do                                  
                  do j = 1, Nrank - m1 + 1                                     
                    A(i,j+Nc) = fact * A(i,j+Nc) 
                    A(i,j+Nc+NmaxC) = fact * A(i,j+Nc+NmaxC)
                    A(i+NmaxL,j+Nc) = fact * A(i+NmaxL,j+Nc)
                    A(i+NmaxL,j+Nc+NmaxC) = fact * A(i+NmaxL,j+Nc+NmaxC)
                  end do
                end if
                m1l = - m1l
                Nc  =   Nc + Nrank - m1 + 1
              end do                      
            end if
          end do
        end do
      else
        ml = - m 
        Nl =   Nrank + (m - 1) * (2 * Nrank - m + 2)
        do l = 1, 2
          do i = 1, Nrank - m + 1
            n = m + i - 1
            do m1 = 0, Mrank
              if (m1 == 0) then
                m1l = m1
                if (mirror) then
                  do j = 1, Nrank                                   
                    n1 = j 
                    sp = 1._O + (-1._O)**(n + n1 + m + m1)
                    sm = 2._O - sp  
                    A(i+Nl,j) = sp * A(i+Nl,j)
                    A(i+Nl,j+NmaxC) = sm * A(i+Nl,j+NmaxC)
                    A(i+Nl+NmaxL,j) = sm * A(i+Nl+NmaxL,j)
                    A(i+Nl+NmaxL,j+NmaxC) = sp * A(i+Nl+NmaxL,j+NmaxC)
                  end do
                end if
                if (Nazimutsym >= 2) then
                  fact = one
                  do p = 2, Nazimutsym
                    arg  = real((ml + m1l) * (p - 1),O) * dalfa
                    fact = fact + exp(im * arg)
                  end do
                  do j = 1, Nrank                            
                    A(i+Nl,j) = fact * A(i+Nl,j)
                    A(i+Nl,j+NmaxC) = fact * A(i+Nl,j+NmaxC)
                    A(i+Nl+NmaxL,j) = fact * A(i+Nl+NmaxL,j)
                    A(i+Nl+NmaxL,j+NmaxC) = fact * A(i+Nl+NmaxL,j+NmaxC)
                  end do
                end if
              else
                m1l = m1
                Nc  = Nrank + (m1 - 1) * (2 * Nrank - m1 + 2)
                do l1 = 1, 2
                  if (mirror) then
                    do j = 1, Nrank - m1 + 1
                      n1 = m1 + j - 1 
                      sp = 1._O + (-1._O)**(n + n1 + m + m1)
                      sm = 2._O - sp              
                      A(i+Nl,j+Nc) = sp * A(i+Nl,j+Nc)
                      A(i+Nl,j+Nc+NmaxC) = sm * A(i+Nl,j+Nc+NmaxC)
                      A(i+Nl+NmaxL,j+Nc) = sm * A(i+Nl+NmaxL,j+Nc)
                      A(i+Nl+NmaxL,j+Nc+NmaxC) = sp * A(i+Nl+NmaxL,j+Nc+NmaxC) 
                    end do
                  end if
                  if (Nazimutsym >= 2) then
                    fact = one
                    do p = 2, Nazimutsym
                      arg  = real((ml + m1l) * (p - 1),O) * dalfa
                      fact = fact + exp(im * arg)
                    end do
                    do j = 1, Nrank - m1 + 1                                               
                      A(i+Nl,j+Nc) = fact * A(i+Nl,j+Nc)
                      A(i+Nl,j+Nc+NmaxC) = fact * A(i+Nl,j+Nc+NmaxC)
                      A(i+Nl+NmaxL,j+Nc) = fact * A(i+Nl+NmaxL,j+Nc)
                      A(i+Nl+NmaxL,j+Nc+NmaxC) = fact * A(i+Nl+NmaxL,j+Nc+NmaxC) 
                    end do
                  end if
                  m1l = - m1l
                  Nc  =  Nc + Nrank - m1 + 1
                end do
              end if
            end do
          end do
          ml = - ml
          Nl =   Nl + Nrank - m + 1
        end do
      end if  
    end do
  end if            
end subroutine AzimMirSym
!***********************************************************************************
subroutine ConvergenceMatrixElem ( index1, index2, TypeGeom, Nsurf, Nparam, ml, mc, &
           Nrank, NintMax, Nint1, Nint2, Nazimutsym, k, surf, kb, ind_ref, mirror,  &
           perfectcond, chiral )  
  use parameters
  use derived_parameters
  implicit none
  integer    :: index1, index2, TypeGeom, Nsurf, Nparam, ml, mc, Nrank, NintMax,    &
                Nint1, Nint2, Nazimutsym                  
  real(O)    :: k, surf(Nsurf), kb
  complex(O) :: ind_ref
  logical    :: mirror, perfectcond, chiral
!      
  integer    :: NmaxL, NmaxC, ml_minus, i, j, Nint1l, Nint2l, NintAL, iparam,       &
                Nintl, pint, iint, imn, jmn   
  real(O)    :: r, theta, phi, dA, param1, param2, pondere, nuv(3), sign, MaxA(2),  &
                MinA(2)
  complex(O) :: ki, zc, zl, zcl, zcr, fact, f, mvl(3), nvl(3), mvc(3), nvc(3), v1,  &
                v2, mixt_product, A11(2,2), A12(2,2), A21(2,2), A22(2,2), expfctL,  &
                expfctC
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG1(:,:), paramG2(:,:), weightsG(:,:) 
  complex(O),allocatable :: mvl1(:,:), nvl1(:,:), mvr1(:,:), nvr1(:,:) 
  complex(O),allocatable :: mv1(:,:), nv1(:,:), mv3(:,:), nv3(:,:)  
!
  NmaxL = Nrank - ml + 1
  NmaxC = Nrank - mc + 1
  ml_minus = - ml                     
  allocate (mv1(3,NmaxC), nv1(3,NmaxC), mv3(3,NmaxL), nv3(3,NmaxL))        
  if (chiral) allocate (mvl1(3,NmaxC), nvl1(3,NmaxC), mvr1(3,NmaxC), nvr1(3,NmaxC))  
  ki   = k * ind_ref
  sign = 1._O
  if (index1 == 3 .and. index2 == 1) sign = - 1._O
  f = sign * im * k * k / Pi 
  if (ml == 1) then 
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
         Nparam, Nintparam, paramG1, paramG2, weightsG, mirror, Nazimutsym) 
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
        call elem_geom3D (TypeGeom, Nsurf, surf, param1, param2, iparam, r, theta,  &
             phi, dA, nuv, mirror, Nazimutsym)                       
        zc = ki * r
        zl = cmplx(k * r,0.0,O)
        if (chiral) then
          zcl = zc / (1._O - kb)
          zcr = zc / (1._O + kb)
        end if                    
        if (index1 == 3 .and. index2 == 1) then
          call MN (3, zl, theta, ml_minus, Nrank, NmaxL, mv3, nv3)
          if (.not. chiral) then
            call MN (1, zc, theta, mc, Nrank, NmaxC, mv1, nv1)             
          else
            call MN (1, zcl, theta, mc, Nrank, NmaxC, mvl1, nvl1)
            call MN (1, zcr, theta, mc, Nrank, NmaxC, mvr1, nvr1)
            call MN_left_right (NmaxC, mvl1, nvl1, mvr1, nvr1, mv1, nv1)
          end if                    
        else if (index1 == 1 .and. index2 == 1) then
          call MN (1, zl, theta, ml_minus, Nrank, NmaxL, mv3, nv3)
          if (.not. chiral) then
            call MN (1, zc, theta, mc, Nrank, NmaxC, mv1, nv1)       
          else 
            call MN (1, zcl, theta, mc, Nrank, NmaxC, mvl1, nvl1)
            call MN (1, zcr, theta, mc, Nrank, NmaxC, mvr1, nvr1)
            call MN_left_right (NmaxC, mvl1, nvl1, mvr1, nvr1, mv1, nv1)
          end if                
        end if
        fact    = f * dA * pondere       
        expfctL =  exp(im * ml_minus * phi)
        expfctC =  exp(im * mc * phi)
        do i = 1, 2
          if (i == 1) then
            imn = 1
          else
            imn = NmaxL
          end if     
          mvl(1) = mv3(1,imn) * expfctL
          mvl(2) = mv3(2,imn) * expfctL
          mvl(3) = mv3(3,imn) * expfctL
          nvl(1) = nv3(1,imn) * expfctL
          nvl(2) = nv3(2,imn) * expfctL
          nvl(3) = nv3(3,imn) * expfctL
          do j = 1, 2
            if (j == 1) then
              jmn = 1
            else
              jmn = NmaxC
            end if 
            mvc(1) = mv1(1,jmn) * expfctC
            mvc(2) = mv1(2,jmn) * expfctC
            mvc(3) = mv1(3,jmn) * expfctC
            nvc(1) = nv1(1,jmn) * expfctC
            nvc(2) = nv1(2,jmn) * expfctC
            nvc(3) = nv1(3,jmn) * expfctC
            if (.not. chiral) then
              v1 = mixt_product (nuv, mvc, nvl)
              v2 = mixt_product (nuv, nvc, mvl)
              if (.not. perfectcond) then
                A11(i,j) = A11(i,j) + (v1 + ind_ref * v2) * fact          
              else 
                A11(i,j) = A11(i,j) + v2 * fact          
              end if
            else 
              v1 = mixt_product (nuv, mvc, nvl)
              v2 = mixt_product (nuv, mvc, mvl)             
              A11(i,j) = A11(i,j) + (v1 + ind_ref * v2) * fact                  
            end if                                                                                
            if (.not. chiral) then             
              v1 = mixt_product (nuv, nvc, nvl)
              v2 = mixt_product (nuv, mvc, mvl)
              if (.not. perfectcond) then
                A12(i,j) = A12(i,j) + (v1 + ind_ref * v2) * fact                        
              else 
                A12(i,j) = A12(i,j) + v2 * fact                 
              end if
              v1 = mixt_product (nuv, mvc, mvl)
              v2 = mixt_product (nuv, nvc, nvl)
              if (.not. perfectcond) then
                A21(i,j) = A21(i,j) + (v1 + ind_ref * v2) * fact                        
              else 
                A21(i,j) = A21(i,j) + v2 * fact                 
              end if              
            else 
              v1 = mixt_product (nuv, nvc, nvl)
              v2 = mixt_product (nuv, nvc, mvl)                   
              A12(i,j) = A12(i,j) + (v1 - ind_ref * v2) * fact
              v1 = mixt_product (nuv, mvc, mvl)
              v2 = mixt_product (nuv, mvc, nvl)                   
              A21(i,j) = A21(i,j) + (v1 + ind_ref * v2) * fact 
            end if
            if (.not. chiral) then 
              v1 = mixt_product (nuv, nvc, mvl)
              v2 = mixt_product (nuv, mvc, nvl)
              if (.not. perfectcond) then
                A22(i,j) = A22(i,j) + (v1 + ind_ref * v2) * fact          
              else 
                A22(i,j) = A22(i,j) + v2 * fact              
              end if
            else 
              v1 = mixt_product (nuv, nvc, mvl)
              v2 = mixt_product (nuv, nvc, nvl)             
              A22(i,j) = A22(i,j) + (v1 - ind_ref * v2) * fact  
            end if                                
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
  deallocate (mv1, nv1, mv3, nv3)
  if (chiral) deallocate (mvl1, nvl1, mvr1, nvr1) 
end subroutine ConvergenceMatrixElem            
   
  
