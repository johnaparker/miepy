! **********************************************************************************
! *         ROUTINES FOR COMPUTING THE Q MATRIX AND THE INCIDENT MATRIX            * 
! *         FOR A COMPOSITE PARTICLE                                               *
! *    -----------------------------------------------------------------------     *
! *    Partial list of subroutines (including the subroutines for a layered        *
! *    particle)                                                                   *
! *      matrix_Q_COMP,                     incident_matrix_COMP,                  *
! *      matrix_Q_DS_COMP,                  incident_matrix_DS_COMP,               *
! *      matrix_Q31_LAY,                    auxmatrices_Q31_LAY,                   *
! *      matrix_Q1_LAY,                     incident_matrix_LAY,                   *
! *      matrix_Q31_DS_LAY,                 auxmatrices_Q31_DS_LAY,                *
! *      matrix_Q1_DS_LAY,                  incident_matrix_DS_LAY                 *
! **********************************************************************************
subroutine matrix_Q_COMP (TypeGeom, index1, index2, k, ind_ref, Nsurfmax, surf, m,  &
           Npart, zpart, Nrankp, NmaxC, NrankL, NmaxL, Nint, Nparammax, Nparam,     &
           Nintparam, paramG, weightsG, A, nap, map)
!------------------------------------------------------------------------------------
! The routine computes the Q matrix of a composite particle for the azimuthal       !
! mode m. Localized sources are used for calculation.                               ! 
!------------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: TypeGeom, index1, index2, m, Nsurfmax, Npart, Nrankp(Npart),        &
                NmaxC, NrankL, NmaxL, Nint, Nparammax, Nintparam(Npart,Nparammax),  &
                Nparam(Npart), nap, map
  real(O)    :: k, zpart(Npart), surf(Npart,Nsurfmax),                              &
                paramG(Npart,Nparammax,Nint), weightsG(Npart,Nparammax,Nint)
  complex(O) :: ind_ref(Npart), A(2*nap,2*map)
!  
  integer    :: i, j, pint, m_minus, iparam, Nintl, ipart, Nparaml
  real(O)    :: r, theta, phi, dA, param, pondere, r1, theta1, n(3), nn(3), tamp,   &
                sign, ct, st
  complex(O) :: kc, ki, zl, fact, f
  complex(O),allocatable :: m1(:,:), n1(:,:), m3(:,:), n3(:,:)
!
  if (index1 == 3 .and. index2 == 1 .and. NmaxL /= NmaxC) then
    print "(/,2x,'Error in subroutine matrix_Q_COMP in module Proces2:')"
    print "(  2x,'the relation NmaxL = NmaxC does not hold;')"
    stop
  end if
  allocate (m1(3,NmaxC), n1(3,NmaxC), m3(3,NmaxL), n3(3,NmaxL))
  kc = cmplx(k,0.0,O)
  do i = 1, 2*NmaxL
    do j = 1, 2*NmaxC
      A(i,j) = zero
    end do
  end do     
  m_minus = - m
  sign = 1._O
  if (index1 == 3 .and. index2 == 1) sign = - 1._O
  f = sign * im * 2._O * k * k
  do ipart = 1, Npart
    ki = k * ind_ref(ipart)
    Nparaml = Nparam(ipart)
    do iparam = 1, Nparaml
      Nintl = Nintparam(ipart,iparam)
      do pint = 1, Nintl
        param   = paramG(ipart,iparam,pint)
        pondere = weightsG(ipart,iparam,pint)
        call elem_geomCOMP (TypeGeom, Npart, ipart, Nsurfmax, surf, param, iparam,  &
             r1, theta1, phi, dA, nn)
        ct = cos(theta1)
        r  = sqrt(r1**2 + zpart(ipart)**2 + 2._O * r1 * zpart(ipart) * ct)
        if (r < MachEps) r = MachEps
        theta = acos((zpart(ipart) + r1 * ct) / r)
        ct   = cos(theta1 - theta)
        st   = sin(theta1 - theta)
        n(1) = nn(1) * ct - nn(2) * st
        n(2) = nn(1) * st + nn(2) * ct
        n(3) = 0._O
        tamp = sqrt(n(1)**2 + n(2)**2)
        if (tamp < MachEps) then
          print "(/,2x,'Error in subroutine matrix_Q_COMP in module Proces2:')"
          print "(  2x,'the module of the normal unit vector is zero;')"
          stop
        end if
        n(1) = n(1) / tamp
        n(2) = n(2) / tamp
        zl   = cmplx(k * r,0.0,O)
        if (index1 == 3 .and. index2 == 1) then       
          call MN_poles_COMP1 (1, ipart, ki, r, theta, m, Npart, zpart, Nrankp,     &
               NmaxC, m1, n1)                                         
          call MN_poles_COMP (3, kc, r, theta, m_minus, Npart, zpart, Nrankp,       &
               NmaxC, m3, n3) 
        else if (index1 == 1 .and. index2 == 1) then
          call MN_poles_COMP1 (1, ipart, ki, r, theta, m, Npart, zpart, Nrankp,     &
               NmaxC, m1, n1)                                                                               
          call MN (1, zl, theta, m_minus, NrankL, NmaxL, m3, n3)
        end if                                   
        fact = f * dA * pondere
        call matQ_COMP(m, NmaxC, NmaxL, ind_ref(ipart), fact, m3, n3, m1, n1, n,    &
             A, nap, map)        
      end do  
    end do
  end do
  deallocate (m1, n1, m3, n3)    
end subroutine matrix_Q_COMP
!***********************************************************************************
subroutine incident_matrix_COMP (TypeGeom, k, Nsurfmax, surf, m, Npart, zpart,      &
           Nrankp, NmaxC, NrankL, NmaxL, Nint, Nparammax, Nparam, Nintparam, paramG,&
           weightsG, A, nap, map)
!----------------------------------------------------------------------------------- 
! The routine computes the incident matrix of a composite particle using           !
! localized sources,                                                               !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: TypeGeom, m, Nsurfmax, Npart, Nrankp(Npart), NmaxC, Nint,           &
                Nparammax, Nparam(Npart), NrankL, NmaxL, Nintparam(Npart,Nparammax),&
                nap, map
  real(O)    :: k, surf(Npart,Nsurfmax), paramG(Npart,Nparammax,Nint),              &
                weightsG(Npart,Nparammax,Nint), zpart(Npart)
  complex(O) :: A(2*nap,2*map)
!
  integer    :: i, j, pint, m_minus, iparam, Nintl, ipart, Nparaml
  real(O)    :: r, theta, phi, dA, param, pondere, r1, theta1, n(3), nn(3), tamp,   &
                ct, st
  complex(O) :: kc, fact, zl, f
  complex(O),allocatable :: m1(:,:), n1(:,:), m3(:,:), n3(:,:)
!
  allocate (m1(3,NmaxL), n1(3,NmaxL), m3(3,NmaxC), n3(3,NmaxC))
  kc = cmplx(k,0.0,O)
  do i = 1, 2*NmaxC
    do j = 1, 2*NmaxL
      A(i,j) = zero
    end do
  end do  
  m_minus = - m
  f = - im * 2._O * k * k
  do ipart = 1, Npart
    Nparaml = Nparam(ipart)
    do iparam = 1, Nparaml
      Nintl = Nintparam(ipart,iparam)
      do pint = 1, Nintl
        param   = paramG(ipart,iparam,pint)
        pondere = weightsG(ipart,iparam,pint)
        call elem_geomCOMP (TypeGeom, Npart, ipart, Nsurfmax, surf, param, iparam,  &
             r1, theta1, phi, dA, nn)
        ct = cos(theta1)
        r  = sqrt(r1**2 + zpart(ipart)**2 + 2._O * r1 * zpart(ipart) * ct)
        if (r < MachEps) r = MachEps
        theta = acos((zpart(ipart) + r1 * ct) / r)
        ct   = cos(theta1 - theta)
        st   = sin(theta1 - theta)
        n(1) = nn(1) * ct - nn(2) * st
        n(2) = nn(1) * st + nn(2) * ct
        n(3) = 0._O
        tamp = sqrt(n(1)**2 + n(2)**2)
        if (tamp < MachEps) then
          print "(/,2x,'Error in subroutine incident_matrix_COMP in module Proces2:')"
          print "(  2x,'the module of the normal unit vector is zero;')"
          stop
        end if
        n(1) = n(1) / tamp
        n(2) = n(2) / tamp
        zl   = cmplx(k * r,0.0,O)        
        call MN (1, zl, theta, m, NrankL, NmaxL, m1, n1)                                                                             
        call MN_poles_COMP (3, kc, r, theta, m_minus, Npart, zpart, Nrankp, NmaxC,  &
             m3, n3)      
        fact = f * dA * pondere
        call matQinc_m(m, NmaxL, NmaxC, fact, m3, n3, m1, n1, n, A, nap, map)        
      end do
    end do
  end do                 
  deallocate (m1, n1, m3, n3)      
end subroutine incident_matrix_COMP
!***********************************************************************************
subroutine matrix_Q_DS_COMP (TypeGeom, index1, index2, k, ind_ref, Nsurfmax, surf,  &
           Npart, Nrankpmax, Nrankp, zRe, zIm, zpart, m, NmaxC, NrankL, NmaxL,      &
           Nint, Nparammax, Nparam, Nintparam, paramG, weightsG, A, nap, map)
!------------------------------------------------------------------------------------
! The routine computes the Q matrix of a composite particle for the azimuthal       !
! mode m. Distributed sources are used for calculation.                             !
!------------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: TypeGeom, index1, index2, m, Nsurfmax, Npart, Nrankpmax,            &
                NmaxC, NrankL, NmaxL, Nint, Nparammax, Nintparam(Npart,Nparammax),  &
                Nparam(Npart), Nrankp(Npart), nap, map
  real(O)    :: k, surf(Npart,Nsurfmax), zpart(Npart), zRe(Npart,Nrankpmax),        &
                zIm(Npart,Nrankpmax), paramG(Npart,Nparammax,Nint),                 &
                weightsG(Npart,Nparammax,Nint)  
  complex(O) :: ind_ref(Npart), A(2*nap,2*map)        
!      
  integer    :: i, j, pint, m_minus, iparam, Nintl, ipart, Nparaml
  real(O)    :: r, theta, phi, dA, param, pondere, r1, theta1, n(3), nn(3), tamp,   &
                sign, ct, st
  complex(O) :: kc, ki, zl, fact, f
  complex(O),allocatable :: m1(:,:), n1(:,:), m3(:,:), n3(:,:) 
!
  if (index1 == 3 .and. index2 == 1 .and. NmaxC /= NmaxL) then
    print "(/,2x,'Error in subroutine matrix_Q_DS_COMP in module Proces2:')"
    print "(  2x,'the relation NmaxC = NmaxL does not hold;')"
    stop
  end if
  allocate (m1(3,NmaxC), n1(3,NmaxC), m3(3,NmaxL), n3(3,NmaxL))
  kc = cmplx(k,0.0,O)
  do i = 1, 2*NmaxL
    do j = 1, 2*NmaxC
      A(i,j) = zero
    end do  
  end do       
  m_minus = - m
  sign = 1._O
  if (index1 == 3 .and. index2 == 1) sign = - 1._O
  f = sign * im * 2._O * k * k
  do ipart = 1, Npart
    ki = k * ind_ref(ipart)
    Nparaml= Nparam(ipart)
    do iparam = 1, Nparaml
      Nintl = Nintparam(ipart,iparam)
      do pint = 1, Nintl
        param   = paramG(ipart,iparam,pint)
        pondere = weightsG(ipart,iparam,pint)
        call elem_geomCOMP (TypeGeom, Npart, ipart, Nsurfmax, surf, param, iparam,  &
             r1, theta1, phi, dA, nn)
        ct = cos(theta1)
        r  = sqrt(r1**2 + zpart(ipart)**2 + 2._O * r1 * zpart(ipart) * ct)
        if (r < MachEps) r = MachEps
        theta = acos((zpart(ipart) + r1 * ct) / r)
        ct   = cos(theta1 - theta)
        st   = sin(theta1 - theta)
        n(1) = nn(1) * ct - nn(2) * st
        n(2) = nn(1) * st + nn(2) * ct
        n(3) = 0._O
        tamp = sqrt(n(1)**2 + n(2)**2)
        if (tamp < MachEps) then
          print "(/,2x,'Error in subroutine matrix_Q_DS_COMP in module Proces2:')"
          print "(  2x,'the module of the normal unit vector is zero;')"
          stop
        end if
        n(1) = n(1) / tamp
        n(2) = n(2) / tamp
        zl   = cmplx(k * r,0.0,O)
        if (index1 == 3 .and. index2 == 1) then          
          call MN_DS_COMP1 (1, ipart, ki, r, theta, Npart, Nrankpmax, Nrankp,       &
               zRe, zIm, m, NmaxC, m1, n1)                                                     
          call MN_DS_COMP (3, kc, r, theta, Npart, Nrankpmax, Nrankp, zRe, zIm,     &
               m_minus, NmaxC, m3, n3)
        else if (index1 == 1 .and. index2 == 1) then              
          call MN_DS_COMP1 (1, ipart, ki, r, theta, Npart, Nrankpmax, Nrankp,       &
               zRe, zIm, m, NmaxC, m1, n1)                                                         
          call MN (1, zl, theta, m_minus, NrankL, NmaxL, m3, n3)
        end if           
        fact = f * dA * pondere
        call matQ_COMP(m, NmaxC, NmaxL, ind_ref(ipart), fact, m3, n3, m1, n1, n,    &
             A, nap, map)          
      end do
    end do
  end do                             
  deallocate (m1, n1, m3, n3)
end subroutine matrix_Q_DS_COMP
!***********************************************************************************
subroutine incident_matrix_DS_COMP (TypeGeom, k, Nsurfmax, surf, Npart, Nrankpmax,  &
           Nrankp, zRe, zIm, zpart, m, NmaxC, NrankL, NmaxL, Nint, Nparammax,       &
           Nparam, Nintparam, paramG, weightsG, A, nap, map)
!-----------------------------------------------------------------------------------
! The routine computes the incident matrix of a composite particle using           !
! distrubuted sources.                                                             !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: TypeGeom, m, Nsurfmax, Npart, Nrankpmax, Nrankp(Npart),             &
                NmaxC, NrankL, NmaxL, Nint, Nparammax, Nparam(Npart),               &
                Nintparam(Npart,Nparammax), nap, map
  real(O)    :: k, surf(Npart,Nsurfmax), zpart(Npart), paramG(Npart,Nparammax,Nint),&
                weightsG(Npart,Nparammax,Nint), zRe(Npart,Nrankpmax),               &
                zIm(Npart,Nrankpmax)
  complex(O) :: A(2*nap,2*map)                       
!      
  integer    :: i, j, pint, m_minus, iparam, Nintl, ipart, Nparaml
  real(O)    :: r, theta, phi, dA, param, pondere, r1, theta1, n(3), nn(3), tamp,   &
                ct, st
  complex(O) :: kc, fact, zl, f
  complex(O),allocatable :: m1(:,:), n1(:,:), m3(:,:), n3(:,:) 
!
  allocate (m1(3,NmaxL), n1(3,NmaxL), m3(3,NmaxC), n3(3,NmaxC))
  kc = cmplx(k,0.0,O)
  do i = 1, 2*NmaxC
    do j = 1, 2*NmaxL
      A(i,j) = zero
    end do  
  end do
  m_minus = - m
  f = - im * 2._O * k * k
  do ipart = 1, Npart
    Nparaml = Nparam(ipart)
    do iparam = 1, Nparaml
      Nintl = Nintparam(ipart,iparam)
      do pint = 1, Nintl
        param   = paramG(ipart,iparam,pint)
        pondere = weightsG(ipart,iparam,pint)
        call elem_geomCOMP (TypeGeom, Npart, ipart, Nsurfmax, surf, param, iparam,  &
             r1, theta1, phi, dA, nn)
        ct = cos(theta1)
        r  = sqrt(r1**2 + zpart(ipart)**2 + 2._O * r1 * zpart(ipart) * ct)
        if (r < MachEps) r = MachEps
        theta = acos((zpart(ipart) + r1 * ct) / r)
        ct   = cos(theta1 - theta)
        st   = sin(theta1 - theta)
        n(1) = nn(1) * ct - nn(2) * st
        n(2) = nn(1) * st + nn(2) * ct
        n(3) = 0._O
        tamp = sqrt(n(1)**2 + n(2)**2)
        if (tamp < MachEps) then
          print "(/,2x,'Error in subroutine incident_matrix_DS_COMP in module Proces2:')"
          print "(  2x,'the module of the normal unit vector is zero;')"
          stop
        end if
        n(1) = n(1) / tamp
        n(2) = n(2) / tamp
        zl   = cmplx(k * r,0.0,O)        
        call MN (1, zl, theta, m, NrankL, NmaxL, m1, n1)                                                                  
        call MN_DS_COMP (3, kc, r, theta, Npart, Nrankpmax, Nrankp, zRe, zIm,       &
             m_minus, NmaxC, m3, n3)                      
        fact = f * dA * pondere
        call matQinc_m(m, NmaxL, NmaxC, fact, m3, n3, m1, n1, n, A, nap, map)         
      end do
    end do
  end do                             
  deallocate (m1, n1, m3, n3)   
end subroutine incident_matrix_DS_COMP
! **********************************************************************************
! *          ROUTINES FOR COMPUTING THE Q MATRIX AND THE INCIDENT MATRIX           * 
! *                         FOR A LAYERED PARTICLE                                 *
! **********************************************************************************
subroutine matrix_Q31_LAY (TypeGeom, k, ind_ref, Nsurfmax, surf, Npart, Nrankp,     &
           Nmaxpmax, Nmaxp, zpart, m, Nmax, Nint, Nparammax, Nparam, Nintparam,     &
           paramG, weightsG, AA, nap, map)
!------------------------------------------------------------------------------------
! The routine computes the Q(3,1) matrix of a layered particle for the azimuthal    !
! mode m. Localized sources are used for calculation.                               !
!------------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: TypeGeom, m, Nsurfmax, Npart, Nmaxpmax, Nmax, Nint, Nparammax,      &
                Nintparam(Npart,Nparammax), Nparam(Npart), Nrankp(Npart),           &
                Nmaxp(Npart), nap, map
  real(O)    :: k, surf(Npart,Nsurfmax), zpart(Npart), paramG(Npart,Nparammax,Nint),&
                weightsG(Npart,Nparammax,Nint)
  complex(O) :: ind_ref(Npart), AA(2*nap,2*map)
!       
  integer    :: i, j, ipart, Nrankl, Nrankc, Nl, Nc  
  complex(O) :: kc, kl, index
  complex(O),allocatable :: a(:,:), b(:,:) 
!
  do i = 1, 2*Nmax
    do j = 1, 2*Nmax
      AA(i,j) = zero
    end do               
  end do  
  do ipart = 1, Npart     
    if (ipart == 1) then
      Nrankl = Nmaxp(ipart)
      Nrankc = Nmaxp(ipart)   
      allocate (a(2*Nrankl,2*Nrankc), b(2*Nrankl,2*Nrankc))          
      kl = cmplx(k,0.0,O)
      kc = k * ind_ref(ipart)
      if (abs(kl) < MachEps) kl = cmplx(MachEps,MachEps,O)
      index = kc / kl
      call auxmatrices_Q31_LAY (TypeGeom, ipart, kl, kc, index, Npart+1,            &
           Nsurfmax, surf, Npart, Nrankp, Nmaxpmax, Nmaxp, zpart, m, Nint,          &
           Nparammax, Nparam, Nintparam, paramG, weightsG, Nrankl, Nrankc,          &
           22, a, Nrankl, Nrankc, b, Nrankl, Nrankc) 
      do i = 1, 2*Nrankl
        do j = 1, 2*Nrankc
          AA(i,j) = a(i,j)
          AA(i,j+2*Nrankc) = b(i,j)
        end do
      end do
      deallocate (a, b)
      Nl = 2 * Nmaxp(ipart)
      Nc = 0        
    else    
      Nrankl = Nmaxp(ipart-1)
      Nrankc = Nmaxp(ipart-1) 
      allocate (a(2*Nrankl,2*Nrankc), b(2*Nrankl,2*Nrankc)) 
      do i = 1, 2*Nrankl
        do j = 1, 2*Nrankc
          a(i,j) = zero
          b(i,j) = zero
          if (i == j) b(i,j) = - one
        end do
      end do                
      do i = 1, 2*Nrankl
        do j = 1, 2*Nrankc
          AA(i+Nl,j+Nc) = a(i,j)
          AA(i+Nl,j+Nc+2*Nrankc) = b(i,j)
        end do
      end do
      deallocate (a, b)
      Nc = Nc + 4 * Nmaxp(ipart-1)
!
      Nrankl = Nmaxp(ipart-1)
      Nrankc = Nmaxp(ipart) 
      allocate (a(2*Nrankl,2*Nrankc), b(2*Nrankl,2*Nrankc))            
      kl = k * ind_ref(ipart-1) 
      kc = k * ind_ref(ipart)
      if (abs(kl) < MachEps) kl = cmplx(MachEps,MachEps,O)
      index = kc / kl
      call auxmatrices_Q31_LAY (TypeGeom, ipart, kl, kc, index, Npart,              &
           Nsurfmax, surf, Npart, Nrankp, Nmaxpmax, Nmaxp, zpart, m, Nint,          &
           Nparammax, Nparam, Nintparam, paramG, weightsG, Nrankl, Nrankc,          &
           12, a, Nrankl, Nrankc, b, Nrankl, Nrankc) 
      do i = 1, 2*Nrankl
        do j = 1, 2*Nrankc
          AA(i+Nl,j+Nc) = a(i,j)
          if(ipart < Npart) AA(i+Nl,j+Nc+2*Nrankc) = b(i,j)
        end do
      end do                 
      Nl = Nl + 2 * Nmaxp(ipart-1)
      Nc = Nc - 4 * Nmaxp(ipart-1)
      deallocate (a, b)
!
      Nrankl = Nmaxp(ipart)
      Nrankc = Nmaxp(ipart-1) 
      allocate (a(2*Nrankl,2*Nrankc), b(2*Nrankl,2*Nrankc)) 
      kl = k * ind_ref(ipart-1) 
      kc = k * ind_ref(ipart-1)
      if (abs(kl) < MachEps) kl = cmplx(MachEps,MachEps,O)
      index = kc / kl
      if (zpart(ipart) /= zpart(ipart-1)) then
        call auxmatrices_Q31_LAY (TypeGeom, ipart, kl, kc, index, Npart+1,          &
             Nsurfmax, surf, Npart, Nrankp, Nmaxpmax, Nmaxp, zpart, m, Nint,        &
             Nparammax, Nparam, Nintparam, paramG, weightsG, Nrankl, Nrankc,        &
             21, a, Nrankl, Nrankc, b, Nrankl, Nrankc)	         	     
      else        
        do i = 1, Nrankl
          do j = 1, Nrankc
            a(i,j) = zero
            a(i,j+Nrankc) = zero
            a(i+Nrankl,j) = zero
            a(i+Nrankl,j+Nrankc) = zero	                	    
            if (j == i) then
              a(i,j) = one
              a(i+Nrankl,j+Nrankc) = one
            end if
          end do
        end do
      end if
      do i = 1, 2*Nrankl
        do j = 1, 2*Nrankc
          AA(i+Nl,j+Nc) = a(i,j)
          AA(i+Nl,j+Nc+2*Nrankc) = zero
        end do
      end do
      Nc = Nc + 4 * Nmaxp(ipart-1)
      deallocate (a, b)
!
      Nrankl = Nmaxp(ipart)
      Nrankc = Nmaxp(ipart) 
      allocate (a(2*Nrankl,2*Nrankc), b(2*Nrankl,2*Nrankc))           
      kl = k * ind_ref(ipart-1) 
      kc = k * ind_ref(ipart)
      if (abs(kl) < MachEps) kl = cmplx(MachEps,MachEps,O)
      index = kc / kl
      call auxmatrices_Q31_LAY (TypeGeom, ipart, kl, kc, index, Npart,              &
           Nsurfmax, surf, Npart, Nrankp, Nmaxpmax, Nmaxp, zpart, m, Nint,          &
           Nparammax, Nparam, Nintparam, paramG, weightsG, Nrankl, Nrankc,          &
           22, a, Nrankl, Nrankc, b, Nrankl, Nrankc) 
      do i = 1, 2*Nrankl
        do j = 1, 2*Nrankc
          AA(i+Nl,j+Nc) = a(i,j)
          if (ipart < Npart) AA(i+Nl,j+Nc+2*Nrankc) = b(i,j)
        end do
      end do
      Nl = Nl + 2 * Nmaxp(ipart)
      deallocate (a, b)       
    end if
  end do  
end subroutine matrix_Q31_LAY
!***********************************************************************************
subroutine auxmatrices_Q31_LAY (TypeGeom, ipart, kl, kc, index, Npartlim, Nsurfmax, &
           surf, Npart, Nrankp, Nmaxpmax, Nmaxp, zpart, m, Nint, Nparammax, Nparam, &
           Nintparam, paramG, weightsG, Nrankl, Nrankc, Block, a, nap, map,         &
           b, nbp, mbp)
!-----------------------------------------------------------------------------------
! The routine computes the auxiliary matrices a and b used in the subroutine       !
! matrix_Q31_LAY.                                                                  !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: TypeGeom, ipart, Npartlim, m, Nsurfmax, Npart, Nmaxpmax, Nint,      &
                Nparammax, Nintparam(Npart,Nparammax), Nparam(Npart),               &
                Nrankp(Npart), Nmaxp(Npart), Nrankl, Nrankc, nap, map, nbp, mbp,    &
                Block
  real(O)    :: surf(Npart,Nsurfmax), zpart(Npart), paramG(Npart,Nparammax,Nint),   &
                weightsG(Npart,Nparammax,Nint)
  complex(O) :: kl, kc, index, a(2*nap,2*map), b(2*nbp,2*mbp)
!       
  integer    :: i, j, pint, m_minus, iparam, Nintl, Nparaml, ipartl
  real(O)    :: r, theta, phi, dA, param, pondere, r1, theta1, n(3), nn(3), sign,   &
                tamp, ct, st
  complex(O) :: fact, f
  complex(O),allocatable :: m1(:,:), n1(:,:), m3(:,:), n3(:,:), m2(:,:), n2(:,:) 
!
  allocate (m1(3,Nmaxpmax), n1(3,Nmaxpmax), m2(3,Nmaxpmax), n2(3,Nmaxpmax),         &
            m3(3,Nmaxpmax), n3(3,Nmaxpmax))  
  m_minus = - m  
  do i = 1, 2*Nrankl
    do j = 1, 2*Nrankc
      a(i,j) = zero
      b(i,j) = zero
    end do
  end do  
  if (Block == 22 .or. Block == 12) then
    sign   = 1._O
    ipartl = ipart
  else if (Block == 21) then
    sign   = - 1._O
    ipartl =   ipart - 1
  end if                  
  Nparaml = Nparam(ipartl)
  f = sign * im * 2._O * kl * kl
  do iparam = 1, Nparaml
    Nintl = Nintparam(ipartl,iparam)
    do pint = 1, Nintl
      param   = paramG(ipartl,iparam,pint)
      pondere = weightsG(ipartl,iparam,pint)
      call elem_geomLAY (TypeGeom, Npart, ipartl, Nsurfmax, surf, param, iparam,    &
           r1, theta1, phi, dA, nn)
      ct = cos(theta1)
      r  = sqrt(r1**2 + zpart(ipartl)**2 + 2._O * r1 * zpart(ipartl) * ct)
      if (r < MachEps) r = MachEps
      theta = acos((zpart(ipartl) + r1 * ct) / r)
      ct   = cos(theta1 - theta)
      st   = sin(theta1 - theta)
      n(1) = nn(1) * ct - nn(2) * st
      n(2) = nn(1) * st + nn(2) * ct
      n(3) = 0._O
      tamp = sqrt(n(1)**2 + n(2)**2)
      if (tamp < MachEps) then
        print "(/,2x,'Error in subroutine auxmatrices_Q31_LAY in module Proces2:')"
        print "(  2x,'the module of the normal unit vector is zero;')"
        stop
      end if
      n(1) = n(1) / tamp  
      n(2) = n(2) / tamp  
      if (Block == 21) then  
        call MN_poles_LAY (1, kc, r, theta, m, ipart-1, Npart, Nrankp, Nmaxpmax,    &
             Nmaxp, zpart, m1, n1)           
        !call MN_poles_LAY (3, kc, r, theta, m, ipart-1, Npart, Nrankp, Nmaxpmax,   &
        !     Nmaxp, zpart, m2, n2)               
        call MN_poles_LAY (3, kl, r, theta, m_minus, ipart, Npart, Nrankp, Nmaxpmax,&
             Nmaxp, zpart, m3, n3)
      else if (Block == 22) then
        call MN_poles_LAY (1, kc, r, theta, m, ipart, Npart, Nrankp, Nmaxpmax,      &
             Nmaxp, zpart, m1, n1)           
        call MN_poles_LAY (3, kc, r, theta, m, ipart, Npart, Nrankp, Nmaxpmax,      &
             Nmaxp, zpart, m2, n2)               
        call MN_poles_LAY (3, kl, r, theta, m_minus, ipart, Npart, Nrankp, Nmaxpmax,&
             Nmaxp, zpart, m3, n3)
      else if (Block == 12) then 
        call MN_poles_LAY (1, kc, r, theta, m, ipart, Npart, Nrankp, Nmaxpmax,      &
             Nmaxp, zpart, m1, n1)           
        call MN_poles_LAY (3, kc, r, theta, m, ipart, Npart, Nrankp, Nmaxpmax,      &
             Nmaxp, zpart, m2, n2)
        call MN_poles_LAY (1, kl, r, theta, m_minus, ipart-1, Npart, Nrankp,        &
             Nmaxpmax, Nmaxp, zpart, m3, n3)
      end if 
      fact = f * dA * pondere
      if (Block == 21) then
        call matQ_LAY_a (m, Nrankc, Nrankl, Nmaxpmax, index, fact, m3, n3, m1, n1,  &
             n, a, nap, map)            
      else       
        call matQ_LAY (m, ipart, Npartlim, Nrankc, Nrankl, Nmaxpmax, index, fact,   &
             m3, n3, m2, n2, m1, n1, n, a, nap, map, b, nbp, mbp)      
      end if	     
    end do
  end do
  deallocate (m1, n1, m2, n2, m3, n3)
end subroutine auxmatrices_Q31_LAY 
!***********************************************************************************
subroutine matrix_Q1_LAY (TypeGeom, indexC, k, ind_ref, Nsurfmax, surf, Npart,      &
           Nrankpmax, Nrankp, Nmaxpmax, Nmaxp, zpart, m, Nint, Nparammax, Nparam,   &
           Nintparam, paramG, weightsG, A, nap, map)
!-----------------------------------------------------------------------------------
! The routine computes the Q(1,1) or the Q(1,3) matrices of a layered particle     !
! for the azimuthal mode m. Localized sources are used for calculation.            !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: TypeGeom, m, Npart, Nsurfmax, Nmaxpmax, Nrankp(Npart),              &
                Nmaxp(Npart), Nint, Nparammax, Nintparam(Npart,Nparammax),          &
                Nparam(Npart), Nrankpmax, indexC, nap, map
  real(O)    :: k, surf(Npart,Nsurfmax), zpart(Npart), paramG(Npart,Nparammax,Nint),&
                weightsG(Npart,Nparammax,Nint)
  complex(O) :: ind_ref(Npart), A(2*nap,2*map)
!           
  integer    :: i, j, pint, m_minus, iparam, Nintl, ipart, Nparaml
  real(O)    :: r, theta, phi, dA, param, pondere, r1, theta1, n(3), nn(3), tamp,   &
                ct, st
  complex(O) :: kc, zl, fact, f
  complex(O),allocatable :: m1(:,:), n1(:,:), m3(:,:), n3(:,:)  
!
  allocate (m1(3,Nmaxpmax), n1(3,Nmaxpmax), m3(3,Nmaxpmax), n3(3,Nmaxpmax))
  do i = 1, 2*Nmaxpmax
    do j = 1, 2*Nmaxpmax
      A(i,j) = zero
    end do
  end do
  m_minus = - m 
  ipart   =   1  
  Nparaml =   Nparam(ipart)
  kc = k * ind_ref(ipart)
  f  = im * 2._O * k * k
  do iparam = 1, Nparaml
    Nintl = Nintparam(ipart,iparam)
    do pint = 1, Nintl
      param   = paramG(ipart,iparam,pint)
      pondere = weightsG(ipart,iparam,pint)
      call elem_geomLAY (TypeGeom, Npart, ipart, Nsurfmax, surf, param, iparam,     &
           r1, theta1, phi, dA, nn)
      ct = cos(theta1)
      r  = sqrt(r1**2 + zpart(ipart)**2 + 2._O * r1 * zpart(ipart) * ct)
      if (r < MachEps) r = MachEps
      theta = acos((zpart(ipart) + r1 * ct) / r)
      ct   = cos(theta1 - theta)
      st   = sin(theta1 - theta)
      n(1) = nn(1) * ct - nn(2) * st
      n(2) = nn(1) * st + nn(2) * ct
      n(3) = 0._O
      tamp = sqrt(n(1)**2 + n(2)**2)
      if (tamp < MachEps) then
        print "(/,2x,'Error in subroutine matrix_Q1_LAY in module Proces2:')"
        print "(  2x,'the module of the normal unit vector is zero;')"
        stop
      end if
      n(1) = n(1) / tamp
      n(2) = n(2) / tamp          
      if (indexC == 1) then      
        call MN_poles_LAY (1, kc, r, theta, m, ipart, Npart, Nrankp, Nmaxpmax,      &
             Nmaxp, zpart, m1, n1)        
      else     
        call MN_poles_LAY (3, kc, r, theta, m, ipart, Npart, Nrankp, Nmaxpmax,      &
             Nmaxp, zpart, m1, n1)
      end if
      zl = cmplx(k * r,0.0,O)              
      call MN (1, zl, theta, m_minus, Nrankpmax, Nmaxpmax, m3, n3)
      fact = f * dA * pondere
      call matQ_COMP(m, Nmaxpmax, Nmaxpmax, ind_ref(ipart), fact, m3, n3, m1, n1,   &
           n, A, nap, map)      
    end do
  end do
  deallocate (m1, n1, m3, n3)     
end subroutine matrix_Q1_LAY
!***********************************************************************************
subroutine incident_matrix_LAY (TypeGeom, k, Nsurfmax, surf, Npart, Nrankpmax,      &
           Nrankp, Nmaxpmax, Nmaxp, zpart, m, Nint, Nparammax, Nparam, Nintparam,   &
           paramG, weightsG, A, nap, map)
!-----------------------------------------------------------------------------------
! The routine computes the incident matrix of a layered particle using             !
! localized sources.                                                               !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: TypeGeom, m, Nsurfmax, Npart, Nrankpmax, Nrankp(Npart), Nint,       &
                Nparammax, Nparam(Npart), Nmaxpmax, Nmaxp(Npart),                   &
                Nintparam(Npart,Nparammax), nap, map
  real(O)    :: k, surf(Npart,Nsurfmax), zpart(Npart), paramG(Npart,Nparammax,Nint),&
                weightsG(Npart,Nparammax,Nint)
  complex(O) :: A(2*nap,2*map)
!      
  integer    :: i, j, pint, m_minus, iparam, Nintl, ipart, Nparaml
  real(O)    :: r, theta, phi, dA, param, pondere, r1, theta1, n(3), nn(3), tamp,   &
                ct, st
  complex(O) :: kc, fact, zl, f
  complex(O),allocatable :: m1(:,:), n1(:,:), m3(:,:), n3(:,:) 
!
  allocate (m1(3,Nmaxpmax), n1(3,Nmaxpmax), m3(3,Nmaxpmax), n3(3,Nmaxpmax))
  kc = cmplx(k,0.0,O)
  do i = 1, 2*Nmaxpmax
    do j = 1, 2*Nmaxpmax
      A(i,j) = zero       
    end do
  end do        
  m_minus = - m
  ipart   =   1
  Nparaml =   Nparam(ipart)
  f = im * 2._O * k * k
  do iparam = 1, Nparaml
    Nintl = Nintparam(ipart,iparam)
    do pint = 1, Nintl
      param   = paramG(ipart,iparam,pint)
      pondere = weightsG(ipart,iparam,pint)
      call elem_geomLAY (TypeGeom, Npart, ipart, Nsurfmax, surf, param, iparam,     &
           r1, theta1, phi, dA, nn)
      ct = cos(theta1)
      r  = sqrt(r1**2 + zpart(ipart)**2 + 2._O * r1 * zpart(ipart) * ct)
      if (r < MachEps) r = MachEps
      theta = acos((zpart(ipart) + r1 * ct) / r)
      ct   = cos(theta1 - theta)
      st   = sin(theta1 - theta)
      n(1) = nn(1) * ct - nn(2) * st
      n(2) = nn(1) * st + nn(2) * ct
      n(3) = 0._O
      tamp = sqrt(n(1)**2 + n(2)**2)
      if (tamp < MachEps) then
        print "(/,2x,'Error in subroutine incident_matrix_LAY in module Proces2:')"
        print "(  2x,'the module of the normal unit vector is zero;')"
        stop
      end if
      n(1) = n(1) / tamp
      n(2) = n(2) / tamp
      zl   = cmplx(k * r,0.0,O)      
      call MN (1, zl, theta, m, Nrankpmax, Nmaxpmax, m1, n1)                                           
      call MN_poles_LAY (3, kc, r, theta, m_minus, ipart, Npart, Nrankp, Nmaxpmax,  &
           Nmaxp, zpart, m3, n3)                        
      fact =  f * dA * pondere
      call matQinc_m(m, Nmaxpmax, Nmaxpmax, fact, m3, n3, m1, n1, n, A, nap, map)      
    end do
  end do                   
  deallocate (m1, n1, m3, n3)      
end subroutine incident_matrix_LAY
!***********************************************************************************
subroutine matrix_Q31_DS_LAY (TypeGeom, k, ind_ref, Nsurfmax, surf, Npart,          &
           Nrankpmax, Nrankp, zRe, zIm, zpart, m, Nrank, Nint, Nparammax, Nparam,   &
           Nintparam, paramG, weightsG, AA, nap, map)
!------------------------------------------------------------------------------------
! The routine computes the Q(3,1) matrix of a layered particle for the azimuthal    !
! mode m. Distributed sources are used for calculation.                             !
!------------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: TypeGeom, m, Nsurfmax, Npart, Nrankpmax, Nrank, Nint, Nparammax,    &
                Nintparam(Npart,Nparammax), Nparam(Npart), Nrankp(Npart), nap, map
  real(O)    :: k, surf(Npart,Nsurfmax), zpart(Npart), paramG(Npart,Nparammax,Nint),&
                weightsG(Npart,Nparammax,Nint), zRe(Npart,Nrankpmax),               &
                zIm(Npart,Nrankpmax)
  complex(O) :: ind_ref(Npart), AA(2*nap,2*map)
!      
  integer    :: i, j, ipart, Nrankl, Nrankc, Nl, Nc
  complex(O) :: kc, kl, index
  complex(O),allocatable :: a(:,:), b(:,:) 
!
  do i = 1, 2*Nrank
    do j = 1, 2*Nrank
      AA(i,j) = zero
    end do
  end do  
  do ipart = 1, Npart
    if (ipart == 1) then
      Nrankl = Nrankp(ipart)
      Nrankc = Nrankp(ipart) 
      allocate (a(2*Nrankl,2*Nrankc), b(2*Nrankl,2*Nrankc))         
      kl = cmplx(k,0.0,O)
      kc = k * ind_ref(ipart)
      if (abs(kl) < MachEps) kl = cmplx(MachEps,MachEps,O)
      index = kc / kl
      call auxmatrices_Q31_DS_LAY (TypeGeom, ipart, kl, kc, index, Npart+1,         &
           Nsurfmax, surf, Npart, Nrankpmax, Nrankp, zRe, zIm, zpart, m, Nint,      &
           Nparammax, Nparam, Nintparam, paramG, weightsG, Nrankl, Nrankc,          &
           22, a, Nrankl, Nrankc, b, Nrankl, Nrankc)
      do i = 1, 2*Nrankl
        do j = 1, 2*Nrankc
          AA(i,j) = a(i,j)
          AA(i,j+2*Nrankc) = b(i,j)
        end do
      end do
      deallocate (a, b)
      Nl = 2 * Nrankp(ipart)
      Nc = 0
    else
      Nrankl = Nrankp(ipart-1)
      Nrankc = Nrankp(ipart-1) 
      allocate (a(2*Nrankl,2*Nrankc), b(2*Nrankl,2*Nrankc))           
      kl = k * ind_ref(ipart-1) 
      kc = k * ind_ref(ipart-1)
      if (abs(kl) < MachEps) kl = cmplx(MachEps,MachEps,O)
      index = kc / kl
      call auxmatrices_Q31_DS_LAY (TypeGeom, ipart, kl, kc, index, Npart+1,         &
           Nsurfmax, surf, Npart, Nrankpmax, Nrankp, zRe, zIm, zpart, m, Nint,      &
           Nparammax, Nparam, Nintparam, paramG, weightsG, Nrankl, Nrankc,          &
           11, a, Nrankl, Nrankc, b, Nrankl, Nrankc)	
      do i = 1, 2*Nrankl
        do j = 1, 2*Nrankc
          AA(i+Nl,j+Nc) = zero
          AA(i+Nl,j+Nc+2*Nrankc) = b(i,j)
        end do
      end do
      deallocate (a, b)
      Nc = Nc + 4 * Nrankp(ipart-1)
!
      Nrankl = Nrankp(ipart-1)
      Nrankc = Nrankp(ipart) 
      allocate (a(2*Nrankl,2*Nrankc), b(2*Nrankl,2*Nrankc))      
      kl = k * ind_ref(ipart-1) 
      kc = k * ind_ref(ipart)
      if (abs(kl) < MachEps) kl = cmplx(MachEps,MachEps,O)
      index = kc / kl
      call auxmatrices_Q31_DS_LAY (TypeGeom, ipart, kl, kc, index, Npart,           &
           Nsurfmax, surf, Npart, Nrankpmax, Nrankp, zRe, zIm, zpart, m, Nint,      &
           Nparammax, Nparam, Nintparam, paramG, weightsG, Nrankl, Nrankc,          &
           12, a, Nrankl, Nrankc, b, Nrankl, Nrankc)
      do i = 1, 2*Nrankl
        do j = 1, 2*Nrankc
          AA(i+Nl,j+Nc) = a(i,j)
          if (ipart < Npart) AA(i+Nl,j+Nc+2*Nrankc) = b(i,j)
        end do   
      end do
      Nl = Nl + 2 * Nrankp(ipart-1)
      Nc = Nc - 4 * Nrankp(ipart-1)
      deallocate (a, b)
!
      Nrankl = Nrankp(ipart)
      Nrankc = Nrankp(ipart-1) 
      allocate (a(2*Nrankl,2*Nrankc), b(2*Nrankl,2*Nrankc))          
      kl = k * ind_ref(ipart-1) 
      kc = k * ind_ref(ipart-1)
      if (abs(kl) < MachEps) kl = cmplx(MachEps,MachEps,O)
      index = kc / kl
      call auxmatrices_Q31_DS_LAY (TypeGeom, ipart, kl, kc, index, Npart+1,         &
           Nsurfmax, surf, Npart, Nrankpmax, Nrankp, zRe, zIm, zpart, m, Nint,      &
           Nparammax, Nparam, Nintparam, paramG, weightsG, Nrankl, Nrankc,          &
           21, a, Nrankl, Nrankc, b, Nrankl, Nrankc)	              	     	  	  	       	     	     	     	      	   	   	   
      do i = 1, 2*Nrankl
        do j = 1, 2*Nrankc
          AA(i+Nl,j+Nc) = a(i,j)
          AA(i+Nl,j+Nc+2*Nrankc) = zero
        end do
      end do 
      Nc = Nc + 4 * Nrankp(ipart-1)
      deallocate (a, b)
!
      Nrankl = Nrankp(ipart)
      Nrankc = Nrankp(ipart) 
      allocate (a(2*Nrankl,2*Nrankc), b(2*Nrankl,2*Nrankc))         
      kl = k * ind_ref(ipart-1) 
      kc = k * ind_ref(ipart)
      if (abs(kl) < MachEps) kl = cmplx(MachEps,MachEps,O)
      index = kc / kl
      call auxmatrices_Q31_DS_LAY (TypeGeom, ipart, kl, kc, index, Npart,           &
           Nsurfmax, surf, Npart, Nrankpmax, Nrankp, zRe, zIm, zpart, m, Nint,      &
           Nparammax, Nparam, Nintparam, paramG, weightsG, Nrankl, Nrankc,          &
           22, a, Nrankl, Nrankc, b, Nrankl, Nrankc)
      do i = 1, 2 * Nrankl
        do j = 1, 2 * Nrankc
          AA(i+Nl,j+Nc) = a(i,j)
          if (ipart < Npart) AA(i+Nl,j+Nc+2*Nrankc) = b(i,j)
        end do
      end do  
      Nl = Nl + 2 * Nrankp(ipart)
      deallocate (a, b)       
    end if
  end do       
end subroutine matrix_Q31_DS_LAY
!***********************************************************************************
subroutine auxmatrices_Q31_DS_LAY (TypeGeom, ipart, kl, kc, index, Npartlim,        &
           Nsurfmax, surf, Npart, Nrankpmax, Nrankp, zRe, zIm, zpart, m, Nint,      &
           Nparammax, Nparam, Nintparam, paramG, weightsG, Nrankl, Nrankc, Block,   &
           a, nap, map, b, nbp, mbp)
!------------------------------------------------------------------------------------
! The routine computes the auxiliary matrices a and b used in the subroutine        !
! matrix_Q31_DS_LAY.                                                                !
!------------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: TypeGeom, ipart, Npartlim, m, Nsurfmax, Npart, Nrankpmax, Nint,     &
                Nparammax, Nintparam(Npart,Nparammax), Nparam(Npart), Nrankp(Npart),&
                Nrankl, Nrankc, Block, nap, map, nbp, mbp
  real(O)    :: surf(Npart,Nsurfmax), zpart(Npart), paramG(Npart,Nparammax,Nint),   &
                weightsG(Npart,Nparammax,Nint), zRe(Npart,Nrankpmax),               &
                zIm(Npart,Nrankpmax)
  complex(O) :: kl, kc, index, a(2*nap,2*map), b(2*nbp,2*mbp)
!      
  integer    :: i, j, pint, m_minus, iparam, Nintl, Nparaml, ipartl
  real(O)    :: r, theta, phi, dA, param, pondere, r1, theta1, n(3), nn(3), sign,   &
                tamp, ct, st
  complex(O) :: fact, f
  complex(O),allocatable :: m1(:,:), n1(:,:), m3(:,:), n3(:,:), m2(:,:), n2(:,:)
!
  allocate (m1(3,Nrankpmax), n1(3,Nrankpmax), m2(3,Nrankpmax), n2(3,Nrankpmax),     &
            m3(3,Nrankpmax), n3(3,Nrankpmax))  
  m_minus = - m 
  do i = 1, 2*Nrankl
    do j = 1, 2*Nrankc
      a(i,j) = zero
      b(i,j) = zero
    end do
  end do
  if (Block == 22 .or. Block == 12) then
    sign   = 1._O
    ipartl = ipart
  else if (Block == 11 .or. Block == 21) then
    sign   = - 1._O
    ipartl = ipart - 1
  end if                         
  Nparaml = Nparam(ipartl)
  f = sign * im * 2._O * kl * kl 
  do iparam = 1, Nparaml
    Nintl = Nintparam(ipartl,iparam)
    do pint = 1, Nintl
      param   = paramG(ipartl,iparam,pint)
      pondere = weightsG(ipartl,iparam,pint)
      call elem_geomLAY (TypeGeom, Npart, ipartl, Nsurfmax, surf, param, iparam,    &
           r1, theta1, phi, dA, nn)
      ct = cos(theta1)
      r  = sqrt(r1**2 + zpart(ipartl)**2 + 2._O * r1 * zpart(ipartl) * ct)
      if (r < MachEps) r = MachEps
      theta = acos((zpart(ipartl) + r1 * ct) / r)
      ct   = cos(theta1 - theta)
      st   = sin(theta1 - theta)
      n(1) = nn(1) * ct - nn(2) * st
      n(2) = nn(1) * st + nn(2) * ct
      n(3) = 0._O
      tamp = sqrt(n(1)**2 + n(2)**2)
      if (tamp < MachEps) then
        print "(/,2x,'Error in subroutine auxmatrices_Q31_DS_LAY in module Proces2:')"
        print "(  2x,'the module of the normal unit vector is zero;')"
        stop
      end if
      n(1) = n(1) / tamp 
      n(2) = n(2) / tamp
      if (Block == 22) then           
        call MN_DS_LAY (1, kc, r, theta, ipart, Npart, Nrankpmax, Nrankp, zRe, zIm,  &
             m, m1, n1)   
        call MN_DS_LAY (3, kc, r, theta, ipart, Npart, Nrankpmax, Nrankp, zRe, zIm,  &
             m, m2, n2)   
        call MN_DS_LAY (3, kl, r, theta, ipart, Npart, Nrankpmax, Nrankp, zRe, zIm,  &
             m_minus, m3, n3) 
      else if (Block == 11) then
       !call MN_DS_LAY (1, kc, r, theta, ipart-1, Npart, Nrankpmax, Nrankp, zRe, zIm,&
       !     m, m1, n1)   
        call MN_DS_LAY (3, kc, r, theta, ipart-1, Npart, Nrankpmax, Nrankp, zRe, zIm,&
             m, m2, n2)   
        call MN_DS_LAY (1, kl, r, theta, ipart-1, Npart, Nrankpmax, Nrankp, zRe, zIm,&
             m_minus, m3, n3) 
      else if (Block == 12) then
        call MN_DS_LAY (1, kc, r, theta, ipart, Npart, Nrankpmax, Nrankp, zRe, zIm,  &
             m, m1, n1)   
        call MN_DS_LAY (3, kc, r, theta, ipart, Npart, Nrankpmax, Nrankp, zRe, zIm,  &
             m, m2, n2)   
        call MN_DS_LAY (1, kl, r, theta, ipart-1, Npart, Nrankpmax, Nrankp, zRe, zIm,&
             m_minus, m3, n3) 
      else if (Block == 21) then
        call MN_DS_LAY (1, kc, r, theta, ipart-1, Npart, Nrankpmax, Nrankp, zRe, zIm,&
             m, m1, n1)   
       !call MN_DS_LAY (3, kc, r, theta, ipart-1, Npart, Nrankpmax, Nrankp, zRe, zIm,&
       !     m, m2, n2)   
        call MN_DS_LAY (3, kl, r, theta, ipart, Npart, Nrankpmax, Nrankp, zRe, zIm,  &
             m_minus, m3, n3) 
      end if
      fact = f * dA * pondere
      if (Block == 11) then
        call matQ_LAY_b (m, ipart, Npartlim, Nrankc, Nrankl, Nrankpmax, index, fact, & 
             m3, n3, m2, n2, n, b, nbp, mbp)            
      else if (Block == 21) then
        call matQ_LAY_a (m, Nrankc, Nrankl, Nrankpmax, index, fact, m3, n3, m1, n1,  &
           n, a, nap, map)            
      else 
        call matQ_LAY(m, ipart, Npartlim, Nrankc, Nrankl, Nrankpmax, index, fact,    &
             m3, n3, m2, n2, m1, n1, n, a, nap, map, b, nbp, mbp)       
      end if	     
    end do
  end do
  deallocate (m1, n1, m2, n2, m3, n3)
end subroutine auxmatrices_Q31_DS_LAY 
!***********************************************************************************
subroutine matrix_Q1_DS_LAY (TypeGeom, indexC, k, ind_ref, Nsurfmax, surf, Npart,   &
           Nrankpmax, Nrankp, zRe, zIm, zpart, m, NrankG, NmaxG, Nint, Nparammax,   &
           Nparam, Nintparam, paramG, weightsG, A, nap, map)
!------------------------------------------------------------------------------------
! The routine computes the Q(1,1) or the Q(1,3) matrices of a layered particle      !
! for the azimuthal mode m. Distributed sources are used for calculation.           !
!------------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: TypeGeom, m, Npart, Nsurfmax, Nrankpmax, Nrankp(Npart),             &
                NrankG, NmaxG, Nint, Nparammax, Nintparam(Npart,Nparammax),         &
                Nparam(Npart), indexC, nap, map
  real(O)    :: k, surf(Npart,Nsurfmax), zpart(Npart), paramG(Npart,Nparammax,Nint),&
                weightsG(Npart,Nparammax,Nint), zRe(Npart,Nrankpmax),               &
                zIm(Npart,Nrankpmax) 
  complex(O) :: ind_ref(Npart), A(2*nap,2*map)
!            
  integer    :: i, j, pint, m_minus, iparam, Nintl, ipart, Nparaml
  real(O)    :: r, theta, phi, dA, param, pondere, r1, theta1, n(3), nn(3), tamp,   &
                ct, st
  complex(O) :: kc, zl, fact, f
  complex(O),allocatable :: m1(:,:), n1(:,:), m3(:,:), n3(:,:)  
!
  allocate (m1(3,Nrankpmax), n1(3,Nrankpmax), m3(3,NmaxG), n3(3,NmaxG))
  do i = 1, 2*NmaxG
    do j = 1, 2*Nrankpmax
      A(i,j) = zero
    end do
  end do
  m_minus = - m         
  ipart   =   1  
  Nparaml =   Nparam(ipart)
  kc = k * ind_ref(ipart)
  f  = im * 2._O * k * k 
  do iparam = 1, Nparaml
    Nintl = Nintparam(ipart,iparam)
    do pint = 1, Nintl
      param   = paramG(ipart,iparam,pint)
      pondere = weightsG(ipart,iparam,pint)
      call elem_geomLAY (TypeGeom, Npart, ipart, Nsurfmax, surf, param, iparam, r1, &
           theta1, phi, dA, nn)
      ct = cos(theta1)
      r  = sqrt(r1**2 + zpart(ipart)**2 + 2._O * r1 * zpart(ipart) * ct)
      if (r < MachEps) r = MachEps
      theta = acos((zpart(ipart) + r1 * ct) / r)
      ct   = cos(theta1 - theta)
      st   = sin(theta1 - theta)
      n(1) = nn(1) * ct - nn(2) * st
      n(2) = nn(1) * st + nn(2) * ct
      n(3) = 0._O
      tamp = sqrt(n(1)**2 + n(2)**2)
      if (tamp < MachEps) then
        print "(/,2x,'Error in subroutine matrix_Q1_DS_LAY in module Proces2:')"
        print "(  2x,'the module of the normal unit vector is zero;')"
        stop
      end if
      n(1) = n(1) / tamp
      n(2) = n(2) / tamp          
      if (indexC == 1) then      
        call MN_DS_LAY (1, kc, r, theta, ipart, Npart, Nrankpmax, Nrankp, zRe, zIm, &
             m, m1, n1)
      else     
        call MN_DS_LAY (3, kc, r, theta, ipart, Npart, Nrankpmax, Nrankp, zRe, zIm, &
             m, m1, n1)
      end if
      zl = cmplx(k * r,0.0,O)              
      call MN (1, zl, theta, m_minus, NrankG, NmaxG, m3, n3)
      fact = f * dA * pondere
      call matQ_COMP(m, Nrankpmax, NmaxG, ind_ref(ipart), fact, m3, n3, m1, n1,     &
           n, A, nap, map)       
    end do
  end do                     
  deallocate (m1, n1, m3, n3)     
end subroutine matrix_Q1_DS_LAY
!***********************************************************************************
subroutine incident_matrix_DS_LAY (TypeGeom, k, Nsurfmax, surf, Npart, Nrankpmax,   &
           Nrankp, zRe, zIm, zpart, m, NrankG, NmaxG, Nint, Nparammax, Nparam,      &
           Nintparam, paramG, weightsG, A, nap, map)
!-----------------------------------------------------------------------------------
! The routine computes the incident matrix of a layered particle using             !
! distributed sources.                                                             !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  integer    :: TypeGeom, m, Nsurfmax, Npart, Nrankpmax, Nrankp(Npart), NrankG,     &
                NmaxG, Nint, Nparammax, Nparam(Npart), Nintparam(Npart,Nparammax),  &
                nap, map
  real(O)    :: k, surf(Npart,Nsurfmax), zpart(Npart), paramG(Npart,Nparammax,Nint),&
                weightsG(Npart,Nparammax,Nint), zRe(Npart,Nrankpmax),               &
                zIm(Npart,Nrankpmax)
  complex(O) :: A(2*nap,2*map)
!      
  integer    :: i, j, pint, m_minus, iparam, Nintl, ipart, Nparaml
  real(O)    :: r, theta, phi, dA, param, pondere, r1, theta1, n(3), nn(3), tamp,   &
                ct, st
  complex(O) :: kc, fact, zl, f
  complex(O),allocatable :: m1(:,:), n1(:,:), m3(:,:), n3(:,:) 
!
  allocate (m1(3,NmaxG), n1(3,NmaxG), m3(3,Nrankpmax), n3(3,Nrankpmax))
  kc = cmplx(k,0.0,O)
  do i = 1, 2*Nrankpmax
    do j = 1, 2*NmaxG
      A(i,j) = zero
    end do
  end do 
  m_minus = - m 
  ipart   =   1
  Nparaml =   Nparam(ipart)
  f = im * 2._O * k * k
  do iparam = 1, Nparaml
    Nintl = Nintparam(ipart,iparam)
    do pint = 1, Nintl
      param   = paramG(ipart,iparam,pint)
      pondere = weightsG(ipart,iparam,pint)
      call elem_geomLAY (TypeGeom, Npart, ipart, Nsurfmax, surf, param, iparam,     &
           r1, theta1, phi, dA, nn)
      ct = cos(theta1)
      r  = sqrt(r1**2 + zpart(ipart)**2 + 2._O * r1 * zpart(ipart) * ct)
      if (r < MachEps) r = MachEps 
      theta = acos((zpart(ipart) + r1 * ct) / r)
      ct   = cos(theta1 - theta)
      st   = sin(theta1 - theta)
      n(1) = nn(1) * ct - nn(2) * st
      n(2) = nn(1) * st + nn(2) * ct
      n(3) = 0._O
      tamp = sqrt(n(1)**2 + n(2)**2)
      if (tamp < MachEps) then
        print "(/,2x,'Error in subroutine incident_matrix_DS_LAY in module Proces2:')"
        print "(  2x,'the module of the normal unit vector is zero;')"
        stop
      end if
      n(1) = n(1) / tamp
      n(2) = n(2) / tamp
      zl   = cmplx(k * r,0.0,O)      
      call MN (1, zl, theta, m, NrankG, NmaxG, m1, n1)                                                           
      call MN_DS_LAY (3, kc, r, theta, ipart, Npart, Nrankpmax, Nrankp, zRe, zIm,    &
           m_minus, m3, n3)                   
      fact = f * dA * pondere
      call matQinc_m(m, NmaxG, Nrankpmax, fact, m3, n3, m1, n1, n, A, nap, map)       
    end do
  end do
  deallocate (m1, n1, m3, n3)      
end subroutine incident_matrix_DS_LAY
