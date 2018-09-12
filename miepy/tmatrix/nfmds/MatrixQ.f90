! **********************************************************************************
! *          ROUTINES FOR COMPUTING THE Q MATRIX AND THE INCIDENT MATRIX           * 
! *                         AT AN INTEGRATION POINT                                *
! *    ------------------------------------------------------------------------    *
! *    Partial list of subroutines:                                                *
! *      matQ,        matQ_m,      matQ_COMP,     matQ_LAY,      matQ_LAY_a,       *
! *      matQ_LAY_b,  matQinc_m,   mvnv,          mvnv_m,        mvnvinc_m         *
! *    The routines for performing convergence tests are specified below.          *
! **********************************************************************************
subroutine matQ (NmaxC, NmaxL, chiral, perfectcond, ind_ref, fact, mv3, nv3,        &
           mv1, nv1, nuv, A, nap, map)
  use parameters
  implicit none
  integer    :: NmaxC, NmaxL, nap, map
  real(O)    :: nuv(3)
  complex(O) :: ind_ref, fact, mv1(3,NmaxC), nv1(3,NmaxC), mv3(3,NmaxL),            &
                nv3(3,NmaxL), A(2*nap,2*map)
  logical    :: chiral, perfectcond
!
  integer    :: i, j
  complex(O) :: mvl(3), nvl(3), mvc(3), nvc(3), v1, v2, mixt_product
!
  do i = 1, NmaxL
    mvl(1) = mv3(1,i)
    mvl(2) = mv3(2,i)
    mvl(3) = mv3(3,i)
    nvl(1) = nv3(1,i)
    nvl(2) = nv3(2,i)
    nvl(3) = nv3(3,i)
    do j = 1, NmaxC
      mvc(1) = mv1(1,j)
      mvc(2) = mv1(2,j)
      mvc(3) = mv1(3,j)
      nvc(1) = nv1(1,j)
      nvc(2) = nv1(2,j)
      nvc(3) = nv1(3,j)
      if (.not. chiral) then
        v1 = mixt_product (nuv, mvc, nvl)
        v2 = mixt_product (nuv, nvc, mvl)          
        if (.not. perfectcond) then
          A(i,j) = A(i,j) + (v1 + ind_ref * v2) * fact
        else 
          A(i,j) = A(i,j) + v2 * fact
        end if
      else 
        v1 = mixt_product (nuv, mvc, nvl)
        v2 = mixt_product (nuv, mvc, mvl)                           
        A(i,j) = A(i,j) + (v1 + ind_ref * v2) * fact
      end if
      if (.not. chiral) then
        v1 = mixt_product (nuv, nvc, nvl)
        v2 = mixt_product (nuv, mvc, mvl)
        if (.not. perfectcond) then
          A(i,j+NmaxC) = A(i,j+NmaxC) + (v1 + ind_ref * v2) * fact
        else 
          A(i,j+NmaxC) = A(i,j+NmaxC) + v2 * fact
        end if
      else 
        v1 = mixt_product (nuv, nvc, nvl)
        v2 = mixt_product (nuv, nvc, mvl)                   
        A(i,j+NmaxC) = A(i,j+NmaxC) + (v1 - ind_ref * v2) * fact
      end if 
      if (.not. chiral) then
        v1 = mixt_product (nuv, mvc, mvl)
        v2 = mixt_product (nuv, nvc, nvl)
        if (.not. perfectcond) then
          A(i+NmaxL,j) = A(i+NmaxL,j) + (v1 + ind_ref * v2) * fact
        else 
          A(i+NmaxL,j) = A(i+NmaxL,j) + v2 * fact
        end if
      else 
        v1 = mixt_product (nuv, mvc, mvl)
        v2 = mixt_product (nuv, mvc, nvl)                   
        A(i+NmaxL,j) = A(i+NmaxL,j) + (v1 + ind_ref * v2) * fact
      end if     
      if (.not. chiral) then
        v1 = mixt_product (nuv, nvc, mvl)
        v2 = mixt_product (nuv, mvc, nvl)
        if (.not. perfectcond) then
          A(i+NmaxL,j+NmaxC) = A(i+NmaxL,j+NmaxC) + (v1 + ind_ref * v2) * fact
        else 
          A(i+NmaxL,j+NmaxC) = A(i+NmaxL,j+NmaxC) + v2 * fact
        end if
      else 
        v1 = mixt_product (nuv, nvc, mvl)
        v2 = mixt_product (nuv, nvc, nvl)            
        A(i+NmaxL,j+NmaxC) = A(i+NmaxL,j+NmaxC) + (v1 - ind_ref * v2) * fact 
      end if
    end do        
  end do
end subroutine matQ
! **********************************************************************************
subroutine matQ_m (m, NmaxC, NmaxL, chiral, perfectcond, mirror, ind_ref, fact,     &
           mv3, nv3, mv1, nv1, nuv, A, nap, map)
  use parameters
  implicit none
  integer    :: m, NmaxC, NmaxL, nap, map
  real(O)    :: nuv(3)
  complex(O) :: ind_ref, fact, mv1(3,NmaxC), nv1(3,NmaxC), mv3(3,NmaxL),            &
                nv3(3,NmaxL), A(2*nap,2*map)
  logical    :: chiral, perfectcond, mirror
!
  integer    :: i, j, ni, nj
  real(O)    :: s
  complex(O) :: mvl(3), nvl(3), mvc(3), nvc(3), v1, v2, mixt_product
!
  do i = 1, NmaxL
    if (abs(m) == 0) then
      ni = i
    else
      ni = abs(m) + i - 1
    end if
    mvl(1) = mv3(1,i)
    mvl(2) = mv3(2,i)
    mvl(3) = mv3(3,i)
    nvl(1) = nv3(1,i)
    nvl(2) = nv3(2,i)
    nvl(3) = nv3(3,i)
    do j = 1, NmaxC
      if (abs(m) == 0) then
        nj = j
      else
        nj = abs(m) + j - 1
      end if
      mvc(1) = mv1(1,j)
      mvc(2) = mv1(2,j)
      mvc(3) = mv1(3,j)
      nvc(1) = nv1(1,j)
      nvc(2) = nv1(2,j)
      nvc(3) = nv1(3,j)
      s = (-1._O)**(ni + nj)
      if (.not. chiral) then
        v1 = mixt_product (nuv, mvc, nvl)
        v2 = mixt_product (nuv, nvc, mvl)
        if (.not. perfectcond) then
          A(i,j) = A(i,j) + (v1 + ind_ref * v2) * fact
          if (mirror) A(i,j) = A(i,j)+ s * (v1 + ind_ref * v2) * fact
        else 
          A(i,j) = A(i,j) + v2 * fact
          if (mirror) A(i,j) = A(i,j) + s * v2 * fact
        end if
      else 
        v1 = mixt_product (nuv, mvc, nvl)
        v2 = mixt_product (nuv, mvc, mvl)                   
        A(i,j) = A(i,j) + (v1 + ind_ref * v2) * fact                    
      end if                              
      if (.not. chiral) then
        if (m /= 0) then
          v1 = mixt_product (nuv, nvc, nvl)
          v2 = mixt_product (nuv, mvc, mvl)
          if (.not. perfectcond) then
            A(i,j+NmaxC) = A(i,j+NmaxC) + (v1 + ind_ref * v2) * fact
            if (mirror) A(i,j+NmaxC) = A(i,j+NmaxC) - s * (v1 + ind_ref * v2) * fact
          else 
            A(i,j+NmaxC) = A(i,j+NmaxC) + v2 * fact
            if (mirror) A(i,j+NmaxC) = A(i,j+NmaxC) - s * v2 * fact
          end if
          v1 = mixt_product (nuv, mvc, mvl)
          v2 = mixt_product (nuv, nvc, nvl)
          if (.not. perfectcond) then
            A(i+NmaxL,j) = A(i+NmaxL,j) + (v1 + ind_ref * v2) * fact
            if (mirror) A(i+NmaxL,j) = A(i+NmaxL,j) - s * (v1 + ind_ref * v2) * fact
          else 
            A(i+NmaxL,j) = A(i+NmaxL,j) + v2 * fact
            if (mirror) A(i+NmaxL,j) = A(i+NmaxL,j) - s * v2 * fact
          end if
        end if
      else 
        v1 = mixt_product (nuv, nvc, nvl)
        v2 = mixt_product (nuv, nvc, mvl)                         
        A(i,j+NmaxC) = A(i,j+NmaxC) + (v1 - ind_ref * v2) * fact
        v1 = mixt_product (nuv, mvc, mvl)
        v2 = mixt_product (nuv, mvc, nvl)                         
        A(i+NmaxL,j) = A(i+NmaxL,j) + (v1 + ind_ref * v2) * fact 
      end if
      if (.not. chiral) then 
        v1 = mixt_product (nuv, nvc, mvl)
        v2 = mixt_product (nuv, mvc, nvl)
        if (.not. perfectcond) then
          A(i+NmaxL,j+NmaxC) = A(i+NmaxL,j+NmaxC) + (v1 + ind_ref * v2) * fact
          if (mirror) A(i+NmaxL,j+NmaxC) = A(i+NmaxL,j+NmaxC) +                     &
                                          s * (v1 + ind_ref * v2) * fact
        else 
          A(i+NmaxL,j+NmaxC) = A(i+NmaxL,j+NmaxC) + v2 * fact
          if (mirror) A(i+NmaxL,j+NmaxC) = A(i+NmaxL,j+NmaxC) + s * v2 * fact
        end if
      else 
        v1 = mixt_product (nuv, nvc, mvl)
        v2 = mixt_product (nuv, nvc, nvl)                   
        A(i+NmaxL,j+NmaxC) = A(i+NmaxL,j+NmaxC) + (v1 - ind_ref * v2) * fact    
      end if                              
    end do   
  end do
end subroutine matQ_m
! **********************************************************************************
subroutine matQ_COMP (m, NmaxC, NmaxL, ind_ref, fact, mv3, nv3, mv1, nv1, nuv,      &
           A, nap, map)
  use parameters
  implicit none
  integer    :: m, NmaxC, NmaxL, nap, map
  real(O)    :: nuv(3)
  complex(O) :: ind_ref, fact, mv1(3,NmaxC), nv1(3,NmaxC), mv3(3,NmaxL),            &
                nv3(3,NmaxL), A(2*nap,2*map)
!
  integer    :: i, j
  complex(O) :: mvl(3), nvl(3), mvc(3), nvc(3), v1, v2, mixt_product
!
  do i = 1, NmaxL 
    mvl(1) = mv3(1,i)
    mvl(2) = mv3(2,i)
    mvl(3) = mv3(3,i)
    nvl(1) = nv3(1,i)
    nvl(2) = nv3(2,i)
    nvl(3) = nv3(3,i)
    do j = 1, NmaxC
      mvc(1) = mv1(1,j)
      mvc(2) = mv1(2,j)
      mvc(3) = mv1(3,j)
      nvc(1) = nv1(1,j)
      nvc(2) = nv1(2,j)
      nvc(3) = nv1(3,j)
      v1 = mixt_product (nuv, mvc, nvl)
      v2 = mixt_product (nuv, nvc, mvl)
      A(i,j) = A(i,j) + (v1 + ind_ref  *v2) * fact
      if (m /= 0) then
        v1 = mixt_product (nuv, nvc, nvl)
        v2 = mixt_product (nuv, mvc, mvl)
        A(i,j+NmaxC) = A(i,j+NmaxC) + (v1 + ind_ref * v2) * fact
        v1 = mixt_product (nuv, mvc, mvl)
        v2 = mixt_product (nuv, nvc, nvl)
        A(i+NmaxL,j) = A(i+NmaxL,j) + (v1 + ind_ref * v2) * fact
      end if
      v1 = mixt_product (nuv, nvc, mvl)
      v2 = mixt_product (nuv, mvc, nvl)
      A(i+NmaxL,j+NmaxC) = A(i+NmaxL,j+NmaxC) + (v1 + ind_ref * v2) * fact
    end do
  end do
end subroutine matQ_COMP
! **********************************************************************************
subroutine matQ_LAY (m, ipart, Npartlim, NmaxC, NmaxL, Nmax, index, fact, mv3, nv3, &
           mv2, nv2, mv1, nv1, nuv, a, nap, map, b, nbp, mbp)
  use parameters
  implicit none
  integer    :: m, ipart, Npartlim, NmaxC, NmaxL, Nmax, nap, map, nbp, mbp
  real(O)    :: nuv(3)
  complex(O) :: index, fact, mv1(3,Nmax), nv1(3,Nmax), mv2(3,Nmax),                 &
                nv2(3,Nmax), mv3(3,Nmax), nv3(3,Nmax), a(2*nap,2*map),              &
                b(2*nbp,2*mbp)
!
  integer    :: i, j
  complex(O) :: mvl(3), nvl(3), mvc(3), nvc(3), v1, v2, mixt_product
!
  do i = 1, NmaxL
    mvl(1) = mv3(1,i)
    mvl(2) = mv3(2,i)
    mvl(3) = mv3(3,i)
    nvl(1) = nv3(1,i)
    nvl(2) = nv3(2,i)
    nvl(3) = nv3(3,i)
    do j = 1, NmaxC
      mvc(1) = mv1(1,j)
      mvc(2) = mv1(2,j)
      mvc(3) = mv1(3,j)
      nvc(1) = nv1(1,j)
      nvc(2) = nv1(2,j)
      nvc(3) = nv1(3,j)
      v1 = mixt_product (nuv, mvc, nvl)
      v2 = mixt_product (nuv, nvc, mvl)
      a(i,j) = a(i,j) + (v1 + index * v2) * fact
      if (m /= 0) then
        v1 = mixt_product (nuv, nvc, nvl)
        v2 = mixt_product (nuv, mvc, mvl)
        a(i,j+NmaxC) = a(i,j+NmaxC) + (v1 + index * v2) * fact
        v1 = mixt_product (nuv, mvc, mvl)
        v2 = mixt_product (nuv, nvc, nvl)
        a(i+NmaxL,j) = a(i+NmaxL,j) + (v1 + index * v2) * fact
      end if
      v1 = mixt_product (nuv, nvc, mvl)
      v2 = mixt_product (nuv, mvc, nvl)
      a(i+NmaxL,j+NmaxC) = a(i+NmaxL,j+NmaxC) + (v1 + index * v2) * fact
    end do
    if (ipart < Npartlim) then
      do j = 1, NmaxC
        mvc(1) = mv2(1,j)
        mvc(2) = mv2(2,j)
        mvc(3) = mv2(3,j)
        nvc(1) = nv2(1,j)
        nvc(2) = nv2(2,j)
        nvc(3) = nv2(3,j)
        v1 = mixt_product (nuv, mvc, nvl)
        v2 = mixt_product (nuv, nvc, mvl)
        b(i,j) = b(i,j) + (v1 + index * v2) * fact
        if (m /= 0) then
          v1 = mixt_product (nuv, nvc, nvl)
          v2 = mixt_product (nuv, mvc, mvl)
          b(i,j+NmaxC) = b(i,j+NmaxC) + (v1 + index * v2) * fact
          v1 = mixt_product (nuv, mvc, mvl)
          v2 = mixt_product (nuv, nvc, nvl)
          b(i+NmaxL,j) = b(i+NmaxL,j) + (v1 + index * v2) * fact
        end if
        v1 = mixt_product (nuv, nvc, mvl)
        v2 = mixt_product (nuv, mvc, nvl)
        b(i+NmaxL,j+NmaxC) = b(i+NmaxL,j+NmaxC) + (v1 + index * v2) * fact
      end do
    end if
  end do
end subroutine matQ_LAY
! **********************************************************************************
subroutine matQ_LAY_a (m, NmaxC, NmaxL, Nmax, index, fact, mv3, nv3, mv1, nv1,      &
           nuv, a, nap, map)
  use parameters
  implicit none
  integer    :: m, NmaxC, NmaxL, Nmax, nap, map
  real(O)    :: nuv(3)
  complex(O) :: index, fact, mv1(3,Nmax), nv1(3,Nmax), mv3(3,Nmax), nv3(3,Nmax),    &
                a(2*nap,2*map)                              
!
  integer    :: i, j
  complex(O) :: mvl(3), nvl(3), mvc(3), nvc(3), v1, v2, mixt_product
!
  do i = 1, NmaxL
    mvl(1) = mv3(1,i)
    mvl(2) = mv3(2,i)
    mvl(3) = mv3(3,i)
    nvl(1) = nv3(1,i)
    nvl(2) = nv3(2,i)
    nvl(3) = nv3(3,i)
    do j = 1, NmaxC
      mvc(1) = mv1(1,j)
      mvc(2) = mv1(2,j)
      mvc(3) = mv1(3,j)
      nvc(1) = nv1(1,j)
      nvc(2) = nv1(2,j)
      nvc(3) = nv1(3,j)
      v1 = mixt_product (nuv, mvc, nvl)
      v2 = mixt_product (nuv, nvc, mvl)
      a(i,j) = a(i,j) + (v1 + index * v2) * fact
      if (m /= 0) then
        v1 = mixt_product (nuv, nvc, nvl)
        v2 = mixt_product (nuv, mvc, mvl)
        a(i,j+NmaxC) = a(i,j+NmaxC) + (v1 + index * v2) * fact
        v1 = mixt_product (nuv, mvc, mvl)
        v2 = mixt_product (nuv, nvc, nvl)
        a(i+NmaxL,j) = a(i+NmaxL,j) + (v1 + index * v2) * fact
      end if
      v1 = mixt_product (nuv, nvc, mvl)
      v2 = mixt_product (nuv, mvc, nvl)
      a(i+NmaxL,j+NmaxC) = a(i+NmaxL,j+NmaxC) + (v1 + index * v2) * fact
    end do    
  end do
end subroutine matQ_LAY_a
! **********************************************************************************
subroutine matQ_LAY_b (m, ipart, Npartlim, NmaxC, NmaxL, Nmax, index, fact, mv3,    &
           nv3, mv2, nv2, nuv, b, nbp, mbp)
  use parameters
  implicit none
  integer    :: m, ipart, Npartlim, NmaxC, NmaxL, Nmax, nbp, mbp
  real(O)    :: nuv(3)
  complex(O) :: index, fact, mv2(3,Nmax), nv2(3,Nmax), mv3(3,Nmax), nv3(3,Nmax),    &
                b(2*nbp,2*mbp)
!
  integer    :: i, j
  complex(O) :: mvl(3), nvl(3), mvc(3), nvc(3), v1, v2, mixt_product
!
  do i = 1, NmaxL
    mvl(1) = mv3(1,i)
    mvl(2) = mv3(2,i)
    mvl(3) = mv3(3,i)
    nvl(1) = nv3(1,i)
    nvl(2) = nv3(2,i)
    nvl(3) = nv3(3,i)    
    if (ipart < Npartlim) then
      do j = 1, NmaxC
        mvc(1) = mv2(1,j)
        mvc(2) = mv2(2,j)
        mvc(3) = mv2(3,j)
        nvc(1) = nv2(1,j)
        nvc(2) = nv2(2,j)
        nvc(3) = nv2(3,j)
        v1 = mixt_product (nuv, mvc, nvl)
        v2 = mixt_product (nuv, nvc, mvl)
        b(i,j) = b(i,j) + (v1 + index * v2) * fact
        if (m /= 0) then
          v1 = mixt_product (nuv, nvc, nvl)
          v2 = mixt_product (nuv, mvc, mvl)
          b(i,j+NmaxC) = b(i,j+NmaxC) + (v1 + index * v2) * fact
          v1 = mixt_product (nuv, mvc, mvl)
          v2 = mixt_product (nuv, nvc, nvl)
          b(i+NmaxL,j) = b(i+NmaxL,j) + (v1 + index * v2) * fact
        end if
        v1 = mixt_product (nuv, nvc, mvl)
        v2 = mixt_product (nuv, mvc, nvl)
        b(i+NmaxL,j+NmaxC) = b(i+NmaxL,j+NmaxC) + (v1 + index * v2) * fact
      end do
    end if
  end do
end subroutine matQ_LAY_b
!***********************************************************************************
subroutine matQinc_m (m, NmaxC, NmaxL, fact, mv3, nv3, mv1, nv1, nuv, A, nap, map)
  use parameters
  implicit none
  integer    :: m, NmaxL, NmaxC, nap, map
  real(O)    :: nuv(3)
  complex(O) :: fact, mv1(3,NmaxC), nv1(3,NmaxC), mv3(3,NmaxL), nv3(3,NmaxL),       &
                A(2*nap,2*map)
!
  integer    :: i, j
  complex(O) :: mvl(3), nvl(3), mvc(3), nvc(3), v1, v2, mixt_product
!
  do i = 1, NmaxL
    mvl(1) = mv3(1,i)
    mvl(2) = mv3(2,i)
    mvl(3) = mv3(3,i)
    nvl(1) = nv3(1,i)
    nvl(2) = nv3(2,i)
    nvl(3) = nv3(3,i)
    do j = 1, NmaxC
      mvc(1) = mv1(1,j)
      mvc(2) = mv1(2,j)
      mvc(3) = mv1(3,j)
      nvc(1) = nv1(1,j)
      nvc(2) = nv1(2,j)
      nvc(3) = nv1(3,j)
      v1 = mixt_product (nuv, mvc, nvl)
      v2 = mixt_product (nuv, nvc, mvl)
      A(i,j) = A(i,j) + (v1 + v2) * fact
      if (m /= 0) then
        v1 = mixt_product (nuv, nvc, nvl)
        v2 = mixt_product (nuv, mvc, mvl)
        A(i,j+NmaxC) = A(i,j+NmaxC) + (v1 + v2) * fact
        v1 = mixt_product (nuv, mvc, mvl)
        v2 = mixt_product (nuv, nvc, nvl)
        A(i+NmaxL,j) = A(i+NmaxL,j) + (v1 + v2) * fact
      end if
      v1 = mixt_product (nuv, nvc, mvl)
      v2 = mixt_product (nuv, mvc, nvl)
      A(i+NmaxL,j+NmaxC) = A(i+NmaxL,j+NmaxC) + (v1 + v2) * fact
    end do
  end do           
end subroutine matQinc_m
! ***********************************************************************************
! *     THE FOLLOWING ROUTINES ARE IMPORTANT FOR THE SPEED AND FLEXIBILITY          *
! *     OF THE PROGRAMS. WE SACRIFICE THE SPEED IN FAVOUR OF THE FLEXIBILTY,        *
! *     CLARITY AND SIMPLICITY OF THE COMPUTER CODE. THE SPEED CAN BE INCREASED     *
! *     IF THE VECTOR SPHERICAL WAVE FUNCTIONS ARE EXPLICITELY COMPUTED.            *
! ***********************************************************************************
subroutine mvnv (index1, index2, chiral, zl, zc, zcl, zcr, theta, phi, Mrank,       &
           Nrank, Nmax, NmaxC, NmaxL, mv3, nv3, mv1, nv1)
  use parameters
  implicit none
  integer    :: index1, index2, Mrank, Nrank, Nmax, NmaxC, NmaxL
  real(O)    :: theta, phi
  complex(O) :: zc, zl, zcl, zcr, mv1(3,NmaxC), nv1(3,NmaxC), mv3(3,NmaxL),         &
                nv3(3,NmaxL)
  logical    :: chiral
!  
  complex(O),allocatable :: mvl1(:,:), nvl1(:,:), mvr1(:,:), nvr1(:,:)                     
!                  
  if (chiral) allocate (mvl1(3,NmaxC), nvl1(3,NmaxC), mvr1(3,NmaxC), nvr1(3,NmaxC))                                 
  if (index1 == 3 .and. index2 == 1) then
    call MN_complete (3, zl, theta, phi, Mrank, Nrank, Nmax,.false.,.false.,        &
         mv3, nv3)
    if (.not. chiral) then
      call MN_complete (1, zc, theta, phi, Mrank, Nrank, Nmax,.true.,.false.,       &
           mv1, nv1)
    else 
      call MN_complete (1, zcl, theta, phi, Mrank, Nrank, Nmax,.true.,.false.,      &
           mvl1, nvl1)
      call MN_complete (1, zcr, theta, phi, Mrank, Nrank, Nmax,.true.,.false.,      &
           mvr1, nvr1)
      call MN_left_right (Nmax, mvl1, nvl1, mvr1, nvr1, mv1, nv1)
    end if
  else if (index1 == 3 .and. index2 == 3) then
    call MN_complete (3, zl, theta, phi, Mrank, Nrank, Nmax,.false.,.false.,        &
         mv3, nv3)
    call MN_complete (3, zc, theta, phi, Mrank, Nrank, Nmax,.true.,.false.,         &
         mv1, nv1)
  else if (index1 == 1 .and. index2 == 1) then
    call MN_complete (1, zl, theta, phi, Mrank, Nrank, Nmax,.false.,.false.,        &
         mv3, nv3)
    if (.not. chiral) then
      call MN_complete (1, zc, theta, phi, Mrank, Nrank, Nmax,.true.,.false.,       &
           mv1, nv1)
    else 
      call MN_complete (1, zcl, theta, phi, Mrank, Nrank, Nmax,.true.,.false.,      &
           mvl1, nvl1)
      call MN_complete (1, zcr, theta, phi, Mrank, Nrank, Nmax,.true.,.false.,      &
           mvr1, nvr1)
      call MN_left_right (Nmax, mvl1, nvl1, mvr1, nvr1, mv1, nv1)
    end if 
  else if (index1 == 1 .and. index2 == 3) then
    call MN_complete (1, zl, theta, phi, Mrank, Nrank, Nmax,.false.,.false.,        &
         mv3, nv3)
    call MN_complete (3, zc, theta, phi, Mrank, Nrank, Nmax,.true.,.false.,         &
         mv1, nv1)
  end if    
  if(chiral) deallocate (mvl1, nvl1, mvr1, nvr1)  
end subroutine mvnv 
! **********************************************************************************
subroutine mvnv_m (index1, index2, chiral, DS, zl, zc, zcl, zcr, kc, ki, kil, kir,  &
           r, theta, m, Nrank, Nmax, NmaxC, NmaxL, zRe, zIm, mv3, nv3, mv1, nv1)
  use parameters
  implicit none
  integer    :: index1, index2, m, m_minus, Nrank, Nmax, NmaxC, NmaxL
  real(O)    :: zRe(Nrank), zIm(Nrank), r, theta
  complex(O) :: kc, ki, kil, kir, zc, zl, zcl, zcr, mv1(3,NmaxC), nv1(3,NmaxC),     &
                mv3(3,NmaxL), nv3(3,NmaxL)
  logical    :: DS, chiral
!  
  complex(O),allocatable :: mvl1(:,:), nvl1(:,:), mvr1(:,:), nvr1(:,:)                     
!
  m_minus = - m            
  if (chiral) allocate (mvl1(3,NmaxC), nvl1(3,NmaxC), mvr1(3,NmaxC), nvr1(3,NmaxC))
  if (.not. DS) then      
    if (index1 == 3 .and. index2 == 1) then
      call MN (3, zl, theta, m_minus, Nrank, Nmax, mv3, nv3)
      if (.not. chiral) then
        call MN (1, zc, theta, m, Nrank, Nmax, mv1, nv1)           
      else
        call MN (1, zcl, theta, m, Nrank, Nmax, mvl1, nvl1)
        call MN (1, zcr, theta, m, Nrank, Nmax, mvr1, nvr1)
        call MN_left_right (Nmax, mvl1, nvl1, mvr1, nvr1, mv1, nv1)
      end if                        
    else if (index1 == 3 .and. index2 == 3) then
      call MN (3, zl, theta, m_minus, Nrank, Nmax, mv3, nv3)
      call MN (3, zc, theta, m, Nrank, Nmax, mv1, nv1)               
    else if (index1 == 1 .and. index2 == 1) then
      call MN (1, zl, theta, m_minus, Nrank, Nmax, mv3, nv3)
      if (.not. chiral) then
        call MN (1, zc, theta, m, Nrank, Nmax, mv1, nv1)             
      else 
        call MN (1, zcl, theta, m, Nrank, Nmax, mvl1, nvl1)
        call MN (1, zcr, theta, m, Nrank, Nmax, mvr1, nvr1)
        call MN_left_right (Nmax, mvl1, nvl1, mvr1, nvr1, mv1, nv1)
      end if
    else if (index1 == 1 .and. index2 == 3) then
      call MN (1, zl, theta, m_minus, Nrank, Nmax, mv3, nv3)
      call MN (3, zc, theta, m, Nrank, Nmax, mv1, nv1)       
    end if
  else 
    if (index1 == 3 .and. index2 == 1) then
      call MN_DS (3, kc, r, theta, zRe, zIm, m_minus, Nrank, mv3, nv3)
      if (.not. chiral) then
        call MN_DS (1, ki, r, theta, zRe, zIm, m, Nrank, mv1, nv1)             
      else 
        call MN_DS (1, kil, r, theta, zRe, zIm, m, Nrank, mvl1, nvl1)
        call MN_DS (1, kir, r, theta, zRe, zIm, m, Nrank, mvr1, nvr1)
        call MN_left_right (Nrank, mvl1, nvl1, mvr1, nvr1, mv1, nv1)
      end if                            
    else if (index1 == 1 .and. index2 == 1) then         
      call MN (1, zl, theta, m_minus, Nrank, Nmax, mv3, nv3)
      if (.not. chiral) then
        call MN_DS (1, ki, r, theta, zRe, zIm, m, Nrank, mv1, nv1)       
      else 
        call MN_DS (1, kil, r, theta, zRe, zIm, m, Nrank, mvl1, nvl1)
        call MN_DS (1, kir, r, theta, zRe, zIm, m, Nrank, mvr1, nvr1)
        call MN_left_right (Nrank, mvl1, nvl1, mvr1, nvr1, mv1, nv1)
      end if               
    end if
  end if                            
  if (chiral) deallocate (mvl1, nvl1, mvr1, nvr1)  
end subroutine mvnv_m  
!***********************************************************************************
subroutine mvnvinc_m (zl, kc, r, theta, m, Nrank, Nmax, zRe, zIm, mv3, nv3,        &
           mv1, nv1)
  use parameters
  implicit none
  integer    :: m, m_minus, Nrank, Nmax
  real(O)    :: zRe(Nrank), zIm(Nrank), r, theta
  complex(O) :: kc, zl, mv1(3,Nmax), nv1(3,Nmax), mv3(3,Nrank), nv3(3,Nrank)
!   
  m_minus = - m
  call MN (1, zl, theta, m, Nrank, Nmax, mv1, nv1)
  call MN_DS (3, kc, r, theta, zRe, zIm, m_minus, Nrank, mv3, nv3) 
end subroutine mvnvinc_m
! **********************************************************************************
! *      ROUTINES FOR PERFORMING CONVERGENCE TESTS USING MISHCHENKO'S PROCEDURE    * 
! *    ------------------------------------------------------------------------    *
! *    Partial list of subroutines:                                                *
! *      estimateNrankMishchenko,           estimateNintMishchenko,                *
! *      CscatCextConvLS,                   CscatCextConvDS,                       *
! *      matQConv,                          matQincConv,                           *
! **********************************************************************************
recursive subroutine estimateNrankMishchenko (TypeGeom, k, ind_ref, Nsurf, surf,    &
           zRe, zIm, Nparam, mirror, perfectcond, DS, ComplexPlane, EpsZIm,         &
           x, delta, Ndgs, Nint, Nrank, NrankMax, Cscat1, Cext1)
  use parameters
  use derived_parameters
  implicit none
  integer    :: Nsurf, Nparam, TypeGeom, Nint, Nrank, Ndgs
  real(O)    :: k, surf(Nsurf), zRe(NrankPD), zIm(NrankPD), EpsZIm, x, delta,       &
                Cscat1, Cext1
  complex(O) :: ind_ref
  logical    :: mirror, perfectcond, DS, ComplexPlane 
! 
  integer    :: irank, i, NrankMin, NrankWis, NrankMax
  real(O)    :: Cext, Cscat, DeltaScat, DeltaExt, ErrScat, ErrExt, ErrMax
  logical    :: diverge, repeat
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG(:,:), weightsG(:,:)
!
  delta = 1.e-3_O
  print "(/,2x, a, 1pe10.2, a)",                                                    &
 '- enter the accuracy of the T-matrix calculations as for instance', delta, ';'
  call read_real (delta)  
  print "(/,2x,'- enter the estimated value of  Ndgs, where the number of integration ')"
  print "(  2x,'points for the Nrank convergence test is Nint = Ndgs * Nrank;         ')"  
  call read_integer (Ndgs)         
  NrankWis = int(x + 4.05_O * x**0.33_O)
  NrankMin = max(4,NrankWis)        
  print "(/,2x, a, i3, a)",                                                         &
  'The estimate of the lower bound of Nrank from Wiscombe''s criterion is ',        &
  NrankMin, ';'  
  print "(/,2x,'- enter the lower and upper bounds of Nrank;')"  
  call read_integer2 (NrankMin, NrankMax)       
  print "( )"
  print "(/,2x,'Convergence Test over Nrank')"
  print "(2x,'iteration',1x,'Nrank',3x,'Nint',6x,'EpsScat',8x,'EpsExt')"
  Cscat1  = 0._O
  Cext1   = 0._O
  diverge = .true.
  irank   = NrankMin - 1
  i = 0
  do while (diverge .and. irank < NrankMax)
    i     = i + 1
    irank = irank + 1
    Nrank = irank
    Nint  = Ndgs * Nrank
    allocate (paramG(Nparam,Nint), weightsG(Nparam,Nint), Nintparam(Nparam))
    call interpolation_listAXSYM (TypeGeom, Nsurf, surf, Nint, Nparam,              &
         Nintparam, paramG, weightsG, mirror)
    if (DS) then  
      call check_MaxNrank (Nrank)    
      call zDSAXSYM (TypeGeom, Nsurf, surf, Nrank, ComplexPlane, EpsZIm, zRe, zIm)
      call CscatCextConvDS (TypeGeom, k, ind_ref, Nsurf, surf, zRe, zIm, Nrank,     &
           Nint, Nparam, Nintparam, paramG, weightsG, perfectcond, Cscat, Cext)
    else
      call CscatCextConvLS (TypeGeom, k, ind_ref, Nsurf, surf, Nrank, Nint,         &
           Nparam, Nintparam, paramG, weightsG, mirror, perfectcond, Cscat, Cext)
    end if       
    DeltaScat = abs(Cscat1 - Cscat)
    DeltaExt  = abs(Cext1  - Cext )
    Cscat1    = Cscat
    Cext1     = Cext
    print "(4x,i3,2(5x,i3),2(2x,1pe13.4))", i, Nrank, Nint, DeltaScat, DeltaExt
    if (abs(Cscat) > MachEps .and. abs(Cext) > MachEps) then          
      ErrScat = DeltaScat / abs(Cscat) 
      ErrExt  = DeltaExt  / abs(Cext)
      ErrMax  = max(ErrScat, ErrExt)
      if (ErrMax <= delta) diverge = .false.     
    else
      exit
    end if
    deallocate (paramG, weightsG, Nintparam)                
  end do
  print "( )"
  if (diverge) then
    print "(/,2x,'Convergence over Nrank is not achieved;')"
    print "(  2x,'- repeat the convergence test with new inputs ?:')"
    print "(  1x,'.true. - yes, .false. - no;')"
    call read_logical (repeat)
    if (repeat) then
      call estimateNrankMishchenko (TypeGeom, k, ind_ref, Nsurf, surf,              &
           zRe, zIm, Nparam, mirror, perfectcond, DS, ComplexPlane, EpsZIm,         &
           x, delta, Ndgs, Nint, Nrank, NrankMax, Cscat1, Cext1)                                 
    end if
  end if  
end subroutine estimateNrankMishchenko
! **********************************************************************************
recursive subroutine estimateNintMishchenko (TypeGeom, k, ind_ref, Nsurf, surf,     &
          zRe, zIm, Nparam, mirror, perfectcond, DS, x, delta, Ndgs, Nint, dNint,   &
          Nrank, NrankMax, Cscat1, Cext1)
  use parameters
  use derived_parameters
  implicit none
  integer    :: Nsurf, Nparam, TypeGeom, Nint, dNint, Nrank, NrankMax, Ndgs
  real(O)    :: k, surf(Nsurf), zRe(NrankPD), zIm(NrankPD), x, delta, Cscat1, Cext1
  complex(O) :: ind_ref
  logical    :: mirror, perfectcond, DS
! 
  integer    :: iint, i, NintMax  
  real(O)    :: Cext, Cscat, DeltaScat, DeltaExt, ErrScat, ErrExt
  logical    :: diverge, repeat
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG(:,:), weightsG(:,:)
!                
  print "(/,2x, a, i3, a)",                                                         &
 '- enter an integer value greater than ', Nint + 1, ' for the upper bound of Nint'
  print "(2x,'as for instance ',i4,';')", (Ndgs + 1) * NrankMax
  call read_integer (NintMax)  
  print "( )"
  print "(/,2x,'Convergence Test over Nint')"
  print "(2x,'iteration',1x,'Nrank',3x,'Nint',6x,'EpsScat',8x,'EpsExt')" 
  diverge = .true.
  iint    = Nint
  i = 0
  do while (diverge .and. iint < NintMax)
    i    = i + 1
    iint = iint + dNint
    Nint = iint
    allocate (paramG(Nparam,Nint), weightsG(Nparam,Nint), Nintparam(Nparam))
    call interpolation_listAXSYM (TypeGeom, Nsurf, surf, Nint, Nparam,              &
         Nintparam, paramG, weightsG, mirror)    
    if (DS) then
      call CscatCextConvDS (TypeGeom, k, ind_ref, Nsurf, surf, zRe, zIm, Nrank,     &
           Nint, Nparam, Nintparam, paramG, weightsG, perfectcond, Cscat, Cext)
    else
      call CscatCextConvLS (TypeGeom, k, ind_ref, Nsurf, surf, Nrank, Nint,         &
           Nparam, Nintparam, paramG, weightsG, mirror, perfectcond, Cscat, Cext)
    end if  
    DeltaScat = abs(Cscat1 - Cscat)
    DeltaExt  = abs(Cext1  - Cext )
    Cscat1    = Cscat
    Cext1     = Cext        
    print "(4x,i3,2(5x,i3),2(2x,1pe13.4))", i, Nrank, Nint, DeltaScat, DeltaExt
    if (abs(Cscat) > MachEps .and. abs(Cext) > MachEps) then
      ErrScat = DeltaScat / abs(Cscat) 
      ErrExt  = DeltaExt  / abs(Cext)                 
      if (ErrScat <= delta .and. ErrExt <= delta) diverge = .false.
    else
      exit
    end if      
    deallocate (paramG, weightsG, Nintparam)
  end do  
  print "( )"
  if (diverge) then
    print "(/,2x,'Convergence over Nint is not achieved;')"
    print "(2x,'- repeat the convergence test with new inputs ?:')"
    print "(1x,'.true. - yes, .false. - no;')"
    call read_logical (repeat)
    if (repeat) then
      call estimateNintMishchenko (TypeGeom, k, ind_ref, Nsurf, surf,               &
           zRe, zIm, Nparam, mirror, perfectcond, DS, x, delta, Ndgs, Nint,         &
           dNint, Nrank, NrankMax, Cscat1, Cext1)                            
    end if
  else
    print "(/,2x,'Nint and Nrank estimates:')"
    print "(  2x,'the estimated values of Nint and Nrank from Mishchenko''s convergence')"
    print "(  2x,'criterion are: Nint = ',i4,' and Nrank = ',i3,';')", Nint, Nrank      
  end if  
end subroutine estimateNintMishchenko    
!***********************************************************************************
subroutine CscatCextConvLS  (TypeGeom, k, ind_ref, Nsurf, surf, Nrank, Nint,        &
           Nparam, Nintparam, paramG, weightsG, mirror, perfectcond, Cscat, Cext)
  use parameters
  use derived_parameters
  implicit none
  integer    :: Nrank, Nsurf, Nint, Nparam, Nintparam(Nparam), TypeGeom
  real(O)    :: k, surf(Nsurf), paramG(Nparam,Nint), weightsG(Nparam,Nint),         &
                Cscat, Cext
  complex(O) :: ind_ref
  logical    :: mirror, perfectcond
!      
  integer    :: i, j, pint, iparam, Nintl, n
  real(O)    :: r, theta, phi, dA, param, pondere, nuv(3), ir, QA, QB,              &
                nm, nr, ft, fl 
  complex(O) :: ki, zc, zl, factp, factm, f 
  real(O),allocatable    :: Pnm(:), dPnm(:), pinm(:), taunm(:)
  complex(O),allocatable :: jb(:), jbd(:), jba(:), jbda(:), hb(:), hbd(:),          &
                            Q31A(:,:), Q31B(:,:), Q11A(:,:), Q11B(:,:),             &
                            mv1(:,:), nv1(:,:), mv3(:,:), nv3(:,:),                 &
                            mv1a(:,:), nv1a(:,:)
! 
  allocate (Q31A(Nrank,Nrank), Q31B(Nrank,Nrank), Q11A(Nrank,Nrank),                &
            Q11B(Nrank,Nrank))  
  allocate (mv1(3,Nrank), nv1(3,Nrank), mv3(3,Nrank), nv3(3,Nrank),                 &
            mv1a(3,Nrank), nv1a(3,Nrank))
  allocate (jb(0:Nrank), jbd(0:Nrank), jba(0:Nrank), jbda(0:Nrank),                 &     
            hb(0:Nrank), hbd(0:Nrank)) 
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank)) 
  do i = 1, Nrank
    do j = 1, Nrank
      Q31A(i,j) = zero
      Q31B(i,j) = zero
      Q11A(i,j) = zero
      Q11B(i,j) = zero
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
      if (abs(zc) < MachEps) zc = cmplx(MachEps,MachEps,O) 
      zl = cmplx(k * r,0.0,O)                    
      if (abs(zl) < MachEps) zl = cmplx(MachEps,MachEps,O)                                 
      call leg_normalized (theta, 0, Nrank, Pnm, dPnm, pinm, taunm)
      call besel_j (zc, Nrank, jb, jbd) 
      call besel_j (zl, Nrank, jba, jbda) 
      call besel_h (zl, Nrank, hb, hbd)     
      do n = 1, Nrank                
        nr = real(n * (n + 1),O)
        nm = 1._O / sqrt(2._O * nr)     
        ft = taunm(n) * nm
        fl = nr * Pnm(n) * nm
        mv1(1,n)  =   zero
        mv1(2,n)  =   zero   
        mv1(3,n)  = - jb(n) * ft         
        nv1(1,n)  =   jb(n) * fl / zc
        nv1(2,n)  =   jbd(n) * ft/ zc
        nv1(3,n)  =   zero
        mv1a(1,n) =   zero
        mv1a(2,n) =   zero
        mv1a(3,n) = - jba(n)  * ft         
        nv1a(1,n) =   jba(n)  * fl / zl
        nv1a(2,n) =   jbda(n) * ft / zl
        nv1a(3,n) =   zero
        mv3(1,n)  =   zero
        mv3(2,n)  =   zero
        mv3(3,n)  = - hb(n)  * ft         
        nv3(1,n)  =   hb(n)  * fl / zl
        nv3(2,n)  =   hbd(n) * ft / zl
        nv3(3,n)  =   zero
      end do                                                                                                                
      factm = - f * dA * pondere 
      call matQConv (Nrank, perfectcond, mirror, ind_ref, factm,                     &
           mv3, nv3, mv1, nv1, nuv, Q31A, Q31B, Nrank, Nrank)
      factp = - factm 
      call matQConv (Nrank, perfectcond, mirror, ind_ref, factp,                     &
           mv1a, nv1a, mv1, nv1, nuv, Q11A, Q11B, Nrank, Nrank)
    end do
  end do                  
  call LU_SYSTEM (Q31A, Nrank, Nrank, Q11A, Nrank, Nrank, Nrank)
  call LU_SYSTEM (Q31B, Nrank, Nrank, Q11B, Nrank, Nrank, Nrank)            
  Cscat = 0._O
  Cext  = 0._O
  do n = 1, Nrank
    ir    = real(2 * n + 1,O)
    QA    = abs(Q11A(n,n))
    QB    = abs(Q11B(n,n))
    Cscat = Cscat + ir * (QA * QA + QB * QB)
    QA    = real(Q11A(n,n))
    QB    = real(Q11B(n,n))
    Cext  = Cext + ir * (QA + QB)  
  end do
  deallocate (Q31A, Q31B, Q11A, Q11B) 
  deallocate (mv1, nv1, mv3, nv3, mv1a, nv1a)
  deallocate (jb, jbd, jba, jbda, hb, hbd)
  deallocate (Pnm, dPnm, pinm, taunm)
end subroutine CscatCextConvLS
!***********************************************************************************
subroutine CscatCextConvDS (TypeGeom, k, ind_ref, Nsurf, surf, zRe, zIm, Nrank,     &
           Nint, Nparam, Nintparam, paramG, weightsG, perfectcond, Cscat, Cext)
  use parameters
  use derived_parameters
  implicit none
  integer    :: Nrank, Nsurf, Nint, Nparam, Nintparam(Nparam), TypeGeom
  real(O)    :: k, surf(Nsurf), paramG(Nparam,Nint), weightsG(Nparam,Nint),         &
                zRe(Nrank), zIm(Nrank), Cscat, Cext
  complex(O) :: ind_ref
  logical    :: perfectcond
!      
  integer    :: i, j, pint, iparam, Nintl, n
  real(O)    :: r, theta, phi, dA, param, pondere, nuv(3), ir, QA, QB, ro, z, nm,   &
                nr, cth, sth, ftr, flr, dz
  complex(O) :: kc, ki, zl, factp, factm, f, RR, sint, cost, argJi, argJc,          &
                sinc, cosc, Pmm, taumm, ft, fl, factt, factl, roc, dzc
  real(O),allocatable    :: Pnm(:), dPnm(:), pinm(:), taunm(:)
  complex(O),allocatable :: jb(:), jbd(:), jba(:), jbda(:), hb(:), hbd(:),          &
                            Q31A(:,:), Q31B(:,:), Q11A(:,:), Q11B(:,:),             &
                            QincA(:,:), QincB(:,:), mv1(:,:), nv1(:,:),             &
                            mv3(:,:), nv3(:,:), mv1a(:,:), nv1a(:,:) 
! 
  allocate (Q31A(Nrank,Nrank), Q31B(Nrank,Nrank),  Q11A(Nrank,Nrank),               &
            Q11B(Nrank,Nrank), QincA(Nrank,Nrank), QincB(Nrank,Nrank))
  allocate (mv1(3,Nrank), nv1(3,Nrank), mv3(3,Nrank), nv3(3,Nrank),                 &
            mv1a(3,Nrank), nv1a(3,Nrank))
  allocate (jb(0:2), jbd(0:2), jba(0:Nrank), jbda(0:Nrank),                         &     
            hb(0:2), hbd(0:2)) 
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))                             
  kc = cmplx(k,0.0,O)
  do i = 1, Nrank
    do j = 1, Nrank
      Q31A(i,j)  = zero
      Q31B(i,j)  = zero
      Q11A(i,j)  = zero
      Q11B(i,j)  = zero
      QincA(i,j) = zero 
      QincB(i,j) = zero
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
      zl = cmplx(k * r,0.0,O)
      if (abs(zl) < MachEps) zl = cmplx(MachEps,MachEps,O)                            
      call leg_normalized (theta, 0, Nrank, Pnm, dPnm, pinm, taunm) 
      call besel_j (zl, Nrank, jba, jbda)         
      sth = sin(theta)
      cth = cos(theta)
      ro  = r * sth
      z   = r * cth
      roc = cmplx(ro,0.0,O)
      do n = 1, Nrank    
        dz    = z - zRe(n)
        dzc   = cmplx(dz,0.0,O)
        RR    = sqrt(roc * roc + (dzc - im * zIm(n))**2)
        if (abs(RR) < MachEps) RR = cmplx(MachEps,MachEps,O)
        sint  = ro / RR
        cost  = (dzc - im * zIm(n)) / RR
        argJi = ki * RR
        argJc = kc * RR    
        call besel_j (argJi, 2, jb, jbd)      
        call besel_h (argJc, 2, hb, hbd)
        Pmm   =   sqrt(3._O / 2._O) * cost              
        taumm = - sqrt(3._O / 2._O) * sint
        sinc  = sth * cost - cth * sint
        cosc  = cth * cost + sth * sint     
        ft    = 0.5_O * taumm 
        fl    = Pmm
        mv1(1,n) =   zero
        mv1(2,n) =   zero
        mv1(3,n) = - jb(1)  * ft         
        factl    =   jb(1)  * fl           
        factt    =   jbd(1) * ft
        nv1(1,n) =   (   factl * cosc + factt * sinc) / argJi
        nv1(2,n) =   ( - factl * sinc + factt * cosc) / argJi
        nv1(3,n) =   zero
        mv3(1,n) =   zero
        mv3(2,n) =   zero
        mv3(3,n) = - hb(1)  * ft         
        factl    =   hb(1)  * fl           
        factt    =   hbd(1) * ft
        nv3(1,n) =   (   factl * cosc + factt * sinc) / argJc
        nv3(2,n) =   ( - factl * sinc + factt * cosc) / argJc
        nv3(3,n) =   zero
        nr  = real(n * (n + 1),O)
        nm  = 1._O / sqrt(2._O * nr)            
        ftr = taunm(n) * nm
        flr = nr * Pnm(n) * nm        
        mv1a(1,n) =   zero
        mv1a(2,n) =   zero
        mv1a(3,n) = - jba(n)  * ftr         
        nv1a(1,n) =   jba(n)  * flr / zl
        nv1a(2,n) =   jbda(n) * ftr / zl
        nv1a(3,n) =   zero
      end do
      factm = - f * dA * pondere 
      call matQConv (Nrank, perfectcond, .false., ind_ref, factm,                   &
           mv3, nv3, mv1, nv1, nuv, Q31A, Q31B, Nrank, Nrank)
      factp = - factm 
      call matQConv (Nrank, perfectcond, .false., ind_ref, factp,                   &
           mv1a, nv1a, mv1, nv1, nuv, Q11A, Q11B, Nrank, Nrank) 
      call matQincConv (Nrank, factp, mv3, nv3, mv1a, nv1a, nuv,                    &
           QincA, QincB, Nrank, Nrank)
    end do
  end do                                
  call LU_SYSTEM_DIRECT (Q31A, Nrank, Nrank, QincA, Nrank, Nrank, Nrank, Nrank)
  call LU_SYSTEM_DIRECT (Q31B, Nrank, Nrank, QincB, Nrank, Nrank, Nrank, Nrank)                                           
  call product_matrices (Nrank, Nrank, Nrank, Q11A, Nrank, Nrank, QincA, Nrank,     &
       Nrank) 
  call product_matrices (Nrank, Nrank, Nrank, Q11B, Nrank, Nrank, QincB, Nrank,     &
       Nrank)              
  Cscat = 0._O
  Cext  = 0._O
  do i = 1, Nrank
    ir    = real(2 * i + 1,O)
    QA    = abs(Q11A(i,i))
    QB    = abs(Q11B(i,i))           
    Cscat = Cscat + ir * (QA * QA + QB * QB)
    QA    = real(Q11A(i,i))
    QB    = real(Q11B(i,i))
    Cext  = Cext + ir * (QA + QB)  
  end do  
  deallocate (Q31A, Q31B, Q11A, Q11B, QincA, QincB) 
  deallocate (mv1, nv1, mv3, nv3, mv1a, nv1a)
  deallocate (jb, jbd, jba, jbda, hb, hbd)
  deallocate (Pnm, dPnm, pinm, taunm)
end subroutine CscatCextConvDS
! **********************************************************************************
subroutine matQConv (Nrank, perfectcond, mirror, ind_ref, fact, mv3, nv3, mv1, nv1,  &
           nuv, QA, QB, nap, map)
  use parameters
  implicit none
  integer    :: Nrank, nap, map
  real(O)    :: nuv(3)
  complex(O) :: ind_ref, fact, mv1(3,Nrank), nv1(3,Nrank), mv3(3,Nrank),            &
                nv3(3,Nrank), QA(nap,map), QB(nap,map)
  logical    :: perfectcond, mirror
!
  integer    :: i, j, ni, nj
  real(O)    :: s
  complex(O) :: mvl(3), nvl(3), mvc(3), nvc(3), v1, v2, mixt_product
!
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
        QA(i,j) = QA(i,j) + (v1 + ind_ref * v2) * fact
        if (mirror) QA(i,j) = QA(i,j) + s * (v1 + ind_ref * v2) * fact
      else 
        QA(i,j) = QA(i,j) + v2 * fact
        if (mirror) QA(i,j) = QA(i,j) + s * v2 * fact
      end if                                              
      v1 = mixt_product (nuv, nvc, mvl)
      v2 = mixt_product (nuv, mvc, nvl)
      if (.not. perfectcond) then
        QB(i,j) = QB(i,j) + (v1 + ind_ref * v2) * fact
        if (mirror) QB(i,j) = QB(i,j) + s * (v1 + ind_ref * v2) * fact
      else 
        QB(i,j) = QB(i,j) + v2 * fact
        if (mirror) QB(i,j) = QB(i,j) + s * v2 * fact
      end if                                      
    end do   
  end do
end subroutine matQConv
!***********************************************************************************
subroutine matQincConv (Nrank, fact, mv3, nv3, mv1, nv1, nuv, QA, QB, nap, map)
  use parameters
  implicit none
  integer    :: Nrank, nap, map
  real(O)    :: nuv(3)
  complex(O) :: fact, mv1(3,Nrank), nv1(3,Nrank), mv3(3,Nrank), nv3(3,Nrank),       &
                QA(nap,map), QB(nap,map)
!
  integer    :: i, j
  complex(O) :: mvl(3), nvl(3), mvc(3), nvc(3), v1, v2, mixt_product
!
  do i = 1, Nrank
    mvl(1) = mv3(1,i)
    mvl(2) = mv3(2,i)
    mvl(3) = mv3(3,i)
    nvl(1) = nv3(1,i)
    nvl(2) = nv3(2,i)
    nvl(3) = nv3(3,i)
    do j = 1, Nrank
      mvc(1) = mv1(1,j)
      mvc(2) = mv1(2,j)
      mvc(3) = mv1(3,j)
      nvc(1) = nv1(1,j)
      nvc(2) = nv1(2,j)
      nvc(3) = nv1(3,j)
      v1 = mixt_product (nuv, mvc, nvl)
      v2 = mixt_product (nuv, nvc, mvl)
      QA(i,j) = QA(i,j) + (v1 + v2) * fact        
      v1 = mixt_product (nuv, nvc, mvl)
      v2 = mixt_product (nuv, mvc, nvl)
      QB(i,j) = QB(i,j) + (v1 + v2) * fact
    end do
  end do           
end subroutine matQincConv
