! *********************************************************************************
! *                          GEOMETRIC TRANSFORMATIONS                            *
! *     ---------------------------------------------------------------------     *
! *     Partial list of subroutine:                                               *                                  
! *       scalar_product,                       vector_product,                   *
! *       mixt_product,                         scalar_product_real,              *
! *       vector_product_real,                  T_cartesian_spherical,            *
! *       T_spherical_cartesian,                m_rotation,                       *
! *       m_rotation_inverse,                   T_cartesian_global_local,         *
! *       T_cartesian_local_global,             T_cartesian_local_globalC (unused)*
! *       T_cartesian_global_localC (unused),   T_spherical_global_local,         *    
! *       T_spherical_local_global,             m_unitvct_spherical,              * 
! *       m_unitvct_spherical_global,           angle_unitvct_ItL_ItG,            *     
! *       parameters_coefficients_ab,           parameters_coefficients_Polab     * 
! *********************************************************************************
function scalar_product (n, a) result (na)
  use parameters
  real(O)    :: n(3)
  complex(O) :: a(3), na
!
  na = n(1) * a(1) + n(2) * a(2) + n(3) * a(3)     
end function scalar_product
! **********************************************************************************
subroutine vector_product (a, b, c)
  use parameters
  complex(O) :: a(3), b(3), c(3)
!
  c(1) = a(2) * b(3) - a(3) * b(2)
  c(2) = a(3) * b(1) - a(1) * b(3)
  c(3) = a(1) * b(2) - a(2) * b(1)   
end subroutine vector_product
! **********************************************************************************
function mixt_product (n, a, b) result (nab)
  use parameters
  real(O)    :: n(3)
  complex(O) :: a(3), b(3), nab
!
  complex(O) :: ab(3), scalar_product
!
  call vector_product (a, b, ab)
  nab = scalar_product (n, ab)   
end function mixt_product
! **********************************************************************************
function scalar_product_real (a, b) result (prod)
  use parameters
  implicit none
  real(O) :: a(3), b(3), prod
!
  prod = a(1) * b(1) + a(2) * b(2) + a(3) * b(3)   
end function scalar_product_real 
! **********************************************************************************
subroutine vector_product_real (a, b, c)
  use parameters
  real(O) :: a(3), b(3), c(3)
!
  c(1) = a(2) * b(3) - a(3) * b(2)
  c(2) = a(3) * b(1) - a(1) * b(3)
  c(3) = a(1) * b(2) - a(2) * b(1)   
end subroutine vector_product_real   
!***********************************************************************************
subroutine T_cartesian_spherical (x, y, z, R, theta, phi)
!-----------------------------------------------------------------------------------
! The routine computes the spherical coordinates given the Cartesian coordinates.  !
! For theta = 0 or theta = Pi, the solution is undetermined with respect to phi,   !
! and in these cases we set (by convention) phi = 0.                               !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  real(O) :: x, y, z, R, theta, phi, za, ro, ro2, thetaA
!  
  za  = abs(z)
  ro2 = x * x + y * y
  ro  = sqrt(ro2)
  R   = sqrt(ro2 + z * z)
  if (ro < MachEps .and. za < MachEps) then
    theta = 0._O
    phi   = 0._O
    return
  end if
  if (za < MachEps) then
    theta = Pi / 2._O
  else
    thetaA = atan(ro / za)
    if (z >= MachEps) then
      theta = thetaA
    else
      theta = Pi - thetaA
    end if
  end if
  if (ro < MachEps) then
    phi = 0._O                  !undetermined
  else    
    phi = atan2(y,x)
    if (phi < 0._O) phi = 2._O * Pi + phi     
  end if               
end subroutine T_cartesian_spherical
! **********************************************************************************
subroutine T_spherical_cartesian (R, theta, phi, x, y, z)
!-----------------------------------------------------------------------------------
! The routine computes the Cartesian coordinates given the spherical coordinates.  !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  real(O) :: R, theta, phi, x, y, z
!
  x = R * sin(theta) * cos(phi)
  y = R * sin(theta) * sin(phi)
  z = R * cos(theta)   
end subroutine T_spherical_cartesian
!***********************************************************************************
subroutine m_rotation (alpha, beta, gamma, R)
!-----------------------------------------------------------------------------------      
! The routine computes the rotation matrix of Euler angles alpha, beta and gamma.  !
! The rotation matrix gives the following relations between the local and the      !
! global coordinates:                                                              !
!                          [x,y,z]^T = [R] * [X,Y,Z]^T                             !
! and                                                                              !
!                          [i,j,k]^T = [R] * [I,J,K]^T.                            !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none      
  real(O) :: alpha, beta, gamma, R(3,3)
!
  R(1,1) =   cos(alpha) * cos(beta) * cos(gamma) - sin(alpha) * sin(gamma)
  R(2,1) = - cos(alpha) * cos(beta) * sin(gamma) - sin(alpha) * cos(gamma)
  R(3,1) =   cos(alpha) * sin(beta)
!
  R(1,2) =   sin(alpha) * cos(beta) * cos(gamma) + cos(alpha) * sin(gamma)
  R(2,2) = - sin(alpha) * cos(beta) * sin(gamma) + cos(alpha) * cos(gamma)
  R(3,2) =   sin(alpha) * sin(beta)
!
  R(1,3) = - sin(beta) * cos(gamma)
  R(2,3) =   sin(beta) * sin(gamma)
  R(3,3) =   cos(beta)   
end subroutine m_rotation
!***********************************************************************************
subroutine m_rotation_inverse (alpha, beta, gamma, RT)
!-----------------------------------------------------------------------------------
! The routine computes the inverse of the rotation matrix of Euler angles          !
! alpha, beta and gamma. This matrix gives the following relations between the     !
! global and the local coordinates:                                                !  
!                         [X,Y,Z]^T = [RT] * [x,y,z]^T                             !
! and                                                                              !
!                         [I,J,K]^T = [RT] * [i,j,k]^T.                            !
!----------------------------------------------------------------------------------- 
  use parameters 
  implicit none
  real(O) :: alpha, beta, gamma, RT(3,3)
!
  RT(1,1) =   cos(alpha) * cos(beta) * cos(gamma) - sin(alpha) * sin(gamma)
  RT(1,2) = - cos(alpha) * cos(beta) * sin(gamma) - sin(alpha) * cos(gamma)
  RT(1,3) =   cos(alpha) * sin(beta)
!
  RT(2,1) =   sin(alpha) * cos(beta) * cos(gamma) + cos(alpha) * sin(gamma)
  RT(2,2) = - sin(alpha) * cos(beta) * sin(gamma) + cos(alpha) * cos(gamma)
  RT(2,3) =   sin(alpha) * sin(beta)
!
  RT(3,1) = - sin(beta) * cos(gamma)
  RT(3,2) =   sin(beta) * sin(gamma)
  RT(3,3) =   cos(beta)      
end subroutine  m_rotation_inverse
!***********************************************************************************
subroutine T_cartesian_global_local (xG, yG, zG, alpha, beta, gamma, xL, yL, zL)
!-----------------------------------------------------------------------------------
! The routine computes the local Cartesian coordinates at a rotation with Euler    !
! angles alpha, beta and gamma, given the global Cartesian coordinates.            !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  real(O) :: xG, yG, zG, alpha, beta, gamma, xL, yL, zL, R(3,3)
!
  call m_rotation (alpha, beta, gamma, R)
  xL = R(1,1) * xG + R(1,2) * yG + R(1,3) * zG
  yL = R(2,1) * xG + R(2,2) * yG + R(2,3) * zG
  zL = R(3,1) * xG + R(3,2) * yG + R(3,3) * zG   
end subroutine T_cartesian_global_local
!***********************************************************************************
subroutine T_cartesian_local_global (xL, yL, zL, alpha, beta, gamma, xG, yG, zG)
!-----------------------------------------------------------------------------------
! The routine computes the global Cartesian coordinates at a rotation with Euler   !
! angles alpha, beta and gamma, given the local Cartesian coordinates.             !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  real(O) :: xL, yL, zL, alpha, beta, gamma, xG, yG, zG, RT(3,3)
!
  call m_rotation_inverse (alpha, beta, gamma, RT)
  xG = RT(1,1) * xL + RT(1,2) * yL + RT(1,3) * zL
  yG = RT(2,1) * xL + RT(2,2) * yL + RT(2,3) * zL
  zG = RT(3,1) * xL + RT(3,2) * yL + RT(3,3) * zL   
end subroutine T_cartesian_local_global
!***********************************************************************************
subroutine T_cartesian_global_localC (xG, yG, zG, alpha, beta, gamma, xL, yL, zL)
!-----------------------------------------------------------------------------------
! The routine computes the local Cartesian coordinates at a rotation with Euler    !
! angles alpha, beta and gamma, given the global Cartesian coordinates.            !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  real(O)    :: alpha, beta, gamma, R(3,3)
  complex(O) :: xG, yG, zG,  xL, yL, zL 
!
  call m_rotation (alpha, beta, gamma, R)
  xL = R(1,1) * xG + R(1,2) * yG + R(1,3) * zG
  yL = R(2,1) * xG + R(2,2) * yG + R(2,3) * zG
  zL = R(3,1) * xG + R(3,2) * yG + R(3,3) * zG   
end subroutine T_cartesian_global_localC
!***********************************************************************************
subroutine T_cartesian_local_globalC (xL, yL, zL, alpha, beta, gamma, xG, yG, zG)
  use parameters
  implicit none
  real(O)    :: alpha, beta, gamma, RT(3,3)
  complex(O) :: xL, yL, zL, xG, yG, zG
!
  call m_rotation_inverse (alpha, beta, gamma, RT)
  xG = RT(1,1) * xL + RT(1,2) * yL + RT(1,3) * zL
  yG = RT(2,1) * xL + RT(2,2) * yL + RT(2,3) * zL
  zG = RT(3,1) * xL + RT(3,2) * yL + RT(3,3) * zL   
end subroutine T_cartesian_local_globalC
! **********************************************************************************
subroutine T_spherical_global_local (thetaG, phiG, alpha, beta, gamma, thetaL, phiL)
!-----------------------------------------------------------------------------------
! The routine computes the local spherical coordinates at a rotation with Euler    !
! angles alpha, beta and gamma, given the global spherical coordinates.            !
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none
  real(O) :: thetaG, phiG, alpha, beta, gamma, thetaL, phiL
  real(O) :: raza, xG, yG, zG, xL, yL, zL, theta, phi
!
  raza = 1._O
  call T_spherical_cartesian (raza, thetaG, phiG, xG, yG, zG)
  call T_cartesian_global_local (xG, yG, zG, alpha, beta, gamma, xL, yL, zL)
  call T_cartesian_spherical (xL, yL, zL, raza, theta, phi)
  thetaL = theta
  phiL   = phi 
!..................................................................................!
!          If thetaL = 0 or thetaL = 180 then phiL is undetermined and we          !
!          choose phiL such that: ItG = ItL and IpG = IpL.                         !
!..................................................................................!
  if (abs(thetaL) < MachEps) then    
    if (abs(gamma) < MachEps) then
      phiL = 0._O
    else
      phiL = 2._O * Pi - gamma
    end if
  end if        
  if (abs(thetaL - Pi) < MachEps) then    
    phiL = Pi - gamma
    if (phiL < 0._O) phiL = 2._O * Pi + phiL
  end if                  
end subroutine T_spherical_global_local  
! **********************************************************************************    
subroutine T_spherical_local_global (thetaL, phiL, alpha, beta, gamma, thetaG, phiG)
!-----------------------------------------------------------------------------------
! The routine computes the global spherical coordinates at a rotation with Euler   !
! angles alpha, beta and gamma, given the local spherical coordinates.             !
!-----------------------------------------------------------------------------------    
  use parameters
  use derived_parameters
  implicit none
  real(O) :: thetaL, phiL, alpha, beta, gamma, thetaG, phiG
  real(O) :: raza, xG, yG, zG, xL, yL, zL, theta, phi
!
  raza = 1._O
  call T_spherical_cartesian (raza, thetaL, phiL, xL, yL, zL)
  call T_cartesian_local_global (xL, yL, zL, alpha, beta, gamma, xG, yG, zG)
  call T_cartesian_spherical (xG, yG, zG, raza, theta, phi)
  thetaG = theta
  phiG   = phi 
!..................................................................................!      
!         Íf thetaG = 0 or thetaG = 180 then phiG is undetermined and we           !
!         choose phiG such that: ItG = ItL and IpG = IpL.                          !    
!..................................................................................!
  if (abs(thetaG) < MachEps) then    
    phiG = Pi + alpha
    if (phiG >= 2._O * Pi) phiG = phiG - 2._O * Pi      
  end if        
  if (abs(thetaG - Pi) < MachEps) phiG = alpha          
end subroutine T_spherical_local_global
! **********************************************************************************
subroutine m_unitvct_spherical (theta, phi, S)
!-----------------------------------------------------------------------------------
! The routine computes the matrix [S]: [Ir,It,Ip]^T = [S] * [i,j,k]^T.             !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  real(O) :: theta, phi, S(3,3)
!
  S(1,1) =   sin(theta) * cos(phi)
  S(1,2) =   sin(theta) * sin(phi)
  S(1,3) =   cos(theta)
!
  S(2,1) =   cos(theta) * cos(phi)
  S(2,2) =   cos(theta) * sin(phi)
  S(2,3) = - sin(theta)
!
  S(3,1) = - sin(phi)
  S(3,2) =   cos(phi)
  S(3,3) =   0._O   
end subroutine m_unitvct_spherical
!***********************************************************************************
subroutine m_unitvct_spherical_global (thetaL, phiL, alpha, beta, gamma, SG)
!-----------------------------------------------------------------------------------
! The routine computes the matrix [SG] as:                                         !       
!               [IrL,ItL,IpL]^T = [S] * [i,j,k]^T                                  !
!                               = [S] * [R] * [I,J,K]^T                            !
!                               = [SG] * [I,J,K]^T,                                ! 
! where:                                                                           !
!               [x,y,z]^T = [R] * [X,Y,Z]^T, or equivalently,                      !
!               [i,j,k]^T = [R] * [I,J,K]^T.                                       !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  real(O) :: thetaL, phiL, alpha, beta, gamma, SG(3,3)
!
  integer :: i, j, k
  real(O) :: S(3,3), R(3,3), sum
!
  call m_unitvct_spherical (thetaL, phiL, S)
  call m_rotation (alpha, beta, gamma, R)      
  do i = 1, 3
    do j = 1, 3
      sum = 0._O
      do k = 1, 3
        sum = sum + S(i,k) * R(k,j)
      end do
      SG(i,j) = sum
    end do  
  end do   
end subroutine m_unitvct_spherical_global
!***********************************************************************************
subroutine angle_unitvct_ItL_ItG (thetaG, phiG, thetaL, phiL, alpha, beta, gamma,   &
           cos_angle, sin_angle)
!-----------------------------------------------------------------------------------
! The routine computes cos_angle = ItL * ItG and sin_angle = ItL * IpG.            !
! For this purpose, we compute the matrix [S] as                                   !
!                       [IrG,ItG,IpG]^T = [S] * [I,J,K]^T                          ! 
! and the matrix [SG] as                                                           !
!                       [IrL,ItL,IpL]^T = [SG] * [I,J,K]^T.                        !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  real(O) :: thetaG, phiG, thetaL, phiL, alpha, beta, gamma, cos_angle, sin_angle
!      
  integer :: i
  real(O) :: ItL(3), ItG(3), IpG(3), SG(3,3), S(3,3), scalar_product_real 
!
  call m_unitvct_spherical (thetaG, phiG, S)      
  call m_unitvct_spherical_global (thetaL, phiL, alpha, beta, gamma, SG)
  do i = 1, 3
    ItL(i) = SG(2,i)   
    ItG(i) = S(2,i)    
    IpG(i) = S(3,i)
  end do
  cos_angle = scalar_product_real (ItL, ItG)
  sin_angle = scalar_product_real (ItL, IpG)  
end subroutine angle_unitvct_ItL_ItG  
! **********************************************************************************
subroutine parameters_coefficients_ab (thetaG, phiG, thetaL, phiL, alpha, beta,     &
           gamma, alphap, e0eT, e0eP)
!-----------------------------------------------------------------------------------
! The routine computes e0eT = e0 * ItL and e0eP = e0 * IpL, where                  !
!                      e0 = cos(alphap) * ItG + sin(alphap) * IpG.                 !
! For this purpose, we compute the matrix [S] as                                   !
!                      [IrG,ItG,IpG]^T = [S] * [I,J,K]^T                           !
! and the matrix [SG] as                                                           !
!                      [IrL,ItL,IpL]^T = [SG] * [I,J,K]^T.                         !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  real(O)  :: thetaG, phiG, thetaL, phiL, alpha, beta, gamma, alphap, e0eT, e0eP
!    
  integer  :: i
  real(O)  :: ItL(3), IpL(3), e0(3), SG(3,3), S(3,3), scalar_product_real
!
  call m_unitvct_spherical (thetaG, phiG, S)
  call m_unitvct_spherical_global (thetaL, phiL, alpha, beta, gamma, SG)
  do i = 1, 3
    ItL(i) = SG(2,i)    
    IpL(i) = SG(3,i)
    e0(i)  = cos(alphap) * S(2,i) + sin(alphap) * S(3,i)                          
  end do
  e0eT = scalar_product_real (e0, ItL)
  e0eP = scalar_product_real (e0, IpL)   
end subroutine parameters_coefficients_ab
! **********************************************************************************
subroutine parameters_coefficients_Polab (thetaG, phiG, thetaL, phiL, alpha, beta,  &
           gamma, epol_theta, epol_phi, e0eT, e0eP)
!-----------------------------------------------------------------------------------
! The routine computes e0eT = e0 * ItL and e0eP = e0 * IpL, where                  !
!                      e0 = epol_theta * ItG + epol_phi * IpG,                     !
! and epol_theta and epol_phi are complex numbers. As before we compute the        !
! matrix [S] as                                                                    !
!                      [IrG,ItG,IpG]^T = [S] * [I,J,K]^T                           !
! and the matrix [SG] as                                                           !
!                      [IrL,ItL,IpL]^T = [SG] * [I,J,K]^T.                         !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  real(O)    :: thetaG, phiG, thetaL, phiL, alpha, beta, gamma
  complex(O) :: epol_theta, epol_phi, e0eT, e0eP
!    
  integer    :: i
  real(O)    :: ItL(3), IpL(3), SG(3,3), S(3,3)
  complex(O) :: e0(3), scalar_product
!
  call m_unitvct_spherical (thetaG, phiG, S)
  call m_unitvct_spherical_global (thetaL, phiL, alpha, beta, gamma, SG)
  do i = 1, 3
    ItL(i) = SG(2,i)    
    IpL(i) = SG(3,i)
    e0(i)  = epol_theta * S(2,i) + epol_phi * S(3,i)                          
  end do
  e0eT = scalar_product (ItL, e0) 
  e0eP = scalar_product (IpL, e0)   
end subroutine parameters_coefficients_Polab


























