! ----------------------------------------------------------------------------------
!                       LIBRARY OF PARTICLE GEOMETRIES                             !   
!    -------------------------------------------------------------------------     !
!                           AXISYMMETRIC PARTICLES                                 !             
!                                                                                  !
! The significance of the geometry parameters specified in the input files is      !
! given below.                                                                     !
!                                                                                  !           
!    Particle      TypeGeom   Nsurf   Nparam                surf                   !
!    spheroid         1         2       1         surf(1) - length of the semi-    !
!                                                           axis along the         !
!                                                           symmetry axis          !
!                                                 surf(2) - length of the second   !
!                                                           semi-axis              !
!    cylinder         2         2       3         surf(1) - half-length of         !
!                                                           the cylinder           !
!                                                 surf(2) - cylinder radius        !
!                                                                                  !
!    rounded          3         2       3         surf(1) - half-length of         !
!     oblate                                                the cylinder           ! 
!    cylinder                                     surf(2) - cylinder radius        !
!                                                           including the rounded  !
!                                                           part                   ! 
!                                                                                  !
! 1. In order to explain the significance of the surface parameters we consider    !
! the geometry depicted in Figure 1.                                               !
!                                                                                  !
!      z ^                                                                         !
!        |          a                                                              !
!        |/- - - - - - - - -/:                                                     !
!        |                   :                                                     !
!        |        L_1        :                                                     !
!        |--------->---------.............                                         !
!        |        t_1        |          /                                          !
!        |                   |          ! b                                        !
!        |               t_2 v   x      !                                          !
!      O o - - - - - - - - - | - >....../.                                         !
!        |                   | L_3      !                                          !
!        |       t_3         |          ! c                                        !
!        |--------<----------.........../.                                         !
!        |       L_3                                                               !
!                                                                                  !
!                                   Figure 1                                       !
!                                                                                  !
! The generatrix curve L is characterized by Nsurf = 3 surface parameters:         !
!                                                                                  !
!                   surf(1) = a, surf(2) = b and surf(3) = c,                      !
!                                                                                  !
! and consists of Nparam = 3 smooth curves L_1, L_2 and L_3, i.e.,                 !
!                                                                                  !
!                          L = L_1  U  L_2  U  L_3.                                !
!                                                                                  !
! On each smooth curve L_i, i = 1,2,... Nparam, we consider the curve parameter    !
! t_i (which varies in the interval I_i) and assume the following parametric       !
! representation:                                                                  !
!                                                                                  !
!                  L_i = {(x,z) / x = x_i(t_i), z = z_i(t_i) with t_i in I_i }.    !
!                                                                                  !
! In the routines contained in the geometry library, "iparam" stands for i and for !
! a given index i, "param" stands for t_i. In order to define the geometry         !
! parameters, we consider a generic smooth curve, for which we omit to indicate    !
! the dependency on the index i,                                                   !
!                                                                                  !
!                  L = {(x,z) / x = x(t), z = z(t) with t in I }                   !
!                                                                                  !
! The magnitude of the position vector is given by                                 !
!                                                                                  !
!                 r(t) = sqrt ( x(t)**2 + z(t)**2 ),                               !
!                                                                                  !
! the zenith angle theta is                                                        !
!                                                                                  !
!                 theta(t) = atan ( x(t) / z(t) ),                                 !
!                                                                                  !
! the azimuthal angle phi is phi = 0, the normal unit vector in spherical          !
! coordinates is given by                                                          !
!                                      r(t)                                        !
!                  n_r(t) = ---------------------------,                           ! 
!                           sqrt ( r(t)**2 + r'(t)**2 )                            !
!                                                                                  !
!                                             r'(t)                                !
!                  n_theta(t) = - ----------------------------- ,                  !
!                                  sqrt ( r(t)**2 + r'(t)**2 )                     !
!                                                                                  !
!                  n_phi(t) = 0,                                                   !
! where                                                                            !
!                           dr                                                     ! 
!                  r'(t) =  --(t),                                                 ! 
!                           dt                                                     !
!                                                                                  !
! and the surface element is defined as                                            !
!                                                                                  !
!                 dA(t) = sqrt ( r(t)**2 + r'(t)**2 ) * r(t) * sin[theta(t)].      !
!                                                                                  !
! 2. To compute the integrals over the particle surface we define nodes and        !
! weights on each smooth curve L_i. The number of integration points on the curve  !
! L_i is                                                                           !
!                                                                                  !
!                      Nint_i = int { c_i * Nint },                                !
!                                                                                  !
! where c_i is a multiplicative factor, c_i < 1, and Nint is the global number of  !
! integration points. Note that in the routines contained in the geometry library, !
! Nintparam(i) stands for Nint_i, where i = 1,2,...,Nparam. The code generates     !
! Nint_i nodes and weights in the interval I_i, and store these quantities in the  !
! arrays:                                                                          !
!                                                                                  ! 
!              paramG(i,j), i = 1,2,...,Nparam and j = 1,2,...,Nintparam(i)        !
!                                                                                  !
! and                                                                              !
!                                                                                  !
!              weightsG(i,j), i = 1,2,...,Nparam and j = 1,2,...,Nintparam(i),     !
!                                                                                  !
! respectively. An important parameter of the routines is the logical variable     !
! "mirror ". For particles with a plane of symmetry perpendicular to the axis of   !
! rotation (mirror symmetric particles with mirror  = .true.) with the surface     !
! parametrization                                                                  !
!                                                                                  !
!                       r(theta) = r(Pi - theta),                                  !
!                                                                                  !
! the T matrix can be computed by integrating theta over the interval [0,Pi/2].    !
!                                                                                  !
! 3. The distributed sources can be situated on the axis of symmetry of the       !
! particle (for highly elongated particles) or in the complex plane (for highly    !
! flattened particles). In Figure 2 we show the generatrix of a prolate cylinder   !
! and the image of the generatrix in the complex plane.                            !
!                                                                                  !
!          z ^                           Rez ^                                     !
!            |------                   ------|------                               !
!            |      |                 |      |      |                              !
!            |      |                 |Rez_4 o      |                              !
!            |      |                 |      |      |                              !
!            |      |                 |Rez_2 o      |                              !
!            |      |      x          |      | O    |    Im z                      !
!          O o - - -| - - - >         |Rez_1 o - - -| - - - >                      !
!            |      |                 |      |      |                              !
!            |      |                 |Rez_3 o      |                              !
!            |      |                 |      |      |                              !
!            |      |                 |Rez_5 o      |                              !
!            |      |                 |      |      |                              !
!            |------                   ------|------                               !
!            |                               |                                     !
!        Generatrix of a      Image of the generatrix in the complex plane and     !
!        prolate cylinder     the distribution of 5 sources on the Rez-axis        !
!                                                                                  !
!                                    Figure 2                                      !
!                                                                                  !
! For this configuration, the sources are distributed along the axis of            !
! symmetry (Rez-axis), i.e., Imz_1 = 0, Imz_2 = 0,...,Imz_5 = 0. In Figure 3       !
! we show the generatrix of an oblate cylinder and the image of the generatrix     !
! in the complex plane.                                                            !
!                                                                                  !
!         z  ^                                    Re z ^                           !
!            |----------------         ----------------|----------------           !
!            |                |       |                |                |          !
!            |                |   x   |      Imz_4     | Imz_1   Imz_5  |  Im z    !
!          O o - - - - - - - -| >    -|- - o - o - o - o - o - o - o - -| - - >    !
!            |                |       |  Imz_6   Imz_2 | O   Imz_3      |          !
!            |                |       |                |                |          !
!            |----------------         ----------------|----------------           !
!            |                                         |                           !
!        Generatrix of an     Image of the generatrix in the complex plane and     !
!        oblate cylinder      the distribution of 6 sources on the Imz-axis        !
!                                                                                  !
!                                     Figure 3                                     !   
!                                                                                  ! 
! For this configuration, the sources are distributed along the Imz-axis, i.e.,    ! 
! Rez_1 = 0, Rez_2 = 0,...,Rez_6 = 0.                                              !
!                                                                                  !
! 4. The user can easily modify the routines: "elem_geomAXSYM",                    !
! "interpolation_listAXSYM" and "zDSAXSYM" to generate particles with other        !
! geometries. A new geometry with the index TypeGeom = 4 can be added to the       !
! existing geometries, but the list of parameters must be maintained.              !
! ----------------------------------------------------------------------------------
subroutine elem_geomAXSYM (TypeGeom, Nsurf, surf, param, iparam, r, theta, phi,     &
           dA, n)
!-----------------------------------------------------------------------------------
! The routine computes the geometry parameters of an axisymmetric surface          !
! for a given curve (surface) parameter.                                           !
!                                                                                  !
! Input parameters:                                                                !
! - TypeGeom - integer, index of the particle geometry.                            !
! - Nsurf (integer) - number of surface parameters.                                !
! - surf (real array) - surface parameters specifying the shape of the particle.   !          
! - param (real) - actual curve (surface) parameter.                               !
! - iparam (integer) - index of the actual curve (surface) parameter.              !
!                                                                                  !
! Output parameters:                                                               !
! - r (real) - magnitude of the position vector.                                   !    
! - theta (real) - zenith angle.                                                   !
! - phi (real) - azimuthal angle.                                                  !
! - dA (real) - surface element.                                                   !
! - n (real array) - normal unit vector in spherical coordinates.                  !
!----------------------------------------------------------------------------------- 
  use parameters
  implicit none  
  integer :: TypeGeom, Nsurf, iparam
  real(O) :: surf(Nsurf), param, r, theta, phi, dA, n(3)
!
  real(O) :: a, b, tamp, invtamp, dr, e, e2, invsin, invcos
!
  select case (TypeGeom)
  case (1)
!   --- spheroid ---
    a  = surf(1)
    b  = surf(2)
    e  = a / b
    e2 = e * e
    theta = param   
    phi  = 0._O  
    if (iparam == 1) then          
      tamp    = cos(theta) * cos(theta) + sin(theta) * sin(theta) * e2
      invtamp = 1._O / sqrt(tamp)
      r  =   a * invtamp
      dr = - a * cos(theta) * sin(theta) * (e2 - 1._O) * invtamp / tamp           
    end if
    tamp    = sqrt(r * r + dr * dr)
    invtamp = 1._O / tamp
    dA   =   tamp * r * sin(theta)      
    n(1) =   r * invtamp
    n(2) = - dr * invtamp
    n(3) =   0._O  
  case (2,3) 
!   --- cylinder and rounded oblate cylinder---
    a = surf(1)
    b = surf(2)
    theta  = param  
    phi    = 0._O               
    if (iparam == 1) then
      invcos = 1._O / cos(theta)      
      r  = a * invcos 
      dr = a * sin(theta) * invcos * invcos  
    else if (iparam == 2) then 
      if (TypeGeom == 2) then
        invsin = 1._O / sin(theta)     
        r  =   b * invsin
        dr = - b * cos(theta) * invsin * invsin  
      else
        e  = b - a
        e2 = e * e
        tamp = sqrt(a * a - e2 * cos(theta) * cos(theta))
        r  = e * sin(theta) + tamp
        dr = e * cos(theta) + e2 * sin(theta) * cos(theta) / tamp 
      end if
    else if (iparam == 3) then 
      invcos = 1._O / cos(theta)     
      r  = - a * invcos
      dr = - a * sin(theta) * invcos * invcos 
    end if
    tamp    = sqrt(r * r + dr * dr)
    invtamp = 1._O / tamp
    dA   =   tamp * r * sin(theta)
    n(1) =   r * invtamp
    n(2) = - dr * invtamp
    n(3) =   0._O   
  end select     
end subroutine elem_geomAXSYM     
! **********************************************************************************
subroutine interpolation_listAXSYM (TypeGeom, Nsurf, surf, Nint, Nparam, Nintparam, &
           paramG, weightsG, mirror)
!-----------------------------------------------------------------------------------
! The routine provides the nodes and weights for computing integrals over an       !
! axisymmmetric surface.                                                           !
!                                                                                  !
! Input parameters:                                                                !
! - TypeGeom (integer) - index of the particle geometry.                           !
! - Nsurf (integer) - number of surface parameters.                                !
! - surf (real array) - surface parameters specifying the shape of the particle.   !
! - Nint (integer) - global number of integration points.                          !
! - Nparam (integer) - number of integration curves.                               !
! - mirror  (logical) - if mirror  = t, the particle is mirror symmetric (the      !
!   plane of symmetry or the plane of reflection is perpendicular to the axis of   !
!   rotation).                                                                     !
!                                                                                  ! 
! Output parameters:                                                               !
! - Nintparam (integer array) - number of integration points (nodes) on all        !
!   integration curves.                                                            !
! - paramG (real array) - integration points (nodes) on all integration curves.    !
! - weightsG (real array) - weights on all integration curves.                     !            
!----------------------------------------------------------------------------------- 
  use parameters
  implicit none
  integer :: TypeGeom, Nsurf, Nint, Nparam, Nintparam(Nparam)
  real(O) :: surf(Nsurf), paramG(Nparam,Nint), weightsG(Nparam,Nint)      
  logical :: mirror 
!
  integer :: iparam, p, Nintl, Ninta, Nintb
  real(O) :: a, b, theta0, atheta, btheta, e  
  real(O),allocatable :: xt(:), wt(:)
!
  do iparam = 1, Nparam
    do p = 1, Nint
      paramG  (iparam,p) = 0._O
      weightsG(iparam,p) = 0._O
    end do
  end do
  select case (TypeGeom)
  case (1)
!   --- spheroid ---    
    Nintparam(1) = Nint     
    Nintl = Nint
    allocate (xt(Nintl),wt(Nintl))      
    atheta = 0._O
    if (mirror) then
      btheta = Pi / 2
    else
      btheta = Pi
    end if
    call Gauss_Legendre (atheta, btheta, Nintl, wt, xt)
    do p = 1, Nintl
      paramG  (1,p) = xt(p)
      weightsG(1,p) = wt(p)
    end do 
    deallocate (xt, wt)
  case (2,3)
!   --- cylinder and rounded oblate cylinder---
    a = surf(1)
    b = surf(2)
    if (TypeGeom == 2) then
      theta0 = atan(b / a)
    else
      e = b - a
      theta0 = atan(e / a)
    end if
    Nintb = int(b * Nint / 2 / (a + b))
    Ninta = Nint - 2 * Nintb
    if (Nintb < 20) then
      Nintb = 20
      if (Nintb > Nint) then
        print "(/,2x,'Error in the input file:')"
        print "(2x,'the number of integration points Nint is too low;')"        
        stop
      end if
    end if
    if (Ninta < 20) then
      Ninta = 20
      if (Ninta > Nint) then
        print "(/,2x,'Error in the input file:')"
        print "(2x,'the number of integration points Nint is too low;')"        
        stop
      end if
    end if
    Nintparam(1) = Nintb
    Nintparam(3) = Nintb
    Nintparam(2) = Ninta
    do iparam = 1, 3
      Nintl = Nintparam(iparam)
      allocate (xt(Nintl), wt(Nintl))
      if (iparam == 1) then
        atheta = 0._O
        btheta = theta0
      else if (iparam == 2) then
        atheta = theta0
        btheta = Pi - theta0
      else if (iparam == 3) then
        atheta = Pi - theta0
        btheta = Pi
      end if
      call Gauss_Legendre (atheta, btheta, Nintl, wt, xt)
      do p = 1, Nintl
        paramG  (iparam,p) = xt(p)
        weightsG(iparam,p) = wt(p)
      end do
      deallocate (wt, xt)
    end do      
  end select
end subroutine interpolation_listAXSYM
!***********************************************************************************
subroutine zDSAXSYM (TypeGeom, Nsurf, surf, Nrank, ComplexPlane, EpsZReIm, zRe, zIm)
!-----------------------------------------------------------------------------------
! The routine generates the coordinates of the distributed sources in the complex  !
! plane for an axisymmetric particle.                                              !
!                                                                                  !
! Input parameters:                                                                !
! - TypeGeom (integer) - index of the particle geometry.                           !
! - Nsurf (integer) - number of surface parameters.                                !
! - surf (real array) - surface parameters.                                        !              
! - Nrank (integer) - number of distributed sources (maximum expansion order).     !
! - ComplexPlane (logical) - if ComplexPlane = t, the distributed sources          !
!   are placed in the complex plane.                                               !
! - EpsZReIm (real) - parameter controlling the distribution of the discrete       !
!   sources.                                                                       !
!                                                                                  ! 
! Output parameters:                                                               !
! - zRe, zIm (real arrays) - coordinates of the distributed sources in the         !
!   complex plane.                                                                 ! 
!----------------------------------------------------------------------------------- 
  use parameters
  implicit none
  integer :: TypeGeom, Nrank, Nsurf
  real(O) :: surf(Nsurf), EpsZReIm, zRe(Nrank), zIm(Nrank)  
  logical :: ComplexPlane
!
  integer :: i
  real(O) :: zmin, zmax      
!  
  select case (TypeGeom)        
  case (1,2,3)
! ---- spheroid & cylinder ---        
    if (.not. ComplexPlane) then
      zmin = - EpsZReIm * surf(1)
      zmax =   EpsZReIm * surf(1)
      do i = 1, Nrank       
        zRe(i) = zmin + (i - 1) * 2._O * zmax / (Nrank - 1)
        zIm(i) = 0._O  
      end do
    else 
      zmin = - EpsZReIm * surf(2)
      zmax =   EpsZReIm * surf(2) 
      do i = 1, Nrank       
        zRe(i) = 0._O
        zIm(i) = zmin + (i - 1) * 2._O * zmax / (Nrank - 1)
      end do
    end if
  end select
end subroutine zDSAXSYM
! ----------------------------------------------------------------------------------
!                           NONAXISYMMETRIC PARTICLES                              !
!                                                                                  ! 
! The significance of the geometry parameters specified in the input files is      !
! given below.                                                                     !
!                                                                                  ! 
!    Particle      TypeGeom   Nsurf   Nparam                surf                   !
!    ellipsoid        1         3       1         surf(1) - length of the semi-    !
!                                                           axis along the x-axis  !          
!                                                 surf(2) - length of the semi-    !
!                                                           axis along the y-axis  ! 
!                                                 surf(3) - length of the semi-    !
!                                                           axis along the z-axis  ! 
!                                                                                  !
!    quadratic        2         2       6         surf(1) - half-length of         !
!      prism                                                the prism              !  
!                                                 surf(2) - half-length of the     !
!                                                           square side            !
!                                                                                  !
!   regular N-hedral  3         2       2         surf(1) - half-length of         !
!       prism                       if mirror  = t;         the prism              !
!   (with azimuthal                     3         surf(2) - half-lenght of the     !
!    symmetry and                   if mirror  = f          side of the regular    !
!    optionally with                                        basis                  !
!    mirror symmetry                                                               !
!                                                                                  !
! 1. In order to explain the significance of the surface parameters we consider    !
! the geometry depicted in Figure 1.                                               !
!                                                                                  !
!                                                                                  ! 
!               :       a        :                            v_i /\   u_i         !
!               :/- - - - - - - /:                            \  /  \              !
!               :                :                             \/    \ /           ! 
!               ----------------- .........                    /\     /            !
!             / | z ^          / |      /            z ^      /  \   / \           !
!           /       |        /   |    /  b             |     /    \ /   \          !
!         /     |   |      /     |  /                  |     \ O_i o     \         !
!        ----------------- ......|....                 |      \   / \    /         ! 
!       |       |   |     |      | /                   |       \ /   \  /          !
!       |           |     | y    | !                   |    S_i /     \/           !
!       |       | O o - - |->    | !                   |       / \    /\           !
!       |       _ _/_ _ _ |_ _ _ | !                 O |        y \  /  \          !
!       |         /       |        !  c                o --------> \/              !
!       |     /  / x      |    /   !                  /                            !
!       |   /             |  /     !                 /                             !
!       | /               |/       !                /                              !
!        ----------------- ......../..           x /                               !
!                                                                                  !
!                                   Figure 1                                       !
!                                                                                  !
! The surface S is characterized by Nsurf = 3 surface parameters:                  !
!                                                                                  !
!                 surf(1) = a, surf(2) = b and surf(3) = c,                        !
!                                                                                  !
! and consists of Nparam = 6 smooth surfaces S_1, S_2,...,S_6, i.e.,               !
!                                                                                  !
!                          S = S_1  U  S_2  U ... U  S_6.                          !
!                                                                                  !
! On each smooth surface S_i, i = 1,2,... Nparam, we consider the surface          !
! parameters u_i and v_i (which vary in the intervals Iu_i and Iv_i, respectively),!
! and assume the following parametric representation:                              !
!                                                                                  !
!        S_i = {(x,y,z) / x = x_i(u_i,v_i), y = y_i(u_i,v_i), z = z_i(u_i,v_i)     !
!                with u_i in Iu_i and v_i in Iv_i }.                               !
!                                                                                  !
! In the routines contained in the geometry library, "iparam" stands for i and for !
! a given index i, "param1" stands for u_i and "param2" stands for v_i. In order   !
! to define the geometry parameters, we consider a generic smooth surface, for     !
! which we omit to indicate the dependency on the index i,                         !
!                                                                                  !
!        S = {(x,y,z) / x = x(u,v), y = y(u,v), z = z(u,v)                         !
!              with u in Iu and v in Iv }                                          !
!                                                                                  !
! The magnitude of the position vector is given                                    !
!                                                                                  !
!               r(u,v) = sqrt ( x(u,v)**2 + y(u,v)**2 + z(u,v)**2 ),               !
!                                                                                  !
! and the zenith angle theta and the azimuthal angle phi are the polar angles of   !
! the position vector r. Setting                                                   !
!                          dx                   dx                                 !
!               x_u(u,v) = --(u,v),  x_v(u,v) = --(u,v),                           !
!                          du                   dv                                 ! 
! and similarly for y_u, y_v, z_u and z_v, and introducing the fundamental         !
! parameters of the surface by                                                     !
!                                                                                  !
!        E(u,v) = x_u(u,v)**2 + y_u(u,v)**2 + z_u(u,v)**2,                         !
!                                                                                  !
!        G(u,v) = x_v(u,v)**2 + y_v(u,v)**2 + z_v(u,v)**2,                         !
!                                                                                  !
!        F(u,v) = x_u(u,v) * x_v(u,v) + y_u(u,v) * y_v(u,v) + z_u(u,v) * z_v(u,v), !
!                                                                                  !        
! we define the surface element as                                                 !
!                                                                                  !
!        dA(u,v) = sqrt ( E(u,v) * G(u,v) - F(u,v)**2 ).                           !
!                                                                                  ! 
! The tangential unit vectors to the line coordinates of S tau1 and tau2 are given !
! by (in Cartesian coordinates)                                                    !
!                                                                                  !
!        tau1_x(u,v) = x_u(u,v) / sqrt (x_u(u,v)**2 + y_u(u,v)**2 + z_u(u,v)**2),  !
!                                                                                  !
!        tau1_y(u,v) = y_u(u,v) / sqrt (x_u(u,v)**2 + y_u(u,v)**2 + z_u(u,v)**2),  !
!                                                                                  !
!        tau1_z(u,v) = z_u(u,v) / sqrt (x_u(u,v)**2 + y_u(u,v)**2 + z_u(u,v)**2),  !
!                                                                                  !
! and                                                                              !
!                                                                                  !
!        tau2_x(u,v) = x_v(u,v) / sqrt (x_v(u,v)**2 + y_v(u,v)**2 + z_v(u,v)**2),  !
!                                                                                  !
!        tau2_y(u,v) = y_v(u,v) / sqrt (x_v(u,v)**2 + y_v(u,v)**2 + z_v(u,v)**2),  !
!                                                                                  !
!        tau2_z(u,v) = z_v(u,v) / sqrt (x_v(u,v)**2 + y_v(u,v)**2 + z_v(u,v)**2),  !
!                                                                                  !
! respectively, while the normal unit vector is (in Cartesian coordinates)         !
!                                                                                  !  
!                n(u,v) = tau1(u,v) x tau2(u,v) / |tau1(u,v) x tau2(u,v)|.         !
!                                                                                  !
! Note that in spherical coordinates, we have                                      !
!                                                                                  !
!      n_r(u,v)     =   sin(theta) * cos(phi) * n_x(u,v) +                         !               
!                       sin(theta) * sin(phi) * n_y(u,v) + cos(theta) * n_z(u,v)   !
!                                                                                  !
!      n_theta(u,v) =   cos(theta) * cos(phi) * n_x(u,v) +                         !              
!                       cos(theta) * sin(phi) * n_y(u,v) - sin(theta) * n_z(u,v)   !
!                                                                                  !
!      n_phi(u,v)   = - sin(phi) * n_x(u,v) + cos(phi) * n_y(u,v)                  !
!                                                                                  !
! 2. To compute the integrals over the particle surface we define nodes and weights! 
! on each smooth surface S_i. The numbers of integration points on the surface S_i !
! are                                                                              !
!                                                                                  !
!                      Nint1_i = int { c1_i * Nint1 },                             !
!                                                                                  !
! and                                                                              !
!                                                                                  !
!                      Nint2_i = int { c2_i * Nint2 },                             !
!                                                                                  !                     
! where c1_i and c2_i are multiplicative factors, c1_i < 1 and c2_i < 1, and Nint1 !
! and Nint2 are the global numbers of integration points. In the routines          !
! contained in the geometry library, Nintparam(i) stands for Nint1_i and Nint2_i,  !
! where i = 1,2,...,Nparam. The code generate Nint1_i nodes and weights in the     !
! interval Iu_i, and Nint2_i nodes and weights in the interval Iv_i. The nodes are !
! stored in the arrays                                                             !
!                                                                                  !
!              paramG1(i,j), i = 1,2,...,Nparam and j = 1,2,...,Nintparam(i)       !
!                                                                                  !
! and                                                                              !
!                                                                                  !
!              paramG2(i,j), i = 1,2,...,Nparam and j = 1,2,...,Nintparam(i)       !
!                                                                                  !
! while the weights are stored in the array                                        !
!                                                                                  !
!              weightsG(i,j), i = 1,2,...,Nparam and j = 1,2,...,Nintparam(i),     !
!                                                                                  !
! respectively. In order to explain the occupation of the arrays "paramG1",        !
! "paramG2" and "weightsG", we fix i and consider the discretization depicted in   !
! Figure 2.                                                                        !   
!                                                                                  !          
!             v ^                                                                  !
!          Nint2|...o-----o-----o-----o                                            !
!             : |   |     |     |     |                                            !
!             : |   |     |     |     |                                            !
!             l |...o-----o-----o-----o                                            !
!             : |   |     |     |     |                                            !
!             : |   |     |     |     |                                            !
!             2 o...o-----o-----o-----o                                            !
!               |   |     |     |     |                                            !
!               |   |     |     |     |                                            !
!             1 o...o-----o-----o-----o                                            !
!               |   :     :     :     :                                            !
!             O o --o-----o-----o-----------> u                                    !
!                   1     2 ... k...Nint1                                          ! 
!                                                                                  !
!                                      Figure 2                                    !
!                                                                                  !
! With {u_k}, k = 1,2,...,Nint1, being a discretization of the interval Iu, and    !
! {v_l}, l = 1,2,...,Nint2, being a discretization of the interval Iv, we have     !
!                                                                                  !
!               paramG1(*,(k-1)*Nint2+l) = u_k                                     !
!                                                                                  !
! and                                                                              !
!                                                                                  !
!               paramG1(*,(k-1)*Nint2+l) = v_l.                                    !
!                                                                                  !   
! Further,                                                                         !
!                                                                                  !
!               weightsG(*,(k-1)*Nint2+l) = wu_k * wv_l,                           !
!                                                                                  !
! where {wu_k} are the weights on the interval Iu and {wv_l} are the weights on    !
! the interval Iv.                                                                 !
!                                                                                  !
! Two important parameters of the routines are the logical variable "mirror " and  !
! the integer variable "Nazimutsym". For particles with a plane of symmetry        !
! perpendicular to the axis of rotation (mirror symmetric particles with mirror  = !
! .true.) with the surface parametrization                                         !
!                                                                                  !
!                       r(theta,phi) = r(Pi - theta,phi),                          !
!                                                                                  !
! the T matrix can be computed by integrating theta over the interval [0,Pi/2].    !
! For particles with azimuthal symmetry, i.e., particles with the surface          !       
! parametrization                                                                  !
!                                                                                  !
!                       r(theta,phi) = r(theta,phi + 2Pi / Nazimutsym),            !
!                                                                                  ! 
! where Nazimutsym >= 2, the T matrix can be computed by integrating phi over      !
! the interval [0,2Pi/Nazimutsym].                                                 !
!                                                                                  ! 
! 3. The routine "SurfaceElemNint" computes the following parameters of each       !
! discretized smooth surface i:                                                    !
! - the total area                                                                 !
!                dAtot = SUM {p = 1,2,...,Nintparam(i)} dA(p) * ponderG(i,p),      ! 
! - the mean area of the surface elements                                          !
!                dAmean = dAtot / Nintparam(i),                                    !
! - the maximum area of the surface elements                                       !
!                dAmax = MAX {p = 1,2,...,Nintparam(i)} dA(p) * ponderG(i,p),      !
! - the minimum area of the surface elements                                       !
!                dAmax = MIN {p = 1,2,...,Nintparam(i)} dA(p) * ponderG(i,p),      !
! - and the number of elements per unit surface (number density).                  !
!                Nd = Nintparam(i) / dAtot                                         !
!                                                                                  !
! 4. The user can easily modify the routines: "elem_geom3D", "interpolation_list3D"!
! and "SurfaceElemNint" to generate particles with other geometries. A new geometry!
! with the index TypeGeom = 4 can be added to the existing geometries, but the     !
! list of parameters must be maintained.                                           ! 
! ----------------------------------------------------------------------------------
subroutine elem_geom3D (TypeGeom, Nsurf, surf, param1, param2, iparam, r, theta,    &
           phi, dA, n, mirror , Nazimutsym)
!-----------------------------------------------------------------------------------
! The routine computes the geometry parameters of a nonaxisymmmetric surface       !
! for a given surface parameter.                                                   !
!                                                                                  !
! Input parameters:                                                                !
! - TypeGeom (integer) - index of the particle geometry.                           !
! - Nsurf (integer) - number of surface parameters.                                !
! - surf (real array) - surface parameters specifying the shape of the particle.   !              
! - param1, param2 (real variables) - actual surface parameters.                   !
! - iparam (integer) - index of the actual surface parameters.                     !
! - mirror  (logical) - if mirror  = t, the particle is mirror symmetric (the plane!
!   of symmetry or the plane of reflection is perpendicular to the axis of         !
!   rotation).                                                                     !
! - Nazimutsym (integer) - number of azimuthal symmetric sections. Note that       !
!   Nazimutsym >= 2.                                                               !
!                                                                                  !
! Output parameters:                                                               !
! - r (real) - module of the position vector.                                      !    
! - theta (real) - zenith angle.                                                   !
! - phi (real) - azimuthal angle.                                                  !
! - dA (real) - surface element.                                                   !
! - n (real array) - normal unit vector in spherical coordinates.                  !
!----------------------------------------------------------------------------------- 
  use parameters
  implicit none  
  integer  :: TypeGeom, Nsurf, iparam, Nazimutsym
  real(O)  :: surf(Nsurf), param1, param2, r, theta, phi, dA, n(3)
  logical  :: mirror 
!
  real(O)  :: a, b, c, u, v, x, y, z, xt, yt, zt, xp, yp, zp, E,F , G, tampon,      &
              tau1(3), tau2(3), nc(3), ro, l, alpha
!
  select case (TypeGeom)
  case (1)
!   --- ellipsoid ---  
    a = surf(1)
    b = surf(2)
    c = surf(3)
    if (iparam == 1) then
      u = param1
      v = param2
      x = a * sin(u) * cos(v)
      y = b * sin(u) * sin(v)
      z = c * cos(u)      
      call T_cartesian_spherical (x, y, z, r, theta, phi)
      xt =   a * cos(u) * cos(v)
      yt =   b * cos(u) * sin(v)
      zt = - c * sin(u)
      xp = - a * sin(u) * sin(v)
      yp =   b * sin(u) * cos(v)
      zp =   0._O
      E  =   xt * xt + yt * yt + zt * zt
      G  =   xp * xp + yp * yp + zp * zp
      F  =   xt * xp + yt * yp + zt * zp
      dA =   sqrt(E * G - F * F)
      tampon  = sqrt(xt * xt + yt * yt + zt * zt)
      tau1(1) = xt / tampon
      tau1(2) = yt / tampon
      tau1(3) = zt / tampon
      tampon = sqrt(xp * xp + yp * yp + zp * zp)
      if (tampon == 0._O) then
        tau2(1) = 0._O
        tau2(2) = 1._O
        tau2(3) = 0._O
      else
        tau2(1) = xp / tampon
        tau2(2) = yp / tampon
        tau2(3) = 0._O
      end if
      call vector_product_real (tau1, tau2, nc)
      tampon = sqrt(nc(1)**2 + nc(2)**2 + nc(3)**2)
      nc(1) =  nc(1) / tampon
      nc(2) =  nc(2) / tampon
      nc(3) =  nc(3) / tampon        
      n(1)  =   sin(theta) * cos(phi) * nc(1) +                                     &
                sin(theta) * sin(phi) * nc(2) + cos(theta) * nc(3)
      n(2)  =   cos(theta) * cos(phi) * nc(1) +                                     &
                cos(theta) * sin(phi) * nc(2) - sin(theta) * nc(3)
      n(3)  = - sin(phi) * nc(1) + cos(phi) * nc(2)
    endif
  case (2)
!   --- quadratic prism --- 
    a = surf(1)
    l = surf(2)
    if (iparam == 1) then
      x = a
      y = param2
      z = param1
      nc(1) = 1._O
      nc(2) = 0._O
      nc(3) = 0._O
    else if (iparam == 2) then
      x = param2
      y = a
      z = param1
      nc(1) = 0._O
      nc(2) = 1._O
      nc(3) = 0._O
    else if (iparam == 3) then
      x = - a
      y =   param2
      z =   param1
      nc(1) = - 1._O
      nc(2) =   0._O
      nc(3) =   0._O
    else if (iparam == 4) then
      x =   param2
      y = - a
      z =   param1
      nc(1) =   0._O
      nc(2) = - 1._O
      nc(3) =   0._O
    else if (iparam == 5) then
      x =    param2
      y =    param1
      z =  - l
      nc(1) =   0._O
      nc(2) =   0._O
      nc(3) = - 1._O
    else if (iparam == 6) then
      x = param2
      y = param1
      z = l
      nc(1) = 0._O
      nc(2) = 0._O
      nc(3) = 1._O
    end if 
    call T_cartesian_spherical (x, y, z, r, theta, phi)
    dA = 1._O
    n(1) =   sin(theta) * cos(phi) * nc(1) + sin(theta) * sin(phi) * nc(2) +        &
             cos(theta) * nc(3)
    n(2) =   cos(theta) * cos(phi) * nc(1) + cos(theta) * sin(phi) * nc(2) -        &
             sin(theta) * nc(3)
    n(3) = - sin(phi) * nc(1) + cos(phi) * nc(2)
  case (3)
!  --- regular N-hedral prism --- 
    a = surf(1)
    l = surf(2) 
    alpha = Pi / Nazimutsym
    b = a / tan(alpha)  
    if (iparam == 1) then
      x  = param2
      y  = b
      z  = param1
      dA = 1._O
      nc(1) = 0._O
      nc(2) = 1._O
      nc(3) = 0._O     
    else if (iparam == 2) then
      ro  = param2
      phi = param1
      x  = ro * sin(phi)
      y  = ro * cos(phi)
      z  = l
      dA = ro
      nc(1) = 0._O
      nc(2) = 0._O
      nc(3) = 1._O 
    else if (iparam == 3 .and. .not. mirror ) then
      ro  = param2
      phi = param1
      x  =   ro * sin(phi)
      y  =   ro * cos(phi)
      z  = - l
      dA =   ro
      nc(1) =   0._O
      nc(2) =   0._O
      nc(3) = - 1._O            
    end if
    call T_cartesian_spherical (x, y, z, r, theta, phi) 
    n(1) =   sin(theta) * cos(phi) * nc(1) + sin(theta) * sin(phi) * nc(2) +        &
             cos(theta) * nc(3)
    n(2) =   cos(theta) * cos(phi) * nc(1) + cos(theta) * sin(phi) * nc(2) -        &
             sin(theta) * nc(3)
    n(3) = - sin(phi) * nc(1) + cos(phi) * nc(2)          
  end select
end subroutine elem_geom3D 
! **********************************************************************************
subroutine interpolation_list3D (TypeGeom, Nsurf, surf, Nint1, Nint2, NintAL,       &
           Nparam, Nintparam, paramG1, paramG2, weightsG, mirror , Nazimutsym)
!-----------------------------------------------------------------------------------
! The routine provides the nodes and weights for computing integrals over a        !
! nonaxisymmmetric surface.                                                        !
!                                                                                  !
! Input parameters:                                                                !
! - TypeGeom (integer) - index of the particle geometry.                           !
! - Nsurf (integer) - number of surface parameters.                                !
! - surf (real array) - surface parameters specifying the shape of the particle.   !
! - Nint1, Nint2 (integer variables) - global numbers of integration points.       !
! - NintAL (integer) - specifies the dimension of Nintparam, paramaG1, paramG2     !
!   and weightsG.                                                                  !
! - Nparam (integer) - number of integration surfaces.                             !
! - mirror  (logical) - if mirror  = t, the particle is mirror symmetric (the plane!
!   of symmetry or the plane of reflection is perpendicular to the axis of         !
!   rotation).                                                                     !
! - Nazimutsym (integer) - number of azimuthal symmetric sections. Note that       !
!   Nazimutsym >= 2.                                                               !
!                                                                                  ! 
! Output parameters:                                                               !
! - Nintparam (integer array) - number of integration points (nodes) on all        !
!   integration surfaces.                                                          !
! - paramG1, paramG2 (real arrays) - integration points (nodes) on all             !
!   integration surfaces.                                                          !
! - weightsG (real array) - weights on all integration surfaces.                   !            
!-----------------------------------------------------------------------------------  
  use parameters
  implicit none
  integer :: TypeGeom, Nsurf, Nint1, Nint2, NintAL, Nparam, Nazimutsym,             &
             Nintparam(Nparam)
  real(O) :: surf(Nsurf), paramG1(Nparam,NintAL*NintAL),                            &
             paramG2(Nparam,NintAL*NintAL), weightsG(Nparam,NintAL*NintAL)
  logical :: mirror 
!
  integer :: i, j, N1, N2, N3, k
  real(O) :: a, b, l, alpha, phi, ro
  real(O),allocatable :: xt(:), wt(:), xp(:), wp(:) 
!
  do i = 1, Nparam
    do j = 1, NintAL * NintAL
      paramG1 (i,j) = 0._O
      paramG2 (i,j) = 0._O
      weightsG(i,j) = 0._O
    end do
  end do
  select case (TypeGeom)
  case (1)
!   --- ellipsoid ----  
    if (mirror) then
      N1 = int(Nint1 / 2) + 1  ! alternative setting N1 = Nint1
    else    
      N1 = Nint1
    end if
    N2 = Nint2
    Nintparam(1) = N1 * N2
    allocate (xt(N1), wt(N1), xp(N2), wp(N2))     
    a = 0._O
    if (.not. mirror) then
      b = Pi
    else
      b = Pi / 2._O
    end if
    call Gauss_Legendre (a, b, N1, wt, xt)
    a = 0._O
    b = 2._O * Pi
    call Gauss_Legendre (a, b, N2, wp, xp)      
    do j = 1, N1
      do i = 1, N2
        paramG1 (1,(j-1)*N2+i) = xt(j)
        paramG2 (1,(j-1)*N2+i) = xp(i)
        weightsG(1,(j-1)*N2+i) = wt(j) * wp(i)
      end do
    end do
    deallocate (xt, wt, xp, wp)
  case (2)
!   --- quadratic prism ---
    a  = surf(1)
    l  = surf(2)       
    N1 = Nint1
    N2 = Nint2
    allocate (xt(N1), wt(N1), xp(N2), wp(N2))       
    call Gauss_Legendre (- l, l, N1, wt, xt)
    call Gauss_Legendre (- a, a, N2, wp, xp)  
    do k = 1, 4
      Nintparam(k) = N1 * N2
      do i = 1, N1
        do j = 1, N2
          paramG1 (k,(i-1)*N2+j) = xt(i)
          paramG2 (k,(i-1)*N2+j) = xp(j)
          weightsG(k,(i-1)*N2+j) = wt(i) * wp(j)
        end do
      end do
    end do
    do k = 5, 6
      Nintparam(k) = N2 * N2
      do i = 1, N2
        do j = 1, N2
          paramG1 (k,(i-1)*N2+j) = xp(i)
          paramG2 (k,(i-1)*N2+j) = xp(j)
          weightsG(k,(i-1)*N2+j) = wp(i) * wp(j)
        end do
      end do
    end do
    deallocate (xt, wt, xp, wp)
  case (3)
!   --- regular N-hedral prism --- 
    a  = surf(1)
    l  = surf(2)
    if (mirror) then
      N1 = int(Nint1 / 2) + 1 ! alternative setting: N1 = Nint1  
    else
      N1 = Nint1
    end if
    N2 = Nint2
    N3 = N2       ! alternative setting: N3 = int(Nint2 / 4) + 3                                      
    Nintparam(1) = N1 * N2  
    allocate (xt(N1), wt(N1), xp(N2), wp(N2))       
    if (mirror) then
      call Gauss_Legendre (0._O, l, N1, wt, xt)
    else
      call Gauss_Legendre (- l, l, N1, wt, xt)
    end if
    call Gauss_Legendre (- a, a, N2, wp, xp)                
    do i = 1, N1
      do j = 1, N2
        paramG1 (1,(i-1)*N2+j) = xt(i)
        paramG2 (1,(i-1)*N2+j) = xp(j)
        weightsG(1,(i-1)*N2+j) = wt(i) * wp(j)
      end do
    end do      
    deallocate (xt, wt)
    allocate (xt(N3), wt(N3))
    Nintparam(2) = N2 * N3 
    if (.not. mirror) Nintparam(3) = N2 * N3
    alpha = Pi / Nazimutsym
    b = a / tan(alpha)                      
    call Gauss_Legendre (- alpha, alpha, N2, wp, xp)
    do i = 1, N2
      phi = xp(i)             
      ro  = b / cos(phi)      
      call Gauss_Legendre (0._O, ro, N3, wt, xt)
      do j = 1, N3      
        paramG1 (2,(i-1)*N3+j) = xp(i)
        paramG2 (2,(i-1)*N3+j) = xt(j)
        weightsG(2,(i-1)*N3+j) = wp(i) * wt(j)  
        if (.not. mirror) then
          paramG1 (3,(i-1)*N3+j) = xp(i)
          paramG2 (3,(i-1)*N3+j) = xt(j)
          weightsG(3,(i-1)*N3+j) = wp(i) * wt(j)                
        end if 
      end do
    end do      
    deallocate (xt, wt, xp, wp)      
  end select
end subroutine interpolation_list3D
! **********************************************************************************
subroutine SurfaceElemNint (wavenumber, TypeGeom, Nparam, Nsurf, surf, mirror ,       &
           Nazimutsym, Nint1, Nint2)
! ----------------------------------------------------------------------------------
! The routine computes some parameters of the discretized surface: the mean,       !
! maximum and minimum areas of the surface elements and the number of elements per !
! unit surface (number density).                                                   !
! ----------------------------------------------------------------------------------
  use parameters
  integer     :: TypeGeom, Nparam, Nsurf, Nint1, Nint2, Nazimutsym
  real(O)     :: wavenumber, surf(Nsurf)
  logical     :: mirror 
!
  integer     :: i, j, k, NintAL, N1, N2, N3, iparam, pint, Nintl
  real(O)     :: a, b, c, u, v, x, y, z, xth, yth, zth, xph, yph, zph, E, F, G,     &
                 dA, dAmean, dAtot, dAmax, dAmin, l, alpha, phi, ro, param1,        &
                 param2, r, theta  
  integer,allocatable :: Nintparam(:)
  real(O),allocatable :: paramG1(:,:), paramG2(:,:), weightsG(:,:), xt(:), wt(:),   &
                         xp(:), wp(:)   
!       
  NintAL = max(Nint1,Nint2)
  allocate (paramG1(Nparam,NintAL * NintAL), paramG2(Nparam,NintAL * NintAL),       &
            weightsG(Nparam,NintAL * NintAL))
  allocate (Nintparam(Nparam))
  do i = 1, Nparam
    do j = 1, NintAL * NintAL
      paramG1 (i,j) = 0._O
      paramG2 (i,j) = 0._O 
      weightsG(1,j) = 0._O
    end do
  end do  
  print "(/,2x,'Surface Elements for the Input Values of Nint1 and Nint2:')"
  print "(  2x,'surface    dAmean      dAmax       dAmin       dAtotal    NoDens')"                
  select case (TypeGeom)
  case(1)
!   --- ellipsoid ----
    if (mirror) then
      N1 = int(Nint1 / 2) + 1   ! alternative setting: N1 = Nint1
    else    
      N1 = Nint1
    end if                       
    N2 = Nint2
    Nintparam(1) = N1 * N2
    allocate (xt(N1), wt(N1), xp(N2), wp(N2))
    a = 0._O
    if (.not. mirror) then
      b = Pi
    else
      b = Pi / 2._O
    end if                  
    call Gauss_Legendre (a, b, N1, wt, xt)
    a = 0._O
    b = 2._O * Pi
    call Gauss_Legendre (a, b, N2, wp, xp)      
    do j = 1, N1
      do i = 1, N2
        paramG1 (1,(j-1)*N2+i) = xt(j)
        paramG2 (1,(j-1)*N2+i) = xp(i)   
        weightsG(1,(j-1)*N2+i) = wt(j) * wp(i)      
      end do
    end do      
    deallocate (xt, wt, xp, wp) 
    a = wavenumber * surf(1)
    b = wavenumber * surf(2)
    c = wavenumber * surf(3) 
    do iparam = 1, Nparam    
      Nintl  = Nintparam(iparam)
      dAtot =   0._O
      dAmax = - 1.e+10_O
      dAmin =   1.e+10_O
      do pint = 1, Nintl
        param1 = paramG1(iparam,pint)
        param2 = paramG2(iparam,pint)      
        u = param1
        v = param2
        x = a * sin(u) * cos(v)
        y = b * sin(u) * sin(v)
        z = c * cos(u)        
        call T_cartesian_spherical (x, y, z, r, theta, phi)
        xth =   a * cos(u) * cos(v)
        yth =   b * cos(u) * sin(v)
        zth = - c * sin(u)
        xph = - a * sin(u) * sin(v)
        yph =   b * sin(u) * cos(v)
        zph =   0._O
        E  = xth * xth + yth * yth + zth * zth
        G  = xph * xph + yph * yph + zph * zph
        F  = xth * xph + yth * yph + zth * zph
        dA = sqrt(E * G - F * F) * weightsG(iparam,pint)         
        dAtot = dAtot + dA
        if (dA > dAmax) dAmax = dA
        if (dA < dAmin) dAmin = dA
      end do
      dAmean = dAtot / Nintl
      print "(4x,i2,4x,1pe10.2,2x,1pe10.2,2x,1pe10.2,2x,1pe10.2,1x,1pe10.2)",       &
              iparam, dAmean, dAmax, dAmin, dAtot, Nintl / dAtot                                                
    end do
  case(2)
!   --- quadratic prism ---
    a  = wavenumber * surf(1)
    l  = wavenumber * surf(2)
    N1 = Nint1
    N2 = Nint2
    allocate (xt(N1), wt(N1), xp(N2), wp(N2))
    call Gauss_Legendre (-l, l, N1, wt, xt)
    call Gauss_Legendre (-a, a, N2, wp, xp)
    do k = 1, 4
      Nintparam(k) = N1 * N2
      do i = 1, N1
        do j = 1, N2
          paramG1 (k,(i-1)*N2+j) = xt(i)
          paramG2 (k,(i-1)*N2+j) = xp(j)   
          weightsG(k,(i-1)*N2+j) = wt(i) * wp(j)      
        end do
      end do  
    end do
    do k = 5, 6
      Nintparam(k) = N2 * N2
      do i = 1, N2
        do j = 1, N2
          paramG1 (k,(i-1)*N2+j) = xp(i)
          paramG2 (k,(i-1)*N2+j) = xp(j)   
          weightsG(k,(i-1)*N2+j) = wp(i) * wp(j)      
        end do
      end do  
    end do 
    deallocate (xt, wt, xp, wp)
    do iparam = 1, Nparam    
      Nintl  = Nintparam(iparam)
      dAtot =   0._O
      dAmax = - 1.e+10_O
      dAmin =   1.e+10_O
      do pint = 1, Nintl        
        dA    = weightsG(iparam,pint)    
        dAtot = dAtot + dA
        if (dA > dAmax) dAmax = dA
        if (dA < dAmin) dAmin = dA
      end do
      dAmean = dAtot / Nintl
      print "(4x,i2,4x,1pe10.2,2x,1pe10.2,2x,1pe10.2,2x,1pe10.2,1x,1pe10.2)",       &
              iparam, dAmean, dAmax, dAmin, dAtot, Nintl / dAtot                        
    end do             
  case(3)
!   --- regular N-hedral prism --- 
    a  = wavenumber * surf(1)
    l  = wavenumber * surf(2)
    if (mirror) then
      N1 = int(Nint1 / 2) + 1  ! alternative setting: N1 = Nint1 
    else    
      N1 = Nint1
    end if    
    N2 = Nint2
    N3 = N2       ! alternative setting: N3 = int(Nint2 / 4) + 3 
    Nintparam(1) = N1 * N2
    allocate (xt(N1), wt(N1), xp(N2), wp(N2))
    if (mirror) then
      call Gauss_Legendre (0._O, l, N1, wt, xt)
    else
      call Gauss_Legendre (-l, l, N1, wt, xt)
    end if
    call Gauss_Legendre (-a, a, N2, wp, xp)
    do i = 1, N1
      do j = 1, N2
        paramG1 (1,(i-1)*N2+j) = xt(i)
        paramG2 (1,(i-1)*N2+j) = xp(j)   
        weightsG(1,(i-1)*N2+j) = wt(i) * wp(j)      
      end do
    end do  
    deallocate (xt, wt)
    allocate (xt(N3), wt(N3))
    Nintparam(2) = N2 * N3
    if (.not. mirror) Nintparam(3) = N2 * N3
    alpha = Pi / Nazimutsym
    b = a / tan(alpha)    
    call Gauss_Legendre (-alpha, alpha, N2, wp, xp)
    do i = 1, N2
      phi = xp(i)
      ro  = b / cos(phi)
      call Gauss_Legendre (0._O, ro, N3, wt, xt)
      do j = 1, N3
        paramG1 (2,(i-1)*N3+j) = xp(i)
        paramG2 (2,(i-1)*N3+j) = xt(j)   
        weightsG(2,(i-1)*N3+j) = wp(i) * wt(j)
        if (.not. mirror) then
          paramG1 (3,(i-1)*N2+j) = xp(i)
          paramG2 (3,(i-1)*N3+j) = xt(j)   
          weightsG(3,(i-1)*N3+j) = wp(i) * wt(j)                
        end if
      end do
    end do  
    deallocate (xt, wt, xp, wp)  
    do iparam = 1, Nparam    
      Nintl  = Nintparam(iparam)
      dAtot =   0._O
      dAmax = - 1.e+10_O
      dAmin =   1.e+10_O
      do pint = 1, Nintl   
        param1 = paramG1(iparam,pint)
        param2 = paramG2(iparam,pint)
        if (iparam == 1) then                            
          dA = weightsG(iparam,pint)     
        else if (iparam == 2 .or. (iparam == 3 .and. .not. mirror)) then
          ro = param2
          dA = ro * weightsG(iparam,pint)   
        end if                        
        dAtot = dAtot + dA
        if (dA > dAmax) dAmax = dA
        if (dA < dAmin) dAmin = dA
      end do
      dAmean = dAtot / Nintl
      print "(4x,i2,4x,1pe10.2,2x,1pe10.2,2x,1pe10.2,2x,1pe10.2,1x,1pe10.2)",       &
              iparam, dAmean, dAmax, dAmin, dAtot, Nintl / dAtot                              
    end do                        
  end select  
  deallocate (paramG1, paramG2, Nintparam)
end subroutine SurfaceElemNint
!-----------------------------------------------------------------------------------
!                              COMPOSITE PARTICLES                                 !
!                                                                                  !
! The significance of the geometry parameters specified in the input files is      !
! given below.                                                                     !
!                                                                                  ! 
!    Particle      TypeGeom   Nsurf   Nparam                surf                   !
!  half-spheroids                                                                  !
!  with offset        1         3       2         surf(1) - length of the semi-    ! 
!    origins                                                axis along the         !
!                                                           symmetry axis          !
!                                                 surf(2) - length of the second   !
!                                                           semi-axis              !
!                                                 surf(3) - distance between the   !
!                                                           half-spheroids basis   !
!                                                           and the origin of the  !
!                                                           local coordinate       !
!                                                           system                 ! 
!                                                                                  ! 
!   3 cylinders       2         2       3         surf(1) - half-length of         !
!                                                           the cylinder           ! 
!                                                 surf(2) - cylinder radius        !
!                                                                                  !
! 1. In order to explain the significance of the surface parameters we consider    !
! the geometry depicted in Figure 1.                                               !
!                                                                                  !
!                                                                    curve iparam  !
!                  a1                                       a1     /               !
!            z ^/- - - /:                           z  ^ /- - -/:/                 !       
!              |        :                              |       /:                  !       
!              |--------............                   |--->---- ......            !                         
!              |        |     /                        | param  |    /             !
!              |        |     !                        |        |    !             !
!              |        |     !                        |        |    ! x1          !
!      ../..O1 o        |     !  b1                 O1 o- - - - | - - ->           !
!    z1  !     |        |     !       x                |        |    ! b1          !
!      ../...O o - - - - - - - - - - - >               |        |    !             !
!        !     |        |     !                        |        |    !             !
!    z2  !     |--------....../.....                   |-------- ..../.            !
!        !     |        |     !                        |                           !
!      ../..O2 o        |     !  b2                                                !
!              |        |     !                                                    ! 
!              |--------....../.....                                               !
!              |   a2   :                                                          !
!              |/- - - /:                                                          !
!                                                                                  !
!           Composite particle             Homogeneous region ipart = 1. The       !
!                                          generatrix consists of Nparam(1) = 3    !
!                                          smooth curves. On the smooth curve      !
!                                          iparam, the curve parameter is param    !
!                                                                                  !
!                                   Figure 1                                       !
!                                                                                  !
! The composite particle consists of Npart = 2 homogeneous regions (cylinders),    !
! and the axial positions of the local coordinate systems (corresponding to the    !
! homogeneous regions) with respect to the global coordinate system of the         !
! composite particle are given by the array "zpart", i.e., zpart(1) = z1 and       !
! zpart(2) = z2. The region ipart = 1 is characterized by Nsurf(ipart) = 2 surface !
! parameters:                                                                      !
!                                                                                  !
!                 surf(ipart,1) = a1 and surf(ipart,2) = b1,                       !
!                                                                                  !
! and the region ipart = 2 is also characterized by Nsurf(ipart) = 2 surface       !
! parameters:                                                                      !
!                                                                                  !
!                 surf(ipart,1) = a2 and surf(ipart,2) = b2.                       !
!                                                                                  !
! In addition to Nsurf(ipart),                                                     !
!                                                                                  !
!              Nsurfmax = MAX {ipart = 1,2,...,Npart} Nsurf(ipart).                !
!                                                                                  !
! is an input parameter of the code. The generatrix of the region "ipart" consists !
! of Nparam(ipart) smooth curves, and                                              !
!                                                                                  ! 
!              Nparammax = MAX {ipart = 1,2,...,Npart} Nparam(ipart)               !
!                                                                                  !
! is also an input parameter of the code. On the smooth curve "iparam" (which      !
! belongs to the region "ipart") we consider the curve parameter "param", which    !
! is defined with respect to the local coordinate system of the region "ipart".    !
! Consequently, the magnitude of the position vector r, the zenith angle theta,    !
! the surface element dA and the normal unit vector in spherical coordinates n (on !
! each smooth curve "iparam" belonging to the region "ipart") correspond to the    !
! LOCAL COORDINATE SYSTEM of the region. The significance of these parameters is   !
! as for axisymmetric particles.                                                   !
!                                                                                  ! 
! 2. To compute the integrals over the particle surface we define nodes and weights!
! on each smooth curve "iparam" (which belongs to the region "ipart"). As for      !
! axisymmetric particles, the numbers of integration points is                     !
!   Nintparam(ipart,iparam), ipart = 1,2,.,Npart, iparam = 1,2,.., Nparam(ipart)   !
! while the nodes and weights are stored in the arrays:                            !
!   paramG(ipart,iparam,j), ipart = 1,2,...,Npart, iparam = 1,2,..,Nparam(ipart),  !
!                           j = 1,2,...,Nintparam(ipart,iparam),                   !
! and                                                                              !
!   weightsG(ipart,iparam,j), ipart = 1,2,...,Npart, iparam = 1,2,..,Nparam(ipart),!
!                             j = 1,2,...,Nintparam(ipart,iparam),                 !
! respectively. Because the curve parameter "param" is defined with respect to the !
! local coordinate system of the region "ipart", the nodes correspond to the LOCAL !  
! COORDINATE SYSTEM of the region.                                                 !
!                                                                                  !
! 3. To distribute the discrete sources we follow the same prescriptions as for    !
! axisymmetric particles. However, the coordinates of the discrete sources are     !
! defined with respect to the GLOBAL COORDINATE SYSTEM. If Nranp(ipart) is the     !
! number of discrete sources for the region "ipart", and zRe'(ipart,j) and         !
! zIm'(ipart,j) with j = 1,2,...,Nranp(ipart), are the coordinates of the discrete !
! sources in the local coordinate system of the region, then                       !
!                                                                                  !
!            zRe(ipart,j) = zRe'(ipart,j) + zpart(ipart) and                       !
!                                                                                  ! 
!            zIm(ipart,j) = zIm'(ipart,j) + zpart(ipart).                          !
!                                                                                  !
! 4. The user can easily modify the routines: "elem_geomCOMP",                     !
! "interpolation_listCOMP" and "zDSCOMP" to generate particles with other          !
! geometries. A new geometry with the index TypeGeom = 3 can be added to the       !
! existing geometries, but the list of parameters must be maintained.              !  
!-----------------------------------------------------------------------------------
subroutine elem_geomCOMP (TypeGeom, Npart, ipart, Nsurfmax, surf, param, iparam, r, &
           theta, phi, dA, n)
!-----------------------------------------------------------------------------------
! The routine computes the geometry parameters of a composite particle.            !
!                                                                                  !
! Input parameters:                                                                !
! - TypeGeom (integer) - index of the particle geometry.                           !
! - Npart (integer) - number of homogeneous regions.                               !
! - ipart (integer) - index of the actual homogeneous region.                      !
! - Nsurfmax (integer) - maximum number of surface parameters.                     !
! - surf (real array) - surface parameters of all regions.                         !              
! - param (real) - actual curve (surface) parameter for the region ipart.          !
! - iparam (integer) - index of the actual curve (surface) parameter for the       !
!   region ipart.                                                                  !
!                                                                                  !
! Output parameters:                                                               !
! - r (real) - magnitude of the position vector.                                   !    
! - theta (real) - zenith angle.                                                   !
! - phi (real) - azimuthal angle.                                                  !
! - dA (real) - surface element.                                                   !
! - n (real array) - normal unit vector in spherical coordinates.                  ! 
!----------------------------------------------------------------------------------- 
  use parameters
  implicit none  
  integer :: TypeGeom, Npart, ipart, Nsurfmax, iparam
  real(O) :: surf(Npart,Nsurfmax), param, r, theta, phi, dA, n(3)
!
  real(O) :: a, b, c, e, tamp, dr, tamp1, tamp2, tamp3, tamp4
!
  select case (TypeGeom)
  case (1)
!   --- half-spheroids with offset origins ---
    phi = 0._O
    a = surf(ipart,1)
    b = surf(ipart,2)
    c = surf(ipart,3)
    e = a / b
    theta  = param
    if (ipart == 1) then      
      if (iparam == 1) then     
        tamp1 =   a * a * cos(theta) * cos(theta) +                                 &
                 (a * a - c * c) * e * e * sin(theta) * sin(theta)
        tamp2 = - a * a * cos(theta) * sin(theta) +                                 &
                 (a * a - c * c) * e * e * sin(theta) * cos(theta)
        tamp3 =   cos(theta) * cos(theta) + e * e * sin(theta) * sin(theta)
        tamp4 = - cos(theta) * sin(theta) + e * e * sin(theta) * cos(theta)
        r  = (sqrt(tamp1) - c * cos(theta)) / tamp3
        dr = (tamp2 / sqrt(tamp1) + c * sin(theta)) / tamp3 -                       &
              2._O * (sqrt(tamp1) - c * cos(theta)) * tamp4 / tamp3 / tamp3          
      else if (iparam == 2) then       
        r  = - c / cos(theta)
        dr = - c * sin(theta) / cos(theta) / cos(theta)
      end if
    else if (ipart == 2) then         
      if (iparam == 1) then     
        r  = c / cos(theta)
        dr = c * sin(theta) / cos(theta) / cos(theta) 
      else if (iparam == 2) then        
        tamp1 =   a * a * cos(theta) * cos(theta) +                                 &
                 (a * a - c * c) * e * e * sin(theta) * sin(theta)
        tamp2 = - a * a * cos(theta) * sin(theta) +                                 &
                 (a * a - c * c) * e * e * sin(theta) * cos(theta)
        tamp3 =   cos(theta) * cos(theta) + e * e * sin(theta) * sin(theta)
        tamp4 = - cos(theta) * sin(theta) + e * e * sin(theta) * cos(theta)
        r  = (sqrt(tamp1) + c * cos(theta)) / tamp3
        dr = (tamp2 / sqrt(tamp1) - c * sin(theta)) / tamp3 -                       &
              2._O * (sqrt(tamp1) + c * cos(theta)) * tamp4 / tamp3 / tamp3                     
      end if
    end if
    tamp =   sqrt(r * r + dr * dr)
    dA   =   tamp * r * sin(theta)
    n(1) =   r / tamp
    n(2) = - dr / tamp
    n(3) =   0._O
  case (2)
!   --- 3 cylinders ---
    phi  = 0._O
    a = surf(ipart,1)
    b = surf(ipart,2)    
    theta = param         
    if (iparam == 1) then      
      r  = a / cos(theta)
      dr = a * sin(theta) / cos(theta) / cos(theta)                    
    else if (iparam == 2) then      
      r  =   b / sin(theta)
      dr = - b * cos(theta) / sin(theta) / sin(theta)                             
    else if (iparam == 3) then      
      r  = - a / cos(theta)
      dr = - a * sin(theta) / cos(theta) / cos(theta)                    
    end if
    tamp =   sqrt(r * r + dr * dr)
    dA   =   tamp * r * sin(theta)
    n(1) =   r / tamp
    n(2) = - dr / tamp      
    n(3) =   0._O
  end select     
end subroutine elem_geomCOMP
! **********************************************************************************
subroutine interpolation_listCOMP (TypeGeom, Npart, Nsurfmax, surf, Nint, Nparammax,&
               Nintparam, paramG, weightsG)
!-----------------------------------------------------------------------------------
! The routine provides the nodes and weights for computing integrals over the      !
! surface of a composite particle.                                                 !
!                                                                                  !
! Input parameters:                                                                !
! - TypeGeom (integer) - index of the particle geometry.                           !
! - Npart (integer) - number of homogeneous regions.                               !
! - Nsurfmax (integer) - maximum number of surface parameters.                     !
! - surf (real array) - surface parameters of all regions.                         !              
! - Nint (integer) - global number of integration points.                          !
! - Nparammax (integer) - maximum number of integration curves.                    !
!                                                                                  ! 
! Output parameters:                                                               !
! - Nintparam (integer array) - number of integration points (nodes) on all        !
!   integration curves and regions.                                                !    
! - paramG (real array) - integration points (nodes) on all integration curves     !
!   and regions.                                                                   !
! - weightsG (real array) - weights on all integration curves and regions.         !
!----------------------------------------------------------------------------------- 
  use parameters                    
  implicit none
  integer :: TypeGeom, Npart, Nsurfmax, Nint, Nparammax, Nintparam(Npart,Nparammax)
  real(O) :: surf(Npart,Nsurfmax), paramG(Npart,Nparammax,Nint),                    &
             weightsG(Npart,Nparammax,Nint)      
!      
  integer  :: ipart, iparam, p, Nintl, NintElip, NintLine, Ninta, Nintb
  real(O)  :: a, b, c, theta0, atheta, btheta, LElip, LLine
  real(O),allocatable :: xt(:), wt(:)
!
  do ipart = 1, Npart
    do iparam = 1, Nparammax
      Nintparam(ipart,iparam) = 0
      do p = 1, Nint
        paramG  (ipart,iparam,p) = 0._O
        weightsG(ipart,iparam,p) = 0._O
      end do
    end do
  end do        
  select case (TypeGeom)
  case (1)
!   --- half-spheroids with offset origins ---          
    do ipart = 1, Npart
      a = surf(ipart,1)
      b = surf(ipart,2)
      c = surf(ipart,3)
      theta0 = atan(b / c)
      LElip = 0.25_O * Pi * (3._O * (a + b) - sqrt((a + 3._O * b) * (3._O * a + b))) 
      LLine = b
      NintElip = int( LElip * Nint / (LElip + LLine) )
      NintLine = Nint - NintElip
      if (NintElip < 20) then
        NintElip = 20
        if (NintElip > Nint) then
          print "(/,2x,'Error in the input file:')"
          print "(2x,'the number of integration points Nint is too low;')"      
          stop
        end if
      end if
      if (NintLine < 20) then
        NintLine = 20
        if (NintLine > Nint) then
          print "(/,2x,'Error in the input file:')"
          print "(2x,'the number of integration points Nint is too low;')"      
          stop
        end if
      end if
      if (ipart == 1) then
        Nintparam(ipart,1) = NintElip
        Nintparam(ipart,2) = NintLine
      else if (ipart == 2) then 
        Nintparam(ipart,1) = NintLine
        Nintparam(ipart,2) = NintElip
      end if       
      do iparam = 1, 2
        Nintl = Nintparam(ipart,iparam)
        allocate (xt(Nintl), wt(Nintl))
        if (ipart == 1) then
          if (iparam == 1) then
            atheta = 0._O
            btheta = Pi - theta0
          else if (iparam == 2) then
            atheta = Pi - theta0
            btheta = Pi
          end if
        else if (ipart == 2) then
          if (iparam == 1) then
            atheta = 0._O
            btheta = theta0
          else if (iparam == 2) then
            atheta = theta0
            btheta = Pi  
          end if
         end if          
         call Gauss_Legendre (atheta, btheta, Nintl, wt, xt)
        do p = 1, Nintl
          paramG  (ipart,iparam,p) = xt(p)
          weightsG(ipart,iparam,p) = wt(p)
        end do
        deallocate (xt, wt)
      end do  
    end do        
  case(2)
!   --- 3 cylinders ---    
    do ipart = 1, Npart
      a = surf(ipart,1)
      b = surf(ipart,2)
      theta0 = atan(b / a)
      Nintb  = int(b * Nint / 2 / (a + b))
      Ninta  = Nint - 2 * Nintb
      if (Nintb < 20) then
        Nintb = 20
        if (Nintb > Nint) then
          print "(/,2x,'Error in the input file:')"
          print "(2x,'the number of integration points Nint is too low;')"      
          stop
        end if
      end if
      if (Ninta < 20) then
        Ninta = 20
        if (Ninta > Nint) then
          print "(/,2x,'Error in the input file:')"
          print "(2x,'the number of integration points Nint is too low;')"      
          stop
        end if
      end if
      Nintparam(ipart,1) = Nintb
      Nintparam(ipart,3) = Nintb
      Nintparam(ipart,2) = Ninta
      do iparam = 1, 3
        Nintl = Nintparam(ipart,iparam)
        allocate (xt(Nintl), wt(Nintl))
        if (iparam == 1) then
          atheta = 0._O
          btheta = theta0
        else if (iparam == 2) then
          atheta = theta0
          btheta = Pi - theta0
        else if (iparam == 3) then
          atheta = Pi - theta0
          btheta = Pi
        end if
        call Gauss_Legendre (atheta, btheta, Nintl, wt, xt)
        do p = 1, Nintl
          paramG  (ipart,iparam,p) = xt(p)
          weightsG(ipart,iparam,p) = wt(p)
        end do
        deallocate (wt, xt)
      end do              
    end do
  end select
end subroutine interpolation_listCOMP
! **********************************************************************************
subroutine zDSCOMP (TypeGeom, Npart, Nsurfmax, surf, Nrankpmax, Nrankp, zpart,      &
           ComplexPlane, EpsZReIm, zRe, zIm)
!-----------------------------------------------------------------------------------
! The routine generates the coordinates of the distributed sources in the complex  !
! plane for a composite particle.                                                  !
!                                                                                  !
! Input parameters:                                                                !
! - TypeGeom (integer) - index of the particle geometry.                           !
! - Npart (integer) - number of homogeneous regions.                               !
! - Nsurfmax (integer) - maximum number of surface parameters.                     !
! - surf (real array) - surface parameters of all regions.                         ! 
! - Nrankpmax (integer) - maximum number of distríbuted sources over all regions.  !            
! - Nrankp (integer array) - numbers of distributed sources (maximum expansion     !
!   orders) for all regions.                                                       !
! - zart (real array) - axial positions of the local coordinate systems of the     !
!   regions with respect to the global coordinate system of the composite          !
!   particle.                                                                      !
! - ComplexPlane (logical array) - if ComplexPlane = t, the distributed sources    !
!   for a specific region are placed in the complex plane.                         !
! - EpsZReIm (real array) - parameters controlling the distribution of the         !
!   discrete sources for all regions.                                              !
!                                                                                  ! 
! Output parameters:                                                               !
! - zRe, zIm (real arrays) - coordinates of the distributed sources in the         !
!   complex plane for all regions.                                                 !
!----------------------------------------------------------------------------------- 
  use parameters
  implicit none
  integer :: TypeGeom, Npart, Nrankpmax, Nrankp(Npart), Nsurfmax
  real(O) :: surf(Npart,Nsurfmax), zpart(Npart), zRe(Npart,Nrankpmax),              &
             zIm(Npart,Nrankpmax), EpsZReIm(Npart)
  logical :: ComplexPlane(Npart)
!    
  integer :: i, ipart, p, Nrankpl
  real(O) :: zmin, zmax, a, b, c  
!
  do ipart = 1, Npart
    do p = 1, Nrankpmax
      zRe(ipart,p) = 0._O
      zIm(ipart,p) = 0._O
    end do
  end do
  select case (TypeGeom)
  case (1)      
!   --- half-spheroids with offset origins ---                          
    do ipart = 1, Npart      
      a = surf(ipart,1)
      b = surf(ipart,2)
      c = surf(ipart,3)                   
      if (.not. ComplexPlane(ipart)) then
        if (ipart == 1) then
          zmin = (1._O - EpsZReIm(ipart)) * a
          zmax = EpsZReIm(ipart) * a
        else if (ipart == 2) then
          zmax = - (1._O - EpsZReIm(ipart)) * a
          zmin = - EpsZReIm(ipart) * a
        end if
        Nrankpl = Nrankp(ipart)
        do i = 1, Nrankpl           
          zRe(ipart,i) = zmin + (i - 1) * (zmax - zmin) / (Nrankpl - 1)
          zIm(ipart,i) = 0._O  
        end do
      else                    
        zmin = - EpsZReIm(ipart) * b * sqrt(1 - c * c / a / a)
        zmax =   EpsZReIm(ipart) * b * sqrt(1 - c * c / a / a)                
        Nrankpl = Nrankp(ipart)
        do i = 1, Nrankpl           
          zRe(ipart,i) = zpart(ipart)
          zIm(ipart,i) = zmin + (i - 1) * (zmax - zmin) / (Nrankpl - 1)
        end do                    
      end if
    end do  
  case (2)
! ---- 3 cylinders ---    
    do ipart = 1, Npart              
      a = surf(ipart,1)
      b = surf(ipart,2)      
      if (.not.ComplexPlane(ipart)) then
        zmin = zpart(ipart) - EpsZReIm(ipart) * a
        zmax = zpart(ipart) + EpsZReIm(ipart) * a
        Nrankpl = Nrankp(ipart)
        do i = 1, Nrankpl           
          zRe(ipart,i) = zmin + (i - 1) * (zmax - zmin) / (Nrankpl - 1)
          zIm(ipart,i) = 0._O  
        end do
      else 
        zmin = - EpsZReIm(ipart) * b
        zmax =   EpsZReIm(ipart) * b
        Nrankpl = Nrankp(ipart)
        do i = 1, Nrankpl           
          zRe(ipart,i) = zpart(ipart)
          zIm(ipart,i) = zmin + (i - 1) * (zmax - zmin) / (Nrankpl - 1) 
        end do      
      end if     
    end do   
  end select      
end subroutine zDSCOMP
!-----------------------------------------------------------------------------------
!                              LAYERED PARTICLES                                   !
!                                                                                  !
! The significance of the geometry parameters specified in the input files is      !
! given below.                                                                     !
!                                                                                  ! 
!    Particle      TypeGeom   Nsurf   Nparam                surf                   !
!     layered         1         2        1        surf(1) - length of the semi-    !
!    spheroids                                              axis along the         !
!                                                           symmetry axis          ! 
!                                                 surf(2) - length of the second   !
!                                                           semi-axis              !    
!                                                                                  !
!     layered         2         2        3        surf(1) - half-length of         !
!    cylinders                                              the cylinder           ! 
!                                                 surf(2) - cylinder radius        !
!                                                                                  !
! 1. In order to explain the significance of the surface parameters we consider    !
! the geometry depicted in Figure 1.                                               !
!                                                                                  !
!                   a1                                      a1                     !
!            z ^/ - - -/:                             z ^/ - - -/:                 !
!              |--------:........                       |        :                 !
!              |        |      /                        |--->---- .......          !
!              |        |      !                        |/param  |    /            !
!              |        |< - - ! - -layer 1            /|        |    !            !
!              |        |      !                     /  |        |    !            !
!              |        |      !                   /    |        |    !            !       
!              |  a2    |      !              curve     |        |    !            !
!              |/ - /:  |      !              iparam    |        |    !            !
!              |     :  |      !                        |        |    !            !
!              |----- ..|....  ! b1                     |        |    ! b1   x1    !
!      ../..O1 o     |  |  /   !                     O1 o- - - - | - - - - - >     !
!    z1  !     |     |  |  !   |        x               |        |    !            !
!      ../..O  o - - - - - - - - - - - >                |        |    !            !
!    z2  !     |     |  |  !   !                        |        |    !            !
!       ./..O2 o     |  |  ! b2!                        |        |    !            !
!              |     |  |  !   !                        |        |    !            !
!              |     |<-| -! - !- - layer 2             |        |    !            !
!              |     |  |  !   !                        |        |    !            !
!              |----- ..|../.  !                        |        |    !            !
!              |--------......./..                      |--------...../..          !
!                                                                                  !
!           Layered particle               Layer ipart = 1. The generatrix         !
!                                          consists of Nparam(1) = 3 smooth        !
!                                          curves. On the smooth curve iparam,     !
!                                          the curve parameter is param.           !
!                                                                                  !
!                                   Figure 1                                       !
!                                                                                  !
! The layered particle consists of Npart = 2 layers, and the axial positions of    !
! the local coordinate systems (of the layers) with respect to the global          !
! coordinate system of the layered particle are given by the array "zpart", i.e.,  !
! zpart(1) = z1 and zpart(2) = z2. The layer ipart = 1 is characterized by         !
! Nsurf(ipart) = 2 surface parameters:                                             !
!                                                                                  !
!                 surf(ipart,1) = a1 and surf(ipart,2) = b1,                       !
!                                                                                  !
! and the layer ipart = 2 is also characterized by Nsurf(ipart) = 2 surface        !
! parameters:                                                                      !
!                                                                                  !
!                 surf(ipart,1) = a2 and surf(ipart,2) = b2.                       !
!                                                                                  !
! In addition to Nsurf(ipart),                                                     !
!                                                                                  !
!              Nsurfmax = MAX {ipart = 1,2,...,Npart} Nsurf(ipart).                !
!                                                                                  !
! is an input parameter of the code. The generatrix of the layer "ipart" consists  !
! of Nparam(ipart) smooth curves, and                                              !
!                                                                                  ! 
!              Nparammax = MAX {ipart = 1,2,...,Npart} Nparam(ipart)               !
!                                                                                  !
! is also an imput parameter of the code. On the smooth curve "iparam" (which      !
! belongs to the layer "ipart") we consider the curve parameter "param", which     !
! is defined with respect to the local coordinate system of the layer "ipart".     !
! Consequently, the magnitude of the position vector r, the zenith angle theta,    !
! the surface element dA and the normal unit vector in spherical coordinates n (on !
! each smooth curve "iparam" belonging to the layer "ipart") correspond to the     !
! LOCAL COORDINATE SYSTEM of the layer. The significance of these parameters is    !
! as for axisymmetric particles.                                                   !
!                                                                                  !
! 2. To compute the integrals over the particle surface we define nodes and weights!
! on each smooth curve "iparam" (which belongs to the layer "ipart"). As for       !
! axisymmetric particles, the numbers of integration points is                     !
!   Nintparam(ipart,iparam), ipart = 1,2,.,Npart, iparam = 1,2,.., Nparam(ipart)   !
! while the nodes and weights are stored in the arrays:                            !
!   paramG(ipart,iparam,j), ipart = 1,2,...,Npart, iparam = 1,2,..,Nparam(ipart),  !
!                           j = 1,2,...,Nintparam(ipart,iparam),                   !
! and                                                                              !
!   weightsG(ipart,iparam,j), ipart = 1,2,...,Npart, iparam = 1,2,..,Nparam(ipart),!
!                             j = 1,2,...,Nintparam(ipart,iparam),                 !
! respectively. Because the curve parameter "param" is defined with respect to the !
! local coordinate system of the layer "ipart", the nodes correspond to the LOCAL  !  
! COORDINATE SYSTEM of the layer.                                                  !
!                                                                                  !
! 3. To distribute the discrete sources we follow the same prescriptions as for    !
! axisymmetric particles. However, the coordinates of the discrete sources are     !
! defined with respect to the GLOBAL COORDINATE SYSTEM. If Nranp(ipart) is the     !
! number of discrete sources for the layer "ipart", and zRe'(ipart,j) and          !
! zIm'(ipart,j) with j = 1,2,...,Nranp(ipart), are the coordinates of the discrete !
! sources in the local coordinate system of the layer, then                        !
!                                                                                  !
!            zRe(ipart,j) = zRe'(ipart,j) + zpart(ipart) and                       !
!                                                                                  ! 
!            zIm(ipart,j) = zIm'(ipart,j) + zpart(ipart).                          !
!                                                                                  !
! 4. The user can easily modify the routines: "elem_geomLAY",                      !
! "interpolation_listLAY" and "zDSLAY" to generate particles with other            !
! geometries. A new geometry with the index TypeGeom = 3 can be added to the       !
! existing geometries, but the list of parameters must be maintained               !  
!-----------------------------------------------------------------------------------
subroutine elem_geomLAY (TypeGeom, Npart, ipart, Nsurfmax, surf, param, iparam, r,  &
           theta, phi, dA, n)
!-----------------------------------------------------------------------------------
! The routine computes the geometry parameters of a layered particle.              !
!                                                                                  !
! Input parameters:                                                                !
! - TypeGeom (integer) - index of the particle geometry.                           !
! - Npart (integer) - number of layers.                                            !
! - ipart (integer) - index of the actual layer.                                   !
! - Nsurfmax (integer) - maximum number of surface parameters.                     !
! - surf (real array) - surface parameters of all layers.                          !              
! - param (real) - actual curve (surface) parameter for the layer ipart.           !
! - iparam (integer) - index of the actual curve (surface) parameter for the       !
!   layer ipart.                                                                   !
!                                                                                  !
! Output parameters:                                                               !
! - r (real) - magnitude of the position vector.                                   !    
! - theta (real) - zenith angle.                                                   !
! - phi (real) - azimuthal angle.                                                  !
! - dA (real) - surface element.                                                   !
! - n (real array) - normal unit vector in spherical coordinates.                  !
!----------------------------------------------------------------------------------- 
  use parameters
  implicit none  
  integer :: TypeGeom, Npart, ipart, Nsurfmax, iparam
  real(O) :: surf(Npart,Nsurfmax), param, r, theta, phi, dA, n(3)
!
  real(O) :: a, b, tamp, dr
!
  select case (TypeGeom)
  case (1)
!   --- concentric spheroids ---
    phi  = 0._O
    a = surf(ipart,1)
    b = surf(ipart,2)
    theta = param   
    if (iparam == 1) then     
      tamp =   cos(theta) * cos(theta) + sin(theta) * sin(theta) * a * a / b / b
      r  =   a / sqrt(tamp)
      dr = - a * cos(theta) * sin(theta) * ((a * a) / (b * b) - 1) /                &
            (tamp * sqrt(tamp))    
    end if          
    tamp =   sqrt(r * r + dr * dr)
    dA   =   tamp * r * sin(theta)
    n(1) =   r / tamp
    n(2) = - dr / tamp
    n(3) =   0._O 
  case (2)
!   --- concentric cylinders ---
    phi  = 0._O
    a = surf(ipart,1)
    b = surf(ipart,2)
    theta = param       
    if (iparam == 1) then      
      r  = a / cos(theta)
      dr = a * sin(theta) / cos(theta) / cos(theta)            
    else if (iparam == 2) then      
      r  =   b / sin(theta)
      dr = - b * cos(theta) / sin(theta) / sin(theta)                 
    else if (iparam == 3) then      
      r  = - a / cos(theta)
      dr = - a * sin(theta) / cos(theta) / cos(theta)            
    end if
    tamp =   sqrt(r * r + dr * dr)
    dA   =   tamp * r * sin(theta)
    n(1) =   r / tamp
    n(2) = - dr / tamp
    n(3) =   0._O                       
  end select
end subroutine elem_geomLAY     
!***********************************************************************************
subroutine interpolation_listLAY (TypeGeom, Npart, Nsurfmax, surf, Nint, Nparammax, &
               Nintparam, paramG, weightsG)
!-----------------------------------------------------------------------------------
! The routine provides the nodes and weights for computing integrals over the      !
! surface of a layered particle.                                                   !
!                                                                                  !
! Input parameters:                                                                !
! - TypeGeom (integer) - index of the particle geometry.                           !
! - Npart (integer) - number of layers.                                            !
! - Nsurfmax (integer) - maximum number of surface parameters.                     !
! - surf (real array) - surface parameters of all layers.                          !              
! - Nint (integer) - global number of integration points.                          !
! - Nparammax (integer) - maximum number of integration curves.                    !
!                                                                                  ! 
! Output parameters:                                                               !
! - Nintparam (integer array) - number of integration points (nodes) on all        !
!   integration curves and layers.                                                 !   
! - paramG (real array) - integration points (nodes) on all integration curves     !
!   and layers.                                                                    !
! - weightsG (real array) - weights on all integration curves and layers.          !
!----------------------------------------------------------------------------------- 
  use parameters
  implicit none
  integer :: TypeGeom, Npart, Nsurfmax, Nint, Nparammax, Nintparam(Npart,Nparammax)
  real(O) :: surf(Npart,Nsurfmax), paramG(Npart,Nparammax,Nint),                    &
                 weightsG(Npart,Nparammax,Nint)
!      
  integer :: ipart, iparam, p, Nint0, Nintl, Ninta, Nintb
  real(O) :: a, b, ab, ab0, atheta, btheta, theta0
  real(O),allocatable :: xt(:), wt(:)
!
  do ipart = 1, Npart
    do iparam = 1, Nparammax
      Nintparam(ipart,iparam) = 0
      do p = 1, Nint
        paramG  (ipart,iparam,p) = 0._O
        weightsG(ipart,iparam,p) = 0._O
      end do
    end do
  end do        
  select case (TypeGeom)
  case (1)
!   --- concentric spheroids ---        
    do ipart = 1, Npart
      a  = surf(ipart,1)
      b  = surf(ipart,2) 
      ab = max(a,b)
      if (ipart == 1) then
        ab0   = ab
        Nint0 = Nint                
      else
        Nint0 = int(ab * Nint / ab0)
      end if
      if (Nint0 < 40) then
        Nint0 = 40
        if (Nint0 > Nint) then
          print "(/,2x,'Error in the input file:')"
          print "(2x,'the number of integration points Nint is too low;')"      
          stop
        end if
      end if
      Nintparam(ipart,1) = Nint0
      Nintl = Nintparam(ipart,1)
      allocate (xt(Nintl),wt(Nintl))
      atheta = 0._O
      btheta = Pi
      call Gauss_Legendre (atheta, btheta, Nintl, wt, xt)
      do p = 1, Nintl
        paramG  (ipart,1,p) = xt(p)
        weightsG(ipart,1,p) = wt(p)
      end do
      deallocate (xt, wt) 
    end do  
  case (2)
!   --- concentric cylinders ---
    do ipart = 1, Npart
      a  = surf(ipart,1)
      b  = surf(ipart,2)
      ab = max(a,b)
      theta0 = atan(b / a)    
      if (ipart == 1) then
        ab0   = ab
        Nint0 = Nint
      else
        Nint0 = int(ab * Nint / ab0)
      end if
      if (Nint0 < 40) Nint0 = 40
      Nintb = int(b * Nint0 / 2 / (a + b))
      Ninta = Nint0 - 2 * Nintb
      if (Nintb < 20) then
        Nintb = 20
        if (Nintb > Nint) then
          print "(/,2x,'Error in the input file:')"
          print "(2x,'the number of integration points Nint is too low;')"      
          stop
        end if
      end if
      if (Ninta < 20) then
        Ninta = 20
        if (Ninta > Nint) then
          print "(/,2x,'Error in the input file:')"
          print "(2x,'the number of integration points Nint is too low;')"      
          stop
        end if
      end if
      Nintparam(ipart,1) = Nintb
      Nintparam(ipart,3) = Nintb
      Nintparam(ipart,2) = Ninta
      do iparam = 1, 3
        Nintl = Nintparam(ipart,iparam)
        allocate (xt(Nintl), wt(Nintl))
        if (iparam == 1) then
          atheta = 0._O
          btheta = theta0
        else if (iparam == 2) then
          atheta = theta0
          btheta = Pi - theta0
        else if (iparam == 3) then
          atheta = Pi - theta0
          btheta = Pi
        end if
        call Gauss_Legendre (atheta, btheta, Nintl, wt, xt)
        do p = 1, Nintl
          paramG  (ipart,iparam,p) = xt(p)
          weightsG(ipart,iparam,p) = wt(p)
        end do
        deallocate (wt, xt)
      end do
    end do    
  end select
end subroutine interpolation_listLAY
!***********************************************************************************
subroutine zDSLAY (TypeGeom, Npart, Nsurfmax, surf, Nrankpmax, Nrankp, zpart,       &
           ComplexPlane, EpsZReIm, zRe, zIm)
!-----------------------------------------------------------------------------------
! The routine generates the coordinates of the distributed sources in the complex  !
! plane for a layered particle.                                                    !
!                                                                                  !
! Input parameters:                                                                !
! - TypeGeom (integer) - index of the particle geometry.                           !
! - Npart (integer) - number of layers.                                            !
! - Nsurfmax (integer) - maximum number of surface parameters.                     !
! - surf (real array) - surface parameters of all layers.                          ! 
! - Nrankpmax (integer) - maximum number of distríbuted sources over all layers.   !            
! - Nrankp (integer array) - numbers of distributed sources (maximum expansion     !
!   orders) for all layers.                                                        !
! - zpart (real array) - axial positions of the local coordinate systems of the    !
!   layers with respect to the global coordinate system of the layered             !
!   particle.                                                                      ! 
! - ComplexPlane (logical array) - if ComplexPlane = t, the distributed sources    !
!   for a specific layer are placed in the complex plane.                          !
! - EpsZReIm (real array) - parameters controlling the distribution of the         !
!   discrete sources for all layers.                                               !
!                                                                                  ! 
! Output parameters:                                                               !
! - zRe, zIm (real arrays) - coordinates of the distributed sources in the         !
!   complex plane for all layers.                                                  !
!----------------------------------------------------------------------------------- 
  use parameters
  implicit none
  integer :: TypeGeom, Npart, Nrankpmax, Nrankp(Npart), Nsurfmax
  real(O) :: surf(Npart,Nsurfmax), zpart(Npart), zRe(Npart,Nrankpmax),              &
             zIm(Npart,Nrankpmax), EpsZReIm(Npart)
  logical :: ComplexPlane(Npart)
!   
  integer :: i, ipart, p, Nrankpl
  real(O) :: zmin, zmax, a, b   
!
  select case (TypeGeom)
  case (1,2)
!   ---- concentric spheroids and cylinders ----        
    do ipart = 1, Npart
      do p = 1, Nrankpmax
        zRe(ipart,p) = 0._O
        zIm(ipart,p) = 0._O
      end do
    end do      
    do ipart = 1, Npart       
      a = surf(ipart,1)
      b = surf(ipart,2)                      
      if (.not. ComplexPlane(ipart)) then
        zmin = zpart(ipart) - EpsZReIm(ipart) * a
        zmax = zpart(ipart) + EpsZReIm(ipart) * a
        Nrankpl = Nrankp(ipart)
        do i = 1, Nrankpl           
          zRe(ipart,i) = zmin + (i - 1) * (zmax - zmin) / (Nrankpl - 1)
          zIm(ipart,i) = 0._O  
        end do          
      else 
        zmin = - EpsZReIm(Npart) * b
        zmax =   EpsZReIm(Npart) * b
        Nrankpl = Nrankp(ipart)
        do i = 1, Nrankpl           
          zRe(ipart,i) = zpart(ipart)
          zIm(ipart,i) = zmin + (i - 1) * (zmax - zmin) / (Nrankpl - 1)
        end do          
      end if
    end do  
  end select
end subroutine zDSLAY
