! **********************************************************************************
! *              GEOMETRY OF SURFACES AND MESH MODEL FOR INTEGRATING               *
! *   -------------------------------------------------------------------------    *
! *   -------------------------------------------------------------------------    *
module surface
use parameters

integer, parameter :: MEDIUM_ISOTROP  = 1 ! medium is isotropic (any type of symmetry is possible)
integer, parameter :: MEDIUM_UNIAXIAL = 2 ! medium is uniaxially anisotropic  (spherial symmetry is not possible)
integer, parameter :: MEDIUM_BIAXIAL  = 3 ! medium is biaxially anisotropic  (axial symmetry is not possible)
integer, parameter :: MEDIUM_GENERAL  = 10 ! with out any symmetry

type t_mesh_item
	real(O) :: theta, phi, r    ! theta, phi and radius vector of integration point
	real(O) :: r_theta, r_phi   ! derivation of radius vector (theta, phi) only informativ
	real(O) :: ns_r, ns_theta, ns_phi ! normal vector
	real(O) :: weight           ! [square of surface element] * [weight of surface element in gauss rule]
end type

type t_mesh
	type(t_mesh_item), pointer, dimension(:) :: items => NULL()
	integer	 :: items_count = 0
	logical  :: symmetry_z = .false., symmetry_y  = .false., symmetry_axial = .false.
	integer  :: symmetry_rotN = 1
end type

type t_obj_operation
	real(O) :: scale       = 1._O                 ! scale factor
	real(O) :: origin(3)   = (/0._O, 0._O, 0._O/) ! Position of origin in obj file
	real(O) :: euler_alpha = 0._O ! Euler angles : alpha
	real(O) :: euler_beta  = 0._O ! Euler angles : beta
	real(O) :: euler_gamma = 0._O ! Euler angles : gamma
	integer :: remesh_type = 0  ! 0 - don't divide each surface element (triangle) into subelements, 1 - otherwise
	integer :: remesh_scale= 1  ! 
end type

contains

subroutine mesg_gauss1D(n,z,w)
	integer, intent(in) :: N
    real(O), intent(out), dimension(:) :: Z,W
	!
    real(O) A /1._O/, B /2._O/, C /3._O/
	real(O) :: F, X, check, PB, PC, DJ, PA
	integer :: IND, K, I, J, M, niter
    
	IND=MOD(N,2)
    K=N/2+IND
    F=DFLOAT(N)

	do I = 1,K
		M=N+1-I
		IF (I.EQ.1) X=A-B/((F+A)*F)
		IF (I.EQ.2) X=(Z(N)-A)*4D0+Z(N)
		IF (I.EQ.3) X=(Z(N-1)-Z(N))*1.6D0+Z(N-1)
		IF (I.GT.3) X=(Z(M+1)-Z(M+2))*C+Z(M+3)
		IF (I.EQ.K.AND.IND.EQ.1) X=0D0
		
		niter=0
		check=1d-20
		
		do while (.true.)
			PB=1D0
			niter=niter+1
			if (niter.gt.100) then
				check=check*10d0
			end if
			PC=X
			DJ=A
			do J = 2,N
				DJ=DJ+A
				PA=PB
				PB=PC
				PC=X*PB+(X*PB-PA)*(DJ-A)/DJ
			end do
			PA=A/((PB-X*PC)*F)
			PB=PA*PC*(A-X*X)
			X=X-PB
			if (ABS(PB) .le. check*ABS(X)) exit
		end do

		Z(M)=X
		W(M)=PA*PA*(A-X*X)
		W(M)=B*W(M)
		if (I .ne. K .or. IND.ne.1) then
			Z(I)=-Z(M)
			W(I)=W(M)
		end if
	end do
end subroutine

subroutine mesh_derivation_to_ns(theta, phi, r, r_theta, r_phi, ns_r, ns_theta, ns_phi, weight)
	real(O), intent(in) :: theta, phi, r, r_theta, r_phi
	real(O), intent(out) :: ns_r, ns_theta, ns_phi, weight

    ns_r     = r*r*sin(theta)
    ns_theta = -r*r_theta*sin(theta)
    ns_phi   = -r*r_phi
    
    weight   = sqrt(ns_r*ns_r + ns_theta*ns_theta + ns_phi*ns_phi)

    ns_r     = ns_r / weight
    ns_theta = ns_theta / weight
    ns_phi   = ns_phi / weight
end subroutine

subroutine mesh_item_get(mesh, index, theta, phi, r, ns_r, ns_theta, ns_phi, weight)
	type(t_mesh), intent(in) :: mesh
	integer, intent(in)  :: index
	real(O), intent(out) :: theta, phi, r, ns_r, ns_theta, ns_phi, weight
	!
	if ((index <= 0) .or. (index > mesh%items_count)) then
		theta = 0
		phi   = 0
		r     = 1
		call mesh_derivation_to_ns(theta, phi, r, 0._O, 0._O, ns_r, ns_theta, ns_phi, weight)
		return
	end if
	theta    = mesh%items(index)%theta
	phi      = mesh%items(index)%phi
	r        = mesh%items(index)%r
	ns_r     = mesh%items(index)%ns_r
	ns_theta = mesh%items(index)%ns_theta
	ns_phi   = mesh%items(index)%ns_phi
	weight   = mesh%items(index)%weight
end subroutine

subroutine mesh_item_set_default(mesh, index, theta, phi, weight)
	type(t_mesh), intent(inout) :: mesh
	integer, intent(in) :: index
	real(O), intent(in) :: theta, phi, weight
	!
	real(O) :: ns_r, ns_theta, ns_phi, dA
    
    call mesh_derivation_to_ns(theta, phi, 1._O, 0._O, 0._O, ns_r, ns_theta, ns_phi, dA)
	call mesh_item_set(mesh, index, theta, phi, 1._O, ns_r*dA, ns_theta*dA, ns_phi*dA, weight)
end subroutine

subroutine mesh_item_set(mesh, index, theta, phi, r, ns_r, ns_theta, ns_phi, weight)
	type(t_mesh), intent(inout) :: mesh
	integer, intent(in) :: index
	real(O), intent(in) :: theta, phi, r, ns_r, ns_theta, ns_phi, weight
	!
	if ((index .le. mesh%items_count) .and. (index .gt. 0) ) then
		mesh%items(index)%theta = theta
		mesh%items(index)%phi   = phi
		mesh%items(index)%r       = r
		mesh%items(index)%ns_r    = ns_r
		mesh%items(index)%ns_theta= ns_theta
		mesh%items(index)%ns_phi  = ns_phi
		mesh%items(index)%weight  = weight
	end if
end subroutine

subroutine mesh_clear(mesh)
	type(t_mesh), intent(inout) :: mesh
	
	if (associated(mesh%items)) then
		deallocate(mesh%items)
	end if
	mesh%items_count = 0
end subroutine

subroutine mesh_expand(mesh, count)
	type(t_mesh), intent(inout) :: mesh
	integer, intent(in)      :: count
	!
	type(t_mesh_item), pointer, dimension(:) :: old_items
	integer :: change_count = 200
	integer :: imin = 1
	!
	change_count = MAX(change_count, count)
	if (associated(mesh%items)) then
		if ((mesh%items_count + count) .le. ubound(mesh%items,1)) then
			return
		end if

		old_items => mesh%items
		allocate(mesh%items(1:ubound(old_items,1)+change_count))
		mesh%items(1:ubound(old_items,1)) = old_items(1:ubound(old_items,1))
		imin = ubound(old_items,1)+1
		deallocate(old_items)
	else
		allocate(mesh%items(1:change_count))
	end if

	do i = imin, ubound(mesh%items,1)
		call mesh_item_set_default(mesh, i, 0._O, 0._O, 1._O)
	end do
end subroutine

subroutine mesh_gauss(mesh, n_theta, n_phi, symmetry_y, symmetry_z, symmetry_rotN)
	use parameters

	type(t_mesh), intent(inout) :: mesh
	integer, intent(in) :: n_theta, n_phi, symmetry_rotN
	logical, intent(in) :: symmetry_y, symmetry_z
	!
	real(O) :: Z(500),W(500),Z1(500),W1(500)
	integer :: i_theta, i_phi
	real(O) :: tet1, tet2, phi1, phi2, Phi, Theta, Coof
	!
	call mesg_gauss1D(n_theta,Z,W)
	call mesg_gauss1D(n_phi,Z1,W1)

	mesh%symmetry_y = symmetry_y
	mesh%symmetry_z = symmetry_z
	mesh%symmetry_rotN = symmetry_rotN

	tet1 = 0._O
	if (symmetry_z) then
		tet2 = Pi / 2._O
	else
		tet2 = Pi
	end if

	phi2 = Pi / symmetry_rotN
	if (symmetry_y) then
		phi1 = 0
	else
		phi1 = -phi2
	end if

	call mesh_expand(mesh, n_theta*n_phi)

	do i_theta = 1,n_theta
	do i_phi   = 1,n_phi
		theta = (Z(i_theta)+1._O)/2._O*(tet2-tet1)+tet1
		phi   = (Z1(i_phi)+1._O)/2._O*(phi2-phi1)+phi1
		coof  = W(i_theta)*(tet2-tet1)/2._O  *  W1(i_phi)*(phi2-phi1)/2._O  /  PI
		mesh%items_count = mesh%items_count + 1

		call mesh_item_set_default(mesh, mesh%items_count, theta, phi, coof)
	end do
	end do

end subroutine

subroutine init_ellipsoid(mesh, a, b, c, n_theta, n_phi, medium)
	use parameters

	type(t_mesh), intent(inout) :: mesh
	integer, intent(in) :: n_theta, n_phi
	real(O), intent(in) :: a, b, c
	integer, intent(in), optional :: medium
	!
	integer :: i, i_min, i_max, med
	real(O) :: phi, theta, weight, r1, r, r_theta, r_phi, ns_r, ns_theta, ns_phi, dA
	real(O) :: a2, b2, c2, CST, SNT, CSP, SNP
	!
	if (present(medium)) then
		med = medium
	else
		med = MEDIUM_ISOTROP
	end if
	i_min = mesh%items_count + 1
	if ((a == b) .and. (med < MEDIUM_BIAXIAL)) then
		call mesh_gauss(mesh, n_theta, n_phi, .true., .false., 4)
	elseif (med .ne. MEDIUM_GENERAL) then
		call mesh_gauss(mesh, n_theta, n_phi, .true., .false., 2)
	else
		call mesh_gauss(mesh, n_theta, n_phi, .false., .false., 1)
	end if
	i_max = mesh%items_count

	a2 = a*a
	b2 = b*b
	c2 = c*c

	do i = i_min, i_max
        call mesh_item_set(mesh, i, theta, phi, r, ns_r, ns_theta, ns_phi, weight)

! calculation of radius and normal
		CST = COS(Theta)
		SNT = SIN(Theta)
		CSP = COS(Phi)
		SNP = SIN(Phi)

		R1 = SNT*SNT*(CSP*CSP/A2+SNP*SNP/B2)+CST*CST/C2

		r = 1d0 / SQRT(R1)

		r_theta = 2d0*CST*SNT*(CSP*CSP/A2+SNP*SNP/B2-1/C2)
		r_theta = -0.5d0*r/R1 * r_theta

		r_phi   = 2d0*SNT*SNT*SNP*CSP*(1/B2-1/A2)
		r_phi   = -0.5d0*r/R1 * r_phi
! calculation of radius and normal

        call mesh_derivation_to_ns(theta, phi, r, r_theta, r_phi, ns_r, ns_theta, ns_phi, dA)
        call mesh_item_set(mesh, i, theta, phi, r, ns_r, ns_theta, ns_phi, weight*dA)
	end do
end subroutine

subroutine init_gauss_random_sphere(mesh, a, a_max, n_theta, n_phi, medium)
#ifdef KRASN
!	use Params
	implicit none

	type(t_mesh), intent(inout) :: mesh
	integer, intent(in) :: a_max, n_theta, n_phi
	complex(O), intent(in), dimension(0:a_max,0:a_max) :: a
	integer, intent(in), optional :: medium
	!
	integer :: i, i_min, i_max, m, n, dem
	real(O) :: phi, theta, weight, r, r_theta, r_phi, ns_r, ns_theta, ns_phi, dA
	complex(O) :: aim, c1
	real(O) D0(0:NN),D1(0:NN),D2(0:NN)
	!
	if (present(medium)) then
		med = medium
	else
		med = MEDIUM_GENERAL
	end if
	i_min = mesh%items_count + 1
	call mesh_gauss(mesh, n_theta, n_phi, .false., .false., 1)
	i_max = mesh%items_count

	aim   = (0._O, 1._O)

	do i = i_min, i_max
        call mesh_item_set(mesh, i, theta, phi, r, ns_r, ns_theta, ns_phi, weight)

! calculation of radius and normal
		r       = 0._O
		r_theta = 0._O
		r_phi   = 0._O

		do M = 0, a_max

		C1 = dcos(M*phi)+aim*dsin(M*phi)
		call VIG(dcos(theta),a_max,M,D0,D1,D2)

		do N = max(1,m), a_max
			r       = r       + dreal(a(N,M)*C1)*D0(N)
			r_theta = r_theta + dreal(a(N,M)*C1)*D2(N)
			r_phi   = r_phi   + dreal(a(N,M)*aim*M*C1)*D0(N)
		end do
		end do

		r       = dexp(r)
		r_theta = r*r_theta
		r_phi   = r*r_phi
! calculation of radius and normal

        call mesh_derivation_to_ns(theta, phi, r, r_theta, r_phi, ns_r, ns_theta, ns_phi, dA)
        call mesh_item_set(mesh, i, theta, phi, r, ns_r, ns_theta, ns_phi, weight*dA)
	end do
#endif
end subroutine

subroutine init_obj(mesh, file, group, operation, medium)
    use utils_str
	implicit none
	type(t_mesh), intent(inout) :: mesh
	character(*), intent(in)    :: file
	character(*), intent(in), optional          :: group
	type(t_obj_operation), intent(in), optional :: operation
	integer, intent(in), optional :: medium
	!
	integer        :: med, ios, line_index, vertex_index, face_index, i, j, ii(200), face_set_index, face_set(20)
	character(500) :: buf
	real(O)        :: vertex(1:NfacePD,3), face(1:NfacePD,3), p(3), v1(3), v2(3), ns(3), nuv(3)
	real(O)        :: r, theta, phi, dA
	logical        :: flag_stop, do_cmd, is_group, is_rotate, is_origin, is_scale
	
	if (present(medium)) then
		med = medium
	else
		med = MEDIUM_GENERAL
	end if

	open(unit=100, file=file, status='OLD', iostat=ios)
	if (ios.ne.0) then
        print "(/,2x,'Error of opening geometry file (.obj);')"
        stop
	end if
	
	vertex_index = 1
	line_index   = 1
	flag_stop    = .false.
	if (.not. present(group)) then
	    is_group = .false.
	else
	    is_group = len_c(group) > 0
	end if
	if (present(group)) then
        is_rotate = (abs(operation%euler_alpha)+abs(operation%euler_beta)+abs(operation%euler_gamma)) .ne. 0._O
        is_origin = (abs(operation%origin(1))+abs(operation%origin(2))+abs(operation%origin(3))) .ne. 0._O
        is_scale  = (abs(operation%scale)) .ne. 1._O
	end if
	do_cmd       = .not. is_group
	do while (.not. flag_stop)
	    read (100, '(A)', iostat = ios) buf
        if (ios /= 0) exit
	    if (buf(1:1)=="#") then
	        ! comment line
	    elseif ((buf(1:2)=="v ") .and. do_cmd) then ! vertex
	        read(buf(2:), *, iostat=ios) p(1), p(2), p(3)
	        if (ios.ne.0) then
                print "(/,2x,'Error of reading vertex position from obj-file (line #'(I3)');')", line_index
                stop
	        end if
			! apply operations (scale, rotate, change origin ) to vertex
			if (is_origin) then
				v1(1:3) = p(1:3) - operation%origin(1:3)
			end if
			if (is_rotate) then
				call T_cartesian_global_local(v1(1), v1(2), v1(3), operation%euler_alpha, &
					operation%euler_beta, operation%euler_gamma, p(1), p(2), p(3))
		    end if
			if (is_scale) then
				p(1:3) = p(1:3) * operation%scale
			end if
	        vertex(vertex_index, 1:3) = p(1:3)
	        vertex_index = vertex_index + 1
	    elseif ((buf(1:2)=="f ") .and. do_cmd) then ! face
            ! read info about face
            ii = 0
            buf = trim(buf(2:))
			! TO DO bug for command "f 1/1/1 2/2/2/ 3/3/3"
	        read(buf, *, iostat=ios) (ii(i),i=1,100)
	        if (ios.gt.0) then
                print "(/,2x,'Error of reading face set from obj-file (line #'(I3)');')", line_index
                stop
    	    end if
	        face_index = 1
	        do while (ii(face_index) .ne. 0)
	            i = ii(face_index)
	            if (i >= vertex_index) then
                    print "(/,2x,'Error of reading face set from obj-file (line #'(I3)');')", line_index
                    stop
    	        end if
    	        face(face_index, 1:3) = vertex(i, 1:3)
	            face_index = face_index + 1
	        end do
	        if (face_index .eq. 1) goto 100
	        
			! remesh operation
			face_set_index = 1
			if (present(operation)) then
				if ((operation%remesh_type > 0) .and. (operation%remesh_scale > 1) .and. (face_index .eq. 4)) then
					! TO DO remesh for remesh_scale = 2, 3, 4
					! face_set_index = > 1
				end if
			end if
			if (face_set_index .eq. 1) then ! without remeshing
				face_set_index = 2
				face_set(1) = 1
				face_set(2) = face_index
			end if
			do j = 1, face_set_index-1
				! calculate radius-vector, normal vector and square of face
				p(1:3)  = 0._O  
				do i = face_set(j), face_set(j+1)-1
	                p(1:3) = p(1:3) + face(i, 1:3)
		        end do
			    p(1:3) = p(1:3) / (face_set(j+1)-face_set(j))
				call T_cartesian_spherical(p(1), p(2), p(3), r, theta, phi)
            
	            ! calculate normal vector and square of face
		        ns(1:3) = 0._O
			    do i = face_set(j)+1, face_set(j+1)-2
				    v1(1:3) = face(i, 1:3)   - face(1, 1:3)
					v2(1:3) = face(i+1, 1:3) - face(1, 1:3)

	                ns(1) = ns(1) + v1(2)*v2(3) - v1(3)*v2(2)
		            ns(2) = ns(2) - v1(1)*v2(3) + v1(3)*v2(1)
			        ns(3) = ns(3) + v1(1)*v2(2) - v1(2)*v2(1)
				end do
	            dA = sqrt(ns(1)*ns(1) + ns(2)*ns(2) + ns(3)*ns(3))
		        ns(1:3) = ns(1:3) / dA
			    nuv(1)  = sin(theta) * cos(phi) * ns(1) + sin(theta) * sin(phi) * ns(2) + cos(theta) * ns(3) 
				nuv(2)  = cos(theta) * cos(phi) * ns(1) + cos(theta) * sin(phi) * ns(2) - sin(theta) * ns(3)
	            nuv(3)  =             -sin(phi) * ns(1) +              cos(phi) * ns(2)

		        call mesh_expand(mesh, 1)
			    mesh%items_count = mesh%items_count + 1
				call mesh_item_set(mesh, mesh%items_count, theta, phi, r, nuv(1), nuv(2), nuv(3), dA)
            end do

            100 continue
	    elseif (buf(1:2)=="g ") then ! group
	        buf = trim(buf(3:))
	        if (is_group) then
	            if (do_cmd) then ! current group is complete
	                flag_stop = .true.
	            elseif (trim(buf) == group) then ! begin of requested group
                    do_cmd    = .true.
	            end if
	        end if
	    else
	        ! something others
	    end if
	    line_index = line_index + 1
	end do
	close(unit=100)
end subroutine

end module
