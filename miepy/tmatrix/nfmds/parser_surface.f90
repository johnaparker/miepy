module parser_surface
use surface
use parser
use intepretator

contains

subroutine parser_surface_do_new(mesh, tree, Nint1, Nint2, medium)
    type(t_tree_item), intent(in) :: tree
    type(t_mesh), intent(inout)   :: mesh 
    integer, optional, intent(in) :: medium
    integer                       :: Nint1, Nint2
    !
    type(t_tree_item), pointer :: tmp_item
    character(500) :: buf
	character(80)  :: obj_file, obj_group, geomtype
	type(t_obj_operation) :: obj_oper
	real(O) :: shape_params(100)
	integer :: med, ios, i, Nint1r, Nint2r
	!
    
	if (present(medium)) then
		med = medium
	else
		med = MEDIUM_ISOTROP
	end if

    call tree_find_value(tree, "shape/record" // char(0), buf)
	if (buf == "file_obj") then     ! obj file
		call tree_find_value(tree, "shape/filename" // char(0), obj_file)
		call evaluate(obj_file)
		call tree_find_value(tree, "shape/group" // char(0), obj_group)
		call evaluate(obj_group)
		! origin
		call tree_find_value(tree, "shape/origin" // char(0), buf)
		call evaluate(buf)
		if (trim(buf) .ne. "") then
			read(buf, fmt=*, iostat=ios) obj_oper%origin(1), obj_oper%origin(2), obj_oper%origin(3)
			if (ios .ne. 0) then
				print "(/,2x,'Error of translating (origin <<'(A20)'>> is not correnct);')", trim(buf)
				stop
			end if
		end if
		! rotation
		call tree_find_value(tree, "shape/rotation" // char(0), buf)
		call evaluate(buf)
		if (trim(buf) .ne. "") then
			read(buf, fmt=*, iostat=ios) obj_oper%euler_alpha, obj_oper%euler_beta, obj_oper%euler_gamma
			if (ios .ne. 0) then
				print "(/,2x,'Error of translating (Rotation <<'(A20)'>> is not correnct);')", trim(buf)
				stop
			end if
		end if
		! scale
		call tree_find_value(tree, "shape/scale" // char(0), buf)
		call evaluate(buf)
		if (trim(buf) .ne. "") then
			read(buf, fmt=*, iostat=ios) obj_oper%scale
			if (ios .ne. 0) then
				print "(/,2x,'Error of translating (Scale <<'(A20)'>> is not correnct);')", trim(buf)
				stop
			end if
		end if
		! remesh
		call tree_find_value(tree, "shape/remesh_div" // char(0), buf)
		call evaluate(buf)
		if (trim(buf) .ne. "") then
			read(buf, fmt=*, iostat=ios) obj_oper%remesh_scale
			if (ios .ne. 0) then
				print "(/,2x,'Error of translating (Remesh <<'(A20)'>> is not correnct);')", trim(buf)
				stop
			end if
		end if
		! generate mesh
		call init_obj(mesh, trim(obj_file), trim(obj_group), obj_oper, MEDIUM_GENERAL)
	elseif (buf == "file_fem") then ! fem file
		print "(/,2x,'Record type <<'(A20)'>> have still not implemented);')", trim(buf)
		stop
		! TO DO ...
	elseif (buf == "analytic") then ! analytic shape
		! Nint1
		call tree_find_value(tree, "mesh/nint1" // char(0), buf)
		call evaluate(buf)
		if (len_c(buf) .ne. 0) then
		    read (buf, fmt=*, iostat=ios) Nint1r
			if (ios .ne. 0) then
				print "(/,2x,'Error of reading Nint1 (Number <<'(A20)'>> is not correnct);')", trim(buf)
				stop
			end if
		else
            Nint1r = Nint1
		end if
		! Nint2
		call tree_find_value(tree, "mesh/nint2" // char(0), buf)
		call evaluate(buf)
		if (len_c(buf) .ne. 0) then
		    read (buf, fmt=*, iostat=ios) Nint2r
			if (ios .ne. 0) then
				print "(/,2x,'Error of reading Nint2 (Number <<'(A20)'>> is not correnct);')", trim(buf)
				stop
			end if
		else
		    Nint2r = Nint2
		end if
		! params
		call tree_find_value(tree, "shape/params" // char(0), buf)
		call evaluate(buf)
		if (len_c(buf) .ne. 0) then
		    shape_params = 0._O
		    read (buf, fmt=*, iostat=ios) (shape_params(i), i=1,ubound(shape_params,1))
			if (ios > 0) then
				print "(/,2x,'Error of reading of parameters (String <<'(A20)'>> is not correnct);')", trim(buf)
				stop
			end if
		end if
		! kind of shape
		call tree_find_value(tree, "shape/kind" // char(0), geomtype)
		call evaluate(geomtype)
		if (geomtype .eq. "ellipsoid") then
		    ! generate mesh
		    call init_ellipsoid(mesh, shape_params(1), shape_params(2), shape_params(3), Nint1r, Nint2r, MEDIUM_GENERAL)
		elseif (geomtype .eq. "hex") then
		    ! TO DO ...
		end if
	else
		print "(/,2x,'Error of translating (record type <<'(A20)'>> is unknown);')", trim(buf)
		stop
	end if
end subroutine

subroutine parser_surface_do(mesh, FileGeom, TypeGeom, Nsurf, surf, rp, np, area, &
         Nface, NintAL, Nint1, Nint2, Nparam, miror, Nazimutsym, GeomTree)
    implicit none
    type(t_mesh), intent(inout)   :: mesh
	! analytic
	integer, intent(in) :: Nsurf, NintAL, Nint1, Nint2, Nparam, Nazimutsym, TypeGeom
	real(O), intent(in) :: surf(Nsurf)
	logical, intent(in) :: miror
	! fem
	logical, intent(in) :: FileGeom
	integer, intent(in) :: Nface
	real(O), intent(in) :: rp(3,NfacePD), np(3,NfacePD), area(NfacePD)
	! new
    type(t_tree_item), optional, intent(in) :: GeomTree
	!
	integer :: Nintparam(Nparam), iparam, pint, Nintl
	real(O) :: paramG1(Nparam,NintAL*NintAL), paramG2(Nparam,NintAL*NintAL), weightsG(Nparam,NintAL*NintAL)
	real(O) :: x, y, z, r, theta, phi, param1, param2, pondere, dA, fact, nuv(3)

	if (FileGeom .and. associated(GeomTree%children) ) then
		call parser_surface_do_new(mesh, GeomTree, Nint1, Nint2)
		return
	end if

	call mesh_clear(mesh)

	if (.not. FileGeom) then
		call interpolation_list3D (TypeGeom, Nsurf, surf, Nint1, Nint2, NintAL, Nparam, &
			 Nintparam, paramG1, paramG2, weightsG, miror, Nazimutsym)
		do iparam = 1, Nparam
			Nintl = Nintparam(iparam)
			do pint = 1, Nintl
				param1  = paramG1(iparam,pint)
				param2  = paramG2(iparam,pint)
				pondere = weightsG(iparam,pint)       
				call elem_geom3D (TypeGeom, Nsurf, surf, param1, param2, iparam, r, theta, phi, dA, nuv, miror, Nazimutsym)                       
				fact = dA * pondere
				call mesh_expand(mesh, 1)
				mesh%items_count = mesh%items_count + 1
				call mesh_item_set(mesh, mesh%items_count, theta, phi, r, nuv(1), nuv(2), nuv(3), fact)
			end do
		end do 
		mesh%symmetry_rotN = Nazimutsym
		mesh%symmetry_z    = miror
	else   
		do pint = 1, Nface
			x = rp(1,pint)
			y = rp(2,pint)
			z = rp(3,pint)
			call T_cartesian_spherical (x, y, z, r, theta, phi)
			dA = area(pint)                        
			nuv(1) = sin(theta) * cos(phi) * np(1,pint) +                               &
					 sin(theta) * sin(phi) * np(2,pint) + cos(theta) * np(3,pint) 
			nuv(2) = cos(theta) * cos(phi) * np(1,pint) +                               &
					 cos(theta) * sin(phi) * np(2,pint) - sin(theta) * np(3,pint)
			nuv(3) = - sin(phi) * np(1,pint) + cos(phi) * np(2,pint)
			fact = dA
			call mesh_expand(mesh, 1)
			mesh%items_count = mesh%items_count + 1
			call mesh_item_set(mesh, mesh%items_count, theta, phi, r, nuv(1), nuv(2), nuv(3), fact)
		end do
		mesh%symmetry_rotN = 1
		mesh%symmetry_z    = .false.
	end if
end subroutine

end module