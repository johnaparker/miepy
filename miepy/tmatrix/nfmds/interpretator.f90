! **********************************************************************************
! *         READ PARAMETERS FROM COMMAND LINE AND USE IT IN INPUT FILES            *
! *   -------------------------------------------------------------------------    *
! *   -------------------------------------------------------------------------    *
module intepretator
    use utils_str

	integer, parameter :: MAXLENGTH_PARAM_NAME  = 80
	integer, parameter :: MAXLENGTH_PARAM_VALUE = 200

	type param_value
		character(MAXLENGTH_PARAM_VALUE) :: value
		type(param_value), pointer		 :: next
	end type

	type param_info
		character(MAXLENGTH_PARAM_NAME)  :: name
		integer							 :: value_count = 0
		type(param_value), pointer       :: value
	end type

	type(param_info), pointer, dimension(:) :: parameter_list
	integer :: parameter_list_count = 0

contains

! expand space in array parameter_list for count parameters
subroutine parameter_list_expand(count)
	implicit none
    !
	integer, intent(in) :: count
	!
	type(param_info), pointer, dimension(:) :: list_old
	integer :: change_count = 10
	!
	if ((parameter_list_count + count) .le. ubound(parameter_list,1)) then
		return
	end if
	change_count = MAX(change_count, count)
	if (associated(parameter_list)) then
		list_old => parameter_list
		allocate(parameter_list(1:ubound(list_old,1)+change_count))
		parameter_list(1:ubound(list_old,1)) = list_old(1:ubound(list_old,1))
		parameter_list(ubound(list_old,1)+1:ubound(parameter_list,1))%name = ""
		deallocate(list_old)
	else
		allocate(parameter_list(1:change_count))
		parameter_list(1:ubound(parameter_list,1))%name = ""
	end if
end subroutine

! find index in array for parameter with name NAME
function parameter_list_get(name, name_length) result (res)
	character(*), intent(in) :: name
	integer, intent(in), optional :: name_length
	integer :: res
	!
	integer i, j, k
	!
	if (present(name_length)) then
		k = len_c(name(1:name_length))
	else
		k = len_c(name)
	end if

	i = 1
	j = 0
	do while (i <= parameter_list_count)
		if (trim(parameter_list(i)%name) .eq. name(1:k)) then
			j = i
			exit
		end if
		i = i + 1
	end do
	res = j
end function

! add new parameter with value or add new value to exist parameter
function parameter_list_add(name, value) result (index)
    character(*), intent(in) :: name, value
    integer                  :: index
    !
    integer :: i
    !
	i = parameter_list_get(name)
	if (i <= 0) then
		call parameter_list_expand(1)
		parameter_list_count = parameter_list_count + 1
		parameter_list(parameter_list_count)%name = name(1:len_c(name))
		i = parameter_list_count
	end if
	
	if (.not. associated(parameter_list(i)%value)) then
		allocate(parameter_list(i)%value)
		parameter_list(i)%value%value(1:) = " "
		parameter_list(parameter_list_count)%value_count = 1
	end if

	parameter_list(i)%value%value = value
    index = i
end function

! load parameters from command line
subroutine parameter_load_command_line
	!
	integer :: i, j, k, param_cmd_cnt
	character(MAXLENGTH_PARAM_VALUE) :: buf
	type(param_value), pointer :: cur_value, tmp_value
	!
	param_cmd_cnt = IARGC()
	call parameter_list_expand(10)
	
	i = 1
	do while (i <= param_cmd_cnt)
		! read next value
		buf(1:len(buf)) = " "
		call GETARG(i, buf)

		if (buf(1:1) .eq. "-") then ! new parameter
			call parameter_list_expand(1)
			parameter_list_count = parameter_list_count + 1
			parameter_list(parameter_list_count)%name = "@" // buf(2:len_c(buf)) // "@"
			nullify(cur_value)
		else ! new value for current parameter
			allocate (tmp_value)
			tmp_value%value = buf
			nullify(tmp_value%next)
			if (associated(cur_value)) then
				cur_value%next => tmp_value
				parameter_list(parameter_list_count)%value_count = parameter_list(parameter_list_count)%value_count + 1
			else
				parameter_list(parameter_list_count)%value => tmp_value
				parameter_list(parameter_list_count)%value_count = 1
			end if
			cur_value => tmp_value
			nullify(tmp_value)
		end if
		i = i + 1
	end do
end subroutine

! evaluate string (replace each @variable|def_value@ to its value from parameter_list)
subroutine evaluate(buf, buf_len)
    use parameters
	character(*), intent(inout)   :: buf
	integer, intent(in), optional :: buf_len
	!
	integer i, j, k, len_value, len_from, len_to, ios
	character(MAXLENGTH_PARAM_VALUE) :: str_from, str_to, str_name
	logical is_replace
	complex(O) tmp_complex
	real(O) tmp_real
    !
	if (present(buf_len)) then
	    len_value = buf_len
	else
	    len_value = len_c(buf)
	end if

	i = 1
	j = 0
	do while (i <= len_value)
		if (buf(i:i) == "@") then
			if (j > 0) then
				str_from = trim(buf(j:i)) // char(0)
				len_from = len_c(str_from)

				is_replace = .false.

				if (str_from(2:2) /= "$") then
				    k = 1
				    do while (k <= len_from)
				        if (str_from(k:k) == "|") then
				            exit
				        end if
				        k = k + 1
				    end do
				    if (k <= len_from) then ! with default value
					    str_to   = trim(str_from(k+1:len_from-1)) // char(0)
						len_to   = len_c(str_to)
                        str_name = trim(str_from(1:k-1)) // str_from(len_from:len_from) // char(0)
						is_replace = .true.
				    else ! without default value
					    str_to   = char(0)
					    str_name = str_from
				    end if
					k = parameter_list_get(str_name)
					if (k > 0) then
						if (associated(parameter_list(k)%value)) then
							str_to = trim(parameter_list(k)%value%value) // char(0)
						else
							str_to(1:1) = char(0)
						end if
						len_to = len_c(str_to)
						is_replace = .true.
					end if
				end if
				if (is_replace) then
					call str_replace(buf, str_from, str_to)
					i = i + len_to - len_from
					len_value = len_value + len_to - len_from
					j = 0
				else
					j = i
				end if
			else
				j = i
			end if
		end if
	    i = i + 1
	end do
end subroutine

end module
