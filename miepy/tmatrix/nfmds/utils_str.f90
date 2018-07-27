! **********************************************************************************
! *                  SOME UTILS FOR STRINGS                                        *
! *   -------------------------------------------------------------------------    *
! *  len_c                                                                         *
! *  str_replace                                                                   *
! *   -------------------------------------------------------------------------    *
module utils_str
contains

! length of string with #0 as EOL character or len(trim(buf))
integer function len_c(buf)
	character(*), intent(in) :: buf
	!
	integer :: i
	!
	i = 1
	do while (i <= len(buf))
		if (buf(i:i) == char(0)) then
			len_c = i - 1
			return
		end if
		i = i + 1
	end do
	i = len(buf)
	do while (i >= 1)
	    if (buf(i:i) .ne. " ") then
	        len_c = i
	        return
	    end if
		i = i - 1
	end do
	len_c = 0
end function

! replace in buf str_from to str_to
subroutine str_replace(buf, str_from, str_to, str_from_len_def, str_to_len_def, buf_len_def)
	character(*), intent(inout) :: buf
	character(*), intent(in)    :: str_from, str_to
	integer, intent(in), optional :: str_from_len_def, str_to_len_def, buf_len_def
	!
	integer i, j, str_from_len, str_to_len, buf_len
	!

	if (present(str_from_len_def)) then
		str_from_len  = str_from_len_def
	else
		str_from_len  = len_c(str_from)
	end if
	if (present(str_to_len_def)) then
		str_to_len  = str_to_len_def
	else
		str_to_len  = len_c(str_to)
	end if
	if (present(buf_len_def)) then
		buf_len     =buf_len_def
	else
		buf_len     =len_c(buf)
	end if
	
	i = 1
	j = 1
	do while (i <= buf_len)
		if (buf(i:i) == str_from(j:j)) then
			j = j + 1
			if (j .gt. str_from_len) then
				buf = buf(:i-j+1) // str_to(1:str_to_len) // buf(i+1:)
				i = i + str_to_len - str_from_len
				buf_len = buf_len + str_to_len - str_from_len
				j = 1
			end if
		else
			i = i - j + 1
			j = 1
		end if
		i = i + 1
	end do
end subroutine

! convert array of char to string
subroutine array_to_char(arr, buf)
    character, dimension(:), intent(in) :: arr
	character(*), intent(inout) :: buf
	integer :: i
	do i=1,ubound(arr,1)
	    buf(i:i) = arr(i)
    end do
	buf(ubound(arr,1)+1:) = ' '
end subroutine


end module utils_str