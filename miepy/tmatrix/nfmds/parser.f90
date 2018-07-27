! **********************************************************************************
! *                  PARSE TEXT FILE INTO TREE STRUCTURE                           *
! *   -------------------------------------------------------------------------    *
! *   -------------------------------------------------------------------------    *
module parser
use utils_str
implicit none

! buffer
character, pointer, dimension(:) :: buffer => NULL()
integer	    					 :: buffer_size = 0

public parser_buffer_clear, parser_buffer_expand, parser_load_file, parser_add
public buffer, buffer_size

! syntax
public read_main
private isCharInSet, isCharSpace, isCharName, read_space, read_space_ext, &
        read_comment, read_name, read_field_value, read_value, read_variable

! tree structure
type t_tree_item
	character(80) :: name
	type(t_tree), pointer            :: children => NULL()
	character, pointer, dimension(:) :: value => NULL()
end type

type t_tree
	type(t_tree_item), pointer, dimension(:) :: items => NULL()
	integer	 :: items_count = 0
end type

public  tree_clear, tree_expand, tree_item_expand, tree_item_set_value, &
        tree_item_get_value, tree_get_item, tree_find, tree_find_value

contains

! **********************************************************************************
! *                                  BUFFER                                        *
! **********************************************************************************
! clear buffer
subroutine parser_buffer_clear()
	if (associated(buffer)) then
	    deallocate(buffer)
    end if
	buffer_size = 0
end subroutine

! expand space in buffer for count characters
subroutine parser_buffer_expand(count)
    integer, intent(in) :: count
	!
	character, pointer, dimension(:) :: old_buffer
	integer :: change_count = 1024
	!
	if ((buffer_size + count) .le. ubound(buffer,1)) then
		return
	end if
	change_count = MAX(change_count, count)
	if (associated(buffer)) then
		old_buffer => buffer
		allocate(buffer(1:ubound(old_buffer,1)+change_count))
		buffer(1:ubound(old_buffer,1)) = old_buffer(1:ubound(old_buffer,1))
		old_buffer(ubound(old_buffer,1)+1:ubound(buffer,1)) = " "
		deallocate(old_buffer)
	else
		allocate(buffer(1:change_count))
		buffer(1:ubound(buffer,1)) = " "
	end if
end subroutine

! load in buffer text file to parse
function parser_load_file(file_name) result (ios)
!    use dfport
	!
	character(*), intent(in) :: file_name
	integer					 :: ios
	!
	integer                  :: i, file_size
	integer, dimension(1:13) :: stat_file
	!
	call parser_buffer_clear()

!	ios = stat(file_name, stat_file)
!	if (ios .ne. 0) return
!	file_size = stat_file(8)
	file_size = 1000

	open(unit=99, file=file_name, status='OLD', access='DIRECT', recl=file_size, iostat=ios)
	if (ios .ne. 0) return

    call parser_buffer_expand(file_size)

	read(unit=99, rec=1, iostat=ios) (buffer(i),i=1,file_size)
	if (ios .ne. 0) return
	
	buffer_size = file_size

	close(unit=99)
end function

! append string into buffer
function parser_add(buf, buf_len) result (ios)
	character(*), intent(in) :: buf
	integer, intent(in), optional :: buf_len
	integer	:: ios
	!
	integer :: i, j, count
	!
	if (present(buf_len)) then
		count = buf_len
	else
		count = len_c(buf)
	end if
	call parser_buffer_expand(count)

    do i = 1,count
	    buffer(i+buffer_size) = buf(i:i)
	end do
	buffer_size = buffer_size + count
	ios = 0
end function

! **********************************************************************************
! *                                  SYNTAX                                        *
! **********************************************************************************
function isCharInSet(ch, set, set_length) result (res)
    character, intent(in)    :: ch
	character(*), intent(in) :: set
	integer, intent(in)      :: set_length
	logical :: res
	!
	integer :: i

    do i=1, set_length
	    if (set(i:i) == ch) then
		    res = .true.
			return
        end if
    end do
	res = .false.
end function

function isCharSpace(ch) result (res)
    character, intent(in) :: ch
	logical               :: res
    
    res = isCharInSet(ch," " // char(9),2)
end function

function isCharName(ch) result (res)
    character, intent(in) :: ch
	logical               :: res
	
	res = ((ch >= 'a') .and. (ch <= 'z')) .or. & 
	      ((ch >= 'A') .and. (ch <= 'Z')) .or. &
		  ((ch >= '0') .and. (ch <= '9')) .or. &
		  isCharInSet(ch,"_",1)
end function

function read_space(cur_index) result (index)
    integer, intent(in) :: cur_index
	integer :: index
	
	index = cur_index
	do while ((index <= buffer_size) .and. isCharSpace(buffer(index)))
	    index = index + 1
    end do
end function

function read_space_ext(cur_index) result (index)
    integer, intent(in) :: cur_index
	integer :: index

    index = cur_index
	do while ((index <= buffer_size) .and. (isCharInSet(buffer(index),char(13)//char(10),2).or.isCharSpace(buffer(index))) )
	    index = index + 1
    end do
end function

function read_comment(cur_index) result (index)
    integer, intent(in) :: cur_index
	integer :: index
	
	index = cur_index
	if (buffer(index) == '#') then
	    do while ((index <= buffer_size) .and. (.not. isCharInSet(buffer(index),char(13)//char(10),2)))
		    index = index + 1
        end do
    end if
end function

function read_name(cur_index) result (index)
    integer, intent(in) :: cur_index
	integer :: index
	
	index = cur_index
	do while ((index <= buffer_size) .and. isCharName(buffer(index)))
	    index = index + 1
    end do
end function

function read_field_value(cur_index) result (index)
    integer, intent(in) :: cur_index
	integer :: index
	
	index = cur_index
	if (buffer(index) .ne. '/') then
	    return
    end if

	index = index + 1
	do while ((index <= buffer_size) .and. (.not. isCharInSet(buffer(index),"/"//char(13)//char(10),3)))
	    index = index + 1
    end do
	index = index + 1
end function

function read_value(cur_index) result (index)
    integer, intent(in) :: cur_index
	integer :: index
	
	index = cur_index
	do while ((index <= buffer_size) .and. (.not. isCharInSet(buffer(index),"#"//char(13)//char(10),3)))
	    index = index + 1
    end do
end function

recursive function read_variable(tree, cur_index) result(index)
    integer, intent(in)  :: cur_index
	type(t_tree_item) :: tree
	integer              :: index
	!
	integer		   :: i1,i2,i3, j,k
    character(100) :: buf
	logical		   :: is_record, flag
	
	index = cur_index

    ! field "typed_value"
	is_record = .false.
	i2 = read_space(index)
	if (buffer(i2) == ':') then
	    i1 = read_space(i2+1)
		i2 = read_name(i1)
		do while (i2 .ne. i1)
		    ! name of field
			call array_to_char(buffer(i1:i2-1), buf)
				
			call tree_item_expand(tree)
			call tree_expand(tree%children, 1)
			k = tree%children%items_count + 1
			tree%children%items_count = k
			tree%children%items(k)%name = buf
#ifdef PARSE_DEBUG
			print *, "Field = ", trim(buf)
#endif

			is_record = is_record .or. (trim(buf) == "record")
			! value of field
			i3 = read_field_value(i2)
			if (i2 .ne. i3) then
			    call array_to_char(buffer(i2+1:i3-2), buf)

				call tree_item_set_value(tree%children%items(k), buf)
#ifdef PARSE_DEBUG
				print *, "Value = ", trim(buf)
#endif
				i1 = i3
			else
			    i1 = i2
			end if
            ! next
			i1 = read_space(i1)
			if (buffer(i1) == ',') then
			    i1 = read_space(i1+1)
				i2 = read_name(i1)
			else
			    i2 = i1
            end if
		end do
	end if

	! field "record_content"
	i1 = read_space(i2)
	if ((buffer(i1) == '=') .and. (.not. is_record) ) then
	    i1 = read_space(i1+1)
		i2 = read_value(i1)
		! value of variable
		call array_to_char(buffer(i1:i2-1), buf)

		call tree_item_set_value(tree, buf)
#ifdef PARSE_DEBUG
		print *, "Variable Value = ", trim(buf)
#endif

	else if ((buffer(i1) == '=') .and. (is_record) ) then
	    i1 = i1 + 1
		i1 = read_space_ext(i1)
		i2 = read_name(i1)
#ifdef PARSE_DEBUG
		print *, "RECORD"
#endif
		do while (i1 .ne. i2)
    		! name of item
	    	call array_to_char(buffer(i1:i2-1), buf)

            call tree_item_expand(tree)
    		call tree_expand(tree%children, 1)
	    	k = tree%children%items_count + 1
		    tree%children%items_count = k
    		tree%children%items(k)%name = buf
#ifdef PARSE_DEBUG
		    print *, "Variable Name = ", trim(buf)
#endif
		
	    	if ((trim(buf) .ne. "done").and.(i1 .ne. i2)) then
		        i3 = read_variable(tree%children%items(k), i2)
	    	    do while (i3 <= buffer_size)
	    	        i3 = read_space_ext(i3)
	    	        i1 = read_comment(i3)
	    	        if (i3 .eq. i1) then
	    	            exit
	    	        else
	    	            i3 = i1
	    	        end if
	    	    end do
	    		i2 = read_name(i1)
	        else
    	        i1 = i2
		    end if
	    end do
#ifdef PARSE_DEBUG
	    print *, "DONE"
#endif
    end if
	index = i2		
end function

function read_main(tree, cur_index) result(index)
    integer, intent(in) :: cur_index
	type(t_tree_item), intent(inout) :: tree
	integer :: index
	!
	integer :: i1, i2
	character(100) :: buf
	
	if ((buffer_size == 0) .or. (cur_index > buffer_size)) then
	    return
	end if
	index = cur_index
	i1 = read_space_ext(index)
	i2 = read_name(i1)
	if (i2 == i1) then
	    return
    end if

	! name of main variable
	call array_to_char(buffer(i1:i2-1), buf)
	tree%name = buf

	index = read_variable(tree, i2)
end function
	
! **********************************************************************************
! *                                  TREE STRUCTURE                                *
! **********************************************************************************
! remove all children of tree node
recursive subroutine tree_clear(tree, istep)
	type(t_tree), pointer :: tree
	integer, optional     :: istep
	!
	integer :: i
	
	if (.not. associated(tree)) return
	if (present(istep) .and. (istep > 1000)) return
	
	do i=1,tree%items_count
	    call tree_clear(tree%items(i)%children, istep)
	    if (associated(tree%items(i)%value)) then
	        deallocate (tree%items(i)%value)
	    end if
	end do
	deallocate(tree%items)
	tree%items_count = 0
end subroutine

! expand space for count children of tree node 
subroutine tree_expand(tree, count)
	integer, intent(in) :: count
	type(t_tree), intent(inout) :: tree
	!
	type(t_tree_item), pointer, dimension(:) :: old_items
	integer :: change_count = 5, i
	!
	change_count = MAX(change_count, count)
	if (associated(tree%items)) then
		if ((tree%items_count + count) .le. ubound(tree%items,1)) then
			return
		end if

		old_items => tree%items
		allocate(tree%items(1:ubound(old_items,1)+change_count))
		tree%items(1:ubound(old_items,1)) = old_items(1:ubound(old_items,1))
		tree%items(ubound(old_items,1)+1:ubound(tree%items,1))%name = ""
		do i=ubound(old_items,1)+1,ubound(tree%items,1)
		    nullify(tree%items(i)%value)
		end do
		deallocate(old_items)
	else
		allocate(tree%items(1:change_count))
		tree%items(1:ubound(tree%items,1))%name = ""
		do i=1,ubound(tree%items,1)
		    nullify(tree%items(i)%value)
		end do
	end if
end subroutine

! add ability of tree item to have children
subroutine tree_item_expand(tree_item)
	type(t_tree_item), intent(inout) :: tree_item
	!
	if (.not. associated(tree_item%children)) then
		allocate(tree_item%children)
	end if
end subroutine

! set value of tree node
subroutine tree_item_set_value(tree_item, value)
	type(t_tree_item), intent(inout) :: tree_item
	character(*), intent(in) :: value
	!
	integer :: i, j
	!
	i = len_c(value)
	if (associated(tree_item%value)) then
		if (ubound(tree_item%value,1) <= i) then
			deallocate(tree_item%value)
		end if
	end if
	if (.not.associated(tree_item%value)) then
		allocate(tree_item%value(1:MAX(i,30)))
	end if
	
	do j=1,i
		tree_item%value(j) = value(j:j)
	end do
	tree_item%value(i+1:) = " "
end subroutine

! get value of tree node
subroutine tree_item_get_value(tree_item, value)
	type(t_tree_item), intent(inout) :: tree_item
	character(*) :: value
	!
	integer :: i,j
	
	if (associated(tree_item%value)) then
		i = ubound(tree_item%value,1)
		do j=1,i
			value(j:j) = tree_item%value(j)
		end do
		do j=i+1,len(value)
			value(j:j) = " "
		end do
	else
	    value = repeat(" ", len(value))
	end if
end subroutine

! find children of tree node with fixed name
function tree_get_item(tree, name) result(index)
	type(t_tree), pointer :: tree
	integer :: index
	character(*), intent(in) :: name
	!
	integer :: i,j, name_len
	
	if (associated(tree)) then
		i = 1
		name_len = len_c(name)
		do while (i <= tree%items_count)
			if (trim(tree%items(i)%name) == name(1:name_len)) then
				index = i
				return
			end if
			i = i + 1
		end do
	end if
	index = 0
end function

! parse path to tree node and return it
subroutine tree_find(tree, path, res)
	character(*), intent(in)   :: path
	type(t_tree_item), target  :: tree
	type(t_tree_item), pointer :: res
	!
	integer :: i, i_len, j, k
	character(80) :: name
	!
	res => tree
	i = 1
	i_len = len_c(path)
	do while (i <= i_len)
		! get next name
		j = i
		do while (j <= i_len)
		    if (path(j:j) == '/') then
		        exit
		    end if
			j = j + 1
		end do
		name = path(i:j-1)
		name(j-i+1:) = " "

		! get item
		k = tree_get_item(res%children, name)
		if (k <= 0) then
			res => tree
			return
		end if

		! next item
		res => res%children%items(k)
		i = j + 1
	end do
end subroutine

! parse path to tree node and return its value
subroutine tree_find_value(tree, path, res)
	character(*), intent(in)   :: path
	type(t_tree_item), target  :: tree
	character(*)               :: res
	!
	type(t_tree_item), pointer :: tree_tmp

	call tree_find(tree, path, tree_tmp)
	if (associated(tree_tmp)) then
		call tree_item_get_value(tree_tmp, res)
	else
	    res = repeat(" ", len(res))
	end if
end subroutine

	
end module
