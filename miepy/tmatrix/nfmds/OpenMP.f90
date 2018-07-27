! read number of threads in OpenMP mode
subroutine load_params_command_line
	use parameters
	use intepretator
	use inputoutput

	character(400) buf
	integer i /1/
	
	call parameter_load_command_line

	! num_threads
	i = parameter_list_get("@num_threads@" // char(0))
	if (i>0) then
		if (parameter_list(i)%value_count > 0) then
			buf = parameter_list(i)%value%value
			read (buf,*) OMP_thread_cnt
		else
			OMP_thread_cnt = 1
		end if
	else
		OMP_thread_cnt = 1
	end if

	! input
	InputCurrentUnit = InputConUnit
	i = parameter_list_get("@input@" // char(0))
	if (i>0) then
		if (parameter_list(i)%value_count > 0) then
			buf(1:len(buf)) = " "
			buf = parameter_list(i)%value%value
			call str_replace(buf,"\n" // char(0), char(10) // char(0))
			! pishem vo vremennij  fail 
			open (unit=InputFileUnit, file=trim(PathTEMP) // "input.txt", status='REPLACE')
			write (InputFileUnit, *) trim(buf)
			close(unit=InputFileUnit)
			! otkrivaem konsol iz vremennogo faila
			open (unit=InputFileUnit, file=trim(PathTEMP) // "input.txt", status='OLD')
			InputCurrentUnit = InputFileUnit
            close(InputFileUnit, status='delete')
		end if
	end if
	
	! output path
	i = parameter_list_get("@output@" // char(0))
	if (i>0) then
		if (parameter_list(i)%value_count > 0) then
		    FileOutput = trim(parameter_list(i)%value%value)
		end if
	end if
end subroutine
