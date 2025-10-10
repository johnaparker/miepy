! **********************************************************************************
! *                  ROUTINES FOR READING THE INPUT DATA                           *
! *   -------------------------------------------------------------------------    *
! *   Partial list of subroutines:                                                 *
! *     read_FileFEM,             read_FileFEMAxsym,        read_HeadFileTmat,     *
! *     read_FileTmat,            read_HeadFileTmatVct,     read_FileTmatVct,      *
! *     read_FileSmat(UNUSED),    read_Tmatrix,             read_integer,          *
! *     read_integerbound,        read_integer2,            read_integer3,         *
! *     read_real,                read_logical,             read_char80,           *
! *     LenString,                XFindPar,                 LcString               *
! **********************************************************************************
subroutine read_FileFEM (filein_name, face_num, face_point, face_normal, face_area)
!-----------------------------------------------------------------------------------
! The routine read the surface parameters from the file filein_name.               !
!-----------------------------------------------------------------------------------  
  use parameters
  implicit none
  integer        :: face_num
  real(O)        :: face_point(3,NfacePD), face_normal(3,NfacePD), face_area(NfacePD)
  character(80)  :: filein_name
!
  integer        :: nfaces, NVvr, nelem, nf, nv, ielem, i, ios
  real(O)        :: r(3), n(3), area
!  
  open (unit = iFEM, file = filein_name, status = "old", position = "rewind")      
  read (iFEM, "(i7)", iostat = ios)  nfaces
  if (ios < 0 ) then
    print "(/,2x,'EOF by reading the number of surfaces from the geometry file;')"
    stop
  else if (ios > 0) then
    print "(/,2x,'Error by reading the number of surfaces from the geometry file;')"
    stop
  end if      
  nelem = 0
  do nf = 1, nfaces      
    read (iFEM, "(i7)", iostat = ios) NVvr
    if (ios < 0 ) then
      print "(/,2x,'EOF by reading the number of vertices from the geometry file;')"
      stop
    else if (ios > 0) then
      print "(/,2x,'Error by reading the number of vertices from the geometry file;')"
      stop
    end if                
    do nv = 1, NVvr                   
      nelem = nelem + 1
      call check_MaxNface (nelem)               
      read (iFEM, "(i7,2x,7(e15.7,2x))", iostat = ios) ielem, r(1), r(2), r(3),     &
            n(1), n(2), n(3), area            
      if (ios /= 0) exit
      do i = 1, 3
        face_point (i,nelem) = r(i)
        face_normal(i,nelem) = n(i)
      end do
      face_area(nelem) = area
    end do
    if (ios < 0 ) then
      print "(/,2x,'EOF detected during the reading of the geometry file;')"
      stop
    else if (ios > 0) then
      print "(/,2x,'Error during the reading of the geometry file;')"
      stop
    end if          
  end do
  face_num = nelem                
  close (unit = iFEM)
end subroutine read_FileFEM
! **********************************************************************************
subroutine read_FileFEMAxsym (filein_name, line_num, line_point, line_normal,       &
           line_area)
!-----------------------------------------------------------------------------------
! The routine reads the curve parameters from the file filein_name.                !
!-----------------------------------------------------------------------------------  
  use parameters
  implicit none
  integer        :: line_num
  real(O)        :: line_point(2,NfacePD), line_normal(2,NfacePD), line_area(NfacePD)
  character(80)  :: filein_name
!
  integer        :: nlines, NVvr, nelem, nf, nv, ielem, i, ios
  real(O)        :: r(2), n(2), area
!  
  open (unit = iFEM, file = filein_name, status = "old", position = "rewind")      
  read (iFEM, "(i7)", iostat = ios)  nlines
  if (ios < 0 ) then
    print "(/,2x,'EOF by reading the number of curves from the geometry file;')"
    stop
  else if (ios > 0) then
    print "(/,2x,'Error by reading the number of curves from the geometry file;')"
    stop
  end if
  nelem = 0
  do nf = 1, nlines      
    read  (iFEM, "(i7)", iostat = ios)  NVvr
    if (ios < 0 ) then
      print "(/,2x,'EOF by reading the number of vertices from the geometry file;')"
      stop
    else if (ios > 0) then
      print "(/,2x,'Error by reading the number of vertices from the geometry file;')"
      stop
    end if
    do nv = 1, NVvr                   
      nelem = nelem + 1
      call check_MaxNface (nelem)               
      read (iFEM, "(i7,2x,5(e15.7,2x))", iostat = ios) ielem, r(1), r(2),           &
            n(1), n(2), area          
      if (ios /= 0) exit
      do i = 1, 2
        line_point (i,nelem) = r(i)    
        line_normal(i,nelem) = n(i)                 
      end do
      line_area(nelem) = area                  
    end do
    if (ios < 0 ) then
      print "(/,2x,'EOF detected during the reading of the geometry file;')"
      stop
    else if (ios > 0) then
      print "(/,2x,'Error during the reading of the geometry file;')"
      stop
    end if          
  end do
  line_num = nelem        
  close (unit = iFEM)
end subroutine read_FileFEMAxsym
! **********************************************************************************
subroutine read_HeadFileTmat (nt, mt)
  use parameters
  implicit none
  integer       :: nt, mt, ios
  character(85) :: string
!
  read (iTmat, "(2x,a)", iostat = ios) string
  if (ios < 0 ) then
    print "(/,2x,'EOF detected during the reading of the T - matrix file;')"
    stop
  else if (ios > 0) then
    print "(/,2x,'Error by reading the string: Half - Dimensions of the T Matrix;')"
    stop
  end if  
  read (iTmat, "(2x,i10,2x,i10)", iostat = ios) nt, mt
  if (ios < 0 ) then
    print "(/,2x,'EOF detected during the reading of the T - matrix file;')"
    stop
  else if (ios > 0) then
    print "(/,2x,'Error during the reading of the dimensions of the T matrix;')"
    stop
  end if 
  read (iTmat, "(2x,a)",iostat = ios) string  
  if (ios < 0 ) then
    print "(/,2x,'EOF detected during the reading of the T - matrix file;')"
    stop
  else if (ios > 0) then
    print "(/,2x,'Error during the reading of the string: T Matrix;')"
    stop
  end if 
end subroutine read_HeadFileTmat
! **********************************************************************************
subroutine read_FileTmat (nt, mt, a)
  use parameters
  implicit none
  integer    :: nt, mt, i, j, ios
  complex(O) :: a(2*nt,2*mt)
!    
  i = 0
  do while (i < 2*nt)
    i = i + 1   
    read (iTmat,"(10(2x,1pe24.15,1x,1pe24.15),/)", iostat = ios) (a(i,j), j = 1, 2*mt)    
    if (ios /= 0) exit          
  end do  
  if (ios < 0 ) then
    print "(/,2x,'EOF detected during the reading of the T matrix;')"
    stop
  else if (ios > 0) then
    print "(/,2x,'Error during the reading of the T matrix;')"
    stop
  end if
end subroutine read_FileTmat
! **********************************************************************************
subroutine read_HeadFileTmatVct (nt)
  use parameters
  implicit none
  integer       :: nt, ios
  character(85) :: string
!
  read (iTmat, "(2x,a)",iostat = ios) string
  if (ios < 0 ) then
    print "(/,2x,'EOF detected during the reading of the T - vector file;')"
    stop
  else if (ios > 0) then
    print "(/,2x,'Error by reading the string: Half - Dimensions of the T Vector;')"
    stop
  end if   
  read (iTmat, "(2x,i10)",iostat = ios) nt
  if (ios < 0 ) then
    print "(/,2x,'EOF detected during the reading of the T - vector file;')"
    stop
  else if (ios > 0) then
    print "(/,2x,'Error during the reading of the dimensions of the T vector;')"
    stop
  end if 
  read (iTmat, "(2x,a)",iostat = ios) string
  if (ios < 0 ) then
    print "(/,2x,'EOF detected during the reading of the T - vector file;')"
    stop
  else if (ios > 0) then
    print "(/,2x,'Error during the reading of the string: T Vector;')"
    stop
  end if 
end subroutine read_HeadFileTmatVct
! **********************************************************************************
subroutine read_FileTmatVct (nt, v)
  use parameters
  implicit none
  integer    :: nt, i, ios
  complex(O) :: v(2*nt)
!
  i = 0
  do while (i < 2*nt)
    i = i + 1
    read (iTmat,"(2x,1pe24.15,1x,1pe24.15)",iostat = ios) v(i)
    if (ios /= 0) exit          
  end do  
  if (ios < 0 ) then
    print "(/,2x,'EOF detected during the reading of the T vector;')"
    stop
  else if (ios > 0) then
    print "(/,2x,'Error during the reading of the T vector;')"
    stop
  end if
end subroutine read_FileTmatVct
! **********************************************************************************
subroutine read_FileSmat (N, Nteta, S, Cext, Cscat, Qext, Qscat)
  use parameters
  implicit none
  integer    :: N, Nteta, i, j, ios
  real(O)    :: Cext, Cscat, Qext, Qscat
  complex(O) :: S(N,Nteta)
!    
  i = 0
  do while (i < N)
    i = i + 1   
    read (iSS,"(10(2x,1pe24.15,1x,1pe24.15),/)", iostat = ios) (S(i,j), j = 1, Nteta)    
    if (ios /= 0) exit          
  end do  
  if (ios < 0 ) then
    print "(/,2x,'EOF detected during the reading of the average quantities <SS*>;')"
    stop
  else if (ios > 0) then
    print "(/,2x,'Error during the reading of the average quantities <SS*>;')"
    stop
  end if
  read (iSS,"(4(2x,1pe13.4))", iostat = ios) Cext, Cscat, Qext, Qscat 
  if (ios < 0 ) then
    print "(/,2x,'EOF detected during the reading of cross sections;')"
    stop
  else if (ios > 0) then
    print "(/,2x,'Error during the reading of cross sections;')"
    stop
  end if
end subroutine read_FileSmat
! **********************************************************************************  
subroutine read_Tmatrix (chiral, Nrank, Mrank, Nmaxmax, TG, ntg, mtg)
  use parameters
  implicit none  
  integer       :: Mrank, Nrank, Nmaxmax, ntg, mtg
  complex(O)    :: TG(2*ntg,2*mtg)
  logical       :: chiral
!
  integer       :: m, Nmax, i, j, l, N0, ntl, mtl
  complex(O),allocatable :: TL(:,:)
!  
  call read_HeadFileTmat (ntl, mtl)
  call check_dimensionMat (ntl, mtl, Nrank)             
  allocate (TL(2*ntl, 2*mtl))   
  do m = 0, Mrank
    if (m == 0) then
      Nmax = Nrank
    else
      Nmax = Nrank - m + 1
    end if                                          
    call read_FileTmat (ntl, mtl, TL)
    if (m == 0) then
      do i = 1, Nmax
        do j = 1, Nmax
          TG(i,j)                = TL(i,j)
          TG(i,j+Nmaxmax)        = TL(i,j+Nmax)
          TG(i+Nmaxmax,j)        = TL(i+Nmax,j)
          TG(i+Nmaxmax,j+Nmaxmax)= TL(i+Nmax,j+Nmax)
        end do          
      end do
    else        
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      do l = 1, 2
        do i = 1, Nmax
          do j = 1, Nmax            
            TG(i+N0,j+N0)                 = TL(i,j)
            TG(i+N0,j+N0+Nmaxmax)         = TL(i,j+Nmax)
            TG(i+N0+Nmaxmax,j+N0)         = TL(i+Nmax,j)
            TG(i+N0+Nmaxmax,j+N0+Nmaxmax) = TL(i+Nmax,j+Nmax)
          end do            
        end do
        N0 = N0 + Nmax
        if (.not. chiral) then
          call matrix_m_negativ (Nmax, Nmax, TL, ntl, mtl)
        else              
          call read_FileTmat (ntl, mtl, TL)
        end if
      end do
    end if                                                                                           
  end do
  deallocate (TL)
end subroutine read_Tmatrix
! **********************************************************************************
recursive subroutine read_integer (int)
  implicit none
  integer  :: int, ierr
!
  read (*, *, iostat = ierr) int
  if (ierr /= 0) then
    print "(/,2x,'Input error during the read statement;')"
    print "(  2x,'- enter the integer variable again;')"
    call read_integer (int)
  end if
end subroutine read_integer
! **********************************************************************************
recursive subroutine read_integerbound (int, imin, imax)
  implicit none
  integer  :: int, imin, imax, ierr
!
  read (*, *, iostat = ierr) int
  if (ierr /= 0) then
    print "(/,2x,'Input error during the read statement;')"
    print "(  2x,'- enter the integer variable again;')"
    call read_integer (int)
  end if
  if (int < imin .or. int > imax) then
    print "(/,2x,'Input error: the integer variable is out of bounds;')"
    print "(  2x,'- enter the integer variable again;')"
    call read_integer (int)
  end if
end subroutine read_integerbound
! **********************************************************************************
recursive subroutine read_integer2 (int1, int2)
  implicit none
  integer  :: int1, int2, ierr
!
  read (*, *, iostat = ierr) int1, int2
  if (ierr /= 0) then
    print "(/,2x,'Input error during the read statement;')"
    print "(  2x,'- enter the integer variables again;')"
    call read_integer2 (int1, int2)
  end if
end subroutine read_integer2
! **********************************************************************************
recursive subroutine read_integer3 (int1, int2, int3)
  implicit none
  integer  :: int1, int2, int3, ierr
!
  read (*, *, iostat = ierr) int1, int2, int3
  if (ierr /= 0) then
    print "(/,2x,'Input error during the read statement;')"
    print "(  2x,'- enter the integer variables again;')"
    call read_integer3 (int1, int2, int3)
  end if
end subroutine read_integer3
! **********************************************************************************
recursive subroutine read_real (ral)
  use parameters
  implicit none
  integer  :: ierr
  real(O)  :: ral
!
  read (*, *, iostat = ierr) ral
  if (ierr /= 0) then
    print "(/,2x,'Input error during the read statement;')"
    print "(  2x,'- enter the real variable again;')"
    call read_real (ral)
  end if
end subroutine read_real
! **********************************************************************************
recursive subroutine read_logical (log)
  implicit none
  integer  :: ierr
  logical  :: log
!
  read (*, *, iostat = ierr) log
  if (ierr /= 0) then
    print "(/,2x,'Input error during the read statement;')"
    print "(  2x,'- enter the logical variable again;')"
    call read_logical (log)
  end if
end subroutine read_logical
! **********************************************************************************
recursive subroutine read_char80 (TypePath, char)
  use parameters
  implicit none
  integer        :: TypePath, ierr, LenString
  character(80)  :: char, charend
!
  read (*, *, iostat = ierr) charend
  if (ierr /= 0) then
    print "(/,2x,'Input error during the read statement;')"
    print "(  2x,'- enter the character type variable again;')"
    call read_char80 (TypePath, char)
  end if
  if (TypePath == 1) then
    char = PathOUTPUT(1:LenString(PathOUTPUT)) // charend(1:LenString(charend))
  else if (TypePath == 2) then
    char = PathTEMP(1:LenString(PathTEMP)) // charend(1:LenString(charend))
  else if (TypePath == 3) then
    char = PathGEOM(1:LenString(PathGEOM)) // charend(1:LenString(charend))
  end if
end subroutine read_char80
! **********************************************************************************
integer function LenString (str) 
  implicit none
  character(*) :: str
!
  LenString = len(str)
  do while (LenString > 0 .and. str(LenString:LenString) == ' ')
    LenString = LenString - 1
  end do   
end function LenString 
!***********************************************************************************        
logical function XFindPar ( nunit, parname )        
  implicit none
  integer       :: nunit  
  character(80) :: parname
!    
  integer       :: i, parlen, ios
  logical       :: more, secondsearch
  character(80) :: name, line
!
! --- length of the parameter name ---
  parlen = 0
  do i = len (parname), 1, - 1
    if (parname (i : i) /= ' ' ) then
      parlen = i      
      exit
    end if
  end do
  if (parlen == 0) then
    XFindPar = .true.
  else          
    name = parname (1 : parlen)
    call LcString ( name (1 : parlen) )    
    more = .true.
    secondsearch = .false.
    do while (more)
      line = ' '
      read (nunit,"(a)", iostat = ios) line
      if (ios < 0) then
        if (secondsearch) then
          XFindPar = .false.
          more     = .false.
        else
          secondsearch = .true.
        end if
      else	  	       	   	      			
        call LcString ( line (1 : parlen) )      
        if (line (1 : parlen) == name (1 : parlen)) then
          XFindPar = .true.
          more     = .false.  	    
        end if
      end if        
    end do                         
  end if       
end function XFindPar
!***********************************************************************************                
subroutine LcString ( string )
  implicit none
  character(*) :: string
  integer      :: ccc, i, offset
!
  offset = ichar ( 'a' ) - ichar ( 'a' )
  do i = 1, len (string)
    ccc = ichar ( string (i : i) )
    if (ccc >= ichar ( 'a' ) .and. ccc <= ichar ( 'z' )) then
      string ( i : i ) = char ( ccc + offset )
    end if
  end do       
end subroutine LcString       
! **********************************************************************************
! *                  ROUTINES FOR WRITING THE OUTPUT DATA                          *
! *    ------------------------------------------------------------------------    *
! *    Partial list pf subroutines:                                                *
! *      write_progress,        write_progress_m,        write_progress_low,       *
! *      write_HeadFileTmat,    write_FileTmat,          write_InfoFileTmat        *
! *      write_HeadFileTmatVct, write_FileTmatVct,       write_FileSmat(UNUSED),   *
! *      write_DSCS,            write_4ConvParam,        write_3ConvParam,         *
! *      write_2ConvParam,      write_2ConvParamAxsym,   write_1ConvParam,         *
! *      write_Effic,           write_NrankConvRes,      write_NintConvRes,        *
! *      write_MrankConvRes,    write_TypeConvHead,      write_3ConvParamReg,      *
! *      write_2ConvParamReg,   write_3ConvParamRegCOMP, write_2ConvParamRegCOMP,  *
! *      write_DSCS1,           write_SCAT,              write_ScatMatrix          *
! *      inputDSCS_SCAT,        inputDSCS_SCATAvrg,      inputDSCS_SCATAvrgSPH,    *
! *      write_ScatMatrixAvrg,  inputDSCS_EMF,           write_DSCSPARTSUB,        *
! *      write_DSCSPARTSUB1,    write_EMF                                          *
! **********************************************************************************
subroutine write_progress (start, ip, Np)
  implicit none
  integer       :: ip, Np
  logical       :: start
!
  if (start) print "(/,2x,'progress of main calculation:')"             
  print "(2x,'- ',i3,'  / ',i3,';')", ip, Np
end subroutine write_progress
! **********************************************************************************
subroutine write_progress_m (start, m, ip, Np)
  implicit none
  integer       :: m, ip, Np
  logical       :: start
!
  if (start)                                                                        &
    print "(/,2x,'progress of main calculation for azimuthal mode m =', i3,':')", m           
  print "(2x,'- ',i3,'  / ',i3,';')", ip, Np
end subroutine write_progress_m
! **********************************************************************************
subroutine write_progress_low
!
  print "(/,2x,'Calculations for the lower order configurations;')"  
end subroutine write_progress_low
! **********************************************************************************
subroutine write_HeadFileTmat (nt, mt)
  use parameters
  implicit none
  integer   :: nt, mt
!
  write (iTmat, "(2x,'Half - Dimensions of the T Matrix:')")
  write (iTmat, "(2x,i10,2x,i10)") nt, mt
  write (iTmat, "(2x,'T Matrix:')")
end subroutine write_HeadFileTmat
! **********************************************************************************
subroutine write_FileTmat (nt, mt, a)
  use parameters
  implicit none
  integer    :: nt, mt, i, j
  complex(O) :: a(2*nt,2*mt)
!
  do i = 1, 2*nt
    write (iTmat,"(10(2x,1pe24.15,1x,1pe24.15),/)") (a(i,j), j = 1, 2*mt)
  end do  
end subroutine write_FileTmat
! **********************************************************************************
subroutine write_InfoFileTmat (FileTmat, Mrank, Nrank, axsym, sphere, chiral)
  use parameters
  implicit none
  integer        :: Nrank, Mrank, LenString       
  character(80)  :: FileTmat, FileTmatRoot, FileTmatInfo  
  logical        :: axsym, sphere, chiral
!
  FileTmatRoot = FileTmat(14:LenString(FileTmat))
  FileTmatInfo = '../TMATFILES/Info' // FileTmatRoot(1:LenString(FileTmatRoot))
  open (unit = iTmatInfo, file = FileTmatInfo, status = 'replace')      
  if (.not. sphere) then
    write (iTmatInfo,"(2x,'The T matrix is stored in the file',a55)") FileTmat
  else
    write (iTmatInfo,"(2x,'The T vector is stored in the file',a55)") FileTmat
  end if
  if (axsym) then
    if (sphere) then
      write (iTmatInfo,"(2x,'The scatterer is a spherical particle.')") 
    else
      if (chiral) then
        write (iTmatInfo,"(2x,'The scatterer is an axisymmetric and chiral particle.')") 
      else
        write (iTmatInfo,"(2x,'The scatterer is an axisymmetric particle.')") 
      end if
    end if
  else
    if (chiral) then
      write (iTmatInfo,"(2x,'The scatterer is a nonaxisymmetric and chiral particle.')")
    else
      write (iTmatInfo,"(2x,'The scatterer is a nonaxisymmetric particle.')")
    end if
  end if  
  if (.not. sphere) then            
    write (iTmatInfo,"(2x,'The dimensions of the T matrix are given by:')") 
    write (iTmatInfo,"(2x,'- maximum expansion order,   Nrank = ',i3,',')") Nrank
    write (iTmatInfo,"(2x,'- number of azimuthal modes, Mrank = ',i3,'.')") Mrank
  else
    write (iTmatInfo,"(2x,'The dimension of the T vector is given by:')") 
    write (iTmatInfo,"(2x,'- maximum expansion order,   Nrank = ',i3,'.')") Nrank
  end if
  close (unit = iTmatInfo)
end subroutine write_InfoFileTmat
! **********************************************************************************
subroutine write_HeadFileTmatVct (nt)
  use parameters
  implicit none
  integer   :: nt
!
  write (iTmat, "(2x,'Half - Dimensions of the T Vector:')")
  write (iTmat, "(2x,i10)") nt
  write (iTmat, "(2x,'T Vector:')")
end subroutine write_HeadFileTmatVct
! **********************************************************************************
subroutine write_FileTmatVct (nt, v)
  use parameters
  implicit none
  integer    :: nt, i
  complex(O) :: v(2*nt)
!  
  do i = 1, 2*nt
    write (iTmat,"(2x,1pe24.15,1x,1pe24.15)") v(i)
  end do  
end subroutine write_FileTmatVct
! **********************************************************************************
subroutine write_FileSmat (N, Nteta, S, Cext, Cscat, Qext, Qscat)
  use parameters
  implicit none
  integer    :: N, Nteta, i, j
  real(O)    :: Cext, Cscat, Qext, Qscat
  complex(O) :: S(N,Nteta)
!
  do i = 1, N
    write (iSS,"(10(2x,1pe24.15,1x,1pe24.15),/)") (S(i,j), j = 1, Nteta)
  end do 
  write (iSS,"(4(2x,1pe13.4))") Cext, Cscat, Qext, Qscat 
end subroutine write_FileSmat
! **********************************************************************************
subroutine write_DSCS (Nteta, ExtThetaDom, h, v)
  use parameters
  implicit none      
  integer   :: Nteta, i
  real(O)   :: teta, h(Nteta), v(Nteta)                 
  logical   :: ExtThetaDom
!
  write (iOutput,"(1x, a, //, 2x, a, 9x, a, 8x, a, /)")                             &
 'normalized differential scattering cross section', 'theta', 'parallel',           &
 'perpendicular'        
  do i = 1, Nteta
    if (ExtThetaDom) then         
      teta = real(i - 1,O) * 360._O / real(Nteta - 1,O)                           
    else 
      teta = real(i - 1,O) * 180._O / real(Nteta - 1,O)     
    end if       
    write (iOutput,"(1x,f6.2,5x,1pe13.4,5x,1pe13.4)") teta, h(i), v(i)   
  end do 
  write (iOutput, "(/)")                  
end subroutine write_DSCS
! **********************************************************************************
subroutine  write_4ConvParam (Nint1, Nint2, Nrank, Mrank)
  use parameters
  implicit none
  integer :: Nint1, Nint2, Nrank, Mrank
!
  write (iOutput,"(1x, a, i5, a, 1x, a, i5, a, 1x, a, i3, a, 1x, a, i3, /)")        &
 'Nint1 = ', Nint1, ',', 'Nint2 = ', Nint2, ',', 'Nrank = ', Nrank,                 &
 ',', 'Mrank = ', Mrank
end subroutine write_4ConvParam
! **********************************************************************************
subroutine  write_3ConvParam (Nint, Nrank, m)
  use parameters
  implicit none
  integer :: Nint, Nrank, m
!
  write (iOutput,"(7x,'Nint = ',i5,',',1x,'Nrank = ',i3,',',1x,'Mrank = ',i3,/)")   &
         Nint, Nrank, m
end subroutine write_3ConvParam
! **********************************************************************************
subroutine  write_2ConvParam (Nrank, Mrank)
  use parameters
  implicit none
  integer :: Nrank, Mrank
! 
  write (iOutput, "(11x,'Nrank = ',i3,',',1x,'Mrank = ',i3,/)") Nrank, Mrank
end subroutine write_2ConvParam
! **********************************************************************************
subroutine  write_2ConvParamAxsym (Nint, Nrank)
  use parameters
  implicit none
  integer :: Nint, Nrank
!
  write (iOutput,"(11x,'Nint = ',i5,',',1x,'Nrank = ',i3,/)") Nint, Nrank
end subroutine write_2ConvParamAxsym 
! **********************************************************************************
subroutine  write_1ConvParam (m)
  use parameters
  implicit none
  integer :: m
!
  write (iOutput,"(18x,'m = ',i3,/)") m
end subroutine write_1ConvParam
! **********************************************************************************
subroutine write_Effic (Qscat, Qext)
  use parameters
  implicit none
  real(O) :: Qscat, Qext
!
  write (iOutput,"(6x,'scattering efficiency = ',1pe13.4)")   Qscat
  write (iOutput,"(6x,'extinction efficiency = ',1pe13.4,/)") Qext
end subroutine write_Effic
! **********************************************************************************
subroutine write_NrankConvRes (NthetaConv, Nteta, epsNrank)
  use parameters
  implicit none
  integer :: NthetaConv, Nteta
  real(O) :: epsNrank
!
  write (iOutput,"(1x, a, i2, a, / , 1x, a, 1f5.2, a, /)")                         &
 '--- the solution converges in ', NthetaConv,     ' points ---',                  &
 '--- with an relative error of ', 100 * epsNrank, ' %    ---'
  if (NthetaConv >= int(0.8*Nteta)) then
    print "(/,2x,'Convergence criterion for Nrank is satisfied;')" 
  else
    print "(/,2x,'Convergence criterion for Nrank is not satisfied;')"
  end if 
end subroutine write_NrankConvRes
! **********************************************************************************
subroutine write_NintConvRes (NthetaConv, Nteta, epsNint)
  use parameters
  implicit none
  integer :: NthetaConv, Nteta
  real(O) :: epsNint
!
  write (iOutput,"(1x, a, i2, a, / , 1x, a, 1f5.2, a, /)")                          &
 '--- the solution converges in ', NthetaConv,     ' points ---',                   &
 '--- with an relative error of ', 100 * epsNint,  ' %    ---'
  if (NthetaConv >= int(0.8*Nteta)) then
    print "(/,2x,'Convergence criterion for Nint is satisfied;')" 
  else
    print "(/,2x,'Convergence criterion for Nint is not satisfied;')"
  end if 
end subroutine write_NintConvRes    
! **********************************************************************************
subroutine write_MrankConvRes (NthetaConv, epsMrank)
  use parameters
  implicit none
  integer :: NthetaConv
  real(O) :: epsMrank
!  
  write (iOutput,"(1x, a, i2, a, / , 1x, a, 1f5.2, a, /)")                          &
 '--- the solution converges in ', NthetaConv,     ' points ---',                   &
 '--- with an relative error of ', 100 * epsMrank, ' %    ---'                          
end subroutine write_MrankConvRes
! **********************************************************************************
subroutine write_TypeConvHead (TypeConv)
  use parameters
  implicit none
  integer :: TypeConv
!
  if (TypeConv == 1) then
    write (iOutput,"(5x,'--- Convergence Test over Nint ---',/)")
  else if (TypeConv == 2) then 
    write (iOutput,"(5x,'--- Convergence Test over Nrank ---',/)")
  else if (TypeConv == 3) then
    write (iOutput,"(5x,'--- Convergence Test over Mrank ---',/)")
  else if (TypeConv == 4) then
    write (iOutput,"(2x,'-- Convergence Test over Nrank and Mrank ---',/)")  
  end if
end subroutine write_TypeConvHead
! **********************************************************************************
subroutine  write_3ConvParamRegCOMP (Nint, m, Nrank, Npart, Nrankp, ReduceOrder)
  use parameters
  implicit none
  integer :: Nint, m, Nrank, Npart, Nrankp(Npart), i
  logical :: ReduceOrder
!
  write (iOutput,"(7x,'Nint = ',i5,',',1x,'Nrank = ',i3,',',1x,'Mrank = ',i3,/)")   &
         Nint, Nrank, m
  do i = 1, Npart
    if (.not. ReduceOrder) then
      write (iOutput,"(7x,'Nrank for region ',i2,', Nrank = ',i3,/)")               &
             i, Nrankp(i)
    else
      write (iOutput,"(7x,'Nrank for region ',i2,', Nrank = ',i3,/)")               &
             i, Nrankp(i) - 1
    end if
  end do
end subroutine write_3ConvParamRegCOMP
! **********************************************************************************
subroutine  write_2ConvParamRegCOMP (Nint, Nrank, Npart, Nrankp)
  use parameters
  implicit none
  integer :: Nint, Nrank, Npart, Nrankp(Npart), i
!
  write (iOutput,"(7x,'Nint = ',i5,',',1x,'Nrank = ',i3,/)") Nint, Nrank
  do i = 1, Npart
    write (iOutput,"(7x,'Nrank for region ',i2,', Nrank = ',i3,/)") i, Nrankp(i)
  end do
end subroutine write_2ConvParamRegCOMP
! **********************************************************************************
subroutine  write_3ConvParamReg (Nint, m, Npart, Nrankp, ReduceOrder)
  use parameters
  implicit none
  integer :: Nint, m, Npart, Nrankp(Npart), i
  logical :: ReduceOrder
!
  write (iOutput,"(7x,'Nint = ',i5,',',8x,'Mrank = ',i3,/)") Nint, m
  do i = 1, Npart
    if (.not. ReduceOrder) then
      write (iOutput,"(7x,'Nrank for region ',i2,', Nrank = ',i3,/)")               &
             i, Nrankp(i)
    else
      write (iOutput,"(7x,'Nrank for region ',i2,', Nrank = ',i3,/)")               &
             i, Nrankp(i) - 1
    end if
  end do
end subroutine write_3ConvParamReg
! **********************************************************************************
subroutine  write_2ConvParamReg (Nint, Npart, Nrankp)
  use parameters
  implicit none
  integer :: Nint, Npart, Nrankp(Npart), i
!
  write (iOutput,"(7x,'Nint = ',i5,/)") Nint
  do i = 1, Npart
    write (iOutput,"(7x,'Nrank for region ',i2,', Nrank = ',i3,/)") i, Nrankp(i)
  end do
end subroutine write_2ConvParamReg
! **********************************************************************************                 
subroutine inputDSCS_SCAT (ComputeDSCS, axsym, sphere, chiral, Mrank, Nrank,        &
           phiGS, tetaGI, phiGI, alfamin, alfamax, Nalfa, betamin, betamax, Nbeta,  &
           gamamin, gamamax, Ngama, epol_beta, epol_alpha, alfapGauss, x0, y0, z0,  &
           w0, TypeExcit, wavelength, anorm, normalized, FileTmat)                                                          
  use parameters
  implicit none
  integer       :: Mrank, Nrank, Nalfa, Nbeta, Ngama, LenString, iOut
  real(O)       :: phiGS, tetaGI, phiGI, alfamin, alfamax, betamin, betamax,        &
                   gamamin, gamamax, alfapGauss, x0, y0, z0, w0, wavelength,        &
                   anorm, snorm, k, s, grd 
  complex(O)    :: epol_beta, epol_alpha		   
  character(5)  :: TypeExcit                  
  character(80) :: FileTmat, FileTmatWrite
  logical       :: ComputeDSCS, sphere, axsym, chiral, normalized
!
  if (ComputeDSCS) then
    iOut = iDSCS
  else
    iOut = iSCAT
  end if
  grd = 180._O / Pi
  write (iOut,"(/,2x,'Input Parameters of the Scattering Problem:',/)")
  write (iOut,"(2x, a, 1pe13.4, a)")                                                &
 'wavelength of the ambient medium, wavelength = ', wavelength, ';'
  if (axsym)   write (iOut,"(2x,'axisymmetric particle;')")
  if (sphere)  write (iOut,"(2x,'spherical particle;')")
  if (chiral)  write (iOut,"(2x,'chiral particle;')")
  FileTmatWrite = FileTmat(14:LenString(FileTmat))
  write (iOut,"(2x,'name of the file containing the T matrix, FileTmat = ',a)")     &
         FileTmatWrite
  write (iOut,"(2x,'maximum expansion order, Nrank = ',i3,';')") Nrank
  write (iOut,"(2x,'maximum azimuthal order, Mrank = ',i3,';')") Mrank     
  write (iOut,"(2x,'average over Euler orientation angles:')")
  write (iOut,"(2x,'alphamin = ',f7.2,', alphamax = ',f7.2,', Nalpha = ',i3,',')")  & 
         alfamin * grd, alfamax * grd, Nalfa 
  write (iOut,"(2x,'betamin  = ',f7.2,', betamax  = ',f7.2,', Nbeta  = ',i3,';')")  &
         betamin * grd, betamax * grd, Nbeta
  if (.not. axsym)                                                                  &
  write (iOut,"(2x,'gammamin = ',f7.2,', gammamax = ',f7.2,', Ngamma = ',i3,';')")  &
         gamamin * grd, gamamax * grd, Ngama                                            
  if (TypeExcit(1:1) == 'P') then
    write (iOut,"(2x,'plane wave excitation;')") 
  else if (TypeExcit(1:1) == 'G') then
    write (iOut,"(2x,'Gaussian beam excitation;')")
  end if
  write (iOut,"(2x,'incident direction, thetaGI = ',f7.2,', phiGI = ',f7.2,';')")   &
         tetaGI * grd, phiGI * grd  
  if (ComputeDSCS) then 	 
    if (TypeExcit(1:1) == 'P') then
      write (iOut,"(2x,'beta  polarization vector = (',1pe10.3,',',1pe10.3,')')")   &
             epol_beta
      write (iOut,"(2x,'alpha polarization vector = (',1pe10.3,',',1pe10.3,')')")   &
             epol_alpha 
    else if (TypeExcit(1:1) == 'G') then                                                                        
      write (iOut,"(2x,'polarization angle, alphap = ',f7.2,';')") alfapGauss * grd
    end if              
  end if     
  if (TypeExcit(1:1) == 'G') then         
    write (iOut,"(2x, a, 1pe10.3, a, 1pe10.3, a, 1pe10.3, a)")                      &
   'Cartesian coordinates of the focal point, x0 = ', x0, ', y0 = ', y0,            &
   ', z0 = ', z0, ';'
    write (iOut,"(2x,'waist radius, w0 = ',1pe10.3,';')") w0
    k = 2._O * Pi / wavelength
    s = 1._O / k / w0
    write (iOut,"(2x,'Gaussian beam parameter, s = ',1pe10.3,';')") s  
  end if
  if (ComputeDSCS) write (iOut,"(2x,'scattering plane, phiGS = ',f7.2,';')")        &
                          phiGS * grd  
  write (iOut,"(2x, a, 1pe10.3, a)")                                                &
 'characteristic length of the scatterer, anorm = ', anorm, ';'  
  if (normalized) then
    snorm = Pi * anorm * anorm
    write (iOut,"(2x,'normalization constant, pi * anorm**2 = ',1pe13.4,';')")      &
           snorm
  end if
  write (iOut,*)
end subroutine inputDSCS_SCAT 
! **********************************************************************************
subroutine write_DSCS1 (Nteta, normalized, ExtThetaDom, h, v, Cscat, Cext,           &
           Qscat, Qext)
  use parameters
  implicit none      
  integer   :: Nteta, i
  real(O)   :: teta, h(Nteta), v(Nteta), Cscat, Cext, Qscat, Qext               
  logical   :: normalized, ExtThetaDom
!   
  write (iDSCS,"(/,2x,'Results:',/)") 
  write (iDSCS,"(2x,'Cross Sections and Efficiencies:')")
  write (iDSCS,"(2x, a, 1pe13.4, 5x, a, 1pe13.4)")                                  &
         '<Cscat> = ', Cscat, '<Qscat> = ', Qscat 
  write (iDSCS,"(2x, a, 1pe13.4, 5x, a, 1pe13.4)")                                  &
         '<Cext>  = ', Cext,  '<Qext>  = ', Qext 	              
  write (iDSCS,*)  
  write (iDSCS,"(2x,'Differential Scattering Cross Section:')")     
  if (normalized) then     
    write (iDSCS,"(2x, a, /, 2x, a, 9x, a, 8x, a,/)")                               &
   'normalized average  differential scattering cross section', 'theta',            &
   'parallel', 'perpendicular'                                                            
  else
    write (iDSCS,"(8x, a, /, 2x, a, 9x, a, 8x, a,/)")                               &
   'average differential scattering cross section', 'theta', 'parallel',            &
   'perpendicular'                                                           
  end if
  do i = 1, Nteta
    if (ExtThetaDom) then         
      teta = real(i - 1,O) * 360._O / real(Nteta - 1,O)                           
    else 
      teta = real(i - 1,O) * 180._O / real(Nteta - 1,O)     
    end if                          
    write (iDSCS,"(1x,f6.2,5x,1pe13.4,5x,1pe13.4)") teta, h(i), v(i)     
  end do     
end subroutine write_DSCS1
! **********************************************************************************
subroutine write_SCAT (Nelem, NameElem, Nphi, phi, Ntheta, thetamin, thetamax,      &
           ZTPAv, NphiAL, NthetaAL, KeAv, CscatXAv, CscatYAv, CextXAv, CextYAv,     &
           QscatXAv, QscatYAv, QextXAv, QextYAv, ComputeAsymPar, gX_t, gX_p, gX_r,  &
           gY_t, gY_p, gY_r, TypeExcit)
  use parameters
  implicit none
  integer        :: Nelem, Ntheta(NphiMax), Nphi, NphiAL, NthetaAL, i, j, k,        &
                    itheta, iphi
  real(O)        :: thetamin(NphiMax), thetamax(NphiMax), phi(NphiMax), CscatXAv,   &
                    CscatYAv, CextXAv, CextYAv, QscatXAv, QscatYAv, QextXAv,        &
                    QextYAv, gX_t, gX_p, gX_r, gY_t, gY_p, gY_r, KeAv(4,4),         &
                    ZTPAv(NphiAL,NthetaAL,Nelem), dtheta, phiGS, thetaGS		   
  character(2)   :: NameElem(Nelem)
  character(5)   :: TypeExcit  
  logical        :: ComputeAsymPar    
!
  write (iSCAT,"(/,2x,'Results:',/)")
  write (iSCAT,"(2x,'Scattering Cross Sections and Efficiencies:')")
  write (iSCAT,"(2x, a, 1pe13.4, 5x, a, 1pe13.4)")                                  &
        '<CscatX> = ',CscatXAv,'<QscatX> = ', QscatXAv
  write (iSCAT,"(2x, a, 1pe13.4, 5x, a, 1pe13.4)")                                  &
        '<CscatY> = ',CscatYAv,'<QscatY> = ', QscatYAv  
  write (iSCAT,*)
  write (iSCAT,"(2x,'Extinction Cross Sections and Efficiencies:')")
  write (iSCAT,"(2x, a, 1pe13.4, 5x, a, 1pe13.4)")                                  &
        '<CextX>  = ', CextXAv,'<QextX>  = ', QextXAv
  write (iSCAT,"(2x, a, 1pe13.4, 5x, a, 1pe13.4)")                                  &
        '<CextY>  = ', CextYAv,'<QextY>  = ', QextYAv  	 		
  write (iSCAT,*)
  if (ComputeAsymPar) then	
    write (iSCAT,"(2x,'Mean direction of propagation of the scattered wave:')")
    write (iSCAT,"(2x, a, 1pe13.4, 5x, a, 1pe13.4, 5x, a, 1pe13.4)")                & 
          '<gX_beta> = ', gX_t,'<gX_alpha> = ', gX_p, '<gX_k> = ', gX_r    	
    write (iSCAT,"(2x, a, 1pe13.4, 5x, a, 1pe13.4, 5x, a, 1pe13.4)")                &
          '<gY_beta> = ', gY_t,'<gY_alpha> = ', gY_p, '<gY_k> = ', gY_r  		
    write (iSCAT,*)
  end if
  if (ComputeAsymPar) then    
    write (iSCAT,"(2x, a)")                                                           &
   'Note: <Cscat>, <Cext> and <g> correspond to linearly X- and Y-polarized waves'         	   	   	   	   	   
  else
    write (iSCAT,"(2x, a)")                                                           &
   'Note: <Cscat> and <Cext> correspond to linearly X- and Y-polarized waves'     
  end if     
  write (iSCAT,*)
  if (TypeExcit(1:5) == 'PLANE') then          
    write (iSCAT,"(2x,'Extinction Matrix:')")
    do i = 1, 4
      write (iSCAT,"(4(2x,1pe13.4))") (KeAv(i,j), j = 1, 4)
    end do
    write (iSCAT,*)  
  end if    
  write (iSCAT,"(2x,'Phase Matrix:')")
  write (iSCAT,"(3x,a3,7x,a5,9x,a2,15(13x,a2))") 'phi','theta', NameElem(1),        &
        (NameElem(k), k = 2, Nelem)
  do iphi = 1, Nphi
    phiGS = phi(iphi)
    if (Ntheta(iphi) /= 1) then
      dtheta = (thetamax(iphi) - thetamin(iphi)) / (Ntheta(iphi) - 1)
    else
      dtheta = 0._O
    end if          
    do itheta = 1, Ntheta(iphi)
      thetaGS = thetamin(iphi) + (itheta - 1) * dtheta  
      write (iSCAT,"(1x,f6.2,5x,f6.2,16(5x,1pe10.3))")  phiGS * 180._O / Pi,        &
             thetaGS * 180._O / Pi, (ZTPAv(iphi,itheta,k), k = 1, Nelem) 
    end do
  end do         
end subroutine write_SCAT
! **********************************************************************************
subroutine inputDSCS_SCATAvrg (ComputeDSCS, wavelength, axsym, chiral, FileTmat,    &
           Nrank, Mrank, anorm, normalized, DoNumAvrg, Nalpha, Nbeta, Ngamma,       &
           NthetaGS, epol_beta, epol_alpha)
  use parameters
  implicit none 
  integer       :: Nrank, Mrank, LenString, Nalpha, Nbeta, Ngamma, NthetaGS, iOut
  real(O)       :: wavelength, anorm, snorm                        
  complex(O)    :: epol_beta, epol_alpha
  character(80) :: FileTmat, FileWrite
  logical       :: ComputeDSCS, axsym, chiral, normalized, DoNumAvrg 
!
  if (ComputeDSCS) then
    iOut = iDSCS
  else
    iOut = iSCAT
  end if
  write (iOut,"(/,2x,'Input Parameters of the Scattering Problem:',/)") 
  write (iOut,"(2x,'wavelength of the ambient medium, wavelength = ',1pe13.4,';')") & 
         wavelength
  if (axsym)   write (iOut,"(2x,'axisymmetric particle;')")
  if (chiral)  write (iOut,"(2x,'chiral particle;')")
  FileWrite = FileTmat(14:LenString(FileTmat))
  write (iOut,"(2x,'name of the file containing the T matrix, FileTmat = ',a)")     &
         FileWrite
  write (iOut,"(2x,'maximum expansion order, Nrank = ',i3,';')") Nrank
  write (iOut,"(2x,'maximum azimuthal order, Mrank = ',i3,';')") Mrank   
  write (iOut,"(2x,'plane wave excitation;')")
  write (iOut,"(2x,'incident direction, thetaGI = 0, phiGI = 0;')")
  if (ComputeDSCS) then    
    write (iOut,"(2x,'beta  polarization vector = (',1pe10.3,',',1pe10.3,')')")     &
           epol_beta
    write (iOut,"(2x,'alpha polarization vector = (',1pe10.3,',',1pe10.3,')')")     &
           epol_alpha          
  end if                        
  write (iOut,"(2x,'scattering plane, phiGS = 0.00')")                         
  write (iOut,"(2x,'characteristic length of the particle, anorm = ',1pe10.3,';')") &
         anorm
  if (normalized) then
    snorm = Pi * anorm * anorm
    write (iOut,"(2x,'normalization constant, pi * anorm**2 = ',1pe13.4,';')")      &
           snorm
  end if        
  if (.not. DoNumAvrg) then
    write (iOut,"(2x, a)")                                                          &
   'analytical averaging procedure for computing the quantities <SS*>;'       
  else
    write (iOut,"(2x, a)")                                                          &
   'numerical averaging procedure for computing the quantities <SS*>;'       
    write (iOut,"(2x,'number of integration points for orientational averaging:')")
    if (axsym) then       
      write (iOut,"(2x,'Nalpha = ',i3,', Nbeta = ',i3,';')") Nalpha, Nbeta
    else
      write (iOut,"(2x,'Nalpha = ',i3,', Nbeta = ',i3,', Ngamma = ',i3,';')")       &
             Nalpha, Nbeta, Ngamma
    end if
  end if       
  write (iOut,"(2x, a, i3, a)")                                                     & 
 'the average quantities <SS*> are computed at NthetaGS = ', NthetaGS,              &
 ' scattering angles;'
  write (iOut,*)
end subroutine inputDSCS_SCATAvrg
! **********************************************************************************
subroutine write_ScatMatrix (MirorSym, axsym, FailHovTest, Nfail, FailMessage,      &
           KeAv, teta, Z, i, Nelem, IndI, IndJ, NameElem, Cscat, Cext, Qscat, Qext, &
           CscatV, CextV, QscatV, QextV, AsymPar, AsymParV)
  use parameters
  implicit none
  integer        :: i, Nfail, Nelem, IndI(Nelem), IndJ(Nelem), k, ii, jj, ifail
  real(O)        :: KeAv(4,4), teta, Z(4,4), Cscat, Cext, Qscat, Qext,              &
                    CscatV, CextV, QscatV, QextV, AsymPar, AsymParV, Elem(16)                      
  character(2)   :: NameElem(Nelem)
  logical        :: MirorSym, axsym, FailHovTest
  character(256) :: FailMessage(10)
!
  if (i == 1) then
    write (iSCAT,"(/,2x,'Results:',/)")
    write (iSCAT,"(  2x,'Cross Sections, Efficiencies and Asymmetry Parameter:')") 
    if (.not.MirorSym) then
      write (iSCAT,"(2x, a, 1pe13.4, 5x, a, 1pe13.4)")                              &
            '<Cscat>_I = ', Cscat,  '<Qscat>_I = ', Qscat 
      write (iSCAT,"(2x, a, 1pe13.4, 5x, a, 1pe13.4)")                              &
            '<Cscat>_V = ', CscatV, '<Qscat>_V = ', QscatV 
      write (iSCAT,*)	    
      write (iSCAT,"(2x, a, 1pe13.4, 5x, a, 1pe13.4)")                              &
            '<Cext>_I  = ', Cext,   '<Qext>_I  = ', Qext
      write (iSCAT,"(2x, a, 1pe13.4, 5x, a, 1pe13.4)")                              &
            '<Cext>_V  = ', CextV,  '<Qext>_V  = ', QextV	    
      write (iSCAT,*)
      write (iSCAT,"(2x, a, 1pe13.4)") '< cos >_I = ', AsymPar
      write (iSCAT,"(2x, a, 1pe13.4)") '< cos >_V = ', AsymParV      	    	                          
    else   
      write (iSCAT,"(2x, a, 1pe13.4, 5x, a, 1pe13.4)")                              &
            '<Cscat> = ', Cscat,  '<Qscat> = ', Qscat
      write (iSCAT,"(2x, a, 1pe13.4, 5x, a, 1pe13.4)")                              &
            '<Cext>  = ', Cext,   '<Qext>  = ', Qext 
      write (iSCAT,"(2x, a, 1pe13.4)") '< cos > = ', AsymPar	    	                                        
    end if
    write (iSCAT,*)          
    write (iSCAT,"(2x,'Extinction Matrix:')")    
    do ii = 1, 4
      write (iSCAT,"(4(2x,1pe13.4))") (KeAv(ii,k), k = 1, 4)
    end do                           
    write (iSCAT,*)          
    write (iSCAT,"(2x,'Scattering Matrix:')")  
    write (iSCAT,"(2x,a5,9x,a2,15(13x,a2))") 'theta', NameElem(1),                  &
          (NameElem(k), k = 2, Nelem)
  end if    
  do k = 1, Nelem
    ii = IndI(k)
    jj = IndJ(k)
    Elem(k) = Z(ii,jj)
  end do         
  write (iSCAT,"(1x,f6.2,16(5x,1pe10.3))")  teta * 180._O / Pi,                     &
        (Elem(k), k = 1, Nelem)
  if (axsym .and. FailHovTest) then            
    write (iSCAT,"(2x,'test of Van der Mee and Hovenier is not satisfied;')")
    do ifail = 1, Nfail
      write (iScat,"(2x,a)") FailMessage(ifail)
    end do
  end if
end subroutine write_ScatMatrix
! **********************************************************************************
! *                           SPHERICAL PARTICLES                                  * 
! **********************************************************************************                 
subroutine inputDSCS_SCATAvrgSPH (ComputeDSCSAvrg, amin, amax, TypeDist, Npar, par, &
           Nrank, Nint, alfap, wavelength, anorm, normalized)                            
  use parameters
  implicit none
  integer  :: TypeDist, Npar, Nrank, Nint, iout, i
  real(O)  :: par(Npar), alfap, amin, amax, wavelength, anorm
  logical  :: ComputeDSCSAvrg, normalized
!
  if (ComputeDSCSAvrg) then
    iOut = iDSCS
  else
    iOut = iSCAT
  end if
  write (iOut,"(/,2x,'Input Parameters of the Scattering Problem:',/)")
  write (iOut,"(2x,'wavelength of the free space, wavelength = ',1pe13.4,';')")     &
         wavelength
  if (TypeDist == 1) then
    write (iOut,"(2x,'modified gamma distribution;')")
  else if (TypeDist == 2) then
    write (iOut,"(2x,'log normal distribution;')")
  else if (TypeDist == 3) then
    write (iOut,"(2x,'gamma distribution;')")
  else if (TypeDist == 4) then
    write (iOut,"(2x,'power law distribution;')")
  else if (TypeDist == 5) then
    write (iOut,"(2x,'modified bimodal log normal distribution;')")
  end if
  write (iOut,"(2x,'parameters of the size distribution:')")
  do i = 1, Npar
    write (iOut,"(2x,'par(',i2,') = ',1pe13.4,',')") i, par(i)
  end do
  write (iOut,"(2x, a, 1pe13.4, a, 1pe13.4, a)")                                    &
 'bounds of the particle size distribution, amin = ', amin, ', amax = ', amax, ';'
  write (iOut,"(2x,'number of integration points, Nint = ',i4,';')") Nint
  write (iOut,"(2x,'maximum expansion order, Nrank = ',i3,';')") Nrank
  if (ComputeDSCSAvrg) write (iOut,"(2x, a, f7.2, a)")                              &
 'polarization angle of the incident plane wave, alphap = ', alfap * 180._O / Pi, ';'   
  write (iOut,"(2x, a, 1pe10.3, a)")                                                &
 'characteristic length of the cluster, anorm = ', anorm, ';'    
  if (normalized) write (iOut,"(2x, a, 1pe13.4, a)")                                &
 'normalization constant, pi * anorm**2 = ', Pi * anorm * anorm, ';'
  write (iOut,*)
end subroutine inputDSCS_SCATAvrgSPH 
! **********************************************************************************
subroutine write_ScatMatrixAvrg (teta, Z, i, Cscat, Cext, Qscat, Qext, AsymPar,     &
           AvrgArea, EffRadius, AvrgRadius, AvrgVolume, VWAvrgRadius, EffVar)
  use parameters
  implicit none
  integer        :: i
  real(O)        :: teta, Z(4,4), Cscat, Cext, Qscat, Qext, AsymPar, AvrgArea,      &
                    EffRadius, AvrgRadius, AvrgVolume, VWAvrgRadius, EffVar
!
  if (i == 1) then
    write (iSCAT,"(/,2x,'Results:',/)")
    write (iSCAT,"(2x,'Cross Sections, Efficiencies and Asymmetry Parameter:')") 
    write (iSCAT,"(2x,'average scattering cross section = ',1pe13.4,';')") Cscat
    write (iSCAT,"(2x,'average extinction cross section = ',1pe13.4,';')") Cext                   
    write (iSCAT,*)
    write (iSCAT,"(2x,'average scattering efficiency    = ',1pe13.4,';')") Qscat
    write (iSCAT,"(2x,'average extinction efficiency    = ',1pe13.4,';')") Qext 
    write (iSCAT,*)
    write (iSCAT,"(2x,'average asymmetry parameter      = ',1pe13.4,';')") AsymPar
    write (iSCAT,*) 
    write (iSCAT,"(2x,'Characteristics of the size distribution:')")
    write (iSCAT,"(2x,'average area of the geometric projection  = ',1pe13.4,';')") &
           AvrgArea
    write (iSCAT,"(2x,'effective radius    = ',1pe13.4,';')") EffRadius
    write (iSCAT,"(2x,'effective variance  = ',1pe13.4,';')") EffVar
    write (iSCAT,"(2x,'average radius      = ',1pe13.4,';')") AvrgRadius
    write (iSCAT,"(2x,'average volume      = ',1pe13.4,';')") AvrgVolume
    write (iSCAT,"(2x,'volume - weighted average radius  = ',1pe13.4,';')")         &
           VWAvrgRadius
    write (iSCAT,*)
    write (iSCAT,"(2x,'Scattering Matrix:')")  
    write (iSCAT,"(2x,a5,9x,a2,3(13x,a2))") 'theta', '11', '33', '12', '34'               
  end if        
  write (iSCAT,"(1x,f6.2,16(5x,1pe10.3))") teta * 180._O / Pi,                      &
         Z(1,1), Z(3,3), Z(1,2), Z(3,4)  
end subroutine write_ScatMatrixAvrg
! **********************************************************************************
! *                       PARTICLE ON OR NEAR A PLANE SURFACE                      *
! **********************************************************************************                 
subroutine inputDSCS_EMF (ComputeDSCS, TypeScat, k, ind_ref_s, z0, Mrank, Nrank,    &
           phiGS, beta, alphap, snorm, normalized, FileTmat)                             
  use parameters
  implicit none
  integer       :: TypeScat, Mrank, Nrank, LenString, iOut
  real(O)       :: k, phiGS, beta, alphap, snorm, z0, wavelength, anorm, grd
  complex(O)    :: ind_ref_s                   
  character(80) :: FileTmat, FileTmatWrite
  logical       :: ComputeDSCS, normalized
!
  wavelength = 2._O * Pi / k
  anorm = sqrt(snorm / Pi) / k
  grd   = 180._O / Pi
  if (ComputeDSCS) then
    iOut = iDSCS
  else
    iOut = iSCAT
  end if
  write (iOut,"(/,2x,'Input Parameters of the Scattering Problem:',/)")
  write (iOut,"(2x,'wavelength of the ambient medium, wavelength = ',1pe13.4,';')") & 
         wavelength
  write (iOut,"(2x, a, 1pe10.3, ',', 1pe10.3, a)")                                  &
 'relative refractive index of the substrate, ind_refSUB = (', ind_ref_s, ');' 
  FileTmatWrite = FileTmat(14:LenString(FileTmat))                 
  write (iOut,"(2x,'name of the file containing the T matrix, FileTmat = ',a40)")   &
         FileTmatWrite
  write (iOut,"(2x,'maximum expansion order, Nrank = ',i3,';')") Nrank
  write (iOut,"(2x,'maximum azimuthal order, Mrank = ',i3,';')") Mrank           
  write (iOut,"(2x,'axial position of the particle, z0 = ',1pe10.3,';')") z0  
  if (TypeScat == 1)  then
    write (iOut,"(2x, a)")                                                          &
   'illumination by a plane wave traveling in the ambient medium;'
    write (iOut,"(2x,'incident angle, beta = ',f7.2,';')") beta * grd                          
  else if (TypeScat == 2)  then
    write (iOut,"(2x,'illumination by a plane wave traveling in the substrate;')")
    write (iOut,"(2x,'incident angle, beta = ',f7.2,';')") beta * grd       
    write (iOut,"(2x, a, f7.2, a)")                                                 &
   'propagation angle of the incident wave, Pi - beta = ', 180._O - beta * grd, ';'        
  end if
  write (iOut,"(2x,'polarization angle of the incident wave, alphap = ',f7.2,';')") &
         alphap * grd
  if (ComputeDSCS) write (iOut,"(2x,'scattering plane, phiGS = ', f7.2,';')")       &
                          phiGS * grd  
  write (iOut,"(2x, a, 1pe10.3, a)")                                                & 
 'characteristic length of the scatterer, anorm = ', anorm, ';'   
  if (normalized) write (iOut,"(2x, a, 1pe13.4, a)")                                &
 'normalization constant, pi * anorm**2 = ', Pi * anorm * anorm, ';'
  write (iOut,*)        
end subroutine inputDSCS_EMF
! **********************************************************************************
subroutine write_DSCSPARTSUB (Nteta, h, v)
  use parameters
  implicit none      
  integer   :: Nteta, i
  real(O)   :: teta, h(Nteta), v(Nteta)                 
!
  write (iOutput,"(1x, a, //, 2x, a, 9x, a, 8x, a, /)")                             &
 'normalized differential scattering cross section', 'theta', 'parallel',           &
 'perpendicular'                        
  do i = 1, Nteta
    teta = 90._O + real(i - 1,O) * 180._O / real(Nteta - 1,O)            
    write (iOutput,"(1x,f6.2,5x,1pe13.4,5x,1pe13.4)")  teta, h(i), v(i)  
  end do 
  write (iOutput, "(/)")                  
end subroutine write_DSCSPARTSUB
! **********************************************************************************
subroutine write_DSCSPARTSUB1 (NthetaGS, normalized, h, v)
  use parameters
  implicit none      
  integer   :: NthetaGS, i
  real(O)   :: teta, h(NthetaGS), v(NthetaGS)                   
  logical   :: normalized
!
  write (iDSCS,"(/,2x,'Results:',/)") 
  if (normalized) then
    write (iDSCS,"(2x, a, /, 2x, a, 9x, a, 8x, a, /)")                              &   
   'normalized differential scattering cross section', 'theta', 'parallel',         &
   'perpendicular' 
  else
    write (iDSCS,"(8x, a, /, 2x, a, 9x, a, 8x, a, /)")                              & 
   'differential scattering cross section', 'theta', 'parallel', 'perpendicular'  
  end if
  do i = 1, NthetaGS
    teta = 90._O + real(i - 1,O) * 180._O / real(NthetaGS - 1,O)                                    
    write (iDSCS,  "(1x,f6.2,5x,1pe13.4,5x,1pe13.4)") teta, h(i), v(i)          
  end do
end subroutine write_DSCSPARTSUB1  
! **********************************************************************************
subroutine write_EMF (Nphi, phi, Ntheta, thetamin, thetamax, emf, NphiAL, NthetaAL)
  use parameters
  implicit none
  integer    :: Nphi, Ntheta(NphiMax), NphiAL, NthetaAL, iphi, itheta
  real(O)    :: phi(NphiMax), thetamin(NphiMax), thetamax(NphiMax), phiGS, thetaGS, &
                dtheta, grd
  complex(O) :: emf(3,NphiAL,NthetaAL)
!
  grd = 180._O / Pi
  write (iSCAT,"(/,2x,'Results:',/)") 
  write (iSCAT,"(2x,'Electromagnetic Fields:')")
  write (iSCAT,"(3x, a3, 7x, a5, 12x, a8, 18x, a13)") 'phi','theta', 'parallel',    &
        'perpendicular'          
  do iphi = 1, Nphi
    phiGS = phi(iphi)
    if (Ntheta(iphi) /= 1) then
      dtheta = (thetamax(iphi) - thetamin(iphi)) / (Ntheta(iphi) - 1)
    else
      dtheta = 0._O
    end if          
    do itheta = 1, Ntheta(iphi)
      thetaGS = thetamin(iphi) + (itheta - 1) * dtheta  
      write                                                                                     &
     (iSCAT,"(1x,f6.2,5x,f6.2,5x,'(',1pe10.3,',',1pe10.3,')', 5x,'(',1pe10.3,',',1pe10.3,')')") & 
      phiGS * grd, thetaGS * grd, emf(2,iphi,itheta), emf(3,iphi,itheta)
    end do      
  end do
end subroutine write_EMF         




     
  







