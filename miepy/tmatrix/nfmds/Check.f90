! **********************************************************************************
! *                  ROUTINES FOR CHECKING THE INPUT DATA                          *
! *    ------------------------------------------------------------------------    *
! *    Partial list of subroutines:                                                *
! *      check_MatPropAXSYM,       check_MatPropNONAXSYM,    check_inputAXSYM,     *
! *      check_inputNONAXSYM,      check_inputSCTAVRGSPH,    check_geomAXSYM,      *
! *      check_geomAXSYMOblate,    check_geom3D,             check_geomCOMP,       *      
! *      check_geomLAY,            check_MrankNrank,         check_chirality,      * 
! *      check_ind_ref,            check_ind_ref1,           check_ind_ref2,       *          
! *      check_ind_ref3,           check_Interpolation,      check_MatrixSolver,   *  
! *      check_Integration,        check_anorm,              check_tetaminmax,     *  
! *      check_Nphi,               check_teta_phiminmax,     check_betaminmax,     *  
! *      check_alfaminmax,         check_gamaminmax,         check_TypeExcit,      *        
! *      check_Random,             check_incident_direction, check_incdir_partsub, *    
! *      check_polarization_angle, check_azimuthal_plane,    check_MatrixElem,     *      
! *      check_StoreAvrgMtrSS,     check_NSimpson,                                 *
! *      check_NalphaNbetaNgamma.                                                  *
! *   Subroutines requiring the restarting of the program:                         * 
! *      check_dimensionMat,       check_dimensionVct,       check_circum_radii,   *
! *      check_MaxNrank,           check_MaxNface.                                 *
! **********************************************************************************
recursive subroutine check_MatPropAXSYM (perfectcond, chiral, kb)
  use parameters
  implicit none
  integer  :: ierr
  real(O)  :: kb
  logical  :: perfectcond, chiral, continuare
!  
  if (perfectcond .and. chiral) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the logical variables perfectcond and chiral are both true;')"
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
    print "(2x,'- enter the logical variables perfectcond and chiral;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) perfectcond, chiral
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the logical variables perfectcond and chiral;')"
      end if
    end do
    call check_MatPropAXSYM (perfectcond, chiral, kb)
  end if
  if (chiral .and. kb == 0._O)                                                      &
  print "(/,2x,'Warning: the chirality parameter is zero;')"  
end subroutine check_MatPropAXSYM
!***********************************************************************************
recursive subroutine check_MatPropNONAXSYM (perfectcond, anisotrop, chiral, kb)
  use parameters
  implicit none
  integer       :: ierr
  real(O)       :: kb
  logical       :: perfectcond, chiral, anisotrop, continuare
!  
  if (perfectcond .and. chiral) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the logical variables perfectcond and chiral are both true;')"
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
    print "(2x,'- enter the logical variables perfectcond and chiral;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) perfectcond, chiral
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the logical variables perfectcond and chiral;')"
      end if
    end do
    call check_MatPropNONAXSYM (perfectcond, anisotrop, chiral, kb)
  end if
  if (perfectcond .and. anisotrop) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the logical variables perfectcond and anisotrop are both true;')"
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
      print "(2x,'- enter the logical variables perfectcond and anisotrop;')"
      ierr = 1
      do while (ierr /= 0)
        read (*, *, iostat = ierr) perfectcond, anisotrop
        if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the logical variables perfectcond and anisotrop;')"
      end if
    end do
    call check_MatPropNONAXSYM (perfectcond, anisotrop, chiral, kb)
  end if
  if (chiral .and. anisotrop) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the logical variables chiral and anisotrop are both true;')"
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
      print "(2x,'- enter the logical variables chiral and anisotrop;')"
      ierr = 1
      do while (ierr /= 0)
        read (*, *, iostat = ierr) chiral, anisotrop
        if (ierr /= 0) then 
          print "(/,2x,'Input error during the read statement;')"
          print "(  2x,'- enter the logical variables chiral and anisotrop;')"
      end if
    end do
    call check_MatPropNONAXSYM (perfectcond, anisotrop, chiral, kb)
  end if  
  if (chiral .and. kb == 0._O)                                                      &
  print "(/,2x,'Warning: the chirality parameter is zero;')"  
end subroutine check_MatPropNONAXSYM
! **********************************************************************************
recursive subroutine check_inputAXSYM (miror, chiral, DS)  
  implicit none
  integer  :: ierr
  logical  :: miror, DS, chiral, continuare
!
  if (miror .and. chiral) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'for chiral particles set the logical variable miror to false;')"
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
    print "(2x,'- enter the logical variable miror;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) miror
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the logical variable miror;')"
      end if
    end do
    call check_inputAXSYM (miror, chiral, DS)
  end if
  if (miror .and. DS) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the logical variables miror and DS are both true;')"
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
      print "(2x,'- enter the logical variables miror and DS;')"
      ierr = 1
      do while (ierr /= 0)
        read (*, *, iostat = ierr) miror, DS
        if (ierr /= 0) then 
          print "(/,2x,'Input error during the read statement;')"
          print "(  2x,'- enter the logical variables miror and DS;')"
      end if
    end do
    call check_inputAXSYM (miror, chiral, DS)
  end if    
end subroutine check_inputAXSYM
!***********************************************************************************
recursive subroutine check_inputNONAXSYM (FileGeom, FileFEM, miror, chiral,         &
          anisotrop, Nazimutsym) 
  implicit none
  integer       :: Nazimutsym, ierr
  logical       :: FileGeom, miror, chiral, anisotrop, continuare
  character(80) :: FileFEM
!
  if (FileGeom .and. FileFEM(1:1) == ' ') then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the name of the geometry file is omitted;')"
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
    print "(2x,'- enter the name of the geometry file;')"
    call read_char80 (3, FileFEM)   
    call check_inputNONAXSYM (FileGeom, FileFEM, miror, chiral, anisotrop,          &
         Nazimutsym)
  end if
  if (miror .and. chiral) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'for chiral particles set the logical variable miror to false;')"
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
    print "(2x,'- enter the logical variable miror;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) miror
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the logical variable miror;')"
      end if
    end do
    call check_inputNONAXSYM (FileGeom, FileFEM, miror, chiral, anisotrop,          &
         Nazimutsym)
  end if 
  if (miror .and. anisotrop) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'for anisotropic particles set the logical variable miror to false;')"
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
      print "(2x,'- enter the logical variable miror;')"
      ierr = 1
      do while (ierr /= 0)
        read (*, *, iostat = ierr) miror
        if (ierr /= 0) then 
          print "(/,2x,'Input error during the read statement;')"
          print "(  2x,'- enter the logical variable miror;')"
      end if
    end do
    call check_inputNONAXSYM (FileGeom, FileFEM, miror, chiral, anisotrop,          &
         Nazimutsym)
  end if     
  if (anisotrop .and. Nazimutsym /= 0) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'for anisotropic particles set the integer variable Nazimutsym to 0;')"
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
      print "(2x,'- enter the integer variable Nazimutsym;')"
      ierr = 1
      do while (ierr /= 0)
        read (*, *, iostat = ierr) Nazimutsym
        if (ierr /= 0) then 
          print "(/,2x,'Input error during the read statement;')"
          print "(  2x,'- enter the integer variable Nazimutsym;')"
      end if
    end do
    call check_inputNONAXSYM (FileGeom, FileFEM, miror, chiral, anisotrop,          &
         Nazimutsym)
  end if
end subroutine check_inputNONAXSYM
! **********************************************************************************
recursive subroutine check_inputSCTAVRGSPH (amin, amax, TypeDist, Npar)                          
  use parameters
  implicit none
  integer  :: TypeDist, Npar, ierr
  real(O)  :: amin, amax  
  logical  :: continuare 
!  
  if (TypeDist > 5) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'incorrect value of the variable TypeDist;')"
    print "(  2x,'the permissive values are:')"
    print "(  2x,'- 1 for the modified gamma distribution,')" 
    print "(  2x,'- 2 for the log normal distribution,')"
    print "(  2x,'- 3 for the gamma distribution,')"
    print "(  2x,'- 4 for the power law distribution,')"
    print "(  2x,'- 5 for the modified bimodal log normal distribution;')"
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
    print "(2x,'- enter the variable TypeDist;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) TypeDist
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the variable TypeDist;')"
      end if
    end do
    call check_inputSCTAVRGSPH (amin, amax, TypeDist, Npar)
  end if 
  if (TypeDist == 1 .and. Npar /= 3) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'incorrect number of parameters for TypeDist = 1, set Npar to 3;')"
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    print "(  2x,'  Note: the program can be continued if the input values of the')"
    print "(  2x,'  three parameters of the distribution are correct;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
    print "(2x,'- enter the number of parameters Npar;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) Npar
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the number of parameters Npar;')"
      end if
    end do
    call check_inputSCTAVRGSPH (amin, amax, TypeDist, Npar)
  else if (TypeDist == 2 .and. Npar /= 2) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'incorrect number of parameters for TypeDist = 2, set Npar to 2;')"
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    print "(  2x,'  Note: the program can be continued if the input values of the')"
    print "(  2x,'  two parameters of the distribution are correct;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
    print "(2x,'- enter the number of parameters Npar;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) Npar
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the number of parameters Npar;')"
      end if
    end do
    call check_inputSCTAVRGSPH (amin, amax, TypeDist, Npar) 
  else if (TypeDist == 3 .and. Npar /= 2) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'incorrect number of parameters for TypeDist = 3, set Npar to 2;')"
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    print "(  2x,'  Note: the program can be continued if the input values of the')"
    print "(  2x,'  two parameters of the distribution are correct;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
    print "(2x,'- enter the number of parameters Npar;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) Npar
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the number of parameters Npar;')"
      end if
    end do
    call check_inputSCTAVRGSPH (amin, amax, TypeDist, Npar)
  else if (TypeDist == 4 .and. Npar /= 1) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'incorrect number of parameters for TypeDist = 4, set Npar to 1;')"
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    print "(  2x,'  Note: the program can be continued if the input value of the')"
    print "(  2x,'  parameter of the distribution is correct;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
    print "(2x,'- enter the number of parameters Npar;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) Npar
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the number of parameters Npar;')"
      end if
    end do
    call check_inputSCTAVRGSPH (amin, amax, TypeDist, Npar)
  else if (TypeDist == 5 .and. Npar /= 5) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'incorrect number of parameters for TypeDist = 5, set Npar to 5;')"
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    print "(  2x,'  Note: the program can be continued if the input values of the')"
    print "(  2x,'  five parameters of the distribution are correct;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
    print "(2x,'- enter the number of parameters Npar;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) Npar
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the number of parameters Npar;')"
      end if
    end do
    call check_inputSCTAVRGSPH (amin, amax, TypeDist, Npar)
  end if
  if (amax < amin) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the lower bound of the size distribution exceeds the upper bound;')"
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
      print "(2x,'- enter the bounds amin and amax;')"
      ierr = 1
      do while (ierr /= 0)
        read (*, *, iostat = ierr) amin, amax
        if (ierr /= 0) then 
          print "(/,2x,'Input error during the read statement;')"
          print "(  2x,'- enter the bounds amin and amax;')"
      end if
    end do
    call check_inputSCTAVRGSPH (amin, amax, TypeDist, Npar) 
  end if              
end subroutine check_inputSCTAVRGSPH
! **********************************************************************************
recursive subroutine check_geomAXSYM (TypeGeom, Nsurf, Nparam)
  implicit none
  integer  :: TypeGeom, Nsurf, Nparam, ierr
  logical  :: continuare   
!
  if (TypeGeom > 3) then
    print "(/,2x,'Warning: the geometry is not contained in the library;')"     
  end if 
  select case (TypeGeom)
  case (1)
    if (Nsurf  /= 2)  then
      print "(/,2x,'Error in the input file: set Nsurf to 2;')"
      print "(  2x, a)",                                                            &
     '- enter true to continue the program or false to stop the program; '
      print "(  2x, a)",                                                            &
     '  Note: the program can be continued if the input values of surf(1)'
      print "(  2x,'  and surf(2) are correct;')"
      call read_logical (continuare)
      if (.not. continuare) stop 
      print "(2x,'- enter the number of surface parameters Nsurf;')"
      ierr = 1
      do while (ierr /= 0)
        read (*, *, iostat = ierr) Nsurf
        if (ierr /= 0) then 
          print "(/,2x,'Input error during the read statement;')"
          print "(  2x,'- enter the number of surface parameters Nsurf;')"
        end if
      end do
      call check_geomAXSYM (TypeGeom, Nsurf, Nparam)
    end if
    if (Nparam /= 1)  then
      print "(/,2x,'Error in the input file: set Nparam to 1;')"
      print "(  2x, a)",                                                            &
     '- enter true to continue the program or false to stop the program;'
      call read_logical (continuare)
      if (.not. continuare) stop 
      print "(2x,'- enter the number of parameters Nparam;')"
      ierr = 1
      do while (ierr /= 0)
        read (*, *, iostat = ierr) Nparam
        if (ierr /= 0) then 
          print "(/,2x,'Input error during the read statement;')"
          print "(  2x,'- enter the number of parameters Nparam;')"
        end if
      end do
      call check_geomAXSYM (TypeGeom, Nsurf, Nparam)
    end if
  case (2,3)
    if (Nsurf  /= 2)  then
      print "(/,2x,'Error in the input file: set Nsurf to 2;')"
      print "(  2x, a)",                                                            &
     '- enter true to continue the program or false to stop the program; '
      print "(  2x, a)",                                                            &
     '  Note: the program can be continued if the input values of surf(1)'
      print "(  2x,'  and surf(2) are correct;')"
      call read_logical (continuare)
      if (.not. continuare) stop 
      print "(2x,'- enter the number of surface parameters Nsurf;')"
      ierr = 1
      do while (ierr /= 0)
        read (*, *, iostat = ierr) Nsurf
        if (ierr /= 0) then 
          print "(/,2x,'Input error during the read statement;')"
          print "(  2x,'- enter the number of surface parameters Nsurf;')"
        end if
      end do
      call check_geomAXSYM (TypeGeom, Nsurf, Nparam)
    end if
    if (Nparam /= 3)  then
      print "(/,2x,'Error in the input file: set Nparam to 3;')" 
      print "(  2x, a)",                                                            &
     '- enter true to continue the program or false to stop the program;'
      call read_logical (continuare)
      if (.not. continuare) stop 
      print "(2x,'- enter the number of parameters Nparam;')"
      ierr = 1
      do while (ierr /= 0)
        read (*, *, iostat = ierr) Nparam
        if (ierr /= 0) then 
          print "(/,2x,'Input error during the read statement;')"
          print "(  2x,'- enter the number of parameters Nparam;')"
        end if
      end do
      call check_geomAXSYM (TypeGeom, Nsurf, Nparam)
    end if
  end select
end subroutine check_geomAXSYM
! **********************************************************************************
recursive subroutine check_geomAXSYMOblate (TypeGeom, Nsurf, surf)
  use parameters
  implicit none
  integer  :: TypeGeom, Nsurf, ierr
  real(O)  :: surf(Nsurf)
  logical  :: continuare   
!
  if (TypeGeom == 3 .and. surf(1) >= surf(2)) then
    print "(/,2x,'Error in the input file:')"   
    print "(  2x,'the relation surf(1) < surf(2) is not satisfied;')"
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
    print "(2x,'- enter the surface parameters surf(1) and surf(2);')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) surf(1), surf(2)
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the surface parameters surf(1) and surf(2);')"
      end if
    end do
    call check_geomAXSYMOblate (TypeGeom, Nsurf, surf)
  end if 
end subroutine check_geomAXSYMOblate
! **********************************************************************************
recursive subroutine check_geom3D (TypeGeom, Nsurf, Nparam, anisotropic, miror,     &
                     Nazimutsym)
  implicit none
  integer  :: TypeGeom, Nsurf, Nparam, Nazimutsym, ierr
  logical  :: anisotropic, miror, continuare 
!
  if (TypeGeom > 3) then
    print "(/,2x,'Warning: the geometry is not contained in the library;')"     
  end if 
  select case (TypeGeom)
  case (1)
    if (Nsurf  /= 3) then
      print "(/,2x,'Error in the input file: set Nsurf to 3;')"
      print "(  2x, a)",                                                            &
     '- enter true to continue the program or false to stop the program;'
      print "(  2x, a)",                                                            &
     '  Note: the program can be continued if the input values of surf(1),'
      print "(  2x,'  surf(2) and surf(3) are correct;')"
      call read_logical (continuare)
      if (.not. continuare) stop 
      print "(2x,'- enter the number of surface parameters Nsurf;')"
      ierr = 1
      do while (ierr /= 0)
        read (*, *, iostat = ierr) Nsurf
        if (ierr /= 0) then 
          print "(/,2x,'Input error during the read statement;')"
          print "(  2x,'- enter the number of surface parameters Nsurf;')"
        end if
      end do
      call check_geom3D (TypeGeom, Nsurf, Nparam, anisotropic, miror, Nazimutsym)
    end if
    if (Nparam /= 1)  then
      print "(/,2x,'Error in the input file: set Nparam to 1;')"
      print "(  2x, a)",                                                            &
     '- enter true to continue the program or false to stop the program;'
      call read_logical (continuare)
      if (.not. continuare) stop 
      print "(2x,'- enter the number of parameters Nparam;')"
      ierr = 1
      do while (ierr /= 0)
        read (*, *, iostat = ierr) Nparam
        if (ierr /= 0) then 
          print "(/,2x,'Input error during the read statement;')"
          print "(  2x,'- enter the number of parameters Nparam;')"
        end if
      end do
      call check_geom3D (TypeGeom, Nsurf, Nparam, anisotropic, miror, Nazimutsym)
    end if
    if (miror .and. anisotropic)  then
      print "(/,2x,'Error in the input file: set miror to false;')"
      print "(  2x, a)",                                                            &
     '- enter true to continue the program or false to stop the program;'
      call read_logical (continuare)
      if (.not. continuare) stop 
      print "(2x,'- enter the logical variable miror;')"
      ierr = 1
      do while (ierr /= 0)
        read (*, *, iostat = ierr) miror
        if (ierr /= 0) then 
          print "(/,2x,'Input error during the read statement;')"
          print "(  2x,'- enter the logical variable miror;')"
        end if
      end do
      call check_geom3D (TypeGeom, Nsurf, Nparam, anisotropic, miror, Nazimutsym)
    end if
    if (Nazimutsym /= 0)  then
      print "(/,2x,'Error in the input file: set Nazimutsym to 0;')"
      print "(  2x, a)",                                                            &
     '- enter true to continue the program or false to stop the program;'
      call read_logical (continuare)
      if (.not. continuare) stop 
      print "(2x,'- enter the integer variable Nazimutsym;')"
      ierr = 1
      do while (ierr /= 0)
        read (*, *, iostat = ierr) Nazimutsym
        if (ierr /= 0) then 
          print "(/,2x,'Input error during the read statement;')"
          print "(  2x,'- enter the integer variable Nazimutsym;')"
        end if
      end do
      call check_geom3D (TypeGeom, Nsurf, Nparam, anisotropic, miror, Nazimutsym)
    end if
  case (2)
    if (Nsurf  /= 2)  then
      print "(/,2x,'Error in the input file: set Nsurf to 2;')"
      print "(  2x, a)",                                                            &
     '- enter true to continue the program or false to stop the program;'
      print "(  2x, a)",                                                            &
     '  Note: the program can be continued if the input values of surf(1)'
      print "(  2x,'  and surf(2) are correct;')"
      call read_logical (continuare)
      if (.not. continuare) stop 
        print "(2x,'- enter the number of surface parameters Nsurf;')"
        ierr = 1
        do while (ierr /= 0)
          read (*, *, iostat = ierr) Nsurf
          if (ierr /= 0) then 
            print "(/,2x,'Input error during the read statement;')"
            print "(  2x,'- enter the number of surface parameters Nsurf;')"
        end if
      end do
      call check_geom3D (TypeGeom, Nsurf, Nparam, anisotropic, miror, Nazimutsym)
    end if
    if (Nparam /= 6)  then
      print "(/,2x,'Error in the input file: set Nparam to 6;')"
      print "(  2x, a)",                                                            &
     '- enter true to continue the program or false to stop the program;'
      call read_logical (continuare)
      if (.not. continuare) stop 
      print "(2x,'- enter the number of parameters Nparam;')"
      ierr = 1
      do while (ierr /= 0)
        read (*, *, iostat = ierr) Nparam
        if (ierr /= 0) then 
          print "(/,2x,'Input error during the read statement;')"
          print "(  2x,'- enter the number of parameters Nparam;')"
        end if
      end do
      call check_geom3D (TypeGeom, Nsurf, Nparam, anisotropic, miror, Nazimutsym)
    end if
    if (miror)  then
      print "(/,2x,'Error in the input file: set miror to false;')"
      print "(  2x, a)",                                                            &
     '- enter true to continue the program or false to stop the program;'
      call read_logical (continuare)
      if (.not. continuare) stop 
      print "(2x,'- enter the logical variable miror;')"
      ierr = 1
      do while (ierr /= 0)
        read (*, *, iostat = ierr) miror
        if (ierr /= 0) then 
          print "(/,2x,'Input error during the read statement;')"
          print "(  2x,'- enter the logical variable miror;')"
        end if
      end do
      call check_geom3D (TypeGeom, Nsurf, Nparam, anisotropic, miror, Nazimutsym)
    end if
      if (Nazimutsym /= 0)  then
        print "(/,2x,'Error in the input file: set Nazimutsym to 0;')"
        print "(  2x, a)",                                                          &
     '- enter true to continue the program or false to stop the program;'
        call read_logical (continuare)
        if (.not. continuare) stop 
        print "(2x,'- enter the integer variable Nazimutsym;')"
        ierr = 1
        do while (ierr /= 0)
          read (*, *, iostat = ierr) Nazimutsym
          if (ierr /= 0) then 
            print "(/,2x,'Input error during the read statement;')"
            print "(  2x,'- enter the integer variable Nazimutsym;')"
        end if
      end do
      call check_geom3D (TypeGeom, Nsurf, Nparam, anisotropic, miror, Nazimutsym)
    end if
  case (3)
    if (.not. anisotropic) then
      if (Nsurf  /= 2)  then
        print "(/,2x,'Error in the input file: set Nsurf to 2;')"
        print "(  2x, a)",                                                          &
       '- enter true to continue the program or false to stop the program;'
        print "(  2x, a)",                                                          &
       '  Note: the program can be continued if the input values of surf(1)'
        print "(  2x,'  and surf(2) are correct;')"
        call read_logical (continuare)
        if (.not. continuare) stop 
        print "(2x,'- enter the number of surface parameters Nsurf;')"
        ierr = 1
        do while (ierr /= 0)
          read (*, *, iostat = ierr) Nsurf
          if (ierr /= 0) then 
            print "(/,2x,'Input error during the read statement;')"
            print "(  2x,'- enter the number of surface parameters Nsurf;')"
          end if
        end do
        call check_geom3D (TypeGeom, Nsurf, Nparam, anisotropic, miror, Nazimutsym)
      end if 
      if (miror .and. Nparam /= 2)  then
        print "(/,2x,'Error in the input file: set Nparam to 2;')"
        print "(  2x, a)",                                                          &
       '- enter true to continue the program or false to stop the program;'
        call read_logical (continuare)
        if (.not. continuare) stop 
        print "(2x,'- enter the number of parameters Nparam;')"
        ierr = 1
        do while (ierr /= 0)
          read (*, *, iostat = ierr) Nparam
          if (ierr /= 0) then 
            print "(/,2x,'Input error during the read statement;')"
            print "(  2x,'- enter the number of parameters Nparam;')"
          end if
        end do
        call check_geom3D (TypeGeom, Nsurf, Nparam, anisotropic, miror, Nazimutsym)
      end if
      if (.not. miror .and. Nparam /= 3)  then
        print "(/,2x,'Error in the input file: set Nparam to 3;')"
        print "(  2x, a)",                                                          &
       '- enter true to continue the program or false to stop the program;'
        call read_logical (continuare)
        if (.not. continuare) stop 
        print "(2x,'- enter the number of parameters Nparam;')"
        ierr = 1
        do while (ierr /= 0)
          read (*, *, iostat = ierr) Nparam
          if (ierr /= 0) then 
            print "(/,2x,'Input error during the read statement;')"
            print "(  2x,'- enter the number of parameters Nparam;')"
          end if
        end do
        call check_geom3D (TypeGeom, Nsurf, Nparam, anisotropic, miror, Nazimutsym)
      end if 
      if (Nazimutsym < 2)  then
        print "(/,2x,'Error in the input file: set Nazimutsym > 1;')"
        print "(  2x, a)",                                                          &
       '- enter true to continue the program or false to stop the program;'
        call read_logical (continuare)
        if (.not. continuare) stop 
        print "(2x,'- enter the integer variable Nazimutsym;')"
        ierr = 1
        do while (ierr /= 0)
          read (*, *, iostat = ierr) Nazimutsym
          if (ierr /= 0) then 
            print "(/,2x,'Input error during the read statement;')"
            print "(  2x,'- enter the integer variable Nazimutsym;')"
          end if
        end do
        call check_geom3D (TypeGeom, Nsurf, Nparam, anisotropic, miror, Nazimutsym)
      end if 
    else
      print "(/,2x,'Error in the input file:')"
      print "(  2x,'the integer variable TypeGeom is 3 and for uniaxial particles ')"
      print "(  2x,'no mirror and azimuthal symmetries are considered')"
      stop      
    end if  
  end select
end subroutine check_geom3D
! **********************************************************************************
recursive subroutine check_geomCOMP (TypeGeom, Npart, Nsurf, Nparam)
  implicit none
  integer  :: TypeGeom, Npart, Nsurf(Npart), Nparam(Npart), ipart, ierr
  logical  :: continuare
!
  if (TypeGeom > 2) then
    print "(/,2x,'Warning: the geometry is not contained in the library;')"     
  end if 
  select case (TypeGeom)
  case (1)
    do ipart = 1, Npart
      if (Nsurf(ipart)  /= 3)  then
        print "(/,2x,'Error in the input file: set Nsurf for region',i3,' to 3;')", &
                ipart
        print "(  2x, a)",                                                          &
       '- enter true to continue the program or false to stop the program;'
        print "(  2x, a)",                                                          &
       '  Note: the program can be continued if the input values of the'
        print "(  2x,'  surface parameters for region',i3,' are correct;')", ipart
        call read_logical (continuare)
        if (.not. continuare) stop                           
        print "(2x,'- enter the number of surface parameters Nsurf for this region;')"
        ierr = 1
        do while (ierr /= 0)
          read (*, *, iostat = ierr) Nsurf(ipart)
          if (ierr /= 0) then 
            print "(/,2x,'Input error during the read statement;')"
            print "(  2x, a)",                                                      &
           '- enter the number of surface parameters Nsurf for this region;'
          end if
        end do
        call check_geomCOMP (TypeGeom, Npart, Nsurf, Nparam)
      end if
      if (Nparam(ipart) /= 2)  then
        print "(/,2x,'Error in the input file: set Nparam for region',i3,' to 2;')",&
                ipart
        print "(  2x, a)",                                                          &
       '- enter true to continue the program or false to stop the program;'
        call read_logical (continuare)
        if (.not. continuare) stop 
        print "(2x,'- enter the number of parameters Nparam for this region;')"
        ierr = 1
        do while (ierr /= 0)
          read (*, *, iostat = ierr) Nparam(ipart)
          if (ierr /= 0) then 
            print "(/,2x,'Input error during the read statement;')"
            print "(  2x,'- enter the number of parameters Nparam for this region;')"
          end if
        end do
        call check_geomCOMP (TypeGeom, Npart, Nsurf, Nparam)
      end if
    end do
  case (2)
    do ipart = 1, Npart
      if (Nsurf(ipart)  /= 2)  then
        print "(/,2x,'Error in the input file: set Nsurf for region',i3,' to 2;')", &
                ipart
        print "(  2x, a)",                                                          &
       '- enter true to continue the program or false to stop the program;' 
        print "(  2x, a)",                                                          &
       '  Note: the program can be continued if the input values of the' 
        print "(  2x,'  surface parameters for region',i3,' are correct;')", ipart
        call read_logical (continuare)
        if (.not. continuare) stop 
        print "(2x,'- enter the number of surface parameters Nsurf for this region;')"
        ierr = 1
        do while (ierr /= 0)
          read (*, *, iostat = ierr) Nsurf(ipart)
          if (ierr /= 0) then 
            print "(/,2x,'Input error during the read statement;')"
            print "(  2x, a)",                                                      &
           '- enter the number of surface parameters Nsurf for this region;'
          end if
        end do
        call check_geomCOMP (TypeGeom, Npart, Nsurf, Nparam)
      end if
      if (Nparam(ipart) /= 3)  then
        print "(/,2x,'Error in the input file: set Nparam for region',i3,' to 3;')",&
                ipart
        print "(  2x, a)",                                                          &
       '- enter true to continue the program or false to stop the program;'
        call read_logical (continuare)
        if (.not. continuare) stop 
        print "(2x,'- enter the number of parameters Nparam for this region;')"
        ierr = 1
        do while (ierr /= 0)
          read (*, *, iostat = ierr) Nparam(ipart)
          if (ierr /= 0) then 
            print "(/,2x,'Input error during the read statement;')"
            print "(  2x,'- enter the number of parameters Nparam for this region;')"
          end if
        end do
        call check_geomCOMP (TypeGeom, Npart, Nsurf, Nparam)
      end if
    end do    
  end select
end subroutine check_geomCOMP
! **********************************************************************************
recursive subroutine check_geomLAY (TypeGeom, Npart, Nsurf, Nparam)
  implicit none
  integer  :: TypeGeom, Npart, Nsurf(Npart), Nparam(Npart), ipart, ierr
  logical  :: continuare
!
  if (TypeGeom > 2) then
    print "(/,2x,'Warning: the geometry is not contained in the library;')"     
  end if 
  select case (TypeGeom)
  case (1)
    do ipart = 1, Npart
      if (Nsurf(ipart)  /= 2)  then
        print "(/,2x,'Error in the input file: set Nsurf for region',i3,' to 2;')", &
                ipart
        print "(  2x, a)",                                                          &
       '- enter true to continue the program or false to stop the program;'
        print "(  2x, a)",                                                          &
       '  Note: the program can be continued if the input values of the'
        print "(  2x,'  surface parameters for region',i3,' are correct;')", ipart
        call read_logical (continuare)
        if (.not. continuare) stop 
        print "(2x,'- enter the number of surface parameters Nsurf for this region;')"
        ierr = 1
        do while (ierr /= 0)
          read (*, *, iostat = ierr) Nsurf(ipart)
          if (ierr /= 0) then 
            print "(/,2x,'Input error during the read statement;')"
            print "(  2x, a)",                                                      &
           '- enter the number of surface parameters Nsurf for this region;'  
          end if
        end do
        call check_geomLAY (TypeGeom, Npart, Nsurf, Nparam)
      end if
      if (Nparam(ipart) /= 1)  then
        print "(/,2x,'Error in the input file: set Nparam for region',i3,' to 1;')",&
                ipart
        print "(  2x, a)",                                                          &
       '- enter true to continue the program or false to stop the program;'
        call read_logical (continuare)
        if (.not. continuare) stop 
        print "(2x,'- enter the number of parameters Nparam for this region;')"
        ierr = 1
        do while (ierr /= 0)
          read (*, *, iostat = ierr) Nparam(ipart)
          if (ierr /= 0) then 
            print "(/,2x,'Input error during the read statement;')"
            print "(  2x,'- enter the number of parameters Nparam for this region;')"
          end if
        end do
        call check_geomLAY (TypeGeom, Npart, Nsurf, Nparam)
      end if
    end do
  case (2)
    do ipart = 1, Npart
      if (Nsurf(ipart)  /= 2)  then
        print "(/,2x,'Error in the input file: set Nsurf for region',i3,' to 2;')", &
                ipart
        print "(  2x, a)",                                                          &
       '- enter true to continue the program or false to stop the program;'
        print "(  2x, a)",                                                          &
       '  Note: the program can be continued if the input values of the'
        print "(  2x,'  surface parameters for region',i3,' are correct;')", ipart
        call read_logical (continuare)
        if (.not. continuare) stop 
        print "(2x,'- enter the number of surface parameters Nsurf for this region;')"
        ierr = 1
        do while (ierr /= 0)
          read (*, *, iostat = ierr) Nsurf(ipart)
          if (ierr /= 0) then 
            print "(/,2x,'Input error during the read statement;')"
            print "(  2x, a)",                                                      &
           '- enter the number of surface parameters Nsurf for this region;'
          end if
        end do
        call check_geomLAY (TypeGeom, Npart, Nsurf, Nparam)
      end if
      if (Nparam(ipart) /= 3)  then
        print "(/,2x,'Error in the input file: set Nparam for region',i3,' to 3;')",&
                ipart
        print "(  2x, a)",                                                          &
       '- enter true to continue the program or false to stop the program;'
        call read_logical (continuare)
        if (.not. continuare) stop 
        print "(2x,'- enter the number of parameters Nparam for this region;')"
        ierr = 1
        do while (ierr /= 0)
          read (*, *, iostat = ierr) Nparam(ipart)
          if (ierr /= 0) then 
            print "(/,2x,'Input error during the read statement;')"
            print "(  2x,'- enter the number of parameters Nparam for this region;')"
          end if
        end do
        call check_geomLAY (TypeGeom, Npart, Nsurf, Nparam)
      end if
    end do    
  end select
end subroutine check_geomLAY
! **********************************************************************************
recursive subroutine check_MrankNrank (Mrank, Nrank)
  implicit none
  integer :: Mrank, Nrank, ierr
  logical :: continuare   
!  
  if (Mrank > Nrank) then
    print "(/,2x,'Input error: Mrank exceeds Nrank;')"
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
    print "(2x,'- enter the estimated values of Nrank and Mrank;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) Nrank, Mrank
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the estimated values of Nrank and Mrank;')"
      end if
    end do
    call check_MrankNrank (Mrank, Nrank)
  end if 
end subroutine check_MrankNrank
! **********************************************************************************
recursive subroutine check_chirality (kb)
  use parameters
  use derived_parameters
  implicit none
  real(O)  :: kb
  integer  :: ierr
  logical  :: continuare
!  
  if (abs(kb - 1._O) <= MachEps) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the chirality parameter kb is one;')"
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
    print "(2x,'- enter the chirality parameter kb;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) kb
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the chirality parameter kb;')"
      end if
    end do
    call check_chirality (kb)
  end if 
end subroutine check_chirality
! **********************************************************************************
recursive subroutine check_ind_ref (ind_ref)
  use parameters
  use derived_parameters
  implicit none
  complex(O)  :: ind_ref
  integer     :: ierr
  logical     :: continuare
!  
  if (aimag(ind_ref) < 0._O .or. abs(ind_ref) < MachEps) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the refractive index is zero or the imaginary part is negative')"               
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    call read_logical (continuare)
    if (.not. continuare) stop                              
    print "(2x,'- enter the refractive index;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) ind_ref
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the refractive index;')"
      end if
    end do
    call check_ind_ref (ind_ref)
  end if 
end subroutine check_ind_ref
! **********************************************************************************
recursive subroutine check_ind_ref1 (ipart, ind_ref)
  use parameters
  use derived_parameters
  implicit none
  integer     :: ipart, ierr
  complex(O)  :: ind_ref
  logical     :: continuare
!    
  if (aimag(ind_ref) < 0._O .or. abs(ind_ref) < MachEps) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the refractive index is zero or the imaginary part')"
    print "(  2x,'is negative for the region/particle ',i2,';')", ipart
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
    print "(2x,'- enter the refractive index for this region;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) ind_ref
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the refractive index for this region;')"
      end if
    end do
    call check_ind_ref1 (ipart, ind_ref)
  end if 
end subroutine check_ind_ref1
! **********************************************************************************
recursive subroutine  check_ind_ref2 (ind_ref1, ind_ref2)  
  use parameters
  use derived_parameters
  implicit none  
  complex(O)  :: ind_ref1, ind_ref2
  integer     :: ierr
  logical     :: continuare
!    
  if (aimag(ind_ref1) < 0._O .or. abs(ind_ref1) < MachEps ) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the refractive index of the particle is zero or')"
    print "(  2x,'the imaginary part is negative;')"
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
    print "(2x,'- enter the refractive index of the particle;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) ind_ref1
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the refractive index of the particle;')"
      end if
    end do      
    call check_ind_ref2 (ind_ref1, ind_ref2)
  end if 
  if (aimag(ind_ref2) < 0._O .or. abs(ind_ref2) < MachEps) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the refractive index of the substrate is zero or')"
    print "(  2x,'the imaginary part is negative;')"
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
    print "(2x,'- enter the refractive index of the substrate;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) ind_ref2
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the refractive index of the substrate;')"
      end if
    end do
    call check_ind_ref2 (ind_ref1, ind_ref2)        
  end if   
end subroutine check_ind_ref2
! **********************************************************************************
recursive subroutine  check_ind_ref3 (ind_ref1, ind_ref2)  
  use parameters
  use derived_parameters
  implicit none  
  complex(O)  :: ind_ref1, ind_ref2
  integer     :: ierr
  logical     :: continuare
!    
  if (aimag(ind_ref1) < 0._O .or. abs(ind_ref1) < MachEps) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the first refractive index of the uniaxial particle')"
    print "(  2x,'is zero or the imaginary part is negative;')"
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
    print "(2x,'- enter the first refractive index of the uniaxial particle;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) ind_ref1
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the first refractive index of the uniaxial particle;')"
      end if
    end do
    call check_ind_ref3 (ind_ref1, ind_ref2)        
  end if 
  if (aimag(ind_ref2) < 0._O .or. abs(ind_ref2) < MachEps) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the second refractive index of the uniaxial particle')"
    print "(  2x,'is zero or the imaginary part is negative;')"
    print "(  2x,'- enter true to continue the program or false to stop the program;')"
    call read_logical (continuare)
    if (.not. continuare) stop 
    print "(2x,'- enter the second refractive index of the uniaxial particle;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) ind_ref2
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the second refractive index of the uniaxial particle;')"
      end if
    end do
    call check_ind_ref3 (ind_ref1, ind_ref2)        
  end if 
end subroutine check_ind_ref3
! **********************************************************************************
recursive subroutine check_Interpolation (TypeInterp)
  implicit none
  integer       :: ierr
  character(20) :: TypeInterp
!
  if ( TypeInterp(1:6) /= 'LINEAR' .and.                                            &
       TypeInterp(1:6) /= 'SPLINE' .and.                                            &
       TypeInterp(1:7) /= 'HERMITE' ) then
    print "(/,2x,'Error in the general input file:')"
    print "(  2x,'incorrect value of the variable TypeInterp;')"
    print "(  2x,'the permissive values are: LINEAR, SPLINE and HERMITE;')"
    print "(  2x,'- enter the character type variable TypeInterp;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) TypeInterp
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the character type variable TypeInterp;')"
      end if
    end do
    call check_Interpolation (TypeInterp)
  end if  
end subroutine check_Interpolation
! **********************************************************************************
recursive subroutine check_MatrixSolver (TypeMatrSolv)
  implicit none
  integer       :: ierr
  character(20) :: TypeMatrSolv
!
  if (TypeMatrSolv(1:3) /= 'LU1'   .and.                                            &
      TypeMatrSolv(1:3) /= 'LU2'   .and.                                            &
      TypeMatrSolv(1:3) /= 'LU3'   .and.                                            &
      TypeMatrSolv(1:4) /= 'BICG'        )  then
    print "(/,2x,'Error in the general input file:')"
    print "(  2x,'incorrect value of the variable TypeMatrSolv;')"
    print "(  2x,'the permissive values are: LU1, LU2, LU3 and BICG;')"
    print "(  2x,'- enter the character type variable TypeMatrSolv;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) TypeMatrSolv
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the character type variable TypeMatrSolv;')"
      end if
    end do
    call check_MatrixSolver (TypeMatrSolv)
  end if  
end subroutine check_MatrixSolver
! **********************************************************************************
recursive subroutine check_Integration (TypeIntegr) 
  implicit none
  integer       :: ierr 
  character(20) :: TypeIntegr
!
  if (TypeIntegr(1:4) /= 'MET1' .and. TypeIntegr(1:4) /= 'MET2') then
    print "(/,2x,'Error in the general input file:')"
    print "(  2x,'incorrect value of the variable TypeIntegr;')"
    print "(  2x,'the permissive values are:')"
    print "(  2x,'MET1 - routines from Numerical Recipes,')"
    print "(  2x,'MET2 - modified routines from Slatec library;')"
    print "(  2x,'- enter the character type variable TypeIntegr;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) TypeIntegr
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the character type variable TypeIntegr;')"
      end if
    end do
    call check_Integration (TypeIntegr) 
  end if  
end subroutine check_Integration
! **********************************************************************************
recursive subroutine check_Random (TypeRND) 
  implicit none
  integer       :: ierr 
  character(20) :: TypeRND
!
  if (TypeRND(1:4) /= 'SLAT' .and. TypeRND(1:4) /= 'LPCK'.and.                      &
      TypeRND(1:4) /= 'ZIGG') then
    print "(/,2x,'Error in the general input file:')"
    print "(  2x,'incorrect value of the variable TypeRND;')"
    print "(  2x,'the permissive values are:')"
    print "(  2x,'SLAT - modified routine from Slatec library,')"
    print "(  2x,'LPCK - modified routine from Lapack library;')"
    print "(  2x,'ZIGG - Ziggurat method of Marsaglia and Tsang;')"
    print "(  2x,'- enter the character type variable TypeRND;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) TypeRND
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the character type variable TypeRND;')"
      end if
    end do
    call check_Random (TypeRND) 
  end if  
end subroutine check_Random
! **********************************************************************************
recursive subroutine check_anorm (anorm)
  use parameters
  implicit none
  real(O)  :: anorm
  integer  :: ierr
!  
  if (anorm <= 0._O) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the characteristic length anorm is negative or zero;')"
    print "(  2x,'- enter the characteristic length anorm;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) anorm
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the characteristic length anorm;')"
      end if
    end do
    call check_anorm (anorm)
  end if 
end subroutine check_anorm 
! **********************************************************************************
recursive subroutine check_tetaminmax (thetamin, thetamax, Ntheta)
  use parameters
  use derived_parameters
  implicit none
  integer  :: Ntheta, ierr
  real(O)  :: thetamin, thetamax, min, max
!  
  min = 0._O
  max = 180._O    
  if (thetamin < min .or. thetamin > max .or. thetamax < min .or. thetamax > max .or.&
      thetamin > thetamax) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the relation 0 <= thetamin <= thetamax <= 180 deg is not satisfied;')"
    print "(  2x,'- enter the variables thetamin and thetamax;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) thetamin, thetamax
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the variables thetamin and thetamax;')"
      end if
    end do
    call check_tetaminmax (thetamin, thetamax, Ntheta)                  
  end if   
  if (abs(thetamax - thetamin) <= MachEps .and. Ntheta /= 1) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the number of sample points Ntheta should be 1;')"      
    print "(  2x,'- enter the number of sample points Ntheta;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) Ntheta
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the number of sample points Ntheta;')"
      end if
    end do
    call check_tetaminmax (thetamin, thetamax, Ntheta)
  end if 
end subroutine check_tetaminmax
! **********************************************************************************
recursive subroutine check_Nphi (Nphi)
  use parameters
  implicit none
  integer  :: Nphi, ierr
!      
  if (Nphi > NphiMax) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x, a)",                                                              &
   'the number of scattering planes Nphi exceeds the maximum dimension NphiMax;'                
    print "(  2x,'- enter the number of scattering planes Nphi < ',i3,';')", NphiMax
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) Nphi
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the number of scattering planes Nphi < ',i3,';')",    &
                NphiMax
      end if
    end do
    call check_Nphi (Nphi)
  end if 
end subroutine check_Nphi
! **********************************************************************************
recursive subroutine check_teta_phiminmax (complete, Nphi, phi, Ntheta, thetamin,   &
                     thetamax)
  use parameters
  use derived_parameters
  implicit none
  integer  :: Nphi, Ntheta(NphiMax), iphi, ierr
  real(O)  :: phi(NphiMax), thetamin(NphiMax), thetamax(NphiMax), min, max, grd
  logical  :: complete 
!
  grd = 180._O / Pi
  do iphi = 1, Nphi
    min = 0._O
    max = 360._O    
    if (phi(iphi) < min .or. phi(iphi) > max) then
      print "(/,2x,'Error in the input file:')"
      print "(  2x,'the relation 0 <= phi <= 360 deg is not satisfied')"
      print "(  2x,'for the scattering plane iphi = ',i3,';')", iphi
      print "(  2x, a, i3, a)",                                                     &
     '- enter the azimuthal angle phi of the scattering plane ', iphi, ';'
      ierr = 1
      do while (ierr /= 0)
        read (*, *, iostat = ierr) phi(iphi)
        if (ierr /= 0) then 
          print "(/,2x,'Input error during the read statement;')"
          print "(  2x, a, i3, a)",                                                 &
         '- enter the azimuthal angle phi of the scattering plane ', iphi, ';'
        end if
      end do
      call check_teta_phiminmax (complete, Nphi, phi, Ntheta, thetamin, thetamax)                       
    end if   
    if (complete) then
      min = 0._O
    else
      min = 90._O
      end if
    max = 180._O    
    if (thetamin(iphi) < min .or. thetamin(iphi) > max .or.                         &
        thetamax(iphi) < min .or. thetamax(iphi) > max .or.                         &
        thetamin(iphi) > thetamax(iphi)) then
      print "(/,2x,'Error in the input file:')"
      if (complete) then
        print "(2x, a)",                                                            &
       'the relation 0 <= thetamin <= thetamax <= 180 deg is not satisfied'
      else
        print "(2x, a)",                                                            &
       'the relation 90 <= thetamin <= thetamax <= 180 deg is not satisfied' 
      end if
      print "(2x,'for the scattering plane phi = ',f7.2,';')", phi(iphi) * grd
      print "(2x,'- enter the variables thetamin and thetamax')"
      print "(2x,'for the scattering plane phi = ',f7.2,';')", phi(iphi) * grd        
      ierr = 1
      do while (ierr /= 0)
        read (*, *, iostat = ierr) thetamin(iphi), thetamax(iphi)
        if (ierr /= 0) then 
          print "(/,2x,'Input error during the read statement;')"
          print "(  2x,'- enter the variables thetamin and thetamax')"
          print "(  2x,'for the scattering plane phi = ',f7.2,';')", phi(iphi) * grd              
        end if
      end do
      call check_teta_phiminmax (complete, Nphi, phi, Ntheta, thetamin, thetamax)                               
    end if   
    if (abs(thetamax(iphi) - thetamin(iphi)) <= MachEps .and. Ntheta(iphi) /= 1) then
      print "(/,2x,'Error in the input file:')"
      print "(  2x,'the number of sample points Ntheta for the scattering plane')"
      print "(  2x,'phi = ',f7.2,' should be 1;')", phi(iphi) * grd                                                     
      print "(  2x,'- enter the number of sample points Ntheta')"
      print "(  2x,'for the scattering plane phi = ',f7.2,';')", phi(iphi) * grd
      ierr = 1
      do while (ierr /= 0)
        read (*, *, iostat = ierr) Ntheta(iphi)
        if (ierr /= 0) then 
          print "(/,2x,'Input error during the read statement;')"             
          print "(  2x,'- enter the number of sample points Ntheta')"
          print "(  2x,'for the scattering plane phi = ',f7.2,';')", phi(iphi) * grd
        end if
      end do
      call check_teta_phiminmax (complete, Nphi, phi, Ntheta, thetamin, thetamax)
    end if 
  end do          
end subroutine check_teta_phiminmax
! **********************************************************************************
recursive subroutine check_betaminmax (betamin, betamax, Nbeta)
  use parameters
  use derived_parameters
  implicit none
  integer  :: Nbeta, ierr
  real(O)  :: betamin, betamax, min, max, deltacos
!  
  min = 0._O
  max = 180._O    
  if (betamin < min .or. betamin > max .or. betamax < min .or. betamax > max .or.   &
      betamin > betamax) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the relation 0 <= betamin <= betamax <= 180 deg is not satisfied;')"
    print "(  2x,'- enter the variables betamin and betamax;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) betamin, betamax
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the variables betamin and betamax;')"
      end if
    end do
    call check_betaminmax (betamin, betamax, Nbeta)                     
  end if   
  deltacos = cos(betamin) - cos(betamax)
  if (abs(deltacos) <= MachEps .and. Nbeta /= 1) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the number of quadrature points Nbeta should be 1;')"   
    print "(  2x,'- enter the number of quadrature points Nbeta;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) Nbeta
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the number of quadrature points Nbeta;')"
      end if
    end do
    call check_betaminmax (betamin, betamax, Nbeta)
  end if 
end subroutine check_betaminmax
! **********************************************************************************
recursive subroutine check_alfaminmax (alphamin, alphamax, Nalpha)
  use parameters
  use derived_parameters
  implicit none
  integer  :: Nalpha, ierr
  real(O)  :: alphamin, alphamax, min, max
!
  min = 0._O
  max = 360._O    
  if (alphamin < min .or. alphamin > max .or. alphamax < min .or. alphamax > max .or.&
      alphamin > alphamax) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the relation 0 <= alphamin <= alphamax <= 360 deg is not satisfied;')"
    print "(  2x,'- enter the variables alphamin and alphamax;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) alphamin, alphamax
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the variables alphamin and alphamax;')"
      end if
    end do
    call check_alfaminmax (alphamin, alphamax, Nalpha)                  
  end if   
  if (abs(alphamax - alphamin) <= MachEps .and. Nalpha /= 1) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the number of quadrature points Nalpha should be 1;')"  
    print "(  2x,'- enter the number of quadrature points Nalpha;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) Nalpha
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the number of quadrature points Nalpha;')"
      end if
    end do
    call check_alfaminmax (alphamin, alphamax, Nalpha)
  end if 
end subroutine check_alfaminmax
! **********************************************************************************
recursive subroutine check_gamaminmax (gammamin, gammamax, Ngamma)
  use parameters
  use derived_parameters
  implicit none
  integer  :: Ngamma, ierr
  real(O)  :: gammamin, gammamax, min, max
!
  min = 0._O
  max = 360._O    
  if (gammamin < min .or. gammamin > max .or. gammamax < min .or. gammamax > max .or.&
      gammamin > gammamax) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the relation 0 <= gammamin <= gammamax <= 360 deg is not satisfied;')"
    print "(  2x,'- enter the variables gammamin and gammamax;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) gammamin, gammamax
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the variables gammamin and gammamax;')"
      end if
    end do
    call check_gamaminmax (gammamin, gammamax, Ngamma)                  
  end if  
  if (abs(gammamax - gammamin) <= MachEps .and. Ngamma /= 1) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the number of quadrature points Ngamma should be 1;')"  
    print "(  2x,'- enter the number of quadrature points Ngamma;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) Ngamma
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the number of quadrature points Ngamma;')"
      end if
    end do
    call check_gamaminmax (gammamin, gammamax, Ngamma)
  end if  
end subroutine check_gamaminmax
! **********************************************************************************
recursive subroutine check_TypeExcit (TypeExcit)
  implicit none  
  character(5) :: TypeExcit 
  integer      :: ierr
!
  if ( TypeExcit(1:5) /= 'PLANE' .and. TypeExcit(1:5) /= 'GAUSS' ) then     
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'incorrect value of the variable TypeExcit;')"
    print "(  2x,'the permissive values are: PLANE and GAUSS;')"
    print "(  2x,'- enter the character type variable TypeExcit;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) TypeExcit
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the character type variable TypeExcit;')"
      end if
    end do
    call check_TypeExcit (TypeExcit)
  end if  
end subroutine check_TypeExcit
! **********************************************************************************
recursive subroutine check_incident_direction (thetaGI, phiGI)
  use parameters
  implicit none
  integer  :: ierr
  real(O)  :: thetaGI, phiGI, thetamin, thetamax, phimin, phimax
!  
  thetamin = 0._O
  thetamax = 180._O
  phimin   = 0._O
  phimax   = 360._O
  if (thetaGI < thetamin .or. thetaGI > thetamax) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the incident polar angle thetaGI varies between 0 and 180 deg;')"
    print "(  2x,'- enter the incident polar angle thetaGI;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) thetaGI
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the incident polar angle thetaGI;')"
      end if
    end do
    call check_incident_direction (thetaGI, phiGI)
  end if  
  if (phiGI < phimin .or. phiGI > phimax) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the incident azimuth angle phiGI varies between 0 and 360 deg;')"
    print "(  2x,'- enter the incident azimuth angle phiGI;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) phiGI
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the incident azimuth angle phiGI;')"
      end if
    end do
    call check_incident_direction (thetaGI, phiGI)
  end if  
end subroutine check_incident_direction
! **********************************************************************************
recursive subroutine check_incdir_partsub (beta)
  use parameters
  implicit none
  integer  :: ierr
  real(O)  :: beta, betamin, betamax
!
  betamin =  0._O
  betamax = 90._O
  if (beta < betamin .or. beta > betamax) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the incident angle beta varies between 0 and 90 deg;')"
    print "(  2x,'- enter the incident angle beta;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) beta
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the incident angle beta;')"
      end if
    end do
    call check_incdir_partsub (beta)    
  end if 
end subroutine check_incdir_partsub
! **********************************************************************************
recursive subroutine check_polarization_angle (alphap)
  use parameters
  implicit none
  integer  :: ierr
  real(O)  :: alphap, phimin, phimax
!        
  phimin   = 0._O
  phimax   = 360._O   
  if (alphap < phimin .or. alphap > phimax) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the polarization angle alphap varies between 0 and 360 deg;')"
    print "(  2x,'- enter the polarization angle alphap;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) alphap
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the polarization angle alphap;')"
      end if
    end do
    call check_polarization_angle (alphap)          
  end if   
end subroutine check_polarization_angle
! **********************************************************************************
recursive subroutine check_azimuthal_plane (phiGS)
  use parameters
  implicit none
  integer  :: ierr
  real(O)  :: phiGS, phimin, phimax
!
  phimin   = 0._O
  phimax   = 360._O   
  if (phiGS < phimin .or. phiGS > phimax) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the azimuth angle phiGS varies between 0 and 360 deg;')"
    print "(  2x,'- enter the azimuth angle phiGS;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) phiGS
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the azimuth angle phiGS;')"
      end if
    end do
    call check_azimuthal_plane (phiGS)              
  end if   
end subroutine check_azimuthal_plane
! **********************************************************************************
recursive subroutine check_MatrixElem (Nelem, MatrixElem) 
  implicit none
  integer  :: Nelem, MatrixElem(16), k, ij, i, j, ierr
  logical  :: fail
!  
  if (Nelem < 1 .or. Nelem > 16) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'incorrect value of the number of matrix elements Nelem;')"
    print "(  2x,'the relation 1 <= Nelem <= 16 is not satisfied;')"
    print "(  2x,'- enter the number of matrix elements Nelem;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) Nelem
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the number of matrix elements Nelem;')"
      end if
    end do
    call check_MatrixElem (Nelem, MatrixElem)
  end if
  fail = .false.
  do k = 1, Nelem
    ij = MatrixElem(k)
    if (ij <= 10) then
      fail = .true.
    else
      j = mod(ij,10)
      i = (ij - j) / 10
      if (i < 1 .or. i > 4 .or. j < 1 .or. j > 4) fail = .true.
    end if
  end do
  if (fail) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'incorrect values of the matrix elements in the array MatrixElem;')"
    print "(  2x,'the permissive values are: 11, 12, 13, 14, 21, 22, 23, 24,')"
    print "(  2x,'31, 32, 33, 34, 41, 42, 43, 44')"
    do k = 1, Nelem
      print "(2x,'- enter the matrix element ',i2,';')", k
      ierr = 1
      do while (ierr /= 0)
        read (*, *, iostat = ierr) MatrixElem(k)
        if (ierr /= 0) then 
          print "(/,2x,'Input error during the read statement;')"
          print "(  2x,'- enter the matrix element ',i2,';')", k
        end if
      end do
    end do 
    call check_MatrixElem (Nelem, MatrixElem)
  end if 
end subroutine check_MatrixElem 
! **********************************************************************************
recursive subroutine check_StoreAvrgMtrSS (ComputeAvrgMtrSS, StoreAvrgMtrSS)  
  implicit none
  integer  :: ierr
  logical  :: ComputeAvrgMtrSS, StoreAvrgMtrSS
!  
  if (StoreAvrgMtrSS .and. .not. ComputeAvrgMtrSS) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'ComputeAvrgMtrSS = f and StoreAvrgMtrSS = t;')"
    print "(  2x,'- enter the logical variables ComputeAvrgMtrSS and StoreAvrgMtrSS;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) ComputeAvrgMtrSS, StoreAvrgMtrSS
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x, a)",                                                            &
       '- enter the logical variables ComputeAvrgMtrSS and StoreAvrgMtrSS;' 
      end if
    end do
    call check_StoreAvrgMtrSS (ComputeAvrgMtrSS, StoreAvrgMtrSS)
  end if
end subroutine check_StoreAvrgMtrSS
! **********************************************************************************
recursive subroutine check_NSimpson (N) 
  implicit none
  integer :: N, ierr
!      
  if (mod(N,2) == 0) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x, a, i5, a)",                                                       & 
   'the number of integration points for the Simpson rule is even: N = ', N, ';'
    print "(  2x,'- enter the number of integration points;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) N
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the number of integration points;')"
      end if
    end do
    call check_NSimpson (N)                     
  end if   
end subroutine check_NSimpson
! **********************************************************************************
recursive subroutine check_NthetaGS (N) 
  implicit none
  integer :: N, ierr
!      
  if (mod(N,2) == 0) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x, a, i5, a)",                                                       & 
   'the number of zenith angles for <SS*> calculation is even: N = ', N, ';'
    print "(  2x,'- enter the number of zenith angles;')"
    ierr = 1
    do while (ierr /= 0)
      read (*, *, iostat = ierr) N
      if (ierr /= 0) then 
        print "(/,2x,'Input error during the read statement;')"
        print "(  2x,'- enter the number of zenith angles;')"
      end if
    end do
    call check_NthetaGS (N)                     
  end if   
end subroutine check_NthetaGS
! **********************************************************************************
recursive subroutine check_NalphaNbetaNgamma (index, N) 
  implicit none
  integer :: index, N, ierr
!      
  if (N <= 1) then
    print "(/,2x,'Error in the input file:')"
    if (index == 1) then
      print "(2x,'the number of division points Nalpha is 1;')"
      print "(2x,'- enter the number of division points Nalpha;')"
      ierr = 1
      do while (ierr /= 0)
        read (*, *, iostat = ierr) N
        if (ierr /= 0) then 
          print "(/,2x,'Input error during the read statement;')"
          print "(  2x,'- enter the number of division points Nalpha;')"
        end if
      end do
      call check_NalphaNbetaNgamma (index, N)                   
    else if (index == 2) then
      print "(2x,'the number of division points Nbeta is 1;')"
      print "(2x,'- enter the number of division points Nbeta;')"
      ierr = 1
      do while (ierr /= 0)
        read (*, *, iostat = ierr) N
        if (ierr /= 0) then 
          print "(/,2x,'Input error during the read statement;')"
          print "(  2x,'- enter the number of division points Nbeta;')"
        end if
      end do
      call check_NalphaNbetaNgamma (index, N)
    else if (index == 3) then
      print "(2x,'the number of division points Ngamma is 1;')"
      print "(2x,'- enter the number of division points Ngamma;')"
      ierr = 1
      do while (ierr /= 0)
        read (*, *, iostat = ierr) N
        if (ierr /= 0) then 
          print "(/,2x,'Input error during the read statement;')"
          print "(  2x,'- enter the number of division points Ngamma;')"
        end if
      end do
      call check_NalphaNbetaNgamma (index, N)
    end if
  end if   
end subroutine check_NalphaNbetaNgamma
! **********************************************************************************
! *                 SUBROUTINES REQUIRING RESTARTING OF THE PROGRAM                * 
! **********************************************************************************
subroutine check_dimensionMat (ntg, mtg, Nmax)
  implicit none
  integer  :: ntg, mtg, Nmax
!
  if (ntg < Nmax .or. mtg < Nmax) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the physical dimensions of the T matrix are smaller than the ')"
    print "(  2x,'real  dimensions: Nrank  or/and  Mrank  should be increased;')"
    stop
  end if
end subroutine check_dimensionMat
! **********************************************************************************
subroutine check_dimensionVct (ntg, Nmax)
  implicit none
  integer  :: ntg, Nmax
!
  if (ntg < Nmax) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the physical  dimension of the T vector is smaller')"
    print "(  2x,'than the real dimension: Nrank should be increased;')"
    stop
  end if
end subroutine check_dimensionVct
! **********************************************************************************
subroutine check_circum_radii (Ncs, Rcs)
  use parameters
  implicit none
  integer  :: Ncs, i
  real(O)  :: Rcs(Ncs)
!
  do i = 1, Ncs - 1
    if (Rcs(i) >= Rcs(i+1)) then
      print "(/,2x,'Error in the input file:')"
      print "(  2x,'the radii of the circumscribing spheres are not in increasing order;')"    
      stop
    end if
  end do
end subroutine check_circum_radii 
! **********************************************************************************
subroutine check_MaxNrank (Nrank)
  use parameters
  implicit none
  integer  :: Nrank
!
  if (NrankPD < Nrank) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the number of distributed sources Nrank exceeds')"
    print "(  2x,'the maximum number of distributed sources NrankPD;')"
    stop
  end if
end subroutine check_MaxNrank
! **********************************************************************************
subroutine check_MaxNface (Nface)
  use parameters
  implicit none
  integer  :: Nface
!
  if (NfacePD < Nface) then
    print "(/,2x,'Error in the input file:')"
    print "(  2x,'the number of surface elements Nface exceeds')"
    print "(  2x,'the maximum number of surface elements NfacePD;')"
    stop
  end if
end subroutine check_MaxNface
! **********************************************************************************

