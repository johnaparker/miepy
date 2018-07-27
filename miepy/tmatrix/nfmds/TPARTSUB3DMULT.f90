subroutine TPARTSUB3DMULT 
  use parameters
  use allocation, only: Mrankp, Nrankp, xp, yp, zp, FileTmatp, axsymp
  implicit none                                                                              
  integer       :: TypeScat, Nrank, Mrank, Npart, NXint                                   
  real(O)       :: wavelength, z0, anorm, snorm, epsNrank, epsMrank, k                                       
  complex(O)    :: ind_refSUB
  logical       :: DoConvTest, PrnProgress  
! -----------------------------------------------------------------------------------
!                                 Read the input file                               ! 
! -----------------------------------------------------------------------------------           	   
  call readinputPARTSUB3DMULT ( TypeScat, wavelength, ind_refSUB, Nrank, Mrank,     &
       Npart, z0, anorm, DoConvTest, NXint, epsNrank, epsMrank, PrnProgress,        &
       k, snorm)  	    
! -----------------------------------------------------------------------------------
!                                      Main                                         !
! -----------------------------------------------------------------------------------   
  open (unit = iOutput, file = FileOutput, status = "replace")
  call printinputPARTSUB3DMULT (TypeScat, wavelength, anorm, z0, epsNrank,          &
       epsMrank, ind_refSUB, Npart, xp, yp, zp, FileTmatp, axsymp, Mrankp, Nrankp,  &
	   Mrank, Nrank)
!		
  call convergence_Nrank_MrankPARTSUBMULT (TypeScat, k, ind_refSUB, z0, snorm,     &
       Npart, FileTmatp, axsymp, xp, yp, zp, Mrankp, Nrankp, Mrank, Nrank,          &
       NXint, epsNrank, epsMrank, PrnProgress) 		
  close (unit = iOutput)           
  deallocate (FileTmatp, axsymp, Mrankp, Nrankp, xp, yp, zp)  
end subroutine TPARTSUB3DMULT
!***********************************************************************************
subroutine readinputPARTSUB3DMULT ( TypeScat, wavelength, ind_refSUB, Nrank, Mrank, &
           Npart, z0, anorm, DoConvTest, NXint, epsNrank, epsMrank, PrnProgress,    &
		   k, snorm)     
  use parameters
  use derived_parameters
  use allocation, only: Mrankp, Nrankp, xp, yp, zp, FileTmatp, axsymp
  implicit none 
  integer       :: TypeScat, Mrank, Nrank, Npart, NXint, i, ios, ipart, NrankPart,  &
                   MrankPart                                    
  real(O)       :: wavelength, z0, anorm, snorm, epsNrank, epsMrank, k, xPart,      &
                   yPart, zPart                                       
  complex(O)    :: ind_refSUB
  character(80) :: string, FileTmatPart
  logical       :: axsymPart, DoConvTest, PrnProgress, XFindPar     
! -----------------------------------------------------------------------------------
!                         Read the input file FileInputPARTSUB3DMULT                ! 
! -----------------------------------------------------------------------------------    
  call DrvParameters 
  open (unit = iInputPARTSUB, file = FileInputPARTSUB3DMULT, status = "old",        &
        position = "rewind")   
  wavelength  = 0.1_O * 2._O * Pi  
  ind_refSUB  = (1.3_O,0._O)
  string      = 'OptProp'
  if (XFindPar (iInputPARTSUB, string)) then
    read (iInputPARTSUB, *, iostat = ios) wavelength
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable wavelength;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) ind_refSUB
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refSUB;')"
      stop
    end if      
  else
    print "(/,2x,'Group name OptProp not found;')"
    stop  
  end if  
  k = 2._O * Pi / wavelength 
!        
  TypeScat = 1 
  Npart  = 2     
  Nrank  = 6
  Mrank  = 4
  z0     = 0.1_O
  anorm  = 0.1_O   
  string = 'GenProp'
  if (XFindPar (iInputPARTSUB, string)) then
    read (iInputPARTSUB, *, iostat = ios) TypeScat
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable TypeScat;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) Npart
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Npart;')"
      stop                  
    end if   
    read (iInputPARTSUB, *, iostat = ios) Nrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nrank;')"
      stop                  
    end if    
    read (iInputPARTSUB, *, iostat = ios) Mrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Mrank;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) z0
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable z0;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) anorm
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable anorm;')"
      stop
    end if
  else
    print "(/,2x,'Group name GeomProp not found;')"
    stop  
  end if 
  call check_anorm (anorm)  
  xPart = k * anorm
  snorm = Pi * xPart * xPart
!
  allocate (FileTmatp(Npart), axsymp(Npart), Mrankp(Npart), Nrankp(Npart),          &
            xp(Npart), yp(Npart), zp(Npart) )   
  do ipart = 1, Npart
    FileTmatPart = '../TMATFILES/T.dat'
    axsymPart  = .true.     
    NrankPart  = 6
    MrankPart  = 3
    string     = 'TmatPart'
    if (XFindPar (iInputPARTSUB, string)) then
      read (iInputPARTSUB, *, iostat = ios) FileTmatPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable FileTmatPart')"
        print "(  2x,'for the particle ',i3,';')", ipart
        stop
      end if
      read (iInputPARTSUB, *, iostat = ios) axsymPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable axsymPart;')"
        print "(  2x,'for the particle ',i3,';')", ipart
        stop
      end if
      read (iInputPARTSUB, *, iostat = ios) NrankPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable NrankPart;')"
        print "(  2x,'for the particle ',i3,';')", ipart 
        stop
      end if
      read (iInputPARTSUB, *, iostat = ios) MrankPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable MrankPart;')"
        print "(  2x,'for the particle ',i3,';')", ipart 
        stop
      end if
    else
      print "(/,2x,'Group name TmatPart not found;')"
      stop  
    end if      
    call check_MrankNrank (MrankPart, NrankPart)    
    FileTmatp(ipart) = FileTmatPart
    axsymp(ipart)  = axsymPart     
    Nrankp(ipart)  = NrankPart  
    Mrankp(ipart)  = MrankPart      
!
    xPart = 0.1_O
    yPart = 0.1_O
    zPart = 0.1_O     
    string    = 'GeomPartProp'
    if (XFindPar (iInputPARTSUB, string)) then
      read (iInputPARTSUB, *, iostat = ios) xPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable xPart')"
        print "(  2x,'for the particle ',i3,';')", ipart
        stop
      end if
      read (iInputPARTSUB, *, iostat = ios) yPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable yPart;')"
        print "(  2x,'for the particle ',i3,';')", ipart
        stop
      end if
      read (iInputPARTSUB, *, iostat = ios) zPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable zPart;')"
        print "(  2x,'for the particle ',i3,';')", ipart
        stop
      end if
    else
      print "(/,2x,'Group name GeomPartProp not found;')"
      stop  
    end if            
    xp(ipart) = xPart
    yp(ipart) = yPart
    zp(ipart) = zPart    
  end do
!
  DoConvTest = .true.
  string     = 'ConvTest'
  if (XFindPar (iInputPARTSUB, string)) then
    read (iInputPARTSUB, *, iostat = ios) DoConvTest
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable DoConvTest;')"
      stop
    end if        
  else
    print "(/,2x,'Group name ConvTest not found;')"
    stop  
  end if 
!
  if (DoConvTest) then
    print "(/,2x,'Convergence Test for a Cluster of Particles on a Plane Surface')"
    print "(  2x,'--------------------------------------------------------------')"
  end if   
!   
  NXint  = 60
  string = 'NintX'
  if (XFindPar (iInputPARTSUB, string)) then
    read (iInputPARTSUB, *, iostat = ios) NXint
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable NXint;')"
      stop
    end if        
  else
    print "(/,2x,'Group name NintX not found;')"
    stop  
  end if 
!
  epsNrank = 5.e-2_O
  epsMrank = 5.e-2_O 
  string   = 'Errors'
  if (XFindPar (iInputPARTSUB, string)) then   
    read (iInputPARTSUB, *, iostat = ios) epsNrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsNrank;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) epsMrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsMrank;')"
      stop
    end if    
  else
    print "(/,2x,'Group name Errors not found;')"
    stop  
  end if  
!
  PrnProgress = .true.
  string      = 'PrintProgress'
  if (XFindPar (iInputPARTSUB, string)) then
    read (iInputPARTSUB, *, iostat = ios) PrnProgress
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable PrnProgress;')"
      stop
    end if        
  else
    print "(/,2x,'Group name PrintProgress not found;')"
    stop  
  end if     
  close (unit = iInputPARTSUB)  
end subroutine readinputPARTSUB3DMULT
!***********************************************************************************
subroutine printinputPARTSUB3DMULT (TypeScat, wavelength, anorm, z0, epsNrank,      &
           epsMrank, ind_refSUB, Npart, xp, yp, zp, FileTmatp, axsymp, Mrankp,      &
           Nrankp, Mrank, Nrank)             
  use parameters
  implicit none
  integer    :: TypeScat, Npart, Mrankp(Npart), Nrankp(Npart), Mrank, Nrank,        &
                LenString, i
  real(O)    :: wavelength, anorm, epsNrank, epsMrank, z0, xp(Npart), yp(Npart),    &
                zp(Npart)
  complex(O) :: ind_refSUB
  logical       :: axsymp(Npart)
  character(80) :: FileTmatp(Npart), FileTmatWrite
!
  write (iOutput,"(/,2x,'Input Parameters of the Scattering Problem:',/)")
  write (iOutput,"(2x, a, 1pe13.4, a)")                                             &
 'wavelength of the ambient medium, wavelength = ', wavelength, ';'
  write (iOutput,"(2x, a, 1pe10.3, ',', 1pe10.3, a)")                               &
 'relative refractive index of the substrate, ind_refSUB = (', ind_refSUB, ');'  
  write (iOutput,*)
  do i = 1, Npart    
    write (iOutput,"(2x,'particle: ',i3)") i
           FileTmatWrite = FileTmatp(i)(14:LenString(FileTmatp(i)))
    write (iOutput,"(2x,'name of the file containing the T matrix = ',a)")          &
           FileTmatWrite 
    if (axsymp(i)) then
      write (iOutput,"(2x,'axisymmetric particle;')")              
    else
      write (iOutput,"(2x,'nonaxisymmetric particle;')")             
    end if
    write (iOutput,"(2x,'maximum expansion order, Nrank = ',i3,';')") Nrankp(i)
    write (iOutput,"(2x,'maximum azimuthal order, Mrank = ',i3,';')") Mrankp(i)
    write (iOutput,"(2x, a, 1pe10.3, a, 1pe10.3, a, 1pe10.3, a)")                   &
   'Cartesian coordinates of the particle: x = ', xp(i), ', y = ', yp(i),           &
   ', z = ', zp(i), ';'  
    write (iOutput,*)           
  end do
  write (iOutput,"(2x,'maximum expansion order of incident wave, Nrank = ',i3,';')")& 
         Nrank
  write (iOutput,"(2x,'maximum azimuthal order of incident wave, Mrank = ',i3,';')")& 
         Mrank
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'characteristic length of the particle, anorm = ', anorm, ';'  
  write (iOutput,"(2x,'axial position of the substrate, z0 = ',1pe10.3,';')") z0 
  if (TypeScat == 1)  then
    write (iOutput,"(2x, a)")                                                       &
   'illumination by a plane wave traveling in the ambient medium;'          
  else if (TypeScat == 2)  then
    write (iOutput,"(2x,'illumination by a plane wave traveling in the substrate;')")   
  end if  
  write (iOutput,*)    
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum expansion order tolerance, epsNrank = ', epsNrank, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum azimuthal order tolerance, epsMrank = ', epsMrank, ';'            
  write (iOutput,"(/)")                       
end subroutine printinputPARTSUB3DMULT
!***********************************************************************************
subroutine convergence_Nrank_MrankPARTSUBMULT (TypeScat, k, ind_ref_s, z0, snorm,   &
           Npart, FileTmatp, axsymp, xp, yp, zp, Mrankp, Nrankp, Mrank, Nrank,      &
           NXint, epsNrank, epsMrank, PrnProgress) 
  use parameters                 
  implicit none 
  integer       :: TypeScat, Npart, Mrankp(Npart), Nrankp(Npart), Mrank, Nrank,     &
                   NXint
  real(O)       :: snorm, k, epsNrank, epsMrank, z0, xp(Npart), yp(Npart), zp(Npart)
  complex(O)    :: ind_ref_s
  character(80) :: FileTmatp(Npart)   
  logical       :: axsymp(Npart), PrnProgress 
!
  integer       :: Nmax, ntpg, mtpg, i, j, Nteta, NthetaConvN, NthetaConvM, NmaxAL, &
                   Nmaxmax, ipart, Nl, Nc, Mrankpl, Nrankpl, Nmaxpl, jpart,         &
                   Mrankpl1, Nrankpl1, Nmaxpl1				      
  real(O)       :: beta0, alpha0, alfap, phiGS, x, y, z, zs, zs1, x1, y1, z1,       &
                   dx, dy, dz
  integer,allocatable    :: Nmaxp(:)   
  real(O),allocatable    :: Xint(:), pondereX(:), h(:), v(:), oldh(:), oldv(:),     &
                            oldh0(:), oldv0(:)
  complex(O),allocatable :: aa(:,:), bb(:,:), bb1(:,:), a(:,:), b(:,:), a1(:,:),    &
                            t(:,:), r(:,:), c(:), c1(:), cc(:), em(:,:)
!  
  phiGS  = 0._O
  Nteta  = 11
  beta0  = Pi / 4._O
  alpha0 = 0._O
  alfap  = Pi / 4._O  
  Nmax   = Nrank + Mrank * (2 * Nrank - Mrank + 1)  
!
  allocate (Xint(NXint), pondereX(NXint))
  call Laguerre (NXint, Xint, pondereX)      
!    
  allocate (Nmaxp(Npart))
  NmaxAL  = Nmax
  Nmaxmax = 0
  do ipart = 1, Npart
    Nmaxp(ipart) = Nrankp(ipart) + Mrankp(ipart) * (2 * Nrankp(ipart) -             &
                   Mrankp(ipart) + 1)
    Nmaxmax      = Nmaxmax + Nmaxp(ipart)
    if (Nmaxp(ipart) > NmaxAL) NmaxAL = Nmaxp(ipart)
    if (Nmaxp(ipart) > Nmax) then
      print "(/,2x, a,i3,  a)",                                                     &
     'Warning: the input values of Nrank and Mrank are too low for particle ',ipart,';'        
    end if
  end do   
  allocate (aa(2*Nmaxmax,2*Nmaxmax), bb(2*Nmaxmax,2*NmaxAL))        
  allocate (a(2*NmaxAL,2*NmaxAL), b(2*NmaxAL,2*NmaxAL), a1(2*NmaxAL,2*NmaxAL),      &
            r(2*NmaxAL,2*NmaxAL))			
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta), oldh0(Nteta),             &
            oldv0(Nteta), em(3,Nteta)) 			    
  call identity_matrix (2*Nmaxmax, aa, 2*Nmaxmax, 2*Nmaxmax)
  if (PrnProgress) call write_progress (.true., 1, Npart+4)            
  Nl = 0  
  do ipart = 1, Npart
    Mrankpl = Mrankp(ipart)
    Nrankpl = Nrankp(ipart)
    Nmaxpl  = Nmaxp(ipart)
    open (unit = iTmat, file = FileTmatp(ipart), status = "old", position = "rewind")
    if (.not. axsymp(ipart)) then
      call read_HeadFileTmat (ntpg, mtpg)
      call check_dimensionMat (ntpg, mtpg, Nmaxpl)            
      allocate (t(2*ntpg,2*mtpg))
      call read_FileTmat (ntpg, mtpg, t) 
    else
      ntpg = Nmaxpl
      mtpg = Nmaxpl
      allocate (t(2*ntpg,2*mtpg))
      call read_Tmatrix (.false., Nrankpl, Mrankpl, Nmaxpl, t, ntpg, mtpg)         
    end if
    close (unit = iTmat)                                                                                                       
    x = xp(ipart)
    y = yp(ipart)
    z = zp(ipart)
    zs = abs(z0 - z)	
    call matrixPARTSUB3D (k, ind_ref_s, zs, Mrankpl, Nrankpl, Nmaxpl,               &
         NXint, Xint, pondereX, r, NmaxAL, NmaxAL)        !r(2Nmaxpl,2Nmaxpl)    
    call product_matrices1 (2*Nmaxpl, 2*Nmaxpl, 2*Nmaxpl, t, 2*ntpg, 2*mtpg,        &
         r, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL)    !t*r--->b(2Nmaxpl,2Nmaxpl)
    call extend_matrix5 (ipart, Npart, Nmaxp, Nl, b, NmaxAL, NmaxAL,                &
         aa, Nmaxmax, Nmaxmax)                            ! extend b into aa  		     
    Nc = 0      
    do jpart = 1, Npart            
      if (jpart /= ipart) then
        Mrankpl1 = Mrankp(jpart)
        Nrankpl1 = Nrankp(jpart)
        Nmaxpl1  = Nmaxp(jpart)
        x1  = xp(jpart)
        y1  = yp(jpart)
        z1  = zp(jpart)
		dx  = x1 - x
        dy  = y1 - y
        dz  = z1 - z 
        zs1 = abs(z0 - z1)		 
		call MatTransAB_mn_m1n1 (3, k, dx, dy, dz, Mrankpl, Nrankpl, Nmaxpl,        &
             Mrankpl1, Nrankpl1, Nmaxpl1, a, NmaxAL, NmaxAL)  !T31--->a(2Nmaxpl,2Nmaxpl1)
		call MatTransAB_mn_m1n1 (1, k, dx, dy, dz, Mrankpl, Nrankpl, Nmaxpl,        &
             Mrankpl1, Nrankpl1, Nmaxpl1, a1, NmaxAL, NmaxAL) !T11--->a1(2Nmaxpl,2Nmaxpl1)
        call matrixPARTSUB3D (k, ind_ref_s, zs1, Mrankpl1, Nrankpl1, Nmaxpl1,       &
             NXint, Xint, pondereX, r, NmaxAL, NmaxAL)        !r(2Nmaxpl1,2Nmaxpl1)
        call product_matrices1 (2*Nmaxpl, 2*Nmaxpl1, 2*Nmaxpl1, a1, 2*NmaxAL,       &
             2*NmaxAL, r, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL)  !T11*r--->b(2Nmaxpl,2Nmaxpl1)
        call sum_matrices (2*Nmaxpl, 2*Nmaxpl1, a, 2*NmaxAL, 2*NmaxAL,                &
             b, 2*NmaxAL, 2*NmaxAL)                           !a+b--->a(2Nmaxpl,2Nmaxpl1)			   
        call product_matrices1 (2*Nmaxpl, 2*Nmaxpl, 2*Nmaxpl1, t, 2*ntpg, 2*mtpg,   &
             a, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL)    !t*a--->b(2Nmaxpl,2Nmaxpl1)                            
        call extend_matrix1 (ipart, jpart, Npart, Nmaxp, Nl, Nc, b, NmaxAL, NmaxAL, &
             aa, Nmaxmax, Nmaxmax)                            ! extend b into aa                          
      end if
      Nc = Nc + 2 * Nmaxp(jpart)      
    end do
    call MatTransAB_mn_m1n1 (1, k, -x, -y, -z, Mrankpl, Nrankpl, Nmaxpl, Mrank,     &
         Nrank, Nmax, a, NmaxAL, NmaxAL)                      !T11--->a(2Nmaxpl,2Nmax)	
    call product_matrices1 (2*Nmaxpl, 2*Nmaxpl, 2*Nmax, t, 2*ntpg, 2*mtpg,          &
         a, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL)        !t*a--->b(2Nmaxpl,2Nmax)      
    call extend_matrix2 (ipart, Npart, Nmaxp, Nmax, Nmaxmax, Nl, b, NmaxAL, NmaxAL, &
         bb, Nmaxmax, NmaxAL)                                 !extend b into bb  
    Nl = Nl + 2 * Nmaxp(ipart)
    if (PrnProgress) call write_progress (.false., ipart+1, Npart+4)
    deallocate (t)
  end do
  call LU_SYSTEM_DIRECT (aa, 2*Nmaxmax, 2*Nmaxmax, bb, 2*Nmaxmax, 2*NmaxAL,         &
       2*Nmaxmax, 2*Nmax)                                     !aa^(-1)*bb--->bb(2Nmaxmax,2Nmax)
  if (PrnProgress) call write_progress (.false., Npart+2, Npart+4)
  deallocate(aa, a, b, a1, r)
!
  allocate( c(2*Nmax), c1(2*Nmax), cc(2*Nmaxmax))
  if (TypeScat == 1) then
    call PWcoefficientsPARTSUB3D (beta0, alpha0, alfap, Mrank, Nrank, Nmax, c)
    call PWcoefficientsPARTSUBrefl3D (beta0, alpha0, alfap, z0, k, ind_ref_s,         &
         Mrank, Nrank, Nmax, c1)
    call sum_vectors (c, c1, 2*Nmax)                           !c + c1--->c(2Nmax)
  else    
    call PWcoefficientsPARTSUBtrans3D (beta0, alpha0, alfap, z0, k, ind_ref_s,        &
	     Mrank, Nrank, Nmax, c)                                !c(2Nmax)
  end if   
  call product_matrix_vector (2*Nmaxmax, 2*Nmax, bb, 2*Nmaxmax, 2*NmaxAL, c, cc) !bb*c--->cc(2Nmaxmax)  
  if (PrnProgress) call write_progress (.false., Npart+3, Npart+4)
!
  do j = 1, 3
    do i = 1, Nteta
      em(j,i) = zero
    end do
  end do   
  call DSCSPARTSUBMULT (1, cc, k, z0, ind_ref_s, xp, yp, zp, Mrankp, Nrankp, Nmaxp, &
       Npart, Nmaxmax, Nteta, phiGS, snorm, em,.true., h, v)  
  call DSCSPARTSUBMULT (2, cc, k, z0, ind_ref_s, xp, yp, zp, Mrankp, Nrankp, Nmaxp, &
       Npart, Nmaxmax, Nteta, phiGS, snorm, em,.true., h, v)                   
  do i = 1, Nteta
    oldh(i)  = h(i)
    oldv(i)  = v(i)
    oldh0(i) = h(i)
    oldv0(i) = v(i)    
  end do  
  write (iOutput, "(1x, a, i4, a, i3, a, i3, /)")                                   &
 'NXint = ', NXint, ', Nrank = ', Nrank, ', Mrank = ', Mrank 
  call write_DSCSPARTSUB (Nteta, h, v) 
  if (PrnProgress) call write_progress (.false., Npart+4, Npart+4)
!  
! --- (Nrank - 1) configuration --- 
  allocate(bb1(2*Nmaxmax,2*NmaxAL)) 
  call copy_matrix (2*Nmaxmax, 2*Nmax, bb, 2*Nmaxmax,2*NmaxAL,                      &
       bb1, 2*Nmaxmax,2*NmaxAL)                                        ! bb1 = bb    
  call matrix_Nrank_SUBMULT (Mrank, Nrank, Nmax, Nmaxmax, bb1, Nmaxmax, NmaxAL)	
  call product_matrix_vector (2*Nmaxmax, 2*Nmax, bb1, 2*Nmaxmax, 2*NmaxAL, c, cc) 
  do j = 1, 3
    do i = 1, Nteta
      em(j,i) = zero
    end do
  end do  
  call DSCSPARTSUBMULT (1, cc, k, z0, ind_ref_s, xp, yp, zp, Mrankp, Nrankp, Nmaxp, &
       Npart, Nmaxmax, Nteta, phiGS, snorm, em,.true., h, v)  
  call DSCSPARTSUBMULT (2, cc, k, z0, ind_ref_s, xp, yp, zp, Mrankp, Nrankp, Nmaxp, &
       Npart, Nmaxmax, Nteta, phiGS, snorm, em,.true., h, v)            
  call delta_DSCSPARTSUB (Mrank, Nteta, h, v, oldh, oldv, epsNrank, NthetaConvN)    
  write (iOutput, "(1x, a, i4, a, i3, a, i3, /)")                                   &
 'NXint = ', NXint, ', Nrank = ', Nrank - 1, ', Mrank = ', Mrank 
  call write_DSCSPARTSUB (Nteta, h, v)
  call write_NrankConvRes (NthetaConvN, Nteta - 2, epsNrank)
!  
! --- (Mrank - 1) configuration ---  
  call copy_matrix (2*Nmaxmax, 2*Nmax, bb, 2*Nmaxmax,2*NmaxAL,                      &
       bb1, 2*Nmaxmax,2*NmaxAL)                                        ! bb1 = bb    
  call matrix_Mrank_SUBMULT (Mrank, Nrank, Nmax, Nmaxmax, bb1, Nmaxmax, NmaxAL)	
  call product_matrix_vector (2*Nmaxmax, 2*Nmax, bb1, 2*Nmaxmax, 2*NmaxAL, c, cc) 
  do j = 1, 3
    do i = 1, Nteta
      em(j,i) = zero
    end do
  end do  
  call DSCSPARTSUBMULT (1, cc, k, z0, ind_ref_s, xp, yp, zp, Mrankp, Nrankp, Nmaxp, &
       Npart, Nmaxmax, Nteta, phiGS, snorm, em,.true., h, v)  
  call DSCSPARTSUBMULT (2, cc, k, z0, ind_ref_s, xp, yp, zp, Mrankp, Nrankp, Nmaxp, &
       Npart, Nmaxmax, Nteta, phiGS, snorm, em,.true., h, v)          
  call delta_DSCSPARTSUB (Mrank, Nteta, h, v, oldh0, oldv0, epsMrank, NthetaConvM)    
  write (iOutput, "(1x, a, i4, a, i3, a, i3, /)")                                   &
 'NXint = ', NXint, ', Nrank = ', Nrank, ', Mrank = ', Mrank - 1
  call write_DSCSPARTSUB (Nteta, h, v)
  call write_MrankConvRes (NthetaConvM, epsMrank)   
!
  if (NthetaConvN >= int(0.8*Nteta) .and. NthetaConvM >= int(0.8*Nteta)) then
    print "(/,2x,'Convergence criteria for Nrank and Mrank are satisfied;')"                                
  else
    print "(/,2x,'Convergence criteria for Nrank and Mrank are not satisfied;')"
  end if 
  call DiffScatCrossSectPARTSUB3DMULT (TypeScat, k, ind_ref_s, z0, snorm,           &
       xp, yp, zp, Mrankp, Nrankp, Nmaxp, Nmaxmax, Npart, Nrank, Mrank, Nmax,       &
       bb, NmaxAL )  
  deallocate (Xint, pondereX, Nmaxp, bb, bb1, c, c1, cc, h, v, oldh, oldv, oldh0,   &
              oldv0, em) 
end subroutine convergence_Nrank_MrankPARTSUBMULT      	 
! **********************************************************************************  
subroutine DiffScatCrossSectPARTSUB3DMULT (TypeScat, k, ind_ref_s, z0, snorm,       &
           xp, yp, zp, Mrankp, Nrankp, Nmaxp, Nmaxmax, Npart, Nrank, Mrank, Nmax,   &
           bb, NmaxAL )                           
  use parameters
  implicit none 
  integer       :: TypeScat, Npart, Mrankp(Npart), Nrankp(Npart), Nmaxp(Npart),     &
                   Nrank, Mrank, Nmaxmax, Nmax, NmaxAL
  real(O)       :: snorm, k, z0, xp(Npart), yp(Npart), zp(Npart)
  complex(O)    :: ind_ref_s, bb(2*Nmaxmax,2*NmaxAL) 
!  
  integer       :: NthetaGS, i, itheta
  real(O)       :: beta, alpha, alphap, phiGS
  character(80) :: FileDSCS
  logical       :: normalized, WriteInputInfo
  real(O),allocatable    :: h(:), v(:)
  complex(O),allocatable :: c(:), c1(:), cc(:), em(:,:)
! -----------------------------------------------------------------------------------
!                                 Read the input file                               !
! -----------------------------------------------------------------------------------               
  call readinputPARTSUB3DMULT1 ( beta, alpha, alphap, NthetaGS, phiGS, normalized,  &
       FileDSCS, WriteInputInfo ) 
! -----------------------------------------------------------------------------------
!                                       Main                                        !
! ----------------------------------------------------------------------------------- 
  print "(/,2x,'Scattering Characteristics of a Cluster of Particles on a Plane Surface')"  
  print "(  2x,'-----------------------------------------------------------------------')"  
!
  allocate( c(2*Nmax), c1(2*Nmax), cc(2*Nmaxmax))
  if (TypeScat == 1) then
    call PWcoefficientsPARTSUB3D (beta, alpha, alphap, Mrank, Nrank, Nmax, c)
    call PWcoefficientsPARTSUBrefl3D (beta, alpha, alphap, z0, k, ind_ref_s,        &
         Mrank, Nrank, Nmax, c1)
    call sum_vectors (c, c1, 2*Nmax)                           
  else    
    call PWcoefficientsPARTSUBtrans3D (beta, alpha, alphap, z0, k, ind_ref_s,       &
	     Mrank, Nrank, Nmax, c)                               
  end if   
  call product_matrix_vector (2*Nmaxmax, 2*Nmax, bb, 2*Nmaxmax, 2*NmaxAL, c, cc) 
!  
  open (unit = iDSCS, file = FileDSCS, status = "replace")
  allocate (h(NthetaGS), v(NthetaGS), em(3,NthetaGS)) 
  do i = 1, 3
    do itheta = 1, NthetaGS
      em(i,itheta) = zero
    end do
  end do 
  call DSCSPARTSUBMULT (1, cc, k, z0, ind_ref_s, xp, yp, zp, Mrankp, Nrankp, Nmaxp, &
       Npart, Nmaxmax, NthetaGS, phiGS, snorm, em, normalized, h, v)  
  call DSCSPARTSUBMULT (2, cc, k, z0, ind_ref_s, xp, yp, zp, Mrankp, Nrankp, Nmaxp, &
       Npart, Nmaxmax, NthetaGS, phiGS, snorm, em,normalized, h, v)                         
  if (WriteInputInfo) call inputDSCS3DMULT (TypeScat, k, ind_ref_s, z0, Mrank,      &
                           Nrank, phiGS, beta, alpha, alphap, snorm, normalized )
  call write_DSCSPARTSUB1 (NthetaGS, normalized, h, v)                    
  close (unit = iDSCS) 
  deallocate (h, v, em, c, c1, cc)
end subroutine DiffScatCrossSectPARTSUB3DMULT
! **********************************************************************************
subroutine inputDSCS3DMULT (TypeScat, k, ind_ref_s, z0, Mrank, Nrank, phiGS, beta,  &
           alpha, alphap, snorm, normalized)                             
  use parameters
  implicit none
  integer       :: TypeScat, Mrank, Nrank
  real(O)       :: k, phiGS, beta, alpha, alphap, snorm, z0, wavelength, anorm, grd
  complex(O)    :: ind_ref_s                   
  logical       :: normalized
!
  wavelength = 2._O * Pi / k
  anorm = sqrt(snorm / Pi) / k
  grd   = 180._O / Pi
  write (iDSCS,"(/,2x,'Input Parameters of the Scattering Problem:',/)")
  write (iDSCS,"(2x,'wavelength of the ambient medium, wavelength = ',1pe13.4,';')")& 
         wavelength
  write (iDSCS,"(2x, a, 1pe10.3, ',', 1pe10.3, a)")                                 &
 'relative refractive index of the substrate, ind_refSUB = (', ind_ref_s, ');' 
  write (iDSCS,"(2x,'maximum expansion order of incident wave, Nrank = ',i3,';')")  &
         Nrank
  write (iDSCS,"(2x,'maximum azimuthal order of incident wave, Mrank = ',i3,';')")  &
         Mrank           
  write (iDSCS,"(2x,'axial position of the substrate, z0 = ',1pe10.3,';')") z0  
  if (TypeScat == 1)  then
    write (iDSCS,"(2x, a)")                                                         &
   'illumination by a plane wave traveling in the ambient medium;'
    write (iDSCS,"(2x,'polar incident angle, beta = ',f7.2,';')") beta * grd                          
  else if (TypeScat == 2)  then
    write (iDSCS,"(2x,'illumination by a plane wave traveling in the substrate;')")
    write (iDSCS,"(2x,'polar incident angle, beta = ',f7.2,';')") beta * grd       
    write (iDSCS,"(2x, a, f7.2, a)")                                                &
   'polar propagation angle of the incident wave, Pi - beta = ', 180._O - beta * grd, ';'        
  end if
  write (iDSCS,"(2x,'azimuthal incident angle, alpha = ',f7.2,';')") alpha * grd 
  write (iDSCS,"(2x,'polarization angle of the incident wave, alphap = ',f7.2,';')")&
         alphap * grd
  write (iDSCS,"(2x,'scattering plane, phiGS = ', f7.2,';')") phiGS * grd  
  write (iDSCS,"(2x, a, 1pe10.3, a)")                                               & 
 'characteristic length of the scatterer, anorm = ', anorm, ';'   
  if (normalized) write (iDSCS,"(2x, a, 1pe13.4, a)")                               &
 'normalization constant, pi * anorm**2 = ', Pi * anorm * anorm, ';'
  write (iDSCS,*)        
end subroutine inputDSCS3DMULT 
! **********************************************************************************
subroutine readinputPARTSUB3DMULT1 ( beta, alpha, alphap, NthetaGS, phiGS,          &
           normalized, FileDSCS, WriteInputInfo ) 
  use parameters
  implicit none     
  integer       :: NthetaGS, ios
  real(O)       :: beta, alpha, alphap, phiGS, deg
  character(80) :: FileDSCS, string
  logical       :: normalized, WriteInputInfo, XFindPar  
! -----------------------------------------------------------------------------------
!                        Read the input file FileInputPARTSUB3DMULT                 !
! -----------------------------------------------------------------------------------             
  open (unit = iInputPARTSUB, file = FileInputPARTSUB3DMULT, status = "old",        &
        position = "rewind")                
  beta   = 0._O
  alpha  = 0._O  
  alphap = 0._O
  string = 'IncWave'
  if (XFindPar (iInputPARTSUB, string)) then
    read (iInputPARTSUB, *, iostat = ios) beta
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable beta;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) alpha
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable alpha;')"
      stop
    end if    
    read (iInputPARTSUB, *, iostat = ios) alphap
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable alphap;')"
      stop
    end if         
  else
    print "(/,2x,'Group name IncWave not found;')"
    stop  
  end if 
  call check_incdir_partsub (beta)
  call check_polarization_angle (alphap) 
  deg    = Pi / 180._O
  beta   = beta   * deg
  alpha  = alpha  * deg   
  alphap = alphap * deg
! 
  NthetaGS = 181 
  phiGS    = 0._O  
  normalized = .true.
  FileDSCS   = '../OUTPUTFILES/DSCSpartsub.dat'
  WriteInputInfo = .true.
  string     = 'DSCS'
  if (XFindPar (iInputPARTSUB, string)) then
    read (iInputPARTSUB, *, iostat = ios) NthetaGS
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable NthetaGS;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) phiGS
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable phiGS;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) normalized
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable normalized;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) FileDSCS
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable FileDSCS;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) WriteInputInfo
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable WriteInputInfo;')"
      stop
    end if
  else
    print "(/,2x,'Group name DSCSEM not found;')"
    stop  
  end if  
  call check_azimuthal_plane (phiGS)
  phiGS = phiGS * deg                 
  close (unit = iInputPARTSUB)
end subroutine readinputPARTSUB3DMULT1          
! **********************************************************************************
subroutine DSCSPARTSUBMULT (TypeField, cc, wavenumber, z0, ind_ref, xp, yp, zp,     &
           Mrankp, Nrankp, Nmaxp, Npart, Nmaxmax, NthetaGS, phiAZIMUT, snorm, em,   &
           normalized, h, v)
!------------------------------------------------------------------------------------     
! The routine computes the differential scattering cross sections of a cluster of   !
! particles on a plane substrate in the azimuthal plane phiAZIMUTH.                 !
! -----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: TypeField, Npart, Nmaxmax, Mrankp(Npart), Nrankp(Npart),            &
                Nmaxp(Npart), NthetaGS
  real(O)    :: phiAZIMUT, snorm, h(NthetaGS), v(NthetaGS), wavenumber, z0,         &
                xp(Npart), yp(Npart), zp(Npart)      
  complex(O) :: ind_ref, em(3,NthetaGS), cc(2*Nmaxmax)
  logical    :: normalized
!
  integer    :: itheta, m, k, N0, l, Mrankpl, Nrankpl, Nmaxpl, Nl, ipart
  real(O)    :: thetaGS, phiGS, fact, x, y, z, zm, arg
  complex(O) :: sum(3), phase
  complex(O),allocatable :: c(:), Minf(:,:), Ninf(:,:)      
!
  if (normalized) then
    fact = snorm
  else 
    fact = wavenumber * wavenumber
  end if
  fact = 1._O / fact  
  do itheta = 1, NthetaGS  
    thetaGS = Pi / 2._O + real(itheta - 1,O) * Pi / real(NthetaGS - 1,O)
    phiGS   = phiAZIMUT
    if (thetaGS > Pi) then
      thetaGS = 2._O * Pi - thetaGS
      phiGS   = phiAZIMUT + Pi
    end if                                  
    sum(1) = zero
    sum(2) = zero
    sum(3) = zero
!
	Nl = 0
    do ipart = 1, Npart
      Mrankpl = Mrankp(ipart)
      Nrankpl = Nrankp(ipart)
      Nmaxpl  = Nmaxp(ipart)                                                                                                      
      x  = xp(ipart)
      y  = yp(ipart)
      z  = zp(ipart)
	  zm = 2._O * abs(z0 - z) + z
!
      allocate (Minf(3,Nmaxpl), Ninf(3,Nmaxpl)) 
	  if (TypeField == 1) then               
        call MN_infinit_complete (thetaGS, phiGS, Mrankpl, Nrankpl, Nmaxpl, .true., &
             Minf, Ninf)
      else if (TypeField == 2) then   
        call MN_infinit_reflection_complete1 (wavenumber, ind_ref, thetaGS, phiGS,  &
             Mrankpl, Nrankpl, Nmaxpl, Minf, Ninf)	    
      end if 
!
      allocate( c(2*Nmaxpl) )
      do k = 1, 2*Nmaxpl
        c(k) = cc(Nl+k)
      end do
!
      if (TypeField == 1) then               
        arg = sin(thetaGS) * cos(phiGS) * x + sin(thetaGS) * sin(phiGS) * y +       &
		      cos(thetaGS) *  z    
      else if (TypeField == 2) then            
        arg = sin(thetaGS) * cos(phiGS) * x + sin(thetaGS) * sin(phiGS) * y +       &
		      cos(thetaGS) *  zm    
      end if       	  		  
      phase = exp( - im * arg * wavenumber )
!
      do m = 0, Mrankpl
        if (m == 0) then
          do k = 1, Nrankpl
            sum(1) = sum(1) + phase * (Minf(1,k) * c(k) + Ninf(1,k) * c(Nmaxpl+k))
            sum(2) = sum(2) + phase * (Minf(2,k) * c(k) + Ninf(2,k) * c(Nmaxpl+k))
            sum(3) = sum(3) + phase * (Minf(3,k) * c(k) + Ninf(3,k) * c(Nmaxpl+k))
          end do
        else          
          N0 = Nrankpl + (m - 1) * (2 * Nrankpl - m + 2)
          do l = 1, 2
            do k = 1, Nrankpl - m + 1
              sum(1) = sum(1) + phase * (Minf(1,N0+k) * c(N0+k) +                   &
			                             Ninf(1,N0+k) * c(Nmaxpl+N0+k))
              sum(2) = sum(2) + phase * (Minf(2,N0+k) * c(N0+k) +                   &
			                             Ninf(2,N0+k) * c(Nmaxpl+N0+k))
              sum(3) = sum(3) + phase * (Minf(3,N0+k) * c(N0+k) +                   &
			                             Ninf(3,N0+k) * c(Nmaxpl+N0+k))
            end do           
            N0 = N0 + Nrankpl - m + 1
          end do
        end if
      end do                          
      Nl = Nl + 2 * Nmaxp(ipart)
	  deallocate(c, Minf, Ninf)
    end do
    em(1,itheta) = em(1,itheta) + sum(1)
    em(2,itheta) = em(2,itheta) + sum(2)
    em(3,itheta) = em(3,itheta) + sum(3)                        
    h(itheta) = abs(em(2,itheta))**2 * fact
    v(itheta) = abs(em(3,itheta))**2 * fact      
  end do       
end subroutine DSCSPARTSUBMULT
!***********************************************************************************
subroutine MN_infinit_reflection_complete1 (wavenumber, ind_ref, theta, phi, Mrank, &
           Nrank, Nmax, Minf, Ninf)
!-----------------------------------------------------------------------------------
! The routine computes the reflected localized vector spherical wave functions     !
! at infinity (excepting the factor exp(j*k*R)/(k*R)) and the azimuthal modes      !
! m = 0,+1,-1,...,+Mrank,-Mrank.                                                   !
!                                                                                  !
! Input parameters:                                                                !
! - wavenumber (real) - wave number in the ambient medium.                         !
! - ind_ref (complex) - relative refractive index of the substrate.                !
! - theta, phi (real variables) - polar angles.                                    !
! - Mrank (integer) - number of azimuthal modes.                                   !
! - Nrank (integer) - maximum expansion order.                                     !
! - Nmax (integer) - specifies the dimension of the vectors.                       !
!                                                                                  !
! Output parameters:                                                               !
! - Minf, Ninf (complex arrays) - vector spherical wave functions.                 !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax, k, m, N0, ml, l, n
  real(O)    :: theta, theta0, phi, wavenumber, arga, mlr, nm
  complex(O) :: ind_ref, Rpar, Rperp, Tpar, Tperp, cost, sint, cost0, fact,         &
                factc, factt, factp, Minf(3,Nmax), Ninf(3,Nmax)
  real(O),allocatable :: Pnm(:), dPnm(:), pinm(:), taunm(:)       
!
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank))
  theta0 = Pi - theta
  cost0 = cmplx(cos(theta0),0.0,O)
  call Fresnel_aer_sub (cost0, ind_ref, Rpar, Rperp, Tpar, Tperp, cost, sint)       
  do m = 0, Mrank                  
    call leg_normalized (theta0, m, Nrank, Pnm, dPnm, pinm, taunm)
    if (m == 0) then
      do k = 1, Nrank
        n     = k               
        nm    = real(2 * n * (n + 1),O)         
        nm    = 1._O / sqrt(nm)             
        factc = (-im)**(n + 1) * nm
        factt = factc * taunm(n)                        
        Minf(1,k) =   zero
        Minf(2,k) =   zero
        Minf(3,k) = - factt * Rperp
        Ninf(1,k) =   zero
        Ninf(2,k) =   im * factt * Rpar
        Ninf(3,k) =   zero
      end do    
    else
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      ml = m
      do l = 1, 2
        mlr   = real(ml,O)
        arga  = mlr * phi
        fact  = exp(im * arga) 
        do k = 1, Nrank - m + 1
          n     = m + k - 1
          nm    = real(2 * n * (n + 1),O)               
          nm    = 1._O / sqrt(nm)            
          factc = fact * (-im)**(n + 1) * nm    
          factp = factc * mlr * pinm(n)
          factt = factc * taunm(n)
          Minf(1,N0+k) =   zero
          Minf(2,N0+k) =   im * factp * Rpar
          Minf(3,N0+k) = - factt * Rperp
          Ninf(1,N0+k) =   zero
          Ninf(2,N0+k) =   im * factt * Rpar
          Ninf(3,N0+k) = - factp * Rperp         
        end do
        N0 =   N0 + Nrank - m + 1
        ml = - ml       
      end do
    end if                         
  end do
  deallocate (Pnm, dPnm, pinm, taunm)      
end subroutine MN_infinit_reflection_complete1
!***********************************************************************************
subroutine matrix_Nrank_SUBMULT (Mrank, Nrank, Nmax, Nmaxmax, a, nap, map)
  use parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax, Nmaxmax, nap, map, m, i, l, N0
  complex(O) :: a(2*nap,2*map)
!
  do m = 0, Mrank
    if (m == 0) then
      do i = 1, 2*Nmaxmax
        a(i,Nrank)      = zero
        a(i,Nrank+Nmax) = zero
      end do
    else          
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      do l = 1, 2
        do i = 1, 2*Nmaxmax
          a(i,N0+Nrank-m+1)      = zero          
          a(i,N0+Nrank-m+1+Nmax) = zero
        end do
        N0 = N0 + Nrank - m + 1
      end do              
    end if
  end do   
end subroutine matrix_Nrank_SUBMULT
!***********************************************************************************
subroutine matrix_Mrank_SUBMULT (Mrank, Nrank, Nmax, Nmaxmax, a, nap, map)
  use parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax, Nmaxmax, nap, map, m, i, j, l, N0
  complex(O) :: a(2*nap,2*map)
!
  m  = Mrank
  N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
  do l = 1, 2           
    do i = 1, 2*Nmaxmax
      do j = 1, Nrank - m + 1
        a(i,j+N0)           = zero
        a(i,j+N0+Nmax)      = zero
      end do
    end do
    N0 = N0 + Nrank - m + 1
  end do                   
end subroutine matrix_Mrank_SUBMULT
!*********************************************************************************** 
subroutine extend_matrix5 (ipart, Npart, Nmaxp, Nl, a, nap, map, aa, naap, maap)
!-----------------------------------------------------------------------------------      
! Extend the matrix a(2*Nmaxp(ipart),2*Nmaxp(ipart)) into the global matrix        !
! aa(2*Nmaxmax,2*Nmaxmax).                                                         !
!-----------------------------------------------------------------------------------      
  use parameters
  implicit none 
  integer    :: ipart, Npart, Nmaxp(Npart), Nl, nap, map, naap, maap, i, j
  complex(O) :: a(2*nap,2*map), aa(2*naap,2*maap)  
!
  do i = 1, 2*Nmaxp(ipart)
    do j = 1, 2*Nmaxp(ipart)
      aa(i+Nl,j+Nl) = aa(i+Nl,j+Nl) - a(i,j)
    end do        
  end do            
end subroutine extend_matrix5