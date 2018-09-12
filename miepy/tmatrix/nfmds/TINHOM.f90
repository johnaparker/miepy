subroutine TINHOM
!------------------------------------------------------------------------------------
! 1. General Considerations                                                         !
! --------------------------                                                        !
! TINHOM is a routine for computing the T matrix and the scattering characteristics !
! of an inhomogeneous, dielectric, axisymmetric particle with an arbitrarily shaped !
! inclusion The program supports calculation for dielectric and axisymmetric host   !
! particles with REAL refractive indices. The main feature of this routine is that  !
! the T matrix of the inclusion is provided as input parameter in the file FileTmat.!
! The inclusion can be a homogeneous, axisymmetric or nonaxisymmetric particle, a   !
! composite or a layered scatterer and an aggregate. This choice enables the        !
! analysis of particles with complex structure and enhances the flexibility of the  !
! program. For computing the T matrix of the inclusion, the refractive index of the !
! ambient medium (the parameter ind_refMed in the corresponding T-matrix program)   !
! must be identical with the relative refractive index of host particle (the        !
! parameter indrefRel in the present T-matrix program). As an example we consider   !
! the scattering problem of an inhomogeneous, axisymmetric particle with an         !
! axisymmetric inclusion. The refractive index of the ambient medium is n_medium,   !
! the relative refractive index of the host particle with respect to the ambient    !
! medium is nRel_host (REAL), while the RELATIVE REFRACTIVE INDEX OF THE INCLUSION  !
! WITH RESPECT TO THE HOST PARTICLE is nRel_inclusion. The T matrix of the          !
! inclusion will be computed with the TAXSYM routine by setting                     !
! ind_refMed = nRel_host and ind_refRel = nRel_inclusion. Note that the T matrix of !
! a spherical inclusion must be computed with the routine TAXSYM and not with the   !
! program TSPHERE (because the T matrix of the inclusion must be specified as a     !
! two-dimensional array and the routine TSPHERE provides an one-dimensional array). ! 
!                                                                                   !
! The T matrix of an inhomogeneous particle can be expressed in terms of the T      !
! matrices of the host particle and the inclusion as                                !
!                                                                                   !
!            T =  { T_host - Q13 * [TR(O1,O2) * T_incl * TR(O2,O1)] * Q31^(-1) }    !
!               * { I - Q33 * [TR(O1,O2) * T_incl * TR(O2,O1)] * Q31^(-1) }^(-1),   !
!                                                                                   !
! where T_host is the transition matrix of the host particle, T_incl is the         !
! transition matrix of the inclusion and TR(O1,O2) is the matrix transforming the   !
! vector spherical wave functions from the origin O1 of the host particle to the    !
! origin O2 of the inclusion (through rotations and translations).                  ! 
!                                                                                   !
! For axisymmetric host particles, the surface integrals reduce to line integrals   !
! over the generatrix. The integrals are evaluated by means of Gauss-Legendre       !
! quadrature method, and note that for host particles with a plane of symmetry      !
! perpendicular to the axis of rotation, the integration can be performed along the !
! half-generatrix curve. The numerical stability and accuracy of T-matrix           !
! computations for an nonsmooth generatrix can be enhanced by dividing the          !
! generatrix into piecewise smooth curves and applying separate quadrature formulas !
! to each curve.                                                                    !
!                                                                                   !
! Host particle geometries currently supported include: spheroids, cylinders and    !
! rounded oblate cylinders. As for axisymmetric particles, the routines             !
! "interpolation_listAXSYM" and "elem_geomAXSYM" provide the required geometry      !
! parameters. For more details, consult the TAXSYM routine.                         !
!                                                                                   !
! 2. Convergence Test                                                               !
! --------------------                                                              !
! Convergence tests precede the scattering characteristics calculation. The         !
! parameters of the convergence procedure correspond to the host particle and they  !
! are:                                                                              !
! - the number of integration points for computing integrals over the generatrix,   !
!   Nint,                                                                           !
! - the maximum expansion order Nrank, and                                          !
! - the maximum azimuthal order Mrank.                                              !
! The significance of the number of integration points Nint is given in the routine !
! "TAXSYM.f90".                                                                     ! 
!                                                                                   !
! For convergence tests, the incident wave is assumed to be a plane wave traveling  !
! along the Z-axis of the global coordinate system and the scattering               !
! characteristics are computed in the azimuthal plane  phi = 0°. The Euler          !
! orientation angles of the host particle are alpha = beta = 45° and gamma = 0°,    !
! and interactive convergence tests over Nint, Nrank and Mrank are performed.       !
!                                                                                   !  
! Nrank and Mrank are assumed to not depend on Nrank1 and Mrank1, which are the     !
! maximum expansion and azimuthal orders of the inclusion. The convergence tests    !
! are carried out over Nint, Nrank and Mrank, while Nrank1 and Mrank1 are kept      !
! constant. For the integration test, the scattering problem is solved for          !
! Nint and Nint + dNint, while Nrank and Mrank are constant. For the convergence    !
! test over the expansion order, the scattering problem is solved for the pairs     !
! (Nrank, Mrank) and (Nrank - 1, Mrank). Finally, for the azimuthal order test,     !
! the cases (Nrank, Mrank) and (Nrank, Mrank - 1) are considered. The normalized    !
! differential scattering cross section will be checked at 20° increments for       !
! convergence within epsX (epsNint,epsNrank or epsMrank) tolerance. If the          !
! calculated results are converged within this tolerance at 80% of the scattering   !
! angles, then convergence is achieved. The T matrix is stored for later use by     !
! other programs, and the values of Nrank and Mrank are printed to the screen and   !
! to the T-matrix information file (see "Description.txt"). These values together   !
! with the T matrix serve as INPUT PARAMETERS for other programs.                   !
!                                                                                   !
! 3. Estimates of Nint, Nrank and Mrank                                             !
! --------------------------------------                                            !
! The above convergence tests require estimates of Nint, Nrank and Mrank.           ! 
! These estimates must be supplied by the user, and for this purpose, Wiscombe's    !
! truncation limit criterion [W. J. Wiscombe, Improved Mie scattering algorithms,   !
! Applied Optics, 19, 1505-1509, 1980] can be used.                                 !
!                                                                                   !
! The truncation limit criterion proposed by Wiscombe provides the estimate         !
!                                                                                   !
!                NrankW = int(x + 4.05 * x**0.33 + 2),                              !
!                                                                                   !
! where x is the size parameter, x = k * Rcirc, k is the wave number and Rcirc is   !
! the radius of the smallest circumscribing sphere. For almost spherical particles, !
! Mrank can be chosen as Nrank - 2,..,Nrank, while for less spherical particles     !
! Mrank can be smaller than Nrank - 2.                                              !
!                                                                                   ! 
! The value of Nint depends on the size parameter, the shape and the relative       !
! refractive index of the host particle. A conservative estimate of Nint is         !
!                                                                                   !
!                        Nint = Ndgs * Nrank,                                       !
!                                                                                   ! 
! where Ndgs = 10, 12, 14, ....                                                     !
!                                                                                   !
! 4. Strategies for Performing Convergence Tests                                    !
! -----------------------------------------------                                   !
! The following strategy for performing convergence tests can be employed:          !
! 1. choose a value of Nrank close to the value predicted by Wiscombe's truncation  !
!    limit criterion (or smaller) and set Mrank = Nrank - 2 (or smaller);           !
! 2. set Nint = Ndgs * Nrank, where Ndgs = 10, 12, ...;                             !
! 3. perform the convergence tests over Nrank and Mrank;                            ! 
! 4. if convergence is not achieved, set Nrank = Nrank + 1, maintain the relations  !
!    Nint = Ndgs * Nrank and Mrank = Nrank - 2, and go to Step 2;                   !
! 5. if convergence is achieved, perform several integration tests with the         !
!    initial Nint-value, Nint = Ndgs * Nrank, until the DSCS converges.             !
! For particles which are too extreme in size and/or aspect ratio, the errors of    !
! the extinction and scattering cross sections (instead of the differential         !
! scattering cross section) can be analyzed.                                        ! 
!                                                                                   !
! 5. Additional Comments                                                            !
! -----------------------                                                           !
! The convergence tests can be switched off by setting the logical variable         !
! DoConvTest to false. In this case, the values of Nint, Nrank and Mrank must be    !
! specified in the input file.                                                      !    
!                                                                                   !
! 6. Input Parameters                                                               !
! --------------------                                                              !
! The parameters specified in the input file "/INPUTFILES/InputINHOM.dat" are       !
! listed below.                                                                     !
!                                                                                   !
! - wavelength (real) - wavelength of the incident light in vacuo.                  !
!                                                                                   !
! - ind_refMed (real) - refractive index of the ambient medium (nonabsorbing        !
!   dielectric medium).                                                             !
!                                                                                   !
! - ind_refRel (REAL) - relative refractive index of the host particle with respect !
!   to the ambient medium.                                                          !
!                                                                                   !
! - TypeGeom (integer) - parameter specifying the type of the host particle         !
!   geometry.                                                                       !
!                                                                                   !
! - Nsurf (integer) - number of surface parameters for the host particle.           !
!                                                                                   !
! - surf (real array: surf(1), surf(2),...,surf(Nsurf)) - surface parameters        !
!   specifying the shape of the host particle. The dimension of the array surf is   !
!   NsurfPD. The integer parameter NsurfPD is specified in the routine              !
!   "Parameters.f90" and has the value NsurfPD = 10. If Nsurf > NsurfPD, the        !
!   execution is automatically terminated.                                          !
!                                                                                   !
! - Nparam (integer) - number of smooth curves forming the generatrix curve of      !
!   the host particle.                                                              !
!                                                                                   !
!   The permissive values of the surface parameters (for the supported geometries)  !
!   are summarized below.                                                           !
!                                                                                   !
!   Host particle TypeGeom   Nsurf   Nparam                surf                     !
!   spheroid         1         2       1         surf(1) - length of the semi-      !
!                                                          axis along the           !
!                                                          symmetry axis            !
!                                                surf(2) - length of the second     ! 
!                                                          semi-axis                !
!                                                                                   !
!   cylinder         2         2       3         surf(1) - half-length of           !
!                                                          the cylinder             !
!                                                surf(2) - cylinder radius          !
!                                                                                   !
!   rounded          3         2       3         surf(1) - half-length of           !
!    oblate                                                the cylinder             !
!   cylinder                                     surf(2) - cylinder radius          ! 
!                                                          including the rounded    !
!                                                          part                     !
!                                                                                   !
! - anorm (real) - characteristic length of the host particle which is used to      !
!   normalize the differential scattering cross sections.                           !
!                                                                                   !
! - Rcirc (real) - characteristic length of the host particle (usually the radius   !
!   of the smallest circumscribing sphere) which is used to compute an estimate     !
!   of the maximum expansion order by using Wiscombe's truncation limit criterion   !
!  (the size parameter is x = k * Rcirc, where k is the wave number in the ambient  !
!   medium). Alternatively, Rcirc can be chosen as the equal-volume sphere radius   !
!   or the surface-equivalent-sphere radius. This parameter must be specified if    !
!   the interactive convergence tests are performed (DoConvTest = t).               !
!                                                                                   !
! - miror (logical) - if miror = t, the host particle is mirror symmetric (the      !
!   plane of symmetry or the plane of reflection is perpendicular to the axis of    !
!   rotation).                                                                      !
!                                                                                   !
! - FileTmat (character(80)) - name of the file containing the T matrix of the      !
!   inclusion.                                                                      !
!                                                                                   !
! - axsym1 (logical) - if axsym1 = t, the inclusion is a rotationally symmetric     !
!   particle (axisymmetric particle).                                               !
!                                                                                   !
! - chiral1 (logical) - if chiral1 = t, the inclusion is an optical active particle  !
!   (chiral particle).                                                              !
!                                                                                   !
! - Nrank1, Mrank1 (integer variables) - the maximum expansion and azimuthal        !
!   orders of the inclusion.                                                        !
!                                                                                   !   
! - x1, y1, z1 (real variables) - Cartesian coordinates specifying the position     !
!   of the inclusion with respect to the coordinate system of the host particle.    !
!                                                                                   ! 
! - alpha1, beta1, gamma1 (real variables) - Euler angles specifying the            !
!   orientation of the coordinate system of the inclusion with respect to the       !
!   coordinate system of the host particle.                                         !
!                                                                                   !
! - DoConvTest (logical) - if DoConvTest = t, the interactive convergence tests     !
!   over Nint, Nrank and Mrank are invoked. An estimates of Nrank is given by       !
!   Wiscombe's truncation limit criterion. If DoConvTest = f, the values of Nint,   !
!   Nrank and Mrank must be supplied in the input file.                             !
!                                                                                   !
! - ExtThetaDom (logical) - if ExtThetaDom = t the DSCS is computed for             !
!   scattering angles ranging from 0° to 180° in the azimuthal plane phiGS = 0°, and!
!   from 180° to 0° in the azimuthal plane phiGS = 180°. The total number of        !
!   scattering angles is 10. If ExtThetaDom = f the DSCS is computed for            !
!   scattering angles ranging from 0 to 180° in the azimuthal plane phiGS = 0.      !
!                                                                                   !
! - Nint (integer) - number of integration points in computing integrals over       !
!   the generatrix curve of the host particle. This parameter is used if the        !
!   convergence tests are not performed (DoConvTest = f).                           !
!                                                                                   !
! - Nrank (integer) - maximum expansion order for the host particle. This           !
!   parameter is used if the convergence tests are not performed (DoConvTest = f).  !
!                                                                                   !
! - Mrank (integer) - maximum azimuthal order for the host particle. This           !
!   parameter is used if the convergence tests are not performed (DoConvTest = f).  !
!                                                                                   !
!   The following table explain the significance of the variables Nint, Nrank and   !
!   Mrank.                                                                          !
!                                                                                   !
!     DoConvTest       Nrank and Mrank      Nint1 and Nint2     TypeConvTest        !
!         t             console input        console input         1 or 2           !
!         f              file input           file input             0              !
!                                                                                   !  
!   The significance of the variable TypeConvTest is summarized below.              !
!                                                                                   !
!           TypeConvTest                    Significance                            !
!               0                        no convergence test                        !  
!               1                      convergence test over Nint                   ! 
!               2                 convergence test over Nrank and Mrank             !
!                                                                                   !
! - epsNint (real) - error tolerance for the integration test.                      !
!                                                                                   !
! - epsNrank (real) - error tolerance for the expansion order test.                 !
!                                                                                   !
! - epsMrank (real) - error tolerance for the azimuthal order test.                 !
!                                                                                   !
! - dNint (integer) - number of division points for the integration test. Note      !
!   that the scattering problem is solved for Nint and Nint + dNint.                ! 
!                                                                                   !
! - FileTmatG (character(80)) - name of the file to which the T matrix of the       !
!   inhomogeneous particle is written.                                              !
!                                                                                   !
! - PrnProgress (logical) - if PrnProgress = t, the progress of calculation         !
!   is printed.                                                                     !
!                                                                                   !
! Note that all input parameters specifying lengths must be provided in the same    !
! units of lengths.                                                                 !
!                                                                                   !
! 7. Logical Scheme                                                                 !
! ------------------                                                                !  
! The organization of the code is as follows:                                       !
! 1. < read input data; the input file provides                                     !
!    - the name of the file containing the T matrix of the inclusion,               !
!    - the maximum expansion and azimuthal orders Nrank1 and Mrank1,                ! 
!    - the Cartesian coordinates specifying the position of the inclusion, and      !
!    - the Euler angles specifying the orientation of the inclusion >               !  
!                                                                                   !
! 2. < T-matrix calculation >                                                       !
!                                                                                   !
!    if ( do convergence test ) then                                                !
!                                                                                   !
!       < the code computes an estimate of Nrank by using                           !
!         Wiscombe's truncation criterion >                                         !
!       < the code prompts for the estimated values of Nint, Nrank and Mrank >      !
!       < the code prompts for the type of convergence test:                        !
!          1 - over Nint, 2 - over Nrank and Mrank >                                !
!                                                                                   !
!    else if ( .not. do convergence test ) then                                     !
!                                                                                   !
!       < the input file provides the values of Nint, Nrank and Mrank >             !
!       type of convergence test = 0 (no convergence test)                          !
!                                                                                   !
!    end if                                                                         !
!    if ( do integration test, i.e., type of convergence test = 1 ) then            !
!       < the code computes the DSCS for Nint and Nint + dNint                      !
!         and write the results to the file "Output.dat" >                          !
!       < the execution is terminated >                                             !
!    end if                                                                         !
!    if ( do expansion and azimuthal order tests,                                   !
!         i.e., type of convergence test = 2 ) then                                 !
!       < the code computes the DSCS for (Nrank,Mrank), (Nrank - 1,Mrank) and       ! 
!         (Nrank, Mrank - 1), and write the results to the file "Output.dat" >      !
!       < the T matrix is stored in the file FileTmatG and the T-matrix             !
!         information file is created >                                             !
!       < Mrank and Nrank are written to the screen and to the T-matrix             !
!         information file >                                                        !      
!       < go to Step 3 >                                                            !
!    end if                                                                         !
!    if ( .not. do convergence test, i.e., type of convergence test = 0 ) then      !
!       < the T matrix is computed and stored in the file FileTmatG >               !
!       < the T-matrix information file is created >                                !
!       < Mrank and Nrank are written to the screen and to the T-matrix             !
!         information file >                                                        !
!       < go to Step 3 >                                                            !
!    end if                                                                         !
! 3. < scattering characteristics calculation >                                     !   
!------------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none 
  integer       :: TypeGeom, Nsurf, Nparam, TypeConvTest, Mrank, Nrank, Nint,       &
                   Mrank1, Nrank1, dNint                    
  real(O)       :: ks, ind_refMed, ind_refRel, wavelength, anorm, surf(NsurfPD),    &
                   snorm, epsNint, epsNrank, epsMrank, x1, y1, z1, alpha1, beta1,   &
                   gamma1, Rcirc 
  character(80) :: FileTmat, FileTmatG 
  logical       :: miror, axsym1, chiral1, DoConvTest, ExtThetaDom, PrnProgress
! -----------------------------------------------------------------------------------
!                              Read the input file                                  ! 
! -----------------------------------------------------------------------------------
  call readinputINHOM ( wavelength, ind_refMed, ind_refRel, TypeGeom,               &
       Nsurf, surf, Nparam, anorm, Rcirc, miror, FileTmat, axsym1, chiral1,         &
       Nrank1, Mrank1, x1, y1, z1, alpha1, beta1, gamma1, DoConvTest,               &
       ExtThetaDom, Nint, Nrank, Mrank, epsNint, epsNrank, epsMrank, dNint,         &
       FileTmatG, PrnProgress, TypeConvTest, ks, snorm )  
! -----------------------------------------------------------------------------------
!                                      Main                                         !
! -----------------------------------------------------------------------------------           
  if (DoConvTest) then
    open (unit = iOutput, file = FileOutput, status = "replace")
    call printinputINHOM (TypeGeom, Nsurf, surf, Nparam, miror, FileTmat, axsym1,   &
         Mrank1, Nrank1, x1, y1, z1, alpha1, beta1, gamma1, wavelength, anorm,      &
         Rcirc, ind_refMed, ind_refRel, dNint, epsNint, epsNrank, epsMrank)
    if (TypeConvTest == 1) then
      call convergence_NintINHOM (TypeGeom, ks, ind_refRel, snorm, Nsurf, surf,     &
           Nparam, Mrank, Nrank, Nint, dNint, miror, x1, y1, z1, alpha1, beta1,     &
           gamma1, Mrank1, Nrank1, FileTmat, axsym1, chiral1, epsNint, ExtThetaDom, &
           PrnProgress)
    else if (TypeConvTest == 2) then          
      call convergence_Nrank_MrankINHOM (TypeGeom, ks, ind_refRel, snorm, Nsurf,    &
           surf, Nparam, Mrank, Nrank, Nint, miror, x1, y1, z1, alpha1, beta1,      &
           gamma1, Mrank1, Nrank1, FileTmat, axsym1, chiral1, epsNrank, epsMrank,   &
           ExtThetaDom, FileTmatG, PrnProgress)
      
    end if
    close (unit = iOutput)  
  else
    call TMatrix_Nrank_MrankINHOM (TypeGeom, ks, ind_refRel, Nsurf, surf, Nparam,   &
         Mrank, Nrank, Nint, miror, x1, y1, z1, alpha1, beta1, gamma1, Mrank1,      &
         Nrank1, FileTmat, axsym1, chiral1, FileTmatG, PrnProgress)    
  end if   
end subroutine TINHOM
!***********************************************************************************
subroutine readinputINHOM ( wavelength, ind_refMed, ind_refRel, TypeGeom,           &
           Nsurf, surf, Nparam, anorm, Rcirc, miror, FileTmat, axsym1, chiral1,     &
           Nrank1, Mrank1, x1, y1, z1, alpha1, beta1, gamma1, DoConvTest,           &
           ExtThetaDom, Nint, Nrank, Mrank, epsNint, epsNrank, epsMrank, dNint,     &
           FileTmatG, PrnProgress, TypeConvTest, ks, snorm )
  use parameters
  use derived_parameters
  implicit none 
  integer       :: TypeGeom, Nsurf, Nparam, TypeConvTest, Mrank, Nrank, Nint,       &
                   Mrank1, Nrank1, dNint, i, ios, NrankW                     
  real(O)       :: ks, ind_refMed, ind_refRel, wavelength, anorm, surf(NsurfPD),    &
                   xpart, snorm, epsNint, epsNrank, epsMrank, x1, y1, z1, alpha1,   &
                   beta1, gamma1, Rcirc, x, grd 
  character(80) :: FileTmat, FileTmatG, string 
  logical       :: miror, axsym1, chiral1, DoConvTest, ExtThetaDom, PrnProgress,    &
                   XFindPar
! -----------------------------------------------------------------------------------
!                         Read the input file FileInputINHOM                        ! 
! ----------------------------------------------------------------------------------- 
  call DrvParameters 
  open (unit = iInputINHOM, file = FileInputINHOM, status = "old",                  &
        position = "rewind")   
  wavelength = 0.1_O * 2._O * Pi
  ind_refMed = 1._O
  ind_refRel = 1.2_O
  string     = 'OptProp'
  if (XFindPar (iInputINHOM, string)) then
    read (iInputINHOM, *, iostat = ios) wavelength
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable wavelength;')"
      stop
    end if
    read (iInputINHOM, *, iostat = ios) ind_refMed
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refMed;')"
      stop
    end if
    read (iInputINHOM, *, iostat = ios) ind_refRel
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refRel;')"
      stop
    end if      
  else
    print "(/,2x,'Group name OptProp not found;')"
    stop  
  end if
  ks = 2._O * Pi * ind_refMed / wavelength   
!
  TypeGeom = 1
  Nsurf = 2
  do i = 1, NsurfPD
    surf(i) = 1._O
  end do
  Nparam = 1
  anorm  = 1._O
  Rcirc  = 1._O
  miror  = .true. 
  string = 'GeomPropHost'
  if (XFindPar (iInputINHOM, string)) then    
    read (iInputINHOM, *, iostat = ios) TypeGeom
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable TypeGeom;')"
      stop
    end if         
    read (iInputINHOM, *, iostat = ios) Nsurf
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nsurf;')"
      stop
    end if
    if (Nsurf > NsurfPD) then
      print "(/,2x,'Input error: Nsurf exceeds NsurfPD;')"                                    
      stop
    end if    
    do i = 1, Nsurf
      read (iInputINHOM, *, iostat = ios) surf(i)
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable surf;')"
        stop
      end if
    end do    
    read (iInputINHOM, *, iostat = ios) Nparam
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nparam;')"
      stop
    end if
    read (iInputINHOM, *, iostat = ios) anorm
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable anorm;')"
      stop
    end if
    read (iInputINHOM, *, iostat = ios) Rcirc
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Rcirc;')"
      stop
    end if
    read (iInputINHOM, *, iostat = ios) miror
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable miror;')"
      stop
    end if                                                
  else
    print "(/,2x,'Group name GeomPropHost not found;')"
    stop  
  end if      
  call check_geomAXSYM (TypeGeom, Nsurf, Nparam)
  call check_geomAXSYMOblate (TypeGeom, Nsurf, surf) 
  call check_anorm (anorm)  
  xpart = ks * anorm
  snorm = Pi * xpart * xpart     
!
  FileTmat = '../TMATFILES/T.dat'
  axsym1   = .true.
  chiral1  = .false.
  Nrank1   = 6
  Mrank1   = 4  
  string = 'TmatIncl'
  if (XFindPar (iInputINHOM, string)) then    
    read (iInputINHOM, *, iostat = ios) FileTmat
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable FileTmat;')"
      stop
    end if         
    read (iInputINHOM, *, iostat = ios) axsym1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable axsym1;')"
      stop
    end if       
    read (iInputINHOM, *, iostat = ios) chiral1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable chiral1;')"
      stop
    end if
    read (iInputINHOM, *, iostat = ios) Nrank1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nrank1;')"
      stop
    end if
    read (iInputINHOM, *, iostat = ios) Mrank1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Mrank1;')"
      stop
    end if                                                   
  else
    print "(/,2x,'Group name TmatIncl not found;')"
    stop  
  end if                              
  call check_MrankNrank (Mrank1, Nrank1)   
!
  x1 = 0.1_O
  y1 = 0.1_O
  z1 = 0.1_O
  alpha1 = 45._O
  beta1  = 45._O
  gamma1 = 0._O
  string = 'GeomPropIncl'
  if (XFindPar (iInputINHOM, string)) then    
    read (iInputINHOM, *, iostat = ios) x1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable x1;')"
      stop
    end if         
    read (iInputINHOM, *, iostat = ios) y1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable y1;')"
      stop
    end if       
    read (iInputINHOM, *, iostat = ios) z1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable z1;')"
      stop
    end if
    read (iInputINHOM, *, iostat = ios) alpha1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable alpha1;')"
      stop
    end if
    read (iInputINHOM, *, iostat = ios) beta1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable beta1;')"
      stop
    end if 
    read (iInputINHOM, *, iostat = ios) gamma1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable gamma1;')"
      stop
    end if
  else
    print "(/,2x,'Group name GeomPropIncl not found;')"
    stop  
  end if    
  grd    = Pi / 180._O
  alpha1 = alpha1 * grd
  beta1  = beta1  * grd
  gamma1 = gamma1 * grd
!
  DoConvTest   = .true.
  ExtThetaDom  = .true.
  string       = 'ConvTest'
  if (XFindPar (iInputINHOM, string)) then
    read (iInputINHOM, *, iostat = ios) DoConvTest
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable DoConvTest;')"
      stop
    end if
    read (iInputINHOM, *, iostat = ios) ExtThetaDom
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ExtThetaDom;')"
      stop
    end if        
  else
    print "(/,2x,'Group name ConvTest not found;')"
    stop  
  end if      
!  
  if (DoConvTest) then    
    print "(/,2x,'Convergence Test for an Inhomogeneous Particle')"
    print "(  2x,'----------------------------------------------')"          
  else
    print "(/,2x,'T-Matrix Computation for an Inhomogeneous Particle')"
    print "(  2x,'--------------------------------------------------')"
  end if
!
  x = ks * Rcirc          
  NrankW = int(x + 4.05_O * x**0.33_O + 2._O)
  if (.not. DoConvTest) then    
    Nint   = 100
    Nrank  =  16
    Mrank  =   8
    string = 'NintNrankMrankHost'
    if (XFindPar (iInputINHOM, string)) then
      read (iInputINHOM, *, iostat = ios) Nint
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Nint;')"
        stop
      end if
      read (iInputINHOM, *, iostat = ios) Nrank
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Nrank;')"
        stop
      end if 
      read (iInputINHOM, *, iostat = ios) Mrank
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Mrank;')"
        stop
      end if
    else
      print "(/,2x,'Group name NintNrankMrankHost not found;')"
      stop  
    end if
    print "(/,2x,'Input values:')"
    print "(  2x,'the input values of  Nint,  Nrank and  Mrank are ',i4,',',i4)",   &
              Nint, Nrank
    print "(  2x,'and ',i4,', respectively, while the estimated value of Nrank')",  &
              Mrank 
    print "(  2x,'from Wiscombe''s criterion is ', i3,';')", NrankW              
  else  
    print "(/,2x,'Nrank estimate:')"                                                    
    print "(  2x, a, i3, a)",                                                       &
   'the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';'
    print "(/,2x,'- enter the estimated values of Nint, Nrank and Mrank, where')"
    print "(  2x,'  Nint = Ndgs * Nrank, Ndgs = 10,12,... , and for almost spherical')"
    print "(  2x,'  particles Mrank = Nrank - 2,...,Nrank, while for less spherical;')"
    print "(  2x,'  Mrank can be smaller than Nrank - 2;')"    
    call read_integer3 (Nint, Nrank, Mrank)
  end if
  call check_MrankNrank (Mrank, Nrank)
  if (Nrank < Nrank1) print "(/,2x,'Warning: Nrank is too low;')"
  if (Mrank < Mrank1) print "(/,2x,'Warning: Mrank is too low;')"
!  
  if (DoConvTest) then      
    print "(/,2x, a)",                                                              &
   '- enter the type of convergence test: 1 - Nint, 2 - Nrank and Mrank;'
    call read_integerbound (TypeConvTest, 1, 2)                      
  else
    TypeConvTest = 0
  end if
!
  epsNint  = 5.e-2_O
  epsNrank = 5.e-2_O
  epsMrank = 5.e-2_O
  dNint    = 4
  string   = 'Errors'
  if (XFindPar (iInputINHOM, string)) then
    read (iInputINHOM, *, iostat = ios) epsNint
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsNint;')"
      stop
    end if
    read (iInputINHOM, *, iostat = ios) epsNrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsNrank;')"
      stop
    end if
    read (iInputINHOM, *, iostat = ios) epsMrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsMrank;')"
      stop
    end if
    read (iInputINHOM, *, iostat = ios) dNint
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable dNint;')"
      stop
    end if
  else
    print "(/,2x,'Group name Errors not found;')"
    stop  
  end if
!     
  FileTmatG = '../TMATFILES/TG.dat'
  string    = 'Tmat'
  if (XFindPar (iInputINHOM, string)) then
    read (iInputINHOM, *, iostat = ios) FileTmatG
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable FileTmatG;')"
      stop
    end if           
  else
    print "(/,2x,'Group name Tmat not found;')"
    stop  
  end if    
!
  PrnProgress = .true.
  string      = 'PrintProgress'
  if (XFindPar (iInputINHOM, string)) then
    read (iInputINHOM, *, iostat = ios) PrnProgress
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable PrnProgress;')"
      stop
    end if           
  else
    print "(/,2x,'Group name PrintProgress not found;')"
    stop  
  end if           
  close (unit = iInputINHOM)
end subroutine readinputINHOM
!***********************************************************************************
subroutine printinputINHOM (TypeGeom, Nsurf, surf, Nparam, miror, FileTmat, axsym1, &
           Mrank1, Nrank1, x1, y1, z1, alfa1, beta1, gama1, wavelength, anorm,      &
           Rcirc, ind_refMed, ind_refRel, dNint, epsNint, epsNrank, epsMrank)
  use parameters
  implicit none
  integer       :: TypeGeom, Nsurf, Nparam, dNint, Mrank1, Nrank1, i, LenString
  real(O)       :: wavelength, anorm, ind_refMed, ind_refRel, surf(Nsurf), x1, y1,  &
                   z1, alfa1, beta1, gama1, epsNint, epsNrank, epsMrank, Rcirc
  character(80) :: FileTmat, FileTmatWrite
  logical       :: miror, axsym1
!
  write (iOutput,"(/,2x,'Input Parameters of the Scattering Problem:',/)")
  write (iOutput,"(2x,'wavelength of the free space, wavelength = ',1pe13.4,';')")  &
         wavelength
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'refractive index of the ambient medium, ind_refMed = ', ind_refMed, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'relative refractive index of the host particle, ind_refRel = ', ind_refRel, ';'
  write (iOutput,*)
  write (iOutput,"(2x,'index of the host particle geometry, TypeGeom = ',i2,';')")  &
         TypeGeom
  if (TypeGeom == 1) then
    write (iOutput,"(2x,'spheroid;')")  
  else if (TypeGeom == 2) then
    write (iOutput,"(2x,'cylinder;')")
  else if (TypeGeom == 3) then
    write (iOutput,"(2x,'rounded oblate cylinder;')")           
  end if
  write (iOutput,"(2x,'number of surface parameters, Nsurf = ',i2,';')") Nsurf
  write (iOutput,"(2x,'surface parameters:')")
  do i = 1, Nsurf   
    write (iOutput,"(2x,'surf(',i2,') = ',1pe10.3,',')") i, surf(i)    
  end do        
  if (TypeGeom == 1) then
    write (iOutput,"(2x,'where surf(1) is the half-axis along the symmetry axis')")
    write (iOutput,"(2x,'and   surf(2) is the second half-axis;')")
  else if (TypeGeom == 2 .or. TypeGeom == 3) then
    write (iOutput,"(2x,'where surf(1) is the half-length of the cylinder')")
    write (iOutput,"(2x,'and   surf(2) is the cylinder radius;')")
  end if     
  write (iOutput,"(2x,'number of integration surfaces, Nparam = ',i2,';')") Nparam
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'characteristic length of the host particle, anorm = ', anorm, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'radius of the circumscribing sphere, Rcirc = ', Rcirc, ';'
  if (miror) write (iOutput,"(2x,'mirror symmetric host particle;')")  
  write (iOutput,*)
  FileTmatWrite = FileTmat(14:LenString(FileTmat))
  write (iOutput,"(2x, a, a)")                                                      &
 'name of the file containing the T matrix of the inclusion = ', FileTmatWrite
  if (axsym1) then
    write (iOutput,"(2x,'axisymmetric inclusion;')")
    write (iOutput,"(2x,'maximum expansion order, Nrank = ',i3,';')") Nrank1
    write (iOutput,"(2x,'maximum azimuthal order, Mrank = ',i3,';')") Mrank1    
  else
    write (iOutput,"(2x,'nonaxisymmetric inclusion;')")
    write (iOutput,"(2x,'maximum expansion order, Nrank =',i3,';')") Nrank1
    write (iOutput,"(2x,'maximum azimuthal order, Mrank =',i3,';')") Mrank1     
  end if
  write (iOutput,"(2x, a, 1pe10.3, a, 1pe10.3, a, 1pe10.3, a)")                     &
 'Cartesian coordinates of the inclusion: x = ', x1, ', y = ', y1, ', z = ', z1, ';'             
  write (iOutput,"(2x, a, 1pe10.3, a, 1pe10.3, a, 1pe10.3, a)")                     &
 'Euler angles: alpha = ', alfa1 * 180._O / Pi, ', beta = ', beta1 * 180._O / Pi,   &
 ', gamma = ', gama1 * 180._O / Pi, ';'  
  write (iOutput,*)
  write (iOutput,"(2x,'integration tolerance, epsNint = ',1pe10.3,';')") epsNint
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum expansion order tolerance, epsNrank = ', epsNrank, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum azimuthal order tolerance, epsMrank = ', epsMrank, ';'          
  write (iOutput,"(2x,'integration step, dNint = ',i4,'.')") dNint
  write (iOutput,"(/)")                       
end subroutine printinputINHOM
!***********************************************************************************
subroutine convergence_NintINHOM (TypeGeom, ks, ind_ref, snorm, Nsurf, surf, Nparam,&
           Mrank, Nrank, Nint, dNint, miror, x1, y1, z1, alfa1, beta1, gama1,       &
           Mrank1, Nrank1, FileTmat, axsym1, chiral1, epsNint, ExtThetaDom,         &
           PrnProgress)
  use parameters
  implicit none
  integer       :: TypeGeom, Nsurf, Nparam, Mrank, Nrank, Nint, dNint, Mrank1,      &
                   Nrank1
  real(O)       :: ks, ind_ref, snorm, surf(Nsurf), x1, y1, z1, alfa1, beta1,       &
                   gama1, epsNint
  character(80) :: FileTmat
  logical       :: miror, axsym1, chiral1, ExtThetaDom, PrnProgress
!
  integer       :: Nmax, Nmax1, Nteta, NthetaConv, iNint, NmaxAL, nt1g, mt1g, i
  real(O)       :: k, alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,  &
                   Cext, Qext 
  complex(O)    :: ind_refC
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG(:,:), weightsG(:,:), h(:), v(:), oldh(:), oldv(:)
  complex(O),allocatable :: a(:,:), b(:,:), d(:,:), t1(:,:), s(:,:), c(:), c1(:)  
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = Pi / 4._O
  beta   = Pi / 4._O
  gama   = 0._O
  alfap  = Pi / 4._O
  Nmax   = Nrank + Mrank * (2 * Nrank - Mrank + 1) 
  Nmax1  = Nrank1 + Mrank1 * (2 * Nrank1 - Mrank1 + 1)
  NmaxAL = max(Nmax,Nmax1) 
  k = ks * ind_ref
  ind_refC = cmplx(ind_ref,0.0,O)
  call write_TypeConvHead (1)
  allocate (a(2*NmaxAL,2*NmaxAL), b(2*NmaxAL,2*NmaxAL), d(2*NmaxAL,2*NmaxAL),       &
            s(2*NmaxAL,2*NmaxAL), c(2*NmaxAL),c1(2*NmaxAL))  
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do
  open (unit = iTmat, file = FileTmat, status = "old", position = "rewind")
  if (.not. axsym1) then
    call read_HeadFileTmat (nt1g, mt1g)
    call check_dimensionMat (nt1g, mt1g, Nmax1)       
    allocate (t1(2*nt1g,2*mt1g))
    call read_FileTmat (nt1g, mt1g, t1) 
  else
    nt1g = Nmax1
    mt1g = Nmax1
    allocate (t1(2*nt1g,2*mt1g))
    call read_Tmatrix (chiral1, Nrank1, Mrank1, Nmax1, t1, nt1g, mt1g)            
  end if
  close (unit = iTmat)
  if (PrnProgress) call write_progress (.true., 1, 9)     
  call MatTransRot_TR_mn_m1n1 (1, k, x1, y1, z1, alfa1, beta1, gama1, Mrank,        &
       Nrank, Nmax, Mrank1, Nrank1, Nmax1, s, NmaxAL, NmaxAL)
  call product_matrices (2*Nmax, 2*Nmax1, 2*Nmax1, s, 2*NmaxAL, 2*NmaxAL, t1,       &
       2*nt1g, 2*mt1g)
  if (PrnProgress) call write_progress (.false., 2, 9)
  call MatRotTrans_RT_mn_m1n1 (1, k, x1, y1, z1, alfa1, beta1, gama1, Mrank1,       &
       Nrank1, Nmax1, Mrank, Nrank, Nmax, a, NmaxAL, NmaxAL)
  call product_matrices (2*Nmax, 2*Nmax1, 2*Nmax, s, 2*NmaxAL, 2*NmaxAL, a,         &
       2*NmaxAL, 2*NmaxAL)
  if (PrnProgress) call write_progress (.false., 3, 9)           
  do iNint = 1, 2
    allocate (paramG(Nparam,Nint), weightsG(Nparam,Nint), Nintparam(Nparam))
    call interpolation_listAXSYM (TypeGeom, Nsurf, surf, Nint, Nparam,              &
         Nintparam, paramG, weightsG, miror)
    call matrix_Q_sym (TypeGeom, 3, 3, ks, ind_refC, Nsurf, surf, Mrank, Nrank,     &
         Nmax, Nint, Nparam, Nintparam, paramG, weightsG, miror,.false., a,         &
         NmaxAL, NmaxAL)    
    call product_matrices (2*Nmax, 2*Nmax, 2*Nmax, a, 2*NmaxAL, 2*NmaxAL, s,        &
         2*NmaxAL, 2*NmaxAL)
    call matrix_Q_sym (TypeGeom, 3, 1, ks, ind_refC, Nsurf, surf, Mrank, Nrank,     &
         Nmax, Nint, Nparam, Nintparam, paramG, weightsG, miror,.false., d,         &
         NmaxAL, NmaxAL)
    call sum_matrices (2*Nmax, 2*Nmax, a, 2*NmaxAL, 2*NmaxAL, d, 2*NmaxAL, 2*NmaxAL)  
    if (PrnProgress) call write_progress (.false., 4+3*(iNint-1), 9)     
    call matrix_Q_sym (TypeGeom, 1, 3, ks, ind_refC, Nsurf, surf, Mrank, Nrank,     &
         Nmax, Nint, Nparam, Nintparam, paramG, weightsG, miror,.false., b,         &
         NmaxAL, NmaxAL)
    call product_matrices (2*Nmax, 2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL, s,        &
         2*NmaxAL, 2*NmaxAL)  
    call matrix_Q_sym (TypeGeom, 1, 1, ks, ind_refC, Nsurf, surf, Mrank, Nrank,     &
         Nmax, Nint, Nparam, Nintparam, paramG, weightsG, miror,.false., d,         &
         NmaxAL, NmaxAL) 
    call sum_matrices (2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL, d, 2*NmaxAL, 2*NmaxAL)  
    if (PrnProgress) call write_progress (.false., 5+3*(iNint-1), 9)     
    call LU_SYSTEM (a, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL, 2*Nmax)
    if (PrnProgress) call write_progress (.false., 6+3*(iNint-1), 9)     
    call minus_matrix (2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL)
    call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrank, Nrank,   &
         Nmax, c)    
    call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL, c, c1)
    call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfa, beta, gama, ks, snorm,   &
         ExtThetaDom,.true., h, v)
    call CQscat (c1, Mrank, Nrank, Nmax, ks, snorm, Cscat, Qscat)
    call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfa, beta, gama,            &
         alfap, ks, snorm, Cext, Qext)
    call delta_DSCS (Nteta, h, v, oldh, oldv, epsNint, NthetaConv)       
    call write_3ConvParam (Nint, Nrank, Mrank)                  
    call write_DSCS (Nteta, ExtThetaDom, h, v)
    call write_Effic (Qscat, Qext)
    Nint = Nint + dNint
    deallocate (paramG, weightsG, Nintparam)
  end do
  call write_NintConvRes (NthetaConv, Nteta, epsNint)
  deallocate (a, b, d, c, c1, t1, s, h, v, oldh, oldv)  
end subroutine convergence_NintINHOM
!***********************************************************************************
subroutine convergence_Nrank_MrankINHOM (TypeGeom, ks, ind_ref, snorm, Nsurf, surf, &
           Nparam, Mrank, Nrank, Nint, miror, x1, y1, z1, alfa1, beta1, gama1,      &
           Mrank1, Nrank1, FileTmat, axsym1, chiral1, epsNrank, epsMrank,           &
           ExtThetaDom, FileTmatG, PrnProgress)
  use parameters
  implicit none
  integer       :: TypeGeom, Nsurf, Nparam, Mrank, Nrank, Nint, Mrank1, Nrank1
  real(O)       :: ks, ind_ref, snorm, surf(Nsurf), x1, y1, z1, alfa1, beta1,       &
                   gama1, epsNrank, epsMrank
  character(80) :: FileTmat, FileTmatG
  logical       :: miror, axsym1, chiral1, ExtThetaDom, PrnProgress
!
  integer       :: Nmax, Nmax1, Nteta, i, NthetaConvN, NthetaConvM, NmaxAL,         &
                   nt1g, mt1g
  real(O)       :: k, alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,  &
                   Cext, Qext, r1, sumang
  complex(O)    :: ind_refC
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG(:,:), weightsG(:,:), h(:), v(:),                 &
                            oldh(:), oldv(:), oldh0(:), oldv0(:)
  complex(O),allocatable :: a(:,:), b(:,:), d(:,:), c(:), c1(:),                    &
                            t1(:,:), s(:,:), q11(:,:), q13(:,:),                    &
                            q31(:,:), q33(:,:)     
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = Pi / 4._O
  beta   = Pi / 4._O
  gama   = 0._O
  alfap  = Pi / 4._O
  Nmax   = Nrank + Mrank * (2 * Nrank - Mrank + 1) 
  Nmax1  = Nrank1 + Mrank1 * (2 * Nrank1 - Mrank1 + 1) 
  NmaxAL = max(Nmax,Nmax1)
  k = ks * ind_ref
  ind_refC = cmplx(ind_ref,0.0,O)               
  call write_TypeConvHead (4)
  allocate (a(2*NmaxAL,2*NmaxAL), b(2*NmaxAL,2*NmaxAL), d(2*NmaxAL,2*NmaxAL),       &
            c(2*NmaxAL), c1(2*NmaxAL))
  allocate (s(2*NmaxAL,2*NmaxAL), q11(2*NmaxAL,2*NmaxAL), q13(2*NmaxAL,2*NmaxAL),   &
            q31(2*NmaxAL,2*NmaxAL), q33(2*NmaxAL,2*NmaxAL))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta), oldh0(Nteta),             &
            oldv0(Nteta))                        
  allocate (paramG(Nparam,Nint), weightsG(Nparam,Nint), Nintparam(Nparam))
  call interpolation_listAXSYM (TypeGeom, Nsurf, surf, Nint, Nparam,                &
       Nintparam, paramG, weightsG, miror)         
! --- Nrank configuration --- 
  open (unit = iTmat, file = FileTmat, status = "old", position = "rewind")
  if (.not. axsym1) then
    call read_HeadFileTmat (nt1g, mt1g)
    call check_dimensionMat (nt1g, mt1g, Nmax1)       
    allocate (t1(2*nt1g,2*mt1g))
    call read_FileTmat (nt1g, mt1g, t1) 
  else
    nt1g = Nmax1
    mt1g = Nmax1
    allocate (t1(2*nt1g,2*mt1g))
    call read_Tmatrix (chiral1, Nrank1, Mrank1, Nmax1, t1, nt1g, mt1g)            
  end if
  close (unit = iTmat)
  open  (unit = iTmat, file = FileTmatG, status = 'replace')
  call  write_HeadFileTmat (NmaxAL, NmaxAL) 
  if (PrnProgress) call  write_progress (.true., 1, 6)
  r1 = sqrt(x1**2 + y1**2 + z1**2)
  sumang = abs(alfa1) + abs(beta1) + abs(gama1)
  if (r1 /= 0 .and. sumang /= 0) then
    call MatTransRot_TR_mn_m1n1 (1, k, x1, y1, z1, alfa1, beta1, gama1, Mrank,      &
         Nrank, Nmax, Mrank1, Nrank1, Nmax1, s, NmaxAL, NmaxAL)
  else if (r1 /= 0 .and. sumang == 0) then 
    call MatTransAB_mn_m1n1 (1, k, x1, y1, z1, Mrank, Nrank, Nmax, Mrank1, Nrank1,  &
         Nmax1, s, NmaxAL, NmaxAL)
  else if (r1 == 0 .and. sumang /= 0) then    
    call MatRot_R_mn_m1n1 (alfa1, beta1, gama1, Mrank, Nrank, Nmax, Mrank1, Nrank1, &
         Nmax1, s, NmaxAL, NmaxAL, .true.)
  else if (r1 == 0 .and. sumang == 0) then
    call identity_transformation (Mrank, Nrank, Nmax, Mrank1, Nrank1, Nmax1, s,     &
         NmaxAL, NmaxAL)    
  end if  
  call product_matrices (2*Nmax, 2*Nmax1, 2*Nmax1, s, 2*NmaxAL, 2*NmaxAL, t1,       &
       2*nt1g, 2*mt1g)
  if (PrnProgress) call write_progress (.false., 2, 6)
  if (r1 /=0 .and. sumang /= 0) then
    call MatRotTrans_RT_mn_m1n1 (1, k, x1, y1, z1, alfa1, beta1, gama1, Mrank1,     &
         Nrank1, Nmax1, Mrank, Nrank, Nmax, a, NmaxAL, NmaxAL)
  else if (r1 /=0 .and. sumang == 0) then 
    call MatTransAB_mn_m1n1 (1, k, -x1, -y1, -z1, Mrank1, Nrank1, Nmax1, Mrank,     &
         Nrank, Nmax, a, NmaxAL, NmaxAL)
  else if (r1 == 0 .and. sumang /= 0) then    
    call MatRot_R_mn_m1n1 (-gama1, -beta1, -alfa1, Mrank1, Nrank1, Nmax1, Mrank,    &
         Nrank, Nmax, a, NmaxAL, NmaxAL, .true.)
  else if (r1 == 0 .and. sumang == 0) then
    call identity_transformation (Mrank1, Nrank1, Nmax1, Mrank, Nrank, Nmax, a,     &
         NmaxAL, NmaxAL)    
  end if   
  call product_matrices (2*Nmax, 2*Nmax1, 2*Nmax, s, 2*NmaxAL, 2*NmaxAL, a,         &
       2*NmaxAL, 2*NmaxAL) 
  if (PrnProgress) call write_progress (.false., 3, 6)          
  call matrix_Q_sym (TypeGeom, 3, 3, ks, ind_refC, Nsurf, surf, Mrank, Nrank, Nmax, &
       Nint, Nparam, Nintparam, paramG, weightsG, miror,.false., a, NmaxAL, NmaxAL)
  call copy_matrix (2*Nmax, 2*Nmax, a, 2*NmaxAL, 2*NmaxAL, q33, 2*NmaxAL, 2*NmaxAL)
  call product_matrices (2*Nmax, 2*Nmax, 2*Nmax, a, 2*NmaxAL, 2*NmaxAL, s,          &
       2*NmaxAL, 2*NmaxAL)
  call matrix_Q_sym (TypeGeom, 3, 1, ks, ind_refC, Nsurf, surf, Mrank, Nrank, Nmax, &
       Nint, Nparam, Nintparam, paramG, weightsG, miror,.false., q31, NmaxAL, NmaxAL)
  call sum_matrices (2*Nmax, 2*Nmax, a, 2*NmaxAL, 2*NmaxAL, q31, 2*NmaxAL, 2*NmaxAL)  
  if (PrnProgress) call write_progress (.false., 4, 6)
  call matrix_Q_sym (TypeGeom, 1, 3, ks, ind_refC, Nsurf, surf, Mrank, Nrank, Nmax, &
       Nint, Nparam, Nintparam, paramG, weightsG, miror,.false., b, NmaxAL, NmaxAL)
  call copy_matrix (2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL, q13, 2*NmaxAL, 2*NmaxAL) 
  call product_matrices (2*Nmax, 2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL, s,          &
       2*NmaxAL, 2*NmaxAL)  
  call matrix_Q_sym (TypeGeom, 1, 1, ks, ind_refC, Nsurf, surf, Mrank, Nrank, Nmax, &
       Nint, Nparam, Nintparam, paramG, weightsG, miror,.false., q11, NmaxAL, NmaxAL) 
  call sum_matrices (2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL, q11, 2*NmaxAL, 2*NmaxAL)  
  if (PrnProgress) call write_progress (.false., 5, 6)
  call LU_SYSTEM (a, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL, 2*Nmax)
  if (PrnProgress) call write_progress (.false., 6, 6)
  call minus_matrix (2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL)
  call write_FileTmat (NmaxAL, NmaxAL, b)       
  call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrank, Nrank,     &
       Nmax,c)    
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL, c, c1)
  call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfa, beta, gama, ks, snorm,     &
       ExtThetaDom,.true., h, v)
  call CQscat (c1, Mrank, Nrank, Nmax, ks, snorm, Cscat, Qscat)
  call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfa, beta, gama,              &
       alfap, ks, snorm, Cext, Qext)       
  call write_3ConvParam (Nint, Nrank, Mrank)    
  call write_DSCS (Nteta, ExtThetaDom, h, v)
  call write_Effic (Qscat, Qext)
  do i = 1, Nteta
    oldh(i)  = h(i)
    oldv(i)  = v(i)
    oldh0(i) = h(i)
    oldv0(i) = v(i)
  end do 
  close (unit = iTmat)  
! --- (Nrank - 1) configuration ---
  if (PrnProgress) call write_progress_low
  call copy_matrix (2*Nmax, 2*Nmax, q33, 2*NmaxAL, 2*NmaxAL, a, 2*NmaxAL, 2*NmaxAL)
  call matrix_Nrank_1_right (Mrank, Nrank, Nmax, a, NmaxAL, NmaxAL)
  call product_matrices (2*Nmax, 2*Nmax, 2*Nmax, a, 2*NmaxAL, 2*NmaxAL, s,          &
       2*NmaxAL, 2*NmaxAL)
  call copy_matrix (2*Nmax, 2*Nmax, q31, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL)
  call matrix_Nrank_1_right (Mrank, Nrank, Nmax, b, NmaxAL, NmaxAL)
  call sum_matrices (2*Nmax, 2*Nmax, a, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL)
  call matrix_Nrank_1_left (Mrank, Nrank, Nmax, a, NmaxAL, NmaxAL)  
  call copy_matrix (2*Nmax, 2*Nmax, q13, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL)
  call matrix_Nrank_1_right (Mrank, Nrank, Nmax, b, NmaxAL, NmaxAL)
  call product_matrices (2*Nmax, 2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL, s,          &
       2*NmaxAL, 2*NmaxAL)
  call copy_matrix (2*Nmax, 2*Nmax, q11, 2*NmaxAL, 2*NmaxAL, d, 2*NmaxAL, 2*NmaxAL)
  call matrix_Nrank_1_right (Mrank, Nrank, Nmax, d, NmaxAL, NmaxAL)
  call sum_matrices (2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL, d, 2*NmaxAL, 2*NmaxAL)
  call matrix_Nrank_1_right (Mrank, Nrank, Nmax, b, NmaxAL, NmaxAL)  
  call LU_SYSTEM (a, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL, 2*Nmax)
  call minus_matrix (2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL) 
  call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrank, Nrank,     &
       Nmax, c)    
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL, c, c1)
  call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfa, beta, gama, ks, snorm,     &
       ExtThetaDom,.true., h, v)
  call CQscat (c1, Mrank, Nrank, Nmax, ks, snorm, Cscat, Qscat)
  call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfa, beta, gama,              &
       alfap, ks, snorm, Cext, Qext)
  call delta_DSCS (Nteta, h, v, oldh, oldv, epsNrank, NthetaConvN)
  call write_3ConvParam (Nint, Nrank - 1, Mrank)
  call write_DSCS (Nteta, ExtThetaDom, h, v)
  call write_Effic (Qscat, Qext)
  call write_NrankConvRes (NthetaConvN, Nteta, epsNrank)
! --- (Mrank - 1) configuration ---
  call copy_matrix (2*Nmax, 2*Nmax, q33, 2*NmaxAL, 2*NmaxAL, a, 2*NmaxAL, 2*NmaxAL)
  call matrix_Mrank_1_right (Mrank, Nrank, Nmax, a, NmaxAL, NmaxAL)
  call product_matrices (2*Nmax, 2*Nmax, 2*Nmax, a, 2*NmaxAL, 2*NmaxAL, s,          &
       2*NmaxAL, 2*NmaxAL)
  call copy_matrix (2*Nmax, 2*Nmax, q31, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL)
  call matrix_Mrank_1_right (Mrank, Nrank, Nmax, b, NmaxAL, NmaxAL)
  call sum_matrices (2*Nmax, 2*Nmax, a, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL)
  call matrix_Mrank_1_left (Mrank, Nrank, Nmax, a, NmaxAL, NmaxAL)
  call copy_matrix (2*Nmax, 2*Nmax, q13, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL)
  call matrix_Mrank_1_right (Mrank, Nrank, Nmax, b, NmaxAL, NmaxAL)
  call product_matrices (2*Nmax, 2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL, s,          &
       2*NmaxAL, 2*NmaxAL)
  call copy_matrix (2*Nmax, 2*Nmax, q11, 2*NmaxAL, 2*NmaxAL, d, 2*NmaxAL, 2*NmaxAL)
  call matrix_Mrank_1_right (Mrank, Nrank, Nmax, d, NmaxAL, NmaxAL)
  call sum_matrices (2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL, d, 2*NmaxAL, 2*NmaxAL)
  call matrix_Mrank_1_right (Mrank, Nrank, Nmax, b, NmaxAL, NmaxAL)
  call LU_SYSTEM (a, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL, 2*Nmax)
  call minus_matrix (2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL)
  call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrank, Nrank,     &
       Nmax, c)    
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL, c, c1)
  call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfa, beta, gama, ks, snorm,     &
       ExtThetaDom,.true., h, v)
  call CQscat (c1, Mrank, Nrank, Nmax, ks, snorm, Cscat, Qscat)
  call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfa, beta, gama,              &
       alfap, ks, snorm, Cext, Qext)
  call delta_DSCS (Nteta, h, v, oldh0, oldv0, epsMrank, NthetaConvM)
  call write_3ConvParam (Nint, Nrank, Mrank - 1)
  call write_DSCS (Nteta, ExtThetaDom, h, v)
  call write_Effic (Qscat, Qext)
  call write_MrankConvRes (NthetaConvM, epsMrank)
  if (NthetaConvN >= int(0.8*Nteta) .and. NthetaConvM >= int(0.8*Nteta)) then
    print "(/,2x,'Convergence criteria for Nrank and Mrank are satisfied;')"                                      
  else
    print "(/,2x,'Convergence criteria for Nrank and Mrank are not satisfied;')"
  end if 
  call write_InfoFileTmat (FileTmatG, Mrank, Nrank, .false., .false., .false.)
  call ScatCharact (ks, FileTmatG, Mrank, Nrank, .false., .false., .false.) 
  print "(/,2x,'T matrix is stored in ',a50)", FileTmatG
  print "(  2x,'The dimensions of the T matrix are given by:')" 
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", Nrank    
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank  
  deallocate (paramG, weightsG, Nintparam) 
  deallocate (a, b, d, c, c1, t1, s, q11, q13, q31, q33, h, v, oldh, oldv,          &
              oldh0, oldv0)     
end subroutine convergence_Nrank_MrankINHOM
!***********************************************************************************
subroutine TMatrix_Nrank_MrankINHOM (TypeGeom, ks, ind_ref, Nsurf, surf, Nparam,    &
           Mrank, Nrank, Nint, miror, x1, y1, z1, alfa1, beta1, gama1, Mrank1,      &
           Nrank1, FileTmat, axsym1, chiral1, FileTmatG, PrnProgress)
  use parameters
  implicit none
  integer       :: TypeGeom, Nsurf, Nparam, Mrank, Nrank, Nint, Mrank1, Nrank1
  real(O)       :: ks, ind_ref, surf(Nsurf), x1, y1, z1, alfa1, beta1, gama1
  character(80) :: FileTmat, FileTmatG
  logical       :: miror, axsym1, chiral1, PrnProgress
!
  integer       :: Nmax, Nmax1, NmaxAL, nt1g, mt1g
  real(O)       :: k, r1, sumang
  complex(O)    :: ind_refC
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG(:,:), weightsG(:,:)
  complex(O),allocatable :: a(:,:), b(:,:), t1(:,:), s(:,:), q11(:,:), q31(:,:)     
!  
  Nmax   = Nrank + Mrank * (2 * Nrank - Mrank + 1) 
  Nmax1  = Nrank1 + Mrank1 * (2 * Nrank1 - Mrank1 + 1) 
  NmaxAL = max(Nmax,Nmax1)
  k = ks * ind_ref
  ind_refC = cmplx(ind_ref,0.0,O)                 
  allocate (a(2*NmaxAL,2*NmaxAL), b(2*NmaxAL,2*NmaxAL))
  allocate (s(2*NmaxAL,2*NmaxAL), q11(2*NmaxAL,2*NmaxAL), q31(2*NmaxAL,2*NmaxAL))                        
  allocate (paramG(Nparam,Nint), weightsG(Nparam,Nint), Nintparam(Nparam))
  call interpolation_listAXSYM (TypeGeom, Nsurf, surf, Nint, Nparam,                &
       Nintparam, paramG, weightsG, miror)         
  open (unit = iTmat, file = FileTmat, status = "old", position = "rewind")
  if (.not. axsym1) then
    call read_HeadFileTmat (nt1g, mt1g)
    call check_dimensionMat (nt1g, mt1g, Nmax1)       
    allocate (t1(2*nt1g,2*mt1g))
    call read_FileTmat (nt1g, mt1g, t1) 
  else
    nt1g = Nmax1
    mt1g = Nmax1
    allocate (t1(2*nt1g,2*mt1g))
    call read_Tmatrix (chiral1, Nrank1, Mrank1, Nmax1, t1, nt1g, mt1g)            
  end if
  close (unit = iTmat)
  open  (unit = iTmat, file = FileTmatG, status = 'replace')
  call  write_HeadFileTmat (NmaxAL, NmaxAL) 
  if (PrnProgress) call  write_progress (.true., 1, 6)
  r1 = sqrt(x1**2 + y1**2 + z1**2)
  sumang = abs(alfa1) + abs(beta1) + abs(gama1)
  if (r1 /= 0 .and. sumang /= 0) then
    call MatTransRot_TR_mn_m1n1 (1, k, x1, y1, z1, alfa1, beta1, gama1, Mrank,      &
         Nrank, Nmax, Mrank1, Nrank1, Nmax1, s, NmaxAL, NmaxAL)
  else if (r1 /= 0 .and. sumang == 0) then 
    call MatTransAB_mn_m1n1 (1, k, x1, y1, z1, Mrank, Nrank, Nmax, Mrank1, Nrank1,  &
         Nmax1, s, NmaxAL, NmaxAL)
  else if (r1 == 0 .and. sumang /= 0) then    
    call MatRot_R_mn_m1n1 (alfa1, beta1, gama1, Mrank, Nrank, Nmax, Mrank1, Nrank1, &
         Nmax1, s, NmaxAL, NmaxAL, .true.)
  else if (r1 == 0 .and. sumang == 0) then
    call identity_transformation (Mrank, Nrank, Nmax, Mrank1, Nrank1, Nmax1, s,     &
         NmaxAL, NmaxAL)    
  end if  
  call product_matrices (2*Nmax, 2*Nmax1, 2*Nmax1, s, 2*NmaxAL, 2*NmaxAL, t1,       &
       2*nt1g, 2*mt1g)
  if (PrnProgress) call write_progress (.false., 2, 6)
  if (r1 /=0 .and. sumang /= 0) then
    call MatRotTrans_RT_mn_m1n1 (1, k, x1, y1, z1, alfa1, beta1, gama1, Mrank1,     &
         Nrank1, Nmax1, Mrank, Nrank, Nmax, a, NmaxAL, NmaxAL)
  else if (r1 /=0 .and. sumang == 0) then 
    call MatTransAB_mn_m1n1 (1, k, -x1, -y1, -z1, Mrank1, Nrank1, Nmax1, Mrank,     &
         Nrank, Nmax, a, NmaxAL, NmaxAL)
  else if (r1 == 0 .and. sumang /= 0) then    
    call MatRot_R_mn_m1n1 (-gama1, -beta1, -alfa1, Mrank1, Nrank1, Nmax1, Mrank,    &
         Nrank, Nmax, a, NmaxAL, NmaxAL, .true.)
  else if (r1 == 0 .and. sumang == 0) then
    call identity_transformation (Mrank1, Nrank1, Nmax1, Mrank, Nrank, Nmax, a,     &
         NmaxAL, NmaxAL)    
  end if   
  call product_matrices (2*Nmax, 2*Nmax1, 2*Nmax, s, 2*NmaxAL, 2*NmaxAL, a,         &
       2*NmaxAL, 2*NmaxAL) 
  if (PrnProgress) call write_progress (.false., 3, 6)          
  call matrix_Q_sym (TypeGeom, 3, 3, ks, ind_refC, Nsurf, surf, Mrank, Nrank, Nmax, &
       Nint, Nparam, Nintparam, paramG, weightsG, miror,.false., a, NmaxAL, NmaxAL)
  call product_matrices (2*Nmax, 2*Nmax, 2*Nmax, a, 2*NmaxAL, 2*NmaxAL, s,          &
       2*NmaxAL, 2*NmaxAL)
  call matrix_Q_sym (TypeGeom, 3, 1, ks, ind_refC, Nsurf, surf, Mrank, Nrank, Nmax, &
       Nint, Nparam, Nintparam, paramG, weightsG, miror,.false., q31, NmaxAL, NmaxAL)
  call sum_matrices (2*Nmax, 2*Nmax, a, 2*NmaxAL, 2*NmaxAL, q31, 2*NmaxAL, 2*NmaxAL)  
  if (PrnProgress) call write_progress (.false., 4, 6)
  call matrix_Q_sym (TypeGeom, 1, 3, ks, ind_refC, Nsurf, surf, Mrank, Nrank, Nmax, &
       Nint, Nparam, Nintparam, paramG, weightsG, miror,.false., b, NmaxAL, NmaxAL) 
  call product_matrices (2*Nmax, 2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL, s,          &
       2*NmaxAL, 2*NmaxAL)  
  call matrix_Q_sym (TypeGeom, 1, 1, ks, ind_refC, Nsurf, surf, Mrank, Nrank, Nmax, &
       Nint, Nparam, Nintparam, paramG, weightsG, miror,.false., q11, NmaxAL, NmaxAL) 
  call sum_matrices (2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL, q11, 2*NmaxAL, 2*NmaxAL)  
  if (PrnProgress) call write_progress (.false., 5, 6)
  call LU_SYSTEM (a, 2*NmaxAL, 2*NmaxAL, b, 2*NmaxAL, 2*NmaxAL, 2*Nmax)
  if (PrnProgress) call write_progress (.false., 6, 6)
  call minus_matrix (2*Nmax, 2*Nmax, b, 2*NmaxAL, 2*NmaxAL)
  call write_FileTmat (NmaxAL, NmaxAL, b)         
  close (unit = iTmat)  
  call write_InfoFileTmat (FileTmatG, Mrank, Nrank, .false., .false., .false.)
  call ScatCharact (ks, FileTmatG, Mrank, Nrank, .false., .false., .false.) 
  print "(/,2x,'T matrix is stored in ',a50)", FileTmatG
  print "(  2x,'The dimensions of the T matrix are given by:')" 
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", Nrank    
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank  
  deallocate (paramG, weightsG, Nintparam) 
  deallocate (a, b, t1, s, q11, q31)    
end subroutine TMatrix_Nrank_MrankINHOM
