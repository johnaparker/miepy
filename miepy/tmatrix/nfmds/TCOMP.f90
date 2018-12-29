subroutine TCOMP
!------------------------------------------------------------------------------------
! 1. General Considerations                                                         !
! --------------------------                                                        !
! TCOMP is a routine for computing the T matrix and the scattering characteristics  !
! of axisymmetric composite particles. A composite particle consists of several     !
! nonenclosing regions (parts), each characterized by arbitrary but constant        !
! values of electric permittivity and magnetic permeability.                        !
!                                                                                   !
! A composite particle is divided into homogeneous regions. The electric and        !
! magnetic fields of a specific region i can be approximated by using localized     !
! or distributed vector spherical wave functions, and Nrankp(i) will be referred    !
! to as the maximum expansion order of the region i.                                !
!                                                                                   !
! The transition matrix of a composite particle can be obtained by considering the  !
! null-field equation for the total electric field inside each region. Note that    !
! for axisymmetric composite particles, the scattering problem decouples over the   !
! azimuthal modes m and the T matrix can be computed separately for each m.         !
!                                                                                   !
! Particle geometries currently supported include: half-spheroids with offset       ! 
! origins and three cylinders. The following routines (from the file "GeomLib.f90") !
! provide the required geometry parameters:                                         !
! - the routine "interpolation_listCOMP" computes the integration nodes and         !
!   weights for each generatrix curve,                                              !
! - the routine "elem_geomCOMP" provides the geometry parameters (position vector,  !
!   polar and azimuthal angles and normal unit vector in spherical coordinates)     !
!   at each integration node, and                                                   !
! - the routine "zDSCOMP" generates the coordinates of the distributed sources.     !
! The generatrix of each homogeneous region is described with respect to a local    !
! coordinate system and the axial positions of the local coordinate systems are     !
! specified in the global coordinate system of the composite particle. By           !
! convention, the COORDINATES OF THE DISTRIBUTED SOURCES corresponding to a         !
! specific region are provided in the GLOBAL COORDINATE SYSTEM. The user can easily !
! modify the above routines to generate particles with other geometries. Note that  !
! the list of parameters must be maintained and only geometries with an analytical  !
! description of the surface can be implemented.                                    !
!                                                                                   !
! 2. Convergence Test                                                               !
! --------------------                                                              !
! Convergence tests precede the scattering characteristics calculation. The         !
! following parameters control the T-matrix computation:                            !
! - the number of integration points Nint,                                          !
! - the maximum expansion order for each region Nrankp,                             !
! - the maximum expansion order for the composite particle Nrank, and               !
! - the maximum azimuthal order Mrank.                                              !
!                                                                                   !
! Before going any further, we will explain the significance of the parameters      !
! controlling the convergence process. Let Nrankp(i) denotes the maximum expansion  !
! order of the region i, i = 1,2,...,Npart, where Npart is the number of            !
! nonenclosing regions. Nrank determines the dimension of the T matrix (centered    !
! at the origin of the global coordinate system) and, by convention, Nrank will be  !
! referred to as the maximum expansion order of the composite particle.             !
! Analogously, Nint will be referred to as the global number of integration points  !
! for the composite particle and its significance will be described below. The      !
! generatrix curve i (corresponding to the homogeneous region i) consists of        !
! NparamPart(i) piecewise smooth curves. On each smooth curve j, j = 1,2,...,       !
! NparamPart(i), the number of integration points is                                !
!                                                                                   !
!               Nint_region_curve(i,j) = int { c(i,j) * Nint },                     !
!                                                                                   !
! where c(i,j) is a multiplicative factor. The parameters c(i,j), do not appear     !
! as explicit parameters; they are hard coded in the routine                        !
! "interpolation_listCOMP" from the file "GeomLib.f90". For half-spheroids, Nint    !
! is the number of integration points corresponding to each half-spheroid. The       !
! generatrix of a half-spheroid is divided into a quarter-ellipse and a straight    !
! line. Denoting by a the semi-axis of the spheroid along the axis of symmetry and  !
! by b the second semi-axis, we have                                                !
!                                                                                   !
!    Nint_region_curve(half-spheroid,quarter-ellipse) = int {Le * Nint / (Le + b)}  !
!                                                                                   !
! and                                                                               !
!                                                                                   !  
!    Nint_region_curve(half-spheroid,line) = Nint - int {Le * Nint / (Le + b)},     !
!                                                                                   !
! where                                                                             !
!                                                                                   !
!             Le = Pi * { 3(a + b) - sqrt[(a + 3b) * (3a + b)] } / 4                !
!                                                                                   !
! is an approximation of the length of the quarter-ellipse.                         !  
! For three cylinders, Nint is the number of integration points corresponding to    !
! each cylinder. The generatrix of a cylinder is divided into three straight lines. ! 
! The first and third lines are perpendicular to the axis of symmetry, while the    !
! second line is parallel to the axis of symmetry. Denoting by a the half-length    !
! of the cylinder and by b the cylinder radius, we have:                            !
!                                                                                   !
!               Nint_region_curve(cylinder,first line)                              !
!      = Nint_region_curve(cylinder,third line) = int {[b/2/(a + b)] * Nint}        !
!                                                                                   !
! and                                                                               !
!                                                                                   !
!                 Nint_region_curve(cylinder, second line)                          ! 
!                 = Nint - 2 * int {[b/2/(a + b)] * Nint }.                         !
!                                                                                   !
! The above selection rules choose the same number of integration points (Nint)     !
! for each homogeneous regions. Therefore, a conservative parameter choice rule is  !
!                                                                                   !  
!                      Nint = max {Nint_region},                                    !
!                                                                                   !
! where Nint_region(i) is the required number of integration points for the         !
! homogeneous region i.                                                             !
!                                                                                   !
! For convergence tests, the incident wave is assumed to be a plane wave traveling  !
! along the Z-axis of the global coordinate system and the scattering               !
! characteristics are computed in the azimuthal plane phi = 0°.                     !
!                                                                                   !
! For the integration and expansion order test, an axisymmetric orientation of the  !
! particle is considered, i.e., the incident wave propagates along the axis of      !
! symmetry of the particle. In this case, all azimuthal modes are zero except       !
! for m = - 1 and m = 1. For the convergence test over the number of integration    !
! points, the scattering problem is solved for Nint and Nint + dNint, while for     !
! the convergence test over the expansion orders, the scattering problem is         !
! solved for                                                                        !
!                                                                                   !
!                 Nrankp(i), i = 1,2,...,Npart, and Nrank,                          !
!                                                                                   !
! and                                                                               !
!                 Nrankp(i) - 1, i = 1,2,...,Npart, and Nrank - 1.                  !
!                                                                                   !
! By convention, the convergence test over the expansion orders will be referred    !
! to as the convergence test over Nrank. The normalized differential scattering     !
! cross section (DSCS) will be checked at 20° increments for convergence within     !
! epsX (epsNint or epsNrank) tolerance. If the calculated results converge within   !
! this tolerance at 80% of the scattering angles, then convergence is achieved.     !
!                                                                                   ! 
! When the convergence tests over Nrank and Nint are finished, we pass to the       !
! azimuthal order test. The program automatically sets the particle to a more       !
! general orientation, i.e., alpha = beta = 45°, and solves the scattering problem  !
! for increasing m values until convergence of the angular scattering is achieved.  !
! The T matrix is stored for later use by other programs, and the values of Nrank   !
! and Mrank are printed to the screen and to the T-matrix information file.  These  !
! values together with the T matrix serve as INPUT PARAMETERS for other programs.   !
!                                                                                   !
! 3. Estimates of Nint and Nrank                                                    !
! -------------------------------                                                   !
! The above convergence tests require estimates for Nint, Nrankp(i),                ! 
! i = 1,2,...,Npart, and Nrank. These estimates are user-defined.                   !
!                                                                                   !
! Reasonable starting values of Nrankp and Nrank are given by Wiscombe's truncation !
! limit criterion [W. J. Wiscombe, Improved Mie scattering algorithms, Applied      !
! Optics, 19, 1505-1509, 1980]. For Nrankp, we have                                 !
!                                                                                   !                     
!                NrankpW = int(xp + 4.05 * xp**0.33 + 2),                           !
!                                                                                   ! 
! where xp is the size parameter of the region, xp = k * lnormPart, k is the wave   !
! number and lnormPart is the radius of the smallest sphere circumscribing the      !
! region. The estimate of Nrank is                                                  !
!                                                                                   !
!                NrankW = int(x + 4.05 * x**0.33 + 2),                              !
!                                                                                   !
! where x is the size parameter of the composite particle, x = k * Rcirc, and       !
! Rcirc is the radius of the smallest sphere circumscribing the composite particle. !
!                                                                                   !
! The value of Nint depends on the type of discrete sources, the size parameter,    !
! and the optical properties of the composite particle. Denoting by Nrankpmax the   !
! maximum value of Nrankp(i), i = 1,2,...,Npart, i.e.,                              !
!                                                                                   !
!                        Nrankpmax = max {Nrankp(i)},                               !
!                                                                                   !
! we propose the following estimate                                                 !
!                                                                                   ! 
!                        Nint = Ndgs * Nrankpmax,                                   !
!                                                                                   !
! where Ndgs = 15, 20, 25, .... For distributed sources, the value of Nint should   !
! be substantially larger than the value corresponding to localized sources.        !
!                                                                                   !
! 4. Strategies for Performing Convergence Tests                                    !
! -----------------------------------------------                                   !
! Two strategies for performing convergence tests over Nint and Nrank are           !
! recommended. The following technique is an extension of the approach used for     !
! homogeneous axisymmetric particles:                                               !
! 1. for each region i, choose a value of Nrankp(i) close to the value predicted    !
!    by Wiscombe's truncation limit criterion and proceed similarly for Nrank;      !
! 2. for the set Nrankp(i) and Nrank, perform the expansion order test with         !
!    Nint = Ndgs * Nrankpmax, where e.g., Ndgs = 15;                                !
! 3. if convergence is not achieved, set Nrankp(i) = Nrankp(i) + 1 for all regions  !
!    i, set Nrank = Nrank + 1, maintain the relation Nint = Ndgs * Nrankpmax and    !
!    go to Step 2;                                                                  !
! 4. if convergence is achieved, perform several integration tests with the         ! 
!    initial Nint-value, Nint = Ndgs * Nrankpmax, until the DSCS converges.         !
! At Step 3, it is not necessary to increase all Nrankp-values (for all regions).   !
! If the user has the feeling that some Nrankp-values are correct, these values     !
! need not to be changed.                                                           !  
!                                                                                   !
! For large and/or highly aspherical composite particles, T-matrix computations     !
! become poorly covergent or divergent. Frequently, a semi-convergent behaviour is  !
! attended: the errors decrease with increasing the maximum expansion order,        !
! attain a relative constant level and afterwards increase. The region of           !
! stability can be localized by using the following technique:                      !
! 1. for each region i, choose a value of Nrankp(i) smaller than the value          !
!    predicted by Wiscombe's truncation limit criterion and proceed similarly for   !
!    Nrank;                                                                         !
! 2. for the set Nrankp(i) and Nrank, perform several integration tests with the    !
!    initial Nint-value, Nint = Ndgs * Nrankpmax, until the DSCS converges;         !
! 3. if convergence is not achieved, set Nrankp(i) = Nrankp(i) - 1,                 !
!    i = 1,2,...,Npart, and Nrank = Nrank - 1, or increase the tolerance epsNint,   !
!    and go to Step 2,                                                              !
! 4. if convergence is achieved (for some Nint), perform several expansion order    !
!    tests for increasing Nrankp- and Nrank-values and monitor the errors of the    !
!    differential scattering cross sections.                                        ! 
! As before, at Steps 3 and 4, not all Nrankp-values need to be changed.            !
!                                                                                   !
! For particles which are too extreme in size and/or aspect ratio, the errors of    !
! the extinction and scattering cross sections (instead of the differential         !
! scattering cross section) can be analyzed.                                        ! 
!                                                                                   !
! 5. Additional Comments                                                            !
! -----------------------                                                           !
! For the expansion order test and distributed sources, two configurations of       !
! sources are considered. For both configurations, the coordinates of the           !
! distributed sources can be generated automatically or can be specified in the     !
! input file.                                                                       !
!                                                                                   !
! The convergence tests over Nint and Nrank can be switched off by setting the      !
! logical variable DoConvTest to false. In this case, the values of Nint,           !
! Nrankp(i), i = 1,2,...,Npart, and Nrank must be specified in the input file.      !
!                                                                                   !
! Nint and Nrank have been determined for the azimuthal orders m = 1 and m = - 1.   ! 
! For higher azimuthal orders and sources distributed in the complex plane, the     !
! chosen values of Nint and Nrank may lead to false convergence. There are two      !
! strategies for improving the numerical stability.                                 !
! 1. The number of integration points can be increased at each azimuthal mode       !
! calculation with the integration step dNintMrank. The parameter dNintMrank        !
! should be found experimentally by repeating the azimuthal order test for          !
! different input values until the differential scattering cross section converges. !
! A conservative estimate is: dNintMrank = (0.1,...,0.2) * Nint.                    !
! 2. For sources generated automatically, an interactive convergence test over      !
! Nint and Nrank can be performed at each azimuthal mode calculation. The general   !
! strategy is to reduce the number of discrete sources (Nrankp(i), i = 1,2,...,     !
! Npart) and to increase the number of integration points (Nint) for increasing     !
! values of m.  However, in the current version of the program, this test is not    !
! active (see the comments in the subroutine "convergence_MrankDSCOMP").            !
!                                                                                   !
! The number of azimuthal modes is the same for all regions, and for localized      !
! sources, the relation Mrank <= Nrankp(i), for all i = 1,2,...,Npart, should be    !
! satisfied. However, Nrankp(i), i = 1,2,...,Npart, are determined before Mrank     !
! is computed. Therefore, for the azimuthal mode calculation m, with |m| > 1,       !
! we set Nrankp(i) = |m|, if Nrankp(i) < |m| for some i. In this case, the maximum  !
! expansion order of a region is m-dependent and                                    !
!                                                                                   !
!                    n = |m|,|m|+1,...,Nrankp(i,|m|).                               !
!                                                                                   !
! 6. Input Parameters                                                               !
! ---------------------                                                             !
! The parameters specified in the input file "/INPUTFILES/InputCOMP.dat" are listed !
! below                                                                             !
!                                                                                   !
! - wavelength (real) - wavelength of the incident light in vacuo.                  !
!                                                                                   !
! - ind_refMed (real) - refractive index of the ambient medium (nonabsorbing        !
!   dielectric medium).                                                             !
!                                                                                   !
! - TypeGeom (integer) - parameter specifying the type of the particle geometry.    !
!                                                                                   !
! - Npart (integer) - number of homogeneous regions.                                !
!                                                                                   !
! - anorm (real) - characteristic length of the composite particle which is used    !
!   to normalize the differential scattering cross sections.                        !
!                                                                                   !
! - Rcirc (real) - characteristic length of the composite particle (usually the     !
!   radius of the smallest circumscribing sphere) which is used to compute an       !
!   estimate of the maximum expansion order by using Wiscombe's truncation limit    !
!   criterion. Alternatively, Rcirc can be chosen as the equal-volume sphere radius !
!   or the surface-equivalent-sphere radius.                                        !
!                                                                                   !
! - DoConvTest (logical) - if DoConvTest = t, the interactive convergence tests     !
!   over Nint and Nrank are invoked. Estimates of Nrank for all homogeneous         !
!   regions are provided by Wiscombe's truncation limit criterion. If               !
!   DoConvTest = f, the values of Nint and Nrank must be supplied in the input      !
!   file.                                                                           !
!                                                                                   !
! - DS (logical) - if DS = t, distributed sources are used for T-matrix             !
!   calculation, otherwise localized sources are employed.                          !
!                                                                                   !
! - autGenDS (logical) - if autGenDS = t, the routine "zDSCOMP" generates the       !
!   coordinates of the distributed sources. Otherwise, the coordinates of the       !
!   distributed sources must be provided by the user. If DS = f, the logical        !
!   variable autGenDS is ignored.                                                   !
!                                                                                   !
! - Nint (integer) - global number of integration points for the composite          !
!   particle. This parameter is used if the convergence tests are not performed     !
!   (DoConvTest = f).                                                               !
!                                                                                   !
! - Nrank (integer) - maximum expansion order for the composite particle. This      !
!   parameter is used if the convergence tests are not performed (DoConvTest = f).  !
!                                                                                   !
! The next parameters (specified in a sequence of group statements) characterize     !
! each homogeneous region.                                                          !
! - ind_refPartRel (complex) - relative refractive index of the actual region       !
!   with respect to the ambient medium. The imaginary part of the relative          !
!   refractive index must be zero for nonabsorbing particles and positive for       !
!   absorbing particles.                                                            !
!                                                                                   !
! - NsurfPart (integer) - number of surface parameters of the actual region.        !
!                                                                                   !
! - NparamPart (integer) - number of smooth curves forming the generatrix curve of  !
!   the actual region.                                                              !
!                                                                                   !
! - OriginPart (real) - axial position of the local coordinate system of the        !
!   actual region with respect to the global coordinate system of the composite     !
!   particle. This parameter can be positive or negative.                           !
!                                                                                   !
! - surfPart (real array: surfPart(1), surfPart(2),...,surfPart(NsurfPart)) -       !
!   surface parameters specifying the shape of the actual region. The dimension of  !
!   the array is NsurfPD. The integer parameter NsurfPD is specified in the         !
!   routine  "Parameters.f90" and has the value NsurfPD = 10.                       !
!   If NsurfPart > NsurfPD, the execution is automatically terminated.              !
!                                                                                   !
!   The permissive values of the surface parameters (for the supported geometries)  !
!   are listed below.                                                               !
!                                                                                   !
!      Particle      TypeGeom  NsurfPart  NparamPart          surfPart              !      
!   half-spheroids                                                                  !
!    with offset         1         3         2         surfPart(1) - length of the  !
!      origins                                             semi-axis along the      !
!                                                          symmetry axis            !    
!                                                      surfPart(2) - length of the  !
!                                                          second semi-axis         !
!                                                      surfPart(3) - distance       !
!                                                          between the half-        ! 
!                                                          spheroid basis and the   !
!                                                          origin of the local      !
!                                                          coordinate system        ! 
!                                                                                   !   
!    3 cylinders         2         2         3         surfPart(1) - half-length    !
!                                                          of the cylinder          !  
!                                                      surfPart(2) - cylinder       !
!                                                          radius                   !
!                                                                                   !
!   For half-spheroids, surf(3) specifies the position of the local coordinate      !
!   system, and this parameter has been introduced for a general characterization   !
!   of the surface. It is recommended to set surfPart(3) = 0.5 * surfPart(1). For   !
!   touching half-spheroids, surf(3) = |OriginPart|, and THIS SPECIFIC SITUATION    !
!   IS CONSIDERED IN THE CODE. Note that the program DOES NOT CHECK if this         !
!   condition is satisfied.                                                         !
!                                                                                   !
! - lnormPart (real) - characteristic length of the actual region (usually the      !
!   radius of the smallest circumscribing sphere) which is used to compute an       !
!   estimate of the maximum expansion order by using Wiscombe's truncation limit    !
!   criterion. Alternatively, lnormPart can be chosen as the equal-volume sphere    !
!   radius or the surface-equivalent-sphere radius.                                 ! 
!                                                                                   !
! - ComplexPlanePart (logical) - if ComplexPlanePart = t, the distributed sources   !
!   are situated in the complex plane. This parameter is used if distributed sources!
!   are required and the coordinates of the distributed sources are automatically   !
!   generated (DS = t and autGenDS = t).                                            !
!                                                                                   !
! - EpsZReImPart (real) - input parameter of the routine "zDSCOMP" which controls  !
!   the distribution of the discrete sources. This parameter is used if distributed !
!   sources are required and the coordinates of the distributed sources are         !
!   automatically generated (DS = t and autGenDS = t).                              !
!                                                                                   !
! - NrankPart (integer) - maximum expansion order for the actual region. This       !
!   parameter is used if the convergence tests are not performed or the coordinates !
!   of the distributed sources are user-defined. More specifically, NrankPart is    !
!   used if                                                                         !
!   - (DS = f, DoConvTest = f),                                                     !
!   - (DS = t, autGenDS = t, DoConvTest = f), or                                    ! 
!   - (DS = t, autGenDS = f, DoConvTest = t/f).                                     !                                                                  !
!                                                                                   !
! - zRePart, zImPart (real arrays: zXPart(1), zXPart(2),...,zXPart(Nrankp),         !
!   X = Re, Im) - coordinates of the distributed sources for the actual region and  !
!   the expansion order NrankPart. These parameters are used if the coordinates of  !
!   the distributed sources are user-defined (DS = t and autGenDS = f). The         !
!   dimension of the arrays zRePart and zImPart is NrankPD and the inequality       !
!   NrankPart <= NrankPD must hold. The integer parameter NrankPD is specified in   !
!   the routine "Parameters.f90" and has the value NrankPD = 200. If                !
!   NrankPart > NrankPD, the execution is automatically terminated. Note that the   !
!   coordinates of the distributed sources are defined with respect to the GLOBAL   !
!   COORDINATE SYSTEM of the composite particle.                                    !
!                                                                                   !
! - zRePart1, zImPart1 (real arrays: zXPart1(1), zXPart1(2),...,zXPart1(Nrankp-1),  !
!   X = Re, Im) - coordinates of the distributed sources for the actual region and  !
!   the expansion order NrankPart - 1. These parameters are used if the             !
!   coordinates of the distributed sources are user-defined and the expansion order !
!   test is performed (DS = t and autGenDS = f and TypeConvTest = 2). The dimension !
!   of the arrays zRePart1 and zImPart1 is NrankPD. As before, zRePart1 and         !
!   zImPart1 are defined with respect to the GLOBAL COORDINATE SYSTEM of the        !
!   composite particle.                                                             !
!                                                                                   !
! NOTE: THE INPUT ARRAYS zRe, zIm AND zRe1, zIm1 MUST BE SPECIFIED IF (DS = t AND   !
! autGenDS = f), AND MUST BE SEPARATED BY A BLANK LINE. IF THE EXPANSION ORDER TEST !
! IS NOT PERFORMED (TypeConvTest /= 2), THE INPUT ARRAYS zRePart1 AND zImPart1      !
! CAN BE OMITTED.                                                                   !
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
! - dNintMrank (integer) - number of division points for azimuthal mode             !
!   calculation. This parameter is used if distributed sources are employed         !
!   (DS = t) and THERE EXISTS AT LEST ONE LAYER with ComplexPlanePart = t. At       !
!   each azimuthal mode calculation, the number of integration points increases     !
!   with dNintMrank. For regions with sources distributed along the symmetry axis,  !
!   the number of integration points is also increased.                             !
!                                                                                   !
! - FileTmat (character(80)) - name of the file to which the T matrix is written.   !
!                                                                                   !
! - PrnProgress (logical) - if PrnProgress = t, the progress of calculation         !
!   is printed.                                                                     !
!                                                                                   !
! Note that all input parameters specifying lengths must be provided in the same    !
! units of lengths.                                                                 !
!                                                                                   !
! The statements with the group names OptRegProp, GeomRegProp, SourceRegPosAut,     !
! NrankReg and SourceRegPosInp characterize a specific region. THESE STATEMENTS     !
! MUST BE REPEATED FOR ALL Npart REGIONS.                                           !
!                                                                                   !
! 7. Logical Scheme                                                                 !
! ------------------                                                                !  
! The organization of the code is as follows:                                       !
! 1. < read input data >                                                            !
! 2. < T-matrix calculation >                                                       !
!                                                                                   !
!    if ( use localized sources ) then                                              !
!                                                                                   !
!       if ( do convergence test) then                                              !
!          < the code computes an estimate of Nrankp(i) for each region i           !
!            by using Wiscombe's truncation criterion, and prompts for the          !
!            input value of Nrankp(i) >                                             !
!          < the code computes an estimate of Nrank for the composite particle      !
!            by using Wiscombe's truncation criterion >                             !
!          < the code prompts for the estimated values of Nint and Nrank for the    !
!            composite particle >                                                   !
!          < the code prompts for the type of convergence test:                     !
!            1 - over Nint, 2 - over Nrank or 3 - over Mrank >                      !
!       else if ( .not. do convergence test) then                                   !
!          < the input file provides the values of Nrankp(i) for each region i >    !
!          < the input file provides the values of Nint and Nrank for the           !
!            composite particle>                                                    !
!          type of convergence test = 3 (convergence test over Mrank)               !     
!       end if                                                                      !
!                                                                                   !
!    else if ( use distributed sources ) then                                       !
!                                                                                   !
!       if ( use automatic sources generation ) then                                !
!          if ( do convergence test) then                                           !
!             < the code computes an estimate of Nrankp(i) for each region i        !
!               by using Wiscombe's truncation criterion, and prompts for the       !
!               input value of Nrankp(i) >                                          !
!             < the code generates the coordinates of the distributed sources for   !
!               Nrankp(i) and all regions i >                                       !
!             < the code computes an estimate of Nrank for the composite particle   !
!               by using Wiscombe's truncation criterion >                          !
!             < the code prompts for the estimated values of Nint and Nrank for the !
!               composite particle>                                                 !
!             < the code prompts for the type of convergence test:                  !
!               1 - over Nint, 2 - over Nrank or 3 - over Mrank >                   !
!             if ( do expansion order test ) then                                   !
!                < the code generates the coordinates of the distributed sources    !
!                  for Nrankp(i) - 1 and all regions i >                            !
!             end if                                                                !
!          else if ( .not. do convergence test) then                                !
!             < the input file provides the values of Nrankp(i) for each region i > !
!             < the code generates the coordinates of the distributed sources for   !
!               Nrankp(i) and all regions i >                                       !
!             < the input file provides the values of Nint and Nrank for the        !
!               composite particle >                                                !
!             type of convergence test = 3 (convergence test over Mrank)            !
!       else if ( .not. use automatic sources generation ) then                     !
!          if ( do convergence test) then                                           !
!             < the input file provides the values of Nrankp(i) for each region i > !
!             < the input file provides the coordinates of the distributed sources  !
!               for Nrankp(i) and all regions i >                                   !
!             < the code computes an estimate of Nrank for the composite particle   !
!               by using Wiscombe's truncation criterion >                          !
!             < the code prompts for the estimated values of Nint and Nrank for the !
!               composite particle >                                                !
!             < the code prompts for the type of convergence test:                  !
!               1 - over Nint, 2 - over Nrank or 3 - over Mrank >                   !
!             if ( do maximum expansion order test ) then                           !
!                < the input file provides the coordinates of the distributed       !
!                  sources for Nrankp(i) - 1 and all regions i >                    !
!             end if                                                                !
!          else if ( .not. do convergence test) then                                !
!             < the input file provides the values of Nrankp(i) for each region i > !
!             < the input file provides the coordinates of the distributed sources  !
!               for Nrankp(i) and all regions i >                                   !
!             < the input file provides the values of Nint and Nrank for the        !
!               composite particle >                                                !
!             type of convergence test = 3 (convergence test over Mrank)            !  
!          end if                                                                   !
!       end if                                                                      !
!                                                                                   !
!    end if                                                                         !
!    if ( do integration test, i.e., type of convergence test = 1 ) then            !
!       < the code computes the DSCS for Nint and Nint + dNint and                  !
!         write the results to the file "Output.dat" >                              !
!       < the execution is terminated >                                             !
!    end if                                                                         !
!    if ( do expansion order test, i.e., type of convergence test = 2 ) then        !
!       < the code computes the DSCS for Nrankp(i), Nrank and Nrankp(i)-1, Nrank-1  ! 
!         and write the results to the file "Output.dat" >                          !
!       < the execution is terminated >                                             !
!    end if                                                                         !
!    if ( do azimuthal order test, i.e., type of convergence test = 3 ) then        !
!       < the code computes the DSCS for increasing m values and                    !
!         write the results to the file "Output.dat">                               !
!       if ( use distributed sources .and. use automatic sources generation         !
!           .and. sources are distributed in the complex plane ) then               ! 
!          < an interactive convergence test over Nint and Nrank can be performed   !
!            at each azimuthal mode m  >                                            ! 
!        end if                                                                     !
!        Note: to make this test active comment out some lines in the subroutine    !
!        "convergence_MrankDSCOMP".                                                 ! 
!       < the T matrix is stored in the file FileTmat and the T-matrix information  !
!         file is created >                                                         !
!       < Mrank and Nrank are printed to the screen and to the T-matrix             !
!         information file >                                                        !
!       < go to Step 3 >                                                            !
!    end if                                                                         !
! 3. < scattering characteristics calculation >                                     !
!------------------------------------------------------------------------------------
  use parameters  
  use allocation, only: Nsurf, Nparam, Nrankp, Nrankp1, surf, zpart, zRe, zIm,      &
                        zRe1, zIm1, lnorm, ind_ref, ComplexPlane, EpsZReIm    
  implicit none 
  integer       :: TypeGeom, Npart, dNint, TypeConvTest, Nsurfmax, Nparammax,       &
                   Nrankpmax, Nrankpmax1, Nint, Nrank, dNintMrank
  real(O)       :: k, ind_refMed, wavelength, anorm, snorm, epsNint, epsNrank,      &
                   epsMrank                                                 
  logical       :: DoConvTest, DS, autGenDS, PrnProgress 
  character(80) :: FileTmat  
! -----------------------------------------------------------------------------------
!                             Read the input file                                   ! 
! -----------------------------------------------------------------------------------  
  call readinputCOMP ( wavelength, ind_refMed, TypeGeom, Npart, anorm,              &
       DoConvTest, DS, autGenDS, Nint, Nrank, epsNint, epsNrank, epsMrank, dNint,   &
       dNintMrank, FileTmat, PrnProgress, k, snorm, Nsurfmax, Nparammax, Nrankpmax, &
       Nrankpmax1, TypeConvTest )                                
! -----------------------------------------------------------------------------------
!                                      Main                                         !
! ----------------------------------------------------------------------------------- 
  open (unit = iOutput, file = FileOutput, status = "replace")   
  call printinputCOMP (TypeConvTest, TypeGeom, Nsurfmax, Nsurf, surf, Nparam, Npart,&
       lnorm, anorm, Nrankpmax, Nrankp, zpart, zRe, zIm, Nrankpmax1, Nrankp1,       &
       zRe1, zIm1, ind_ref, dNint, wavelength, ind_refMed, epsNint, epsNrank,       &
       epsMrank, DS, autGenDS)
  if (DoConvTest) then         
    if (TypeConvTest == 1) then
      if (.not. DS) then
        call convergence_NintCOMP (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,     &
             Nparammax, Nparam, Npart, Nrankp, zpart, Nrank, Nint, dNint,           &
             epsNint, PrnProgress)
      else 
        call convergence_NintDSCOMP (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,   &
             Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zRe, zIm, zpart, Nrank,   &
             Nint, dNint, epsNint, PrnProgress)
      end if
    else if (TypeConvTest == 2) then
      if (.not. DS) then
        call convergence_NrankCOMP (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,    &
             Nparammax, Nparam, Npart, Nrankp, zpart, Nrank, Nint, epsNrank,        &
             PrnProgress)  
      else 
        call convergence_NrankDSCOMP (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,  &
             Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zRe, zIm, Nrankpmax1,     &
             Nrankp1, zRe1, zIm1, zpart, Nrank, Nint, epsNrank, PrnProgress)
      end if
    else 
      if (.not. DS) then
        call convergence_MrankCOMP (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,    &
             Nparammax, Nparam, Npart, Nrankp, zpart, Nrank, Nint, epsMrank,        &
             FileTmat, PrnProgress)
      else 
        call convergence_MrankDSCOMP (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,  &
             Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zRe, zIm, zpart, Nrank,   &
             Nint, ComplexPlane, EpsZReIm, autGenDS, epsMrank, epsNrank,            &
             dNintMrank, FileTmat, PrnProgress)
      end if
    end if   
  else    
    if (.not. DS) then
      call convergence_MrankCOMP (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,      &
           Nparammax, Nparam, Npart, Nrankp, zpart, Nrank, Nint, epsMrank,          &
           FileTmat, PrnProgress)
    else 
      call convergence_MrankDSCOMP (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,    &
           Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zRe, zIm, zpart, Nrank,     &
           Nint, ComplexPlane, EpsZReIm, autGenDS, epsMrank, epsNrank,              &
           dNintMrank, FileTmat, PrnProgress)
    end if
  end if  
  close (unit = iOutput)          
  deallocate (zRe, zIm, zRe1, zIm1) 
  deallocate (surf, Nsurf, Nparam, Nrankp, Nrankp1, ind_ref, zpart, lnorm)
  deallocate (ComplexPlane, EpsZReIm)  
end subroutine TCOMP
!***********************************************************************************
subroutine readinputCOMP ( wavelength, ind_refMed, TypeGeom, Npart, anorm,          &
           DoConvTest, DS, autGenDS, Nint, Nrank, epsNint, epsNrank, epsMrank,      &
           dNint, dNintMrank, FileTmat, PrnProgress, k, snorm, Nsurfmax, Nparammax, &
           Nrankpmax, Nrankpmax1, TypeConvTest ) 
  use parameters
  use derived_parameters
  use allocation, only: Nsurf, Nparam, Nrankp, Nrankp1, surf, zpart, zRe, zIm,      &
                        zRe1, zIm1, lnorm, EpsZReIm, ind_ref, ComplexPlane  
  implicit none 
  integer       :: TypeGeom, Npart, NsurfPart, NparamPart, dNint, TypeConvTest, i,  &
                   ipart, isurf, Nsurfmax, Nparammax, Nrankpmax, Nrankpmax1, Nint,  &
                   dNintMrank, NrankPart, Nrank, NrankW, ios
  real(O)       :: k, ind_refMed, wavelength, anorm, surfPart(NsurfPD), xpart,      &
                   snorm, epsNint, epsNrank, epsMrank, OriginPart, zRePart(NrankPD),&
                   x, zImPart(NrankPD), zRePart1(NrankPD), zImPart1(NrankPD),       &
                   lnormPart, EpsZReImPart, Rcirc                        
  complex(O)    :: ind_refPartRel
  logical       :: DoConvTest, DS, autGenDS, InputDS, ComplexPlanePart,             &
                   PrnProgress, XFindPar  
  character(80) :: FileTmat, string  
! -----------------------------------------------------------------------------------
!                          Read the input file FileInputCOMP                        ! 
! ----------------------------------------------------------------------------------- 
  call DrvParameters  
  open (unit = iInputCOMP, file = FileInputCOMP, status = "old", position = "rewind")   
  wavelength = 0.1_O * 2._O * Pi
  ind_refMed = 1._O
  string     = 'OptProp'
  if (XFindPar (iInputCOMP, string)) then
    read (iInputCOMP, *, iostat = ios) wavelength                
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable wavelength;')"
      stop
    end if
    read (iInputCOMP, *, iostat = ios) ind_refMed
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refMed;')"
      stop
    end if         
  else
    print "(/,2x,'Group name OptProp not found;')"
    stop  
  end if      
  k = 2._O * Pi * ind_refMed / wavelength
!
  TypeGeom = 1
  Npart    = 2
  anorm    = 1._O
  Rcirc    = 1._O
  string   = 'GeomProp'
  if (XFindPar (iInputCOMP, string)) then
    read (iInputCOMP, *, iostat = ios) TypeGeom
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable TypeGeom;')"
      stop
    end if
    read (iInputCOMP, *, iostat = ios) Npart
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Npart;')"
      stop
    end if
    read (iInputCOMP, *, iostat = ios) anorm
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable anorm;')"
      stop
    end if
    read (iInputCOMP, *, iostat = ios) Rcirc
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Rcirc;')"
      stop
    end if    
  else
    print "(/,2x,'Group name GeomProp not found;')"
    stop  
  end if 
  call check_anorm (anorm)  
  xpart = k * anorm
  snorm = Pi * xpart * xpart     
!
  DoConvTest = .true.
  string     = 'ConvTest'
  if (XFindPar (iInputCOMP, string)) then
    read (iInputCOMP, *, iostat = ios) DoConvTest
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable DoConvTest;')"
      stop
    end if          
  else
    print "(/,2x,'Group name ConvTest not found;')"
    stop  
  end if             
!
  DS = .false.
  autGenDS = .true. 
  string     = 'Sources'
  if (XFindPar (iInputCOMP, string)) then
    read (iInputCOMP, *, iostat = ios) DS
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable DS;')"
      stop
    end if
    read (iInputCOMP, *, iostat = ios) autGenDS
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable autGenDS;')"
      stop
    end if         
  else
    print "(/,2x,'Group name Sources not found;')"
    stop  
  end if    
  InputDS = .false.
  if (DS .and. .not. autGenDS) InputDS = .true. 
!
  allocate (Nsurf(Npart), Nparam(Npart), ind_ref(Npart), zpart(Npart), lnorm(Npart))
  Nsurfmax  = 1
  Nparammax = 1
  do ipart = 1, Npart
    ind_refPartRel = (1.5_O,0._O)
    string = 'OptRegProp'
    if (XFindPar (iInputCOMP, string)) then
      read (iInputCOMP, *, iostat = ios) ind_refPartRel
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable ind_refPartRel;')"
        print "(  2x,'for the region ',i3,';')", ipart
        stop
      end if                                                   
    else
      print "(/,2x,'Group name OptRegProp not found;')"
      print "(  2x,'for the region ',i3,';')", ipart
      stop  
    end if      
    call check_ind_ref1 (ipart, ind_refPartRel)
    ind_ref(ipart) = ind_refPartRel                        
!
    NsurfPart  = 3
    NparamPart = 2
    OriginPart = 0._O
    do isurf = 1, NsurfPD
      surfPart(isurf) = 1._O
    end do
    lnormPart = 1._O
    string    = 'GeomRegProp'
    if (XFindPar (iInputCOMP, string)) then
      read (iInputCOMP, *, iostat = ios) NsurfPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable NsurfPart;')"
        print "(  2x,'for the region ',i3,';')", ipart
        stop
      end if      
      if (NsurfPart > NsurfPD) then
        print "(/,2x,'Input error: NsurfPart exceeds NsurfPD for the region',i3)", &
                ipart                                    
        stop
      end if            
      read (iInputCOMP, *, iostat = ios) NparamPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable NparamPart;')"
        print "(  2x,'for the region ',i3,';')", ipart
        stop
      end if
      read (iInputCOMP, *, iostat = ios) OriginPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable OriginPart;')"
        print "(  2x,'for the region ',i3,';')", ipart
        stop
      end if
      do isurf = 1, NsurfPart
        read (iInputCOMP, *, iostat = ios) surfPart(isurf)
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable surfPart;')"
          print "(  2x,'for the region ',i3,';')", ipart
          stop
        end if
      end do
      read (iInputCOMP, *, iostat = ios) lnormPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable lnormPart;')"
        print "(  2x,'for the region ',i3,';')", ipart
        stop
      end if
    else
      print "(/,2x,'Group name GeomRegProp not found;')"
      print "(  2x,'for the region ',i3,';')", ipart
      stop  
    end if    
    Nsurf(ipart) = NsurfPart
    Nparam(ipart)= NparamPart
    zpart(ipart) = OriginPart
    lnorm(ipart) = lnormPart
    if (Nsurf(ipart)  > Nsurfmax)   Nsurfmax  = Nsurf(ipart)
    if (Nparam(ipart) > Nparammax)  Nparammax = Nparam(ipart)
  end do
  rewind (unit = iInputCOMP)
  allocate (surf(Npart,Nsurfmax))
  do ipart = 1, Npart 
    do isurf = 1, Nsurfmax
      surf(ipart,isurf) = 0._O
    end do 
    NsurfPart  = 3
    NparamPart = 2
    OriginPart = 0._O
    do isurf = 1, NsurfPD
      surfPart(isurf) = 1._O
    end do
    lnormPart = 1._O
    string    = 'GeomRegProp'
    if (XFindPar (iInputCOMP, string)) then
      read (iInputCOMP, *, iostat = ios) NsurfPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable NsurfPart;')"
        print "(  2x,'for the region ',i3,';')", ipart
        stop
      end if
      read (iInputCOMP, *, iostat = ios) NparamPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable NparamPart;')"
        print "(  2x,'for the region ',i3,';')", ipart
        stop
      end if
      read (iInputCOMP, *, iostat = ios) OriginPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable OriginPart;')"
        print "(  2x,'for the region ',i3,';')", ipart
        stop
      end if
      do isurf = 1, NsurfPart
        read (iInputCOMP, *, iostat = ios) surfPart(isurf)
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable surfPart;')"
          print "(  2x,'for the region ',i3,';')", ipart
          stop
        end if
      end do
      read (iInputCOMP, *, iostat = ios) lnormPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable lnormPart;')"
        print "(  2x,'for the region ',i3,';')", ipart
        stop
      end if
    else
      print "(/,2x,'Group name GeomRegProp not found;')"
      print "(  2x,'for the region ',i3,';')", ipart
      stop  
    end if
    do isurf = 1, NsurfPart
      surf(ipart,isurf) = surfPart(isurf)
    end do                    
  end do 
  call check_geomCOMP (TypeGeom, Npart, Nsurf, Nparam)
!
  allocate (ComplexPlane(Npart), EpsZReIm(Npart))
  rewind (unit = iInputCOMP)
  do ipart = 1, Npart
    ComplexPlane(ipart) = .false.
    EpsZReIm(ipart) = 0.95_O
    if (DS .and. autGenDS) then
      ComplexPlanePart = .false.
      EpsZReImPart = 0.95_O
      string       = 'SourceRegPosAut'
      if (XFindPar (iInputCOMP, string)) then
        read (iInputCOMP, *, iostat = ios) ComplexPlanePart
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable ComplexPlanePart;')"
          print "(  2x,'for the region ',i3,';')", ipart
          stop
        end if
        read (iInputCOMP, *, iostat = ios) EpsZReImPart
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable EpsZReImPart;')"
          print "(  2x,'for the region ',i3,';')", ipart
          stop
        end if               
      else
        print "(/,2x,'Group name SourceRegPosAut not found;')"
        print "(  2x,'for the region ',i3,';')", ipart
        stop  
      end if        
      ComplexPlane(ipart) = ComplexPlanePart
      EpsZReIm(ipart) = EpsZReImPart
    end if
  end do
!  
  if (DoConvTest) then
    print "(/,2x,'Convergence Test for a Composite Particle')" 
    print "(  2x,'-----------------------------------------')"       
  else
    print "(/,2x,'Convergence Test for a Composite Particle over Mrank')" 
    print "(  2x,'----------------------------------------------------')"    
  end if                          
  allocate (Nrankp(Npart), Nrankp1(Npart)) 
  rewind (unit = iInputCOMP)
  Nrankpmax  = 1
  Nrankpmax1 = 1
  do ipart = 1, Npart
    x = k * lnorm(ipart)
    NrankW = int(x + 4.05_O * x**0.33_O + 2._O)
    if (.not. DoConvTest .or. InputDS) then
      NrankPart = 17
      string    = 'NrankReg'
      if (XFindPar (iInputCOMP, string)) then
        read (iInputCOMP, *, iostat = ios) NrankPart
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable NrankPart;')"
          print "(  2x,'for the region ',i3,';')", ipart
          stop
        end if
      else
        print "(/,2x,'Group name NrankReg not found;')"
        print "(  2x,'for the region ',i3,';')", ipart
        stop  
      end if       
      Nrankp(ipart) = NrankPart
      if (ipart == 1) print "(/,2x,'Nrank input values:')"
      print "(2x,'the input value of Nrank for the region ',i3,' is ', i4,',')",    &
              ipart, NrankPart
      print "(2x, a, i3, a)",                                                       &
     'while the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';' 
    else    
      if (ipart == 1) print "(/,2x,'Nrank estimates:')"  
      print "(2x,'the estimated value of Nrank from Wiscombe''s criterion')" 
      print "(2x,'for region ',i2,' is ',i3,';')", ipart, NrankW 
      print "(2x,'- enter the estimated value of Nrank for region ',i2)", ipart
      call read_integer (Nrankp(ipart))
    end if
    if (Nrankp(ipart) > Nrankpmax) Nrankpmax = Nrankp(ipart)
    Nrankp1(ipart) = Nrankp(ipart) - 1
    if (Nrankp1(ipart) > Nrankpmax1) Nrankpmax1 = Nrankp1(ipart)
  end do    
!
  rewind (unit = iInputCOMP)
  x = k * Rcirc
  NrankW = int(x + 4.05_O * x**0.33_O + 2._O)    
  if (DoConvTest) then
    print "(/,2x,'Global Nrank estimate:')"                                                            
    print "(  2x, a, i3, a)",                                                       &  
   'the estimated value of the global Nrank from Wiscombe''s criterion is ',        &
    NrankW, ';'
    print "(/,2x, a)",                                                              &
   '- enter the estimated values of Nint and Nrank, where Nint = Ndgs * Nrankpmax'
    print "(  2x,'  Ndgs = 15,20,..., and Nrankpmax = ',i4,';')", Nrankpmax  
    call read_integer2 (Nint, Nrank)                
  else    
    Nint   = 100
    Nrank  =  17    
    string = 'NintNrank'
    if (XFindPar (iInputCOMP, string)) then
      read (iInputCOMP, *, iostat = ios) Nint
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Nint;')"
        stop
      end if
      read (iInputCOMP, *, iostat = ios) Nrank
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Nrank;')"
        stop
      end if                
    else
      print "(/,2x,'Group name NintNrank not found;')"
      stop  
    end if  
    print "(/,2x,'Input values:')"
    print "(  2x, a, i4, a, i4, a)",                                                & 
   'the input values of Nint and Nrank are ', Nint, ' and ', Nrank, ', respectively,'  
    print "(  2x, a, i3, a)",                                                       &
   'while the estimate of the global Nrank from Wiscombe''s criterion is ', NrankW,','
    print "(  2x, a, i3, a)",                                                       &
   'and Nint = Ndgs * Nrankpmax, with Ndgs = 15,20,..., and Nrankpmax = ',          &
    Nrankpmax,';'                                       
  end if
!
  if (DoConvTest) then     
    print "(/,2x, a)",                                                              &
   '- enter the type of convergence test: 1 - Nint, 2 - Nrank, 3 - Mrank;'
    call read_integerbound (TypeConvTest, 1, 3)     
  else    
    TypeConvTest = 3 
  end if
!
  allocate (zRe(Npart,Nrankpmax), zIm(Npart,Nrankpmax), zRe1(Npart,Nrankpmax1),     &
            zIm1(Npart,Nrankpmax1)) 
  rewind (unit = iInputCOMP)
  do ipart = 1, Npart    
    do i = 1, Nrankpmax         
      zRe(ipart,i) = 0._O
      zIm(ipart,i) = 0._O
    end do
    do i = 1, Nrankpmax1            
      zRe1(ipart,i) = 0._O
      zIm1(ipart,i) = 0._O
    end do
    if (DS) then
      call check_MaxNrank (Nrankpmax)
      if (autGenDS) then
        call zDSCOMP (TypeGeom, Npart, Nsurfmax, surf, Nrankpmax, Nrankp, zpart,     &
             ComplexPlane, EpsZReIm, zRe, zIm)
        if (TypeConvTest == 2) call zDSCOMP (TypeGeom, Npart, Nsurfmax, surf,        &
                                    Nrankpmax1, Nrankp1, zpart, ComplexPlane,        &
                                    EpsZReIm, zRe1, zIm1)
      else
        do i = 1, NrankPD
          zRePart(i)  = 0._O
          zImPart(i)  = 0._O
          zRePart1(i) = 0._O
          zImPart1(i) = 0._O
        end do 
        string = 'SourceRegPosInp'
        if (XFindPar (iInputCOMP, string)) then
          do i = 1, Nrankp(ipart)
            read (iInputCOMP, *, iostat = ios) zRePart(i), zImPart(i)
            if (ios /= 0) then
              print "(/,2x, a)",                                                    &
             'Error by reading the input variables zRePart and zImPart;'	      
              print "(  2x,'for the region ',i3,';')", ipart
              stop
            end if 
          end do
          if (TypeConvTest == 2) then	  
            read (iInputCOMP, *)
            do i = 1, Nrankp(ipart) - 1
              read (iInputCOMP, *, iostat = ios) zRePart1(i), zImPart1(i)
              if (ios /= 0) then
                print "(/,2x, a)",                                                  &
               'Error by reading the input variables zRePart1 and zImPart1;'		
                print "(  2x,'for the region ',i3,';')", ipart
                stop
              end if 
            end do
          end if 	    
        else
          print "(/,2x,'Group name SourceRegPosInp not found;')"
          print "(  2x,'for the region ',i3,';')", ipart
          stop  
        end if
        do i = 1, Nrankp(ipart)
          zRe(ipart,i) = zRePart(i)  
          zIm(ipart,i) = zImPart(i)               
        end do
        if (TypeConvTest == 2) then             
          do i = 1, Nrankp1(ipart)
            zRe1(ipart,i) = zRePart1(i) 
            zIm1(ipart,i) = zImPart1(i) 
          end do          
        end if  
      end if                                
    end if
  end do         
!
  epsNint  = 5.e-2_O
  epsNrank = 5.e-2_O
  epsMrank = 5.e-2_O
  dNint = 4  
  dNintMrank = 10
  string     = 'Errors'
  if (XFindPar (iInputCOMP, string)) then
    read (iInputCOMP, *, iostat = ios) epsNint
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsNint;')"
      stop
    end if
    read (iInputCOMP, *, iostat = ios) epsNrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsNrank;')"
      stop
    end if
    read (iInputCOMP, *, iostat = ios) epsMrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsMrank;')"
      stop
    end if 
    read (iInputCOMP, *, iostat = ios) dNint
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable dNint;')"
      stop
    end if  
    read (iInputCOMP, *, iostat = ios) dNintMrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable dNintMrank;')"
      stop
    end if                 
  else
    print "(/,2x,'Group name Errors not found;')"
    stop  
  end if
!
  FileTmat = '../TMATFILES/T.dat'
  string   = 'Tmat' 
  if (XFindPar (iInputCOMP, string)) then
    read (iInputCOMP, *, iostat = ios) FileTmat
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable FileTmat;')"
      stop
    end if             
  else
    print "(/,2x,'Group name Tmat not found;')"
    stop  
  end if   
!
  PrnProgress = .true.
  string   = 'PrintProgress' 
  if (XFindPar (iInputCOMP, string)) then
    read (iInputCOMP, *, iostat = ios) PrnProgress
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable PrnProgress;')"
      stop
    end if             
  else
    print "(/,2x,'Group name PrintProgress not found;')"
    stop  
  end if        
  close (unit = iInputCOMP)            
end subroutine readinputCOMP      
!***********************************************************************************
subroutine printinputCOMP (ic, TypeGeom, Nsurfmax, Nsurf, surf, Nparam, Npart,      &
           lnorm, anorm, Nrankpmax, Nrankp, zpart, zRe, zIm, Nrankpmax1, Nrankp1,   &
           zRe1, zIm1, ind_ref, dNint, wavelength, ind_refMed, epsNint, epsNrank,   &
           epsMrank, DS, autGenDS)
  use parameters
  implicit none  
  integer    :: ic, TypeGeom, Nsurfmax, Npart, Nrankpmax, Nrankpmax1, dNint,        &
                Nrankp(Npart), Nrankp1(Npart), Nparam(Npart), Nsurf(Npart),         &
                NsurfPart, Nrank, i, j
  real(O)    :: wavelength, anorm, surf(Npart,Nsurfmax), zpart(Npart), lnorm(Npart),&
                zRe(Npart,Nrankpmax), zIm(Npart,Nrankpmax),                         &
                zRe1(Npart,Nrankpmax1), zIm1(Npart,Nrankpmax1),                     &
                epsNint, epsNrank, epsMrank, ind_refMed
  complex(O) :: ind_ref(Npart)  
  logical    :: DS,autGenDS
!
  write (iOutput,"(/,2x,'Input Parameters of the Scattering Problem:',/)")
  write (iOutput,"(2x,'wavelength of the free space, wavelength = ',1pe13.4,';')")  &
         wavelength
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'refractive index of the ambient medium, ind_refMed = ', ind_refMed, ';'
  write (iOutput,"(2x,'relative refractive index of each homogeneous region:')")
  do i = 1, Npart         
    if (i /= Npart) then
      write (iOutput,"(2x,'ind_ref(',i2,') = (',1pe13.4,',',1pe13.4,'),')")         &
             i, ind_ref(i)
    else
      write (iOutput,"(2x,'ind_ref(',i2,') = (',1pe13.4,',',1pe13.4,');')")         &
             i, ind_ref(i)	     	     	     
    end if
  end do    
  write (iOutput,*)
  write (iOutput,"(2x,'index of the particle geometry, TypeGeom = ',i2,';')")       &
         TypeGeom
  if (TypeGeom == 1) then
    write (iOutput,"(2x,'half-spheroids;')")    
  else if (TypeGeom == 2) then
    write (iOutput,"(2x,'three cylinders;')")   
  end if
  write (iOutput,"(2x,'number of homogeneous regions, Npart = ',i2,';')") Npart
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'characteristic length of the composite particle, anorm = ', anorm, ';'
  do i = 1,Npart
    write (iOutput,*)
    write (iOutput,"(2x,'region: ',i2)") i
    write (iOutput,"(2x,'number of surface parameters, Nsurf = ',i2,';')") Nsurf(i)
    write (iOutput,"(2x,'surface parameters:')")
    NsurfPart = Nsurf(i)
    do j = 1, NsurfPart
      write (iOutput,"(2x,'surf(',i2,') = ',1pe10.3,',')") j, surf(i,j)      
    end do       
    if (TypeGeom == 1) then
      write (iOutput,"(2x,'where surf(1) is the semi-axis along the symmetry axis,')")
      write (iOutput,"(2x,'surf(2) is the second semi-axis and')")
      write (iOutput,"(2x,'surf(3) is the axial position of the particle origin;')")
    else if (TypeGeom == 2) then
      write (iOutput,"(2x,'where surf(1) is the half-length of the cylinder')")
      write (iOutput,"(2x,'and   surf(2) is the cylinder radius;')")
    end if  
    write (iOutput,"(2x,'number of integration surfaces, Nparam = ',i2,';')")       &
           Nparam(i)
    write (iOutput,"(2x, a, 1pe10.3, a)")                                           &
   'caracteristic length of the region, lnorm = ', lnorm(i), ';'
  end do
  write (iOutput,*)
  if (.not. DS) then
    write (iOutput,"(2x,'computation with localized sources;')")
    write (iOutput,"(2x,'origins of the localized sorces:')")
    do i = 1, Npart
      if (i /= Npart) then
        write (iOutput,"(2x,'zpart(',i2,') = ',1pe10.3,',')") i, zpart(i)
      else
        write (iOutput,"(2x,'zpart(',i2,') = ',1pe10.3,';')") i, zpart(i)
      end if
    end do
  else
    write (iOutput,"(2x,'computation with distributed sources;')")
    if (autGenDS) then
      write (iOutput,"(2x,'the sources are generated automatically;')") 
    else
      write (iOutput,"(2x,'the sources are user-defined;')")
    end if                      
    write (iOutput,*)       
    write (iOutput,"(2x,'z - coordinates of the sources:')")
    write (iOutput,"(2x,'first configuration:')")
    do i = 1, Npart       
      write (iOutput,"(2x,'region: ',i2)") i
      Nrank = Nrankp(i)
      do j = 1, Nrank       
        if (j /= Nrank) then
          write (iOutput,"(2x,'zRe(',i2,') = ',1pe10.3,', zIm(',i2,') = ',1pe10.3,',')")&
                 j, zRe(i,j), j, zIm(i,j)
        else
          write (iOutput,"(2x,'zRe(',i2,') = ',1pe10.3,', zIm(',i2,') = ',1pe10.3,';')")&
                 j, zRe(i,j), j, zIm(i,j)
        end if
      end do
      write (iOutput,*)
    end do
    if (ic == 2) then         
      write (iOutput,"(2x,'second configuration:')")
      do i = 1, Npart       
        write (iOutput,"(2x,'region: ',i2)") i
        Nrank = Nrankp1(i)
        do j = 1, Nrank     
          if (j /= Nrank) then
            write (iOutput,"(2x,'zRe(',i2,') = ',1pe10.3,', zIm(',i2,') = ',1pe10.3,',')")&
                   j, zRe1(i,j), j, zIm1(i,j)
          else
            write (iOutput,"(2x,'zRe(',i2,') = ',1pe10.3,', zIm(',i2,') = ',1pe10.3,';')")&
                   j, zRe1(i,j), j, zIm1(i,j)
          end if
        end do
        write (iOutput,*)
      end do
    end if  
  end if    
  write (iOutput,"(2x,'integration tolerance, epsNint = ',1pe10.3,';')") epsNint
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum expansion order tolerance, epsNrank = ', epsNrank, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum azimuthal order tolerance, epsMrank = ', epsMrank, ';'          
  write (iOutput,"(2x,'integration step, dNint = ',i4,'.')") dNint  
  write (iOutput,"(/)")                       
end subroutine printinputCOMP
! **********************************************************************************
subroutine convergence_NintCOMP (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,       &
           Nparammax, Nparam, Npart, Nrankp, zpart, Nrank, Nint, dNint, epsNint,    &
           PrnProgress)    
  use parameters 
  implicit none
  integer    :: TypeGeom, Nsurfmax, Nparammax, Npart, Nrank, Nint, dNint,           &
                Nrankp(Npart), Nparam(Npart)
  real(O)    :: k, surf(Npart,Nsurfmax), zpart(Npart), snorm, epsNint                
  complex(O) :: ind_ref(Npart)
  logical    :: PrnProgress
!        
  integer    :: Mstart, Mrank, NmaxC, NmaxL, Nmaxmax, Nteta, i, m, iNint,           &
                dimension_column, dimension_row, NthetaConv, ipart, NrankAL
  real(O)    :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,        &
                Cext, Qext  
  integer,allocatable    :: Nintparam(:,:)
  real(O),allocatable    :: paramG(:,:,:), weightsG(:,:,:), h(:), v(:), oldh(:),    &
                            oldv(:)
  complex(O),allocatable :: a(:,:), b(:,:), c(:), c1(:), cc(:)
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = 0._O
  beta   = 0._O
  gama   = 0._O
  alfap  = Pi / 4._O
  Mstart = 1
  Mrank  = 1     
  m = 1
  NrankAL = 0
  do ipart = 1, Npart
    NrankAL = NrankAL + Nrankp(ipart)
  end do
  NrankAL = max(Nrank,NrankAL)
  NmaxC   = dimension_column (m, Npart, Nrankp)
  NmaxL   = dimension_row (m, Nrank)
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  call write_TypeConvHead (1)
  allocate (a(2*NrankAL,2*NrankAL), b(2*NrankAL,2*NrankAL)) 
  allocate (c(2*Nrank), c1(2*Nrank))  
  allocate (cc(2*Nmaxmax)) 
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do
  if (PrnProgress) call write_progress (.true., 1, 9)
  do iNint = 1, 2 
    allocate (paramG(Npart,Nparammax,Nint), weightsG(Npart,Nparammax,Nint),         &
              Nintparam(Npart,Nparammax))             
    call interpolation_listCOMP (TypeGeom, Npart, Nsurfmax, surf, Nint, Nparammax,  &
         Nintparam, paramG, weightsG)       
    call matrix_Q_COMP (TypeGeom, 3, 1, k, ind_ref, Nsurfmax, surf, m, Npart,       &
         zpart, Nrankp, NmaxC, Nrank, NmaxC, Nint, Nparammax, Nparam, Nintparam,    &
         paramG, weightsG, a, NrankAL, NrankAL) 
    if (PrnProgress) call write_progress (.false., 2+4*(iNint-1), 9)            
    call incident_matrix_COMP (TypeGeom, k, Nsurfmax, surf, m, Npart, zpart,        &
         Nrankp, NmaxC, Nrank, NmaxL, Nint, Nparammax, Nparam, Nintparam, paramG,   &
         weightsG, b, NrankAL, NrankAL)      
    if (PrnProgress) call write_progress (.false., 3+4*(iNint-1), 9)                         
    call LU_SYSTEM_DIRECT (a, 2*NrankAL, 2*NrankAL, b, 2*NrankAL, 2*NrankAL,        &
         2*NmaxC, 2*NmaxL)
    if (PrnProgress) call write_progress (.false., 4+4*(iNint-1), 9)
    call matrix_Q_COMP (TypeGeom, 1, 1, k, ind_ref, Nsurfmax, surf, m, Npart,       &
         zpart, Nrankp, NmaxC, Nrank, NmaxL, Nint, Nparammax, Nparam, Nintparam,    &
         paramG, weightsG, a, NrankAL, NrankAL)                          
    call product_matrices (2*NmaxL, 2*NmaxC, 2*NmaxL, a, 2*NrankAL, 2*NrankAL,      &
         b, 2*NrankAL, 2*NrankAL)
    if (PrnProgress) call write_progress (.false., 5+4*(iNint-1), 9)                     
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,     &
         NmaxL, c)
    call product_matrix_vector (2*NmaxL, 2*NmaxL, a, 2*NrankAL, 2*NrankAL, c, c1)
    call extend_vector_positive (c1, cc, m, Mstart, Nrank, NmaxL, Nmaxmax)                               
    call matrix_m_negativ (NmaxL, NmaxL, a, NrankAL, NrankAL)
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,    &
         NmaxL, c)
    call product_matrix_vector (2*NmaxL, 2*NmaxL, a, 2*NrankAL, 2*NrankAL, c, c1)
    call extend_vector_negative (c1, cc, m, Nrank, NmaxL, Nmaxmax)
    call DSCS (cc, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm, &
        .false.,.true., h, v)
    call CQscat (cc, Mrank, Nrank, Nmaxmax, k, snorm, Cscat, Qscat)
    call CQext (cc, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,         &
         alfap, k, snorm, Cext, Qext)
    call delta_DSCS (Nteta, h, v, oldh, oldv, epsNint, NthetaConv)    
    call write_3ConvParamRegCOMP (Nint, m, Nrank, Npart, Nrankp, .false.)   
    call write_DSCS (Nteta,.false., h, v)
    call write_Effic (Qscat, Qext)
    Nint = Nint + dNint
    deallocate (paramG, weightsG, Nintparam)
  end do
  call write_NintConvRes (NthetaConv, Nteta, epsNint)
  deallocate (a, b, c, c1, cc, h, v, oldh, oldv)
end subroutine convergence_NintCOMP
! **********************************************************************************
subroutine convergence_NrankCOMP (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,      &
           Nparammax, Nparam, Npart, Nrankp, zpart, Nrank, Nint, epsNrank,          &
           PrnProgress)    
  use parameters 
  implicit none
  integer    :: TypeGeom, Nsurfmax, Nparammax, Npart, Nrank, Nint, Nrankp(Npart),   &
                Nparam(Npart)
  real(O)    :: k, surf(Npart,Nsurfmax), zpart(Npart), snorm, epsNrank                
  complex(O) :: ind_ref(Npart)
  logical    :: PrnProgress 
!       
  integer    :: Mstart, Mrank, NmaxC, NmaxL, Nmaxmax, Nteta, i, m, NthetaConv,      &
                dimension_column, dimension_row, ipart, NrankAL
  real(O)    :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,        &
                Cext, Qext  
  integer,allocatable    :: Nintparam(:,:)
  real(O),allocatable    :: paramG(:,:,:), weightsG(:,:,:), h(:), v(:), oldh(:),    &
                            oldv(:)
  complex(O),allocatable :: a(:,:), b(:,:), c(:), c1(:), cc(:), a1(:,:), b1(:,:),   &
                            a0(:,:)
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = 0._O
  beta   = 0._O
  gama   = 0._O
  alfap  = Pi / 4._O
  Mstart = 1
  Mrank  = 1   
  m      = 1
  NrankAL = 0
  do ipart = 1, Npart
    NrankAL = NrankAL + Nrankp(ipart)
  end do
  NrankAL = max(Nrank,NrankAL)
  NmaxC   = dimension_column (m, Npart, Nrankp)
  NmaxL   = dimension_row (m, Nrank) 
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  call write_TypeConvHead (2)
  allocate (a(2*NrankAL,2*NrankAL), b(2*NrankAL,2*NrankAL)) 
  allocate (c(2*Nrank), c1(2*Nrank))
  allocate (a1(2*NrankAL,2*NrankAL), b1(2*NrankAL,2*NrankAL), a0(2*NrankAL,2*NrankAL))  
  allocate (cc(2*Nmaxmax)) 
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))      
  allocate (paramG(Npart,Nparammax,Nint), weightsG(Npart,Nparammax,Nint),           &
            Nintparam(Npart,Nparammax))           
  call interpolation_listCOMP (TypeGeom, Npart, Nsurfmax, surf, Nint, Nparammax,    &
       Nintparam, paramG, weightsG)     
  if (PrnProgress) call write_progress (.true., 1, 9)
! --- Nrank configuration ---  
  call matrix_Q_COMP (TypeGeom, 3, 1, k, ind_ref, Nsurfmax, surf, m, Npart, zpart,  &
       Nrankp, NmaxC, Nrank, NmaxC, Nint, Nparammax, Nparam, Nintparam, paramG,     &
       weightsG, a, NrankAL, NrankAL)
  if (PrnProgress) call write_progress (.false., 2, 9)
  call copy_matrix (2*NmaxC, 2*NmaxC, a, 2*NrankAL, 2*NrankAL, a1,                  &
       2*NrankAL, 2*NrankAL)
  call incident_matrix_COMP (TypeGeom, k, Nsurfmax, surf, m, Npart, zpart, Nrankp,  &
       NmaxC, Nrank, NmaxL, Nint, Nparammax, Nparam, Nintparam, paramG, weightsG,   &
       b, NrankAL, NrankAL)
  if (PrnProgress) call write_progress (.false., 3, 9)
  call copy_matrix (2*NmaxC, 2*NmaxL, b, 2*NrankAL, 2*NrankAL, a0,                  &
       2*NrankAL, 2*NrankAL)
  call LU_SYSTEM_DIRECT (a, 2*NrankAL, 2*NrankAL, b, 2*NrankAL, 2*NrankAL,          &
       2*NmaxC, 2*NmaxL)
  if (PrnProgress) call write_progress (.false., 4, 9)
  call matrix_Q_COMP (TypeGeom, 1, 1, k, ind_ref, Nsurfmax, surf, m, Npart, zpart,  &
       Nrankp, NmaxC, Nrank, NmaxL, Nint, Nparammax, Nparam, Nintparam, paramG,     &
       weightsG, a, NrankAL, NrankAL)  
  call copy_matrix (2*NmaxL, 2*NmaxC, a, 2*NrankAL, 2*NrankAL, b1,                  &
       2*NrankAL, 2*NrankAL)
  call product_matrices (2*NmaxL, 2*NmaxC, 2*NmaxL, a, 2*NrankAL, 2*NrankAL, b,     &
       2*NrankAL, 2*NrankAL)
  if (PrnProgress) call write_progress (.false., 5, 9)               
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,       &
       NmaxL, c)
  call product_matrix_vector (2*NmaxL, 2*NmaxL, a, 2*NrankAL, 2*NrankAL, c, c1)
  call extend_vector_positive (c1, cc, m, Mstart, Nrank, NmaxL, Nmaxmax)                                 
  call matrix_m_negativ (NmaxL, NmaxL, a, NrankAL, NrankAL)
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,      &
       NmaxL, c)
  call product_matrix_vector (2*NmaxL, 2*NmaxL, a, 2*NrankAL, 2*NrankAL, c, c1)
  call extend_vector_negative (c1, cc, m, Nrank, NmaxL, Nmaxmax)
  call DSCS (cc, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm,   &
      .false.,.true., h, v)
  call CQscat (cc, Mrank, Nrank, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,           &
       alfap, k, snorm, Cext, Qext)
  do i = 1, Nteta
    oldh(i) = h(i)
    oldv(i) = v(i)
  end do
  call write_3ConvParamRegCOMP (Nint, m, Nrank, Npart, Nrankp, .false.)
  call write_DSCS (Nteta,.false., h, v)
  call write_Effic (Qscat, Qext)
! --- (Nrank - 1) configuration --- 
  call copy_matrix (2*NmaxC, 2*NmaxC, a1, 2*NrankAL, 2*NrankAL, a,                  &
       2*NrankAL, 2*NrankAL)
  call internal_matrix_Nrank_m (m, Npart, Nrankp, NmaxC, a, NrankAL, NrankAL)        
  if (PrnProgress) call write_progress (.false., 6, 9)
  call copy_matrix (2*NmaxC, 2*NmaxL, a0, 2*NrankAL, 2*NrankAL, b,                  &
       2*NrankAL, 2*NrankAL)
  call incident_matrix_Nrank_m (m, Npart, Nrankp, NmaxC, NmaxL, b, NrankAL, NrankAL)                            
  if (PrnProgress) call write_progress (.false., 7, 9)
  call LU_SYSTEM_DIRECT (a, 2*NrankAL, 2*NrankAL, b, 2*NrankAL, 2*NrankAL,          &
       2*NmaxC, 2*NmaxL)    
  if (PrnProgress) call write_progress (.false., 8, 9)
  call copy_matrix (2*NmaxL, 2*NmaxC, b1, 2*NrankAL, 2*NrankAL, a,                  &
       2*NrankAL, 2*NrankAL)
  call difusion_matrix_Nrank_m (m, Npart, Nrankp, NmaxC, NmaxL, a, NrankAL, NrankAL)                        
  call product_matrices (2*NmaxL, 2*NmaxC, 2*NmaxL, a, 2*NrankAL, 2*NrankAL, b,     &
       2*NrankAL,2*NrankAL) 
  if (PrnProgress) call write_progress (.false., 9, 9)          
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,       &
       NmaxL, c)
  call product_matrix_vector (2*NmaxL, 2*NmaxL, a, 2*NrankAL, 2*NrankAL, c, c1)
  call extend_vector_positive (c1, cc, m, Mstart, Nrank, NmaxL, Nmaxmax)                                 
  call matrix_m_negativ (NmaxL, NmaxL, a, NrankAL, NrankAL)
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,      &
       NmaxL, c)
  call product_matrix_vector (2*NmaxL, 2*NmaxL, a, 2*NrankAL, 2*NrankAL, c, c1)
  call extend_vector_negative (c1, cc, m, Nrank, NmaxL, Nmaxmax)
  call DSCS (cc, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm,   &
      .false.,.true., h, v)
  call CQscat (cc, Mrank, Nrank, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,           &
       alfap, k, snorm, Cext, Qext)
  call delta_DSCS (Nteta, h, v, oldh, oldv, epsNrank, NthetaConv)  
  call write_3ConvParamRegCOMP (Nint, m, Nrank - 1, Npart, Nrankp, .true.)
  call write_DSCS (Nteta,.false., h, v)
  call write_Effic (Qscat, Qext)
  call write_NrankConvRes (NthetaConv, Nteta, epsNrank)
  deallocate (a1, b1, a0)
  deallocate (a, b, c, c1, cc, h, v, oldh, oldv, paramG, weightsG, Nintparam)
end subroutine convergence_NrankCOMP
! **********************************************************************************
subroutine convergence_MrankCOMP (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,      &
           Nparammax, Nparam, Npart, Nrankp, zpart, Nrank, Nint, epsMrank,          &
           FileTmat, PrnProgress)          
  use parameters 
  implicit none
  integer       :: TypeGeom, Nsurfmax, Nparammax, Npart, Nrank, Nint,               &
                   Nrankp(Npart), Nparam(Npart)
  real(O)       :: k, surf(Npart,Nsurfmax), zpart(Npart), snorm, epsMrank                
  complex(O)    :: ind_ref(Npart)
  character(80) :: FileTmat
  logical       :: PrnProgress
!             
  integer       :: Mstart, Mrank, NmaxC, NmaxL, Nmaxmax, Nteta, i, m, NthetaConv,   &
                   dimension_column, dimension_row, ipart, NrankAL
  real(O)       :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,     &
                   Cext, Qext
  integer,allocatable    :: Nintparam(:,:)
  real(O),allocatable    :: paramG(:,:,:), weightsG(:,:,:), h(:), v(:), oldh(:),    &
                            oldv(:)
  complex(O),allocatable :: a(:,:), b(:,:), c(:), c1(:), cc(:)
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = Pi / 4._O
  beta   = Pi / 4._O
  gama   = 0._O
  alfap  = Pi / 4._O  
  Mstart = 0
  Mrank  = Nrank   
  NrankAL = 0
  do ipart = 1, Npart
    NrankAL = NrankAL + Nrankp(ipart)
  end do
  NrankAL = max(Nrank,NrankAL)
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  open (unit = iTmat, file = FileTmat, status = 'replace')
  call  write_HeadFileTmat (NrankAL, NrankAL)    
  call  write_TypeConvHead (3)
  call  write_2ConvParamRegCOMP (Nint, Nrank, Npart, Nrankp)
  allocate (a(2*NrankAL,2*NrankAL), b(2*NrankAL,2*NrankAL))
  allocate (c(2*Nrank), c1(2*Nrank))            
  allocate (cc(2*Nmaxmax)) 
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do     
  call zero_matrix (2*NrankAL, 2*NrankAL, a, 2*NrankAL, 2*NrankAL)
  allocate (paramG(Npart,Nparammax,Nint), weightsG(Npart,Nparammax,Nint),           &
            Nintparam(Npart,Nparammax))           
  call interpolation_listCOMP (TypeGeom, Npart, Nsurfmax, surf, Nint, Nparammax,    &
       Nintparam, paramG, weightsG)     
  Mrank = - 1
  do m = Mstart, Nrank    
    call write_1ConvParam (m)
    Mrank = Mrank + 1
    do i = 1, Npart
      if (Nrankp(i) < m) Nrankp(i) = m
    end do      
    NmaxC = dimension_column (m, Npart, Nrankp)
    NmaxL = dimension_row (m, Nrank)
    if (PrnProgress) call write_progress_m (.true., m, 1, 5)  
    call matrix_Q_COMP (TypeGeom, 3, 1, k, ind_ref, Nsurfmax, surf, m, Npart,       &
         zpart, Nrankp, NmaxC, Nrank, NmaxC, Nint, Nparammax, Nparam, Nintparam,    &
         paramG, weightsG, a, NrankAL, NrankAL)
    if (PrnProgress) call write_progress_m (.false., m, 2, 5)  
    call incident_matrix_COMP (TypeGeom, k, Nsurfmax, surf, m, Npart, zpart,        &
         Nrankp, NmaxC, Nrank, NmaxL, Nint, Nparammax, Nparam, Nintparam, paramG,   &
         weightsG, b, NrankAL, NrankAL)
    if (PrnProgress) call write_progress_m (.false., m, 3, 5)                                
    call LU_SYSTEM_DIRECT (a, 2*NrankAL, 2*NrankAL, b, 2*NrankAL, 2*NrankAL,        &
         2*NmaxC, 2*NmaxL)
    if (PrnProgress) call write_progress_m (.false., m, 4, 5) 
    call matrix_Q_COMP (TypeGeom, 1, 1, k, ind_ref, Nsurfmax, surf, m, Npart,       &
         zpart, Nrankp, NmaxC, Nrank, NmaxL, Nint, Nparammax, Nparam, Nintparam,    &
         paramG, weightsG, a, NrankAL, NrankAL)                      
    call product_matrices (2*NmaxL, 2*NmaxC, 2*NmaxL, a, 2*NrankAL, 2*NrankAL, b,       &
         2*NrankAL, 2*NrankAL)
    if (PrnProgress) call write_progress_m (.false., m, 5, 5) 
    call write_FileTmat (NrankAL, NrankAL, a)   
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,     &
         NmaxL, c)
    call product_matrix_vector (2*NmaxL, 2*NmaxL, a, 2*NrankAL, 2*NrankAL, c, c1)
    call extend_vector_positive (c1, cc, m, Mstart, Nrank, NmaxL, Nmaxmax)                               
    if (m /= 0) then
      call matrix_m_negativ (NmaxL, NmaxL, a, NrankAL, NrankAL)
      call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,  &
           NmaxL, c)
      call product_matrix_vector (2*NmaxL, 2*NmaxL, a, 2*NrankAL, 2*NrankAL, c, c1)
      call extend_vector_negative (c1, cc, m, Nrank, NmaxL, Nmaxmax)
    end if  
    call DSCS (cc, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm, &
        .false.,.true., h, v)    
    call delta_DSCS (Nteta, h, v, oldh, oldv, epsMrank, NthetaConv)    
    call write_DSCS (Nteta,.false., h, v)
    if (NthetaConv >= int(0.8*Nteta)) exit
  end do    
  close (unit = iTmat)
  call CQscat (cc, Mrank, Nrank, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,           &
       alfap, k, snorm, Cext, Qext)
  call write_Effic (Qscat, Qext)              
  call write_MrankConvRes (NthetaConv, epsMrank)
  if (NthetaConv >= int(0.8*Nteta)) then
    print "(/,2x,'Convergence criterion for Mrank is satisfied;')"       
  else
    print "(/,2x,'Convergence criterion for Mrank is not satisfied;')"
  end if      
  call write_InfoFileTmat (FileTmat, Mrank, Nrank, .true., .false., .false.)
  call ScatCharact (k, FileTmat, Mrank, Nrank, .true., .false., .false.)    
  print "(/,2x,'T matrix is stored in ',a50)", FileTmat
  print "(  2x,'The dimensions of the T matrix are given by:')"
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", Nrank  
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank  
  deallocate (a, b, c, c1, cc, h, v, oldh, oldv, paramG, weightsG, Nintparam)
end subroutine convergence_MrankCOMP    
! **********************************************************************************
subroutine convergence_NintDSCOMP (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,     &
           Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zRe, zIm, zpart, Nrank,     &
           Nint, dNint, epsNint, PrnProgress)      
  use parameters 
  implicit none
  integer    :: TypeGeom, Nsurfmax, Nparammax, Npart, Nrankpmax, Nrank, Nint,       &
                dNint, Nrankp(Npart), Nparam(Npart)
  real(O)    :: k, surf(Npart,Nsurfmax), zpart(Npart), snorm, epsNint,              &
                zRe(Npart,Nrankpmax), zIm(Npart,Nrankpmax)
  complex(O) :: ind_ref(Npart)
  logical    :: PrnProgress
!        
  integer    :: Mstart, Mrank, Nteta, NmaxC, NmaxL, Nmaxmax, i, m, iNint,           &
                NthetaConv, ipart, dimension_row, NrankAL
  real(O)    :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,        &
                Cext, Qext  
  integer,allocatable    :: Nintparam(:,:)
  real(O),allocatable    :: paramG(:,:,:), weightsG(:,:,:), h(:), v(:), oldh(:),    &
                            oldv(:)
  complex(O),allocatable :: a(:,:), b(:,:), c(:), c1(:), cc(:)
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = 0._O
  beta   = 0._O
  gama   = 0._O
  alfap  = Pi / 4._O
  Mstart = 1
  Mrank  = 1  
  m = 1
  NmaxC = 0
  do ipart = 1, Npart
    NmaxC = NmaxC + Nrankp(ipart)
  end do
  NmaxL   = dimension_row (m, Nrank)
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  NrankAL = max(NmaxC,Nrank)
  call write_TypeConvHead (1)
  allocate (a(2*NrankAL,2*NrankAL), b(2*NrankAL,2*NrankAL))
  allocate (c(2*Nrank), c1(2*Nrank))        
  allocate (cc(2*Nmaxmax)) 
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do
  if (PrnProgress) call write_progress (.true., 1, 9)
  do iNint = 1, 2 
    allocate (paramG(Npart,Nparammax,Nint), weightsG(Npart,Nparammax,Nint),         &
              Nintparam(Npart,Nparammax))             
    call interpolation_listCOMP (TypeGeom, Npart, Nsurfmax, surf, Nint, Nparammax,  &
         Nintparam, paramG, weightsG)       
    call matrix_Q_DS_COMP (TypeGeom, 3, 1, k, ind_ref, Nsurfmax, surf, Npart,       &
         Nrankpmax, Nrankp, zRe, zIm, zpart, m, NmaxC, Nrank, NmaxC, Nint,          &
         Nparammax, Nparam, Nintparam, paramG, weightsG, a, NrankAL, NrankAL)
    if (PrnProgress) call write_progress (.false., 2+4*(iNint-1), 9)     
    call incident_matrix_DS_COMP (TypeGeom, k, Nsurfmax, surf, Npart, Nrankpmax,    &
         Nrankp, zRe, zIm, zpart, m, NmaxC, Nrank, NmaxL, Nint, Nparammax, Nparam,  &
         Nintparam, paramG, weightsG, b, NrankAL, NrankAL)
    if (PrnProgress) call write_progress (.false., 3+4*(iNint-1), 9)     
    call LU_SYSTEM_DIRECT (a, 2*NrankAL, 2*NrankAL, b, 2*NrankAL, 2*NrankAL,        &
         2*NmaxC, 2*NmaxL)             
    if (PrnProgress) call write_progress (.false., 4+4*(iNint-1), 9)     
    call matrix_Q_DS_COMP (TypeGeom, 1, 1, k, ind_ref, Nsurfmax, surf, Npart,       &
         Nrankpmax, Nrankp, zRe, zIm, zpart, m, NmaxC, Nrank, NmaxL, Nint,          &
         Nparammax, Nparam, Nintparam, paramG, weightsG, a, NrankAL, NrankAL)  
    call product_matrices (2*NmaxL, 2*NmaxC, 2*NmaxL, a, 2*NrankAL, 2*NrankAL,      &
         b, 2*NrankAL, 2*NrankAL)   
    if (PrnProgress) call write_progress (.false., 5+4*(iNint-1), 9)     
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,     &
         NmaxL, c)
    call product_matrix_vector (2*NmaxL, 2*NmaxL, a, 2*NrankAL, 2*NrankAL, c, c1)
    call extend_vector_positive (c1, cc, m, Mstart, Nrank, NmaxL, Nmaxmax)                                   
    call matrix_m_negativ (NmaxL, NmaxL, a, NrankAL, NrankAL)
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,    &
         NmaxL, c)
    call product_matrix_vector (2*NmaxL, 2*NmaxL, a, 2*NrankAL, 2*NrankAL, c, c1)
    call extend_vector_negative (c1, cc, m, Nrank, NmaxL, Nmaxmax)
    call DSCS (cc, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm, &
        .false.,.true., h, v)
    call CQscat (cc, Mrank, Nrank, Nmaxmax, k, snorm, Cscat, Qscat)
    call CQext (cc, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,         &
         alfap, k, snorm, Cext, Qext)
    call delta_DSCS (Nteta, h, v, oldh, oldv, epsNint, NthetaConv)
    call write_3ConvParamRegCOMP (Nint, m, Nrank, Npart, Nrankp, .false.)
    call write_DSCS (Nteta,.false., h, v)
    call write_Effic (Qscat, Qext)
    Nint = Nint + dNint
    deallocate (paramG, weightsG, Nintparam)
  end do
  call write_NintConvRes (NthetaConv, Nteta, epsNint)
  deallocate (a, b, c, c1, cc, h, v, oldh, oldv)
end subroutine convergence_NintDSCOMP
! **********************************************************************************   
subroutine convergence_NrankDSCOMP (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,    &
           Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zRe, zIm, Nrankpmax1,       &
           Nrankp1, zRe1, zIm1, zpart, Nrank, Nint, epsNrank, PrnProgress)         
  use parameters 
  implicit none
  integer    :: TypeGeom, Nsurfmax, Nparammax, Npart, Nrankpmax, Nrankpmax1,        &
                Nrank, Nint, Nrankp(Npart), Nrankp1(Npart), Nparam(Npart)
  real(O)    :: k, surf(Npart,Nsurfmax), zpart(Npart), snorm, epsNrank,             &
                zRe(Npart,Nrankpmax), zIm(Npart,Nrankpmax),                         &
                zRe1(Npart,Nrankpmax1), zIm1(Npart,Nrankpmax1)
  complex(O) :: ind_ref(Npart)
  logical    :: PrnProgress
!        
  integer    :: Mstart, Mrank, Nteta, NmaxC, NmaxL, Nmaxmax, i, m, NthetaConv,      &
                ipart, dimension_row, Nrank1, NrankAL
  real(O)    :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,        &
                Cext, Qext  
  integer,allocatable    :: Nintparam(:,:)
  real(O),allocatable    :: paramG(:,:,:), weightsG(:,:,:), h(:), v(:), oldh(:),    &
                            oldv(:)
  complex(O),allocatable :: a(:,:), b(:,:), c(:), c1(:), cc(:)
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = 0._O
  beta   = 0._O
  gama   = 0._O
  alfap  = Pi / 4._O
  Mstart = 1
  Mrank  = 1  
  m = 1
  NmaxC = 0
  do ipart = 1, Npart
    NmaxC = NmaxC + Nrankp(ipart)
  end do
  NmaxL   = dimension_row (m, Nrank)
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  NrankAL = max(NmaxC,Nrank)
  call write_TypeConvHead (2)
  allocate (a(2*NrankAL,2*NrankAL), b(2*NrankAL,2*NrankAL)) 
  allocate (c(2*Nrank), c1(2*Nrank))        
  allocate (cc(2*Nmaxmax)) 
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))    
  allocate (paramG(Npart,Nparammax,Nint), weightsG(Npart,Nparammax,Nint),           &
            Nintparam(Npart,Nparammax))           
  call interpolation_listCOMP (TypeGeom, Npart, Nsurfmax, surf, Nint, Nparammax,    &
       Nintparam, paramG, weightsG)     
  if (PrnProgress) call write_progress (.true., 1, 9)
! --- Nrank configuration ---  
  call matrix_Q_DS_COMP (TypeGeom, 3, 1, k, ind_ref, Nsurfmax, surf, Npart,         &
       Nrankpmax, Nrankp, zRe, zIm, zpart, m, NmaxC, Nrank, NmaxC, Nint, Nparammax, &
       Nparam, Nintparam, paramG, weightsG, a, NrankAL, NrankAL)
  if (PrnProgress) call write_progress (.false., 2, 9)
  call incident_matrix_DS_COMP (TypeGeom, k, Nsurfmax, surf, Npart, Nrankpmax,      &
       Nrankp, zRe, zIm, zpart, m, NmaxC, Nrank, NmaxL, Nint, Nparammax, Nparam,    &
       Nintparam, paramG, weightsG, b, NrankAL, NrankAL)
  if (PrnProgress) call write_progress (.false., 3, 9)
  call LU_SYSTEM_DIRECT (a, 2*NrankAL, 2*NrankAL, b, 2*NrankAL, 2*NrankAL,          &
       2*NmaxC, 2*NmaxL)               
  if (PrnProgress) call write_progress (.false., 4, 9)
  call matrix_Q_DS_COMP (TypeGeom, 1, 1, k, ind_ref, Nsurfmax, surf, Npart,         &
       Nrankpmax, Nrankp, zRe, zIm, zpart, m, NmaxC, Nrank, NmaxL, Nint, Nparammax, &
       Nparam, Nintparam, paramG, weightsG, a, NrankAL, NrankAL)        
  call product_matrices (2*NmaxL, 2*NmaxC, 2*NmaxL, a, 2*NrankAL, 2*NrankAL,        &
       b, 2*NrankAL, 2*NrankAL)
  if (PrnProgress) call write_progress (.false., 5, 9)           
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,       &
       NmaxL, c)
  call product_matrix_vector (2*NmaxL, 2*NmaxL, a, 2*NrankAL, 2*NrankAL, c, c1)
  call extend_vector_positive (c1, cc, m, Mstart, Nrank, NmaxL, Nmaxmax)                                     
  call matrix_m_negativ (NmaxL, NmaxL, a, NrankAL, NrankAL)
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,      &
       NmaxL, c)
  call product_matrix_vector (2*NmaxL, 2*NmaxL, a, 2*NrankAL, 2*NrankAL, c, c1)
  call extend_vector_negative (c1, cc, m, Nrank, NmaxL, Nmaxmax)
  call DSCS (cc, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm,   &
      .false.,.true., h, v)
  call CQscat (cc, Mrank, Nrank, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,           &
       alfap, k, snorm, Cext, Qext)
  do i = 1, Nteta
    oldh(i) = h(i)
    oldv(i) = v(i)
  end do
  call write_3ConvParamRegCOMP (Nint, m, Nrank, Npart, Nrankp, .false.)
  call write_DSCS (Nteta,.false., h, v)
  call write_Effic (Qscat, Qext)
! --- Nrank1 configuration ---
  NmaxC = 0
  do ipart = 1, Npart
    NmaxC = NmaxC + Nrankp1(ipart)
  end do
  Nrank1  = Nrank - 1
  NmaxL   = dimension_row (m, Nrank1)
  Nmaxmax = Nrank1 + Mrank * (2 * Nrank1 - Mrank + 1)
  call matrix_Q_DS_COMP (TypeGeom, 3, 1, k, ind_ref, Nsurfmax, surf, Npart,         &
       Nrankpmax1, Nrankp1, zRe1, zIm1, zpart, m, NmaxC, Nrank1, NmaxC, Nint,       &
       Nparammax, Nparam, Nintparam, paramG, weightsG, a, NrankAL, NrankAL)
  if (PrnProgress) call write_progress (.false., 6, 9)
  call incident_matrix_DS_COMP (TypeGeom, k, Nsurfmax, surf, Npart, Nrankpmax1,     &
       Nrankp1, zRe1, zIm1, zpart, m, NmaxC, Nrank1, NmaxL, Nint, Nparammax,        &
       Nparam, Nintparam, paramG, weightsG, b, NrankAL, NrankAL)
  if (PrnProgress) call write_progress (.false., 7, 9)
  call LU_SYSTEM_DIRECT (a, 2*NrankAL, 2*NrankAL, b, 2*NrankAL, 2*NrankAL,          &
       2*NmaxC, 2*NmaxL)              
  if (PrnProgress) call write_progress (.false., 8, 9)
  call matrix_Q_DS_COMP (TypeGeom, 1, 1, k, ind_ref, Nsurfmax, surf, Npart,         &
       Nrankpmax1, Nrankp1, zRe1, zIm1, zpart, m, NmaxC, Nrank1, NmaxL, Nint,       &
       Nparammax, Nparam, Nintparam, paramG, weightsG, a, NrankAL, NrankAL)        
  call product_matrices (2*NmaxL, 2*NmaxC, 2*NmaxL, a, 2*NrankAL, 2*NrankAL,        &
       b, 2*NrankAL, 2*NrankAL)  
  if (PrnProgress) call write_progress (.false., 9, 9)
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank1,      &
       NmaxL, c)
  call product_matrix_vector (2*NmaxL, 2*NmaxL, a, 2*NrankAL, 2*NrankAL, c, c1)
  call extend_vector_positive (c1, cc, m, Mstart, Nrank1, NmaxL, Nmaxmax)                                    
  call matrix_m_negativ (NmaxL, NmaxL, a, NrankAL, NrankAL)
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank1,     &
       NmaxL, c)
  call product_matrix_vector (2*NmaxL, 2*NmaxL, a, 2*NrankAL, 2*NrankAL, c, c1)
  call extend_vector_negative (c1, cc, m, Nrank1, NmaxL, Nmaxmax)
  call DSCS (cc, Mrank, Nrank1, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm,  &
      .false.,.true., h, v)
  call CQscat (cc, Mrank, Nrank1, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, Nrank1, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,          &
       alfap, k, snorm, Cext, Qext)
  call delta_DSCS (Nteta, h, v, oldh, oldv, epsNrank, NthetaConv)
  call write_3ConvParamRegCOMP (Nint, m, Nrank1, Npart, Nrankp1, .false.)
  call write_DSCS (Nteta,.false., h, v)
  call write_Effic (Qscat, Qext)
  call write_NrankConvRes (NthetaConv, Nteta, epsNrank)
  deallocate (a, b, c, c1, cc, h, v, oldh, oldv, paramG, weightsG, Nintparam)
end subroutine convergence_NrankDSCOMP 
! **********************************************************************************
subroutine convergence_MrankDSCOMP (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,    &
           Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zRe, zIm, zpart, Nrank,     &
           Nint, ComplexPlane, EpsZReIm, autGenDS, epsMrank, epsNrank, dNintMrank,  &
           FileTmat, PrnProgress)      
  use parameters 
  implicit none
  integer       :: TypeGeom, Nsurfmax, Nparammax, Npart, Nrankpmax, Nrank, Nint,    &
                   Nrankp(Npart), Nparam(Npart), dNintMrank
  real(O)       :: k, surf(Npart,Nsurfmax), zpart(Npart), snorm, epsMrank,          &
                   epsNrank, zRe(Npart,Nrankpmax), zIm(Npart,Nrankpmax),            &
                   EpsZReIm(Npart)
  complex(O)    :: ind_ref(Npart)
  character(80) :: FileTmat
  logical       :: autGenDS, PrnProgress, ComplexPlane(Npart)
!             
  integer       :: Mstart, Mrank, Nteta, NmaxC, NmaxL, Nmaxmax, i, m, NthetaConv,   &
                   ipart, dimension_row, NrankAL, NrankpmaxL, NrankpmaxL1
  real(O)       :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,     &
                   Cext, Qext
  logical       :: ComplexDS, ChangeDS, more, continuare, ConvTest
  integer,allocatable    :: Nintparam(:,:), Nrankp1(:)
  real(O),allocatable    :: paramG(:,:,:), weightsG(:,:,:), h(:), v(:), oldh(:),    &
                            oldv(:), zReL(:,:), zImL(:,:), zReL1(:,:), zImL1(:,:)
  complex(O),allocatable :: a(:,:), b(:,:), c(:), c1(:), cc(:)
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = Pi / 4._O
  beta   = Pi / 4._O
  gama   = 0._O
  alfap  = Pi / 4._O
  Mstart = 0
  Mrank  = Nrank  
  NmaxC = 0
  do ipart = 1, Npart
    NmaxC = NmaxC + Nrankp(ipart)
  end do
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  NrankAL = max(NmaxC,Nrank)
  open (unit = iTmat, file = FileTmat, status = 'replace')
  call write_HeadFileTmat (NrankAL, NrankAL) 
  call write_TypeConvHead (3)
  call write_2ConvParamRegCOMP (Nint, Nrank, Npart, Nrankp)
  allocate (a(2*NrankAL,2*NrankAL), b(2*NrankAL,2*NrankAL))
  allocate (c(2*Nrank), c1(2*Nrank))       
  allocate (cc(2*Nmaxmax)) 
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do
  call zero_matrix (2*NrankAL, 2*NrankAL, a, 2*NrankAL, 2*NrankAL)     
  allocate (paramG(Npart,Nparammax,Nint), weightsG(Npart,Nparammax,Nint),           &
            Nintparam(Npart,Nparammax))           
  call interpolation_listCOMP (TypeGeom, Npart, Nsurfmax, surf, Nint, Nparammax,    &
       Nintparam, paramG, weightsG)
  ComplexDS = .false.
  do ipart = 1, Npart           
    do i = 1,Nrankp(ipart)                    
      if (zIm(ipart,i) /= 0._O) ComplexDS = .true.
    end do                      
  end do
  ChangeDS = .false.
! -----------------------------------------------------------------------------------  
! If the sources are distributed in the complex plane and they are automatically    !
! generated, an interactive convergence test over Nint and Nrank can be performed   !
! at each azimuthal mode calculation. To make this test active comment out the next !
! lines.                                                                            !  
! -----------------------------------------------------------------------------------   
! if (ComplexDS .and. autGenDS) then
!   print "(/,2x, a)",                                                              &
!   '- enter true to change the values of Nint and Nrank at each azimuthal mode and'	                                  
!   print "( 2x, 'false otherwise')"        
!   call read_logical (ChangeDS)    
! end if
  Mrank = - 1
  do m = Mstart, Nrank    
    call write_1ConvParam (m)
    Mrank = Mrank + 1
    NmaxL = dimension_row (m, Nrank)                                     
!   --- local convergence test for sources distributed in the complex plane ---    
    NrankpmaxL = Nrankpmax
    allocate (zReL(Npart,NrankpmaxL), zImL(Npart,NrankpmaxL))   
    do ipart = 1, Npart	    
      do i = 1, Nrankp(ipart)
        zReL(ipart,i) = zRe(ipart,i)
        zImL(ipart,i) = zIm(ipart,i)	    	
      end do        
    end do    
    if (ComplexDS .and. m > 1) then
      if (ChangeDS) then
        print "(/,2x,'Azimuthal mode: m = ',i3)", m  
        more = .true.
        do while (more)
          print "(  2x, '- enter the number of integration points Nint')"
          call read_integer (Nint)
          deallocate (paramG, weightsG, Nintparam)      
          allocate (paramG(Npart,Nparammax,Nint), weightsG(Npart,Nparammax,Nint),   &
                    Nintparam(Npart,Nparammax))            
          call interpolation_listCOMP (TypeGeom, Npart, Nsurfmax, surf, Nint,       &
               Nparammax, Nintparam, paramG, weightsG)       	   	   	   	   
          NrankpmaxL = 0 
	  NmaxC      = 0          	   
          do ipart = 1, Npart
            if (ComplexPlane(ipart)) then
              print "(2x,'- enter the estimated value of Nrank for region ',i2)",   &
	              ipart
              call read_integer (Nrankp(ipart))		
	    end if	 
            if (Nrankp(ipart) > NrankpmaxL)  NrankpmaxL  = Nrankp(ipart)
            NmaxC = NmaxC + Nrankp(ipart)	     	  
          end do
	  if (NrankpmaxL > Nrankpmax) then
            print "(/,2x,'Input error:')"
            print "(  2x, a)",                                                      &
           'the number of discrete sources exceeds the number of discrete'	    
            print "(  2x,'sources corresponding to the initial configuration;')"
            stop  
	  end if	  	  	  	  	  	  	  
          deallocate (zReL, zImL)	   	   
          allocate (zReL(Npart,NrankpmaxL), zImL(Npart,NrankpmaxL))
          call zDSCOMP (TypeGeom, Npart, Nsurfmax, surf, NrankpmaxL, Nrankp, zpart, &
               ComplexPlane, EpsZReIm, zReL, zImL)		     
          print "(/,2x, a)",                                                        &
         '- enter true to perform a convergence test over the new values of Nrank and'	                                  
          print "( 2x, 'false otherwise')"             
	  call read_logical (ConvTest)		
	  if (ConvTest) then
	    allocate (Nrankp1(Npart))        
	    NrankpmaxL1 = 0    	   
            do ipart = 1, Npart          
	      Nrankp1(ipart) = Nrankp(ipart) - 1	          
	      if (Nrankp1(ipart) > NrankpmaxL1) NrankpmaxL1 = Nrankp1(ipart) 	   
            end do  	     
	    allocate (zReL1(Npart,NrankpmaxL1), zImL1(Npart,NrankpmaxL1))       
            call zDSCOMP (TypeGeom, Npart, Nsurfmax, surf, NrankpmaxL1, Nrankp1,    &
                 zpart, ComplexPlane, EpsZReIm, zReL1, zImL1)    	      		   
	    call convergence_NrankDSCOMP_m (m, TypeGeom, k, ind_ref, snorm,         &
                 Nsurfmax, surf, Nparammax, Nparam, Npart, NrankpmaxL, Nrankp,      &
                 zReL, zImL, NrankpmaxL1, Nrankp1, zReL1, zImL1, zpart, Nrank,      &
                 Nint, Nintparam, paramG, weightsG, epsNrank) 
            deallocate (zReL1, zImL1, Nrankp1)	     
          end if 	       	  	 	             	  	   	    	   	  	
	  if (ConvTest) then
            print "(2x, a)",                                                        &
           '- enter true for new input values of Nint and Nrank and false to continue;'	
            call read_logical (continuare)
            if (.not. continuare) more = .false.             	      
	  else
            more = .false.
	  end if  
        end do
      else
        deallocate (paramG, weightsG, Nintparam) 
        Nint = Nint + dNintMrank	     
        allocate (paramG(Npart,Nparammax,Nint), weightsG(Npart,Nparammax,Nint),     &
                  Nintparam(Npart,Nparammax))            
        call interpolation_listCOMP (TypeGeom, Npart, Nsurfmax, surf, Nint,         &
             Nparammax, Nintparam, paramG, weightsG)                                      
      end if  		   	   	   	   	   	   	                    
    end if      
!   --- end local convergence test ---                                                                                                                          
    if (PrnProgress) call write_progress_m (.true., m, 1, 5)   
    call matrix_Q_DS_COMP (TypeGeom, 3, 1, k, ind_ref, Nsurfmax, surf, Npart,       &
         NrankpmaxL, Nrankp, zReL, zImL, zpart, m, NmaxC, Nrank, NmaxC, Nint,       &
         Nparammax, Nparam, Nintparam, paramG, weightsG, a, NrankAL, NrankAL)	     	 
    if (PrnProgress) call write_progress_m (.false., m, 2, 5)  
    call incident_matrix_DS_COMP (TypeGeom, k, Nsurfmax, surf, Npart, NrankpmaxL,   &
         Nrankp, zReL, zImL, zpart, m, NmaxC, Nrank, NmaxL, Nint, Nparammax,        &
         Nparam, Nintparam, paramG, weightsG, b, NrankAL, NrankAL)	     	 
    if (PrnProgress) call write_progress_m (.false., m, 3, 5) 
    call LU_SYSTEM_DIRECT (a, 2*NrankAL, 2*NrankAL, b, 2*NrankAL, 2*NrankAL,        &
         2*NmaxC, 2*NmaxL)              	 
    if (PrnProgress) call write_progress_m (.false., m, 4, 5) 
    call matrix_Q_DS_COMP (TypeGeom, 1, 1, k, ind_ref, Nsurfmax, surf, Npart,       &
         NrankpmaxL, Nrankp, zReL, zImL, zpart, m, NmaxC, Nrank, NmaxL, Nint,          &
         Nparammax, Nparam, Nintparam, paramG, weightsG, a, NrankAL, NrankAL)  	     
    call product_matrices (2*NmaxL, 2*NmaxC, 2*NmaxL, a, 2*NrankAL, 2*NrankAL,      &
         b, 2*NrankAL, 2*NrankAL)
    if (PrnProgress) call write_progress_m (.false., m, 5, 5) 
    call write_FileTmat (NrankAL, NrankAL, a)    
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,     &
         NmaxL, c)
    call product_matrix_vector (2*NmaxL, 2*NmaxL, a, 2*NrankAL, 2*NrankAL, c, c1)
    call extend_vector_positive (c1, cc, m, Mstart, Nrank, NmaxL, Nmaxmax)                                   
    if (m /= 0) then
      call matrix_m_negativ (NmaxL, NmaxL, a, NrankAL, NrankAL)
      call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,  &
           NmaxL, c)
      call product_matrix_vector (2*NmaxL, 2*NmaxL, a, 2*NrankAL, 2*NrankAL, c, c1)
      call extend_vector_negative (c1, cc, m, Nrank, NmaxL, Nmaxmax)
    end if  
    call DSCS (cc, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm, &
        .false.,.true., h, v)    
    call delta_DSCS (Nteta, h, v, oldh, oldv, epsMrank, NthetaConv)    
    call write_DSCS (Nteta,.false., h, v)        
    deallocate (zReL, zImL)        
    if (NthetaConv >= int(0.8d0*Nteta)) exit
  end do 
  close (unit = iTmat)
  call CQscat (cc, Mrank, Nrank, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,           &
       alfap, k, snorm, Cext, Qext)
  call write_Effic (Qscat, Qext)              
  call write_MrankConvRes (NthetaConv, epsMrank)
  if (NthetaConv >= int(0.8d0*Nteta)) then
    print "(/,2x,'Convergence criterion for Mrank is satisfied;')"     
  else
    print "(/,2x,'Convergence criterion for Mrank is not satisfied;')"
  end if 
  call write_InfoFileTmat (FileTmat, Mrank, Nrank, .true., .false., .false.)
  call ScatCharact (k, FileTmat, Mrank, Nrank, .true., .false., .false.) 
  print "(/,2x,'T matrix is stored in ',a50)", FileTmat
  print "(  2x,'The dimensions of the T matrix are given by:')"
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", Nrank 
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank       
  deallocate (a, b, c, c1, cc, h, v, oldh, oldv, paramG, weightsG, Nintparam)
end subroutine convergence_MrankDSCOMP
! **********************************************************************************   
subroutine convergence_NrankDSCOMP_m (m, TypeGeom, k, ind_ref, snorm, Nsurfmax,     &
           surf, Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zRe, zIm, Nrankpmax1, &
           Nrankp1, zRe1, zIm1, zpart, Nrank, Nint, Nintparam, paramG, weightsG,    &                                             
           epsNrank)         
  use parameters 
  implicit none
  integer    :: m, TypeGeom, Nsurfmax, Nparammax, Npart, Nrankpmax, Nrankpmax1,     &
                Nrank, Nint, Nrankp(Npart), Nrankp1(Npart), Nparam(Npart),          &
                Nintparam(Npart,Nparammax)                
  real(O)    :: k, surf(Npart,Nsurfmax), zpart(Npart), snorm, epsNrank,             &
                zRe(Npart,Nrankpmax), zIm(Npart,Nrankpmax),                         &
                zRe1(Npart,Nrankpmax1), zIm1(Npart,Nrankpmax1),                     &
                paramG(Npart,Nparammax,Nint), weightsG(Npart,Nparammax,Nint)                
  complex(O) :: ind_ref(Npart)  
!        
  integer    :: Mstart, Mrank, Nteta, NmaxC, NmaxL, Nmaxmax, i, NthetaConv,         &
                ipart, dimension_row, NrankAL
  real(O)    :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,        &
                Cext, Qext, teta  
  real(O),allocatable    :: h(:), v(:), oldh(:), oldv(:)
  complex(O),allocatable :: a(:,:), b(:,:), c(:), c1(:), cc(:)
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = Pi / 4._O
  beta   = Pi / 4._O
  gama   = 0._O
  alfap  = Pi / 4._O
  Mstart = m
  Mrank  = m    
  NmaxC = 0
  do ipart = 1, Npart
    NmaxC = NmaxC + Nrankp(ipart)
  end do
  NmaxL   = dimension_row (m, Nrank)
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  NrankAL = max(NmaxC,Nrank)  
  allocate (a(2*NrankAL,2*NrankAL), b(2*NrankAL,2*NrankAL)) 
  allocate (c(2*Nrank), c1(2*Nrank))        
  allocate (cc(2*Nmaxmax)) 
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))                   
! --- Nrank configuration ---  
  call matrix_Q_DS_COMP (TypeGeom, 3, 1, k, ind_ref, Nsurfmax, surf, Npart,         &
       Nrankpmax, Nrankp, zRe, zIm, zpart, m, NmaxC, Nrank, NmaxC, Nint, Nparammax, &
       Nparam, Nintparam, paramG, weightsG, a, NrankAL, NrankAL) 
  call incident_matrix_DS_COMP (TypeGeom, k, Nsurfmax, surf, Npart, Nrankpmax,      &
       Nrankp, zRe, zIm, zpart, m, NmaxC, Nrank, NmaxL, Nint, Nparammax, Nparam,    &
       Nintparam, paramG, weightsG, b, NrankAL, NrankAL)
  call LU_SYSTEM_DIRECT (a, 2*NrankAL, 2*NrankAL, b, 2*NrankAL, 2*NrankAL,          &
       2*NmaxC, 2*NmaxL)               
  call matrix_Q_DS_COMP (TypeGeom, 1, 1, k, ind_ref, Nsurfmax, surf, Npart,         &
       Nrankpmax, Nrankp, zRe, zIm, zpart, m, NmaxC, Nrank, NmaxL, Nint, Nparammax, &
       Nparam, Nintparam, paramG, weightsG, a, NrankAL, NrankAL)        
  call product_matrices (2*NmaxL, 2*NmaxC, 2*NmaxL, a, 2*NrankAL, 2*NrankAL,        &
       b, 2*NrankAL, 2*NrankAL)
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,       &
       NmaxL, c)
  call product_matrix_vector (2*NmaxL, 2*NmaxL, a, 2*NrankAL, 2*NrankAL, c, c1)
  call extend_vector_positive (c1, cc, m, Mstart, Nrank, NmaxL, Nmaxmax)                                     
  call matrix_m_negativ (NmaxL, NmaxL, a, NrankAL, NrankAL)
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,      &
       NmaxL, c)
  call product_matrix_vector (2*NmaxL, 2*NmaxL, a, 2*NrankAL, 2*NrankAL, c, c1)
  call extend_vector_negative (c1, cc, m, Nrank, NmaxL, Nmaxmax)
  call DSCS (cc, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm,   &
      .false.,.true., h, v)
  call CQscat (cc, Mrank, Nrank, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,           &
       alfap, k, snorm, Cext, Qext)
  do i = 1, Nteta
    oldh(i) = h(i)
    oldv(i) = v(i)
  end do
  print "(1x,'Convergence Test over Nrank')"
  print "(1x,'---------------------------')"  
  print "(1x,'Nint = ',i5,',',1x,'Nrank = ',i3,',',1x,'Mrank = ',i3)", Nint,        &
          Nrank, m  
  do ipart = 1, Npart
    print "(1x,'Nrank for region ',i2,', Nrank = ',i3)", ipart, Nrankp(ipart)  
  end do          
  print "(/,1x, a, /, 2x, a, 9x, a, 8x, a)",                                        &
 'normalized differential scattering cross section', 'theta', 'parallel',           &
 'perpendicular' 
  do i = 1, Nteta    
    teta = real(i - 1,O) * 180._O / real(Nteta - 1,O)                
    print "(1x,f6.2,5x,1pe13.4,5x,1pe13.4)", teta, h(i), v(i)   
  end do  
  print "(1x,'scattering efficiency = ',1pe13.4)",   Qscat
  print "(1x,'extinction efficiency = ',1pe13.4,/)", Qext          
! --- Nrank1 configuration (only the DS configuration is changed) ---
  NmaxC = 0
  do ipart = 1, Npart
    NmaxC = NmaxC + Nrankp1(ipart)
  end do    
  call matrix_Q_DS_COMP (TypeGeom, 3, 1, k, ind_ref, Nsurfmax, surf, Npart,         &
       Nrankpmax1, Nrankp1, zRe1, zIm1, zpart, m, NmaxC, Nrank, NmaxC, Nint,        &
       Nparammax, Nparam, Nintparam, paramG, weightsG, a, NrankAL, NrankAL)
  call incident_matrix_DS_COMP (TypeGeom, k, Nsurfmax, surf, Npart, Nrankpmax1,     &
       Nrankp1, zRe1, zIm1, zpart, m, NmaxC, Nrank, NmaxL, Nint, Nparammax,         &
       Nparam, Nintparam, paramG, weightsG, b, NrankAL, NrankAL)
  call LU_SYSTEM_DIRECT (a, 2*NrankAL, 2*NrankAL, b, 2*NrankAL, 2*NrankAL,          &
       2*NmaxC, 2*NmaxL)              
  call matrix_Q_DS_COMP (TypeGeom, 1, 1, k, ind_ref, Nsurfmax, surf, Npart,         &
       Nrankpmax1, Nrankp1, zRe1, zIm1, zpart, m, NmaxC, Nrank, NmaxL, Nint,        &
       Nparammax, Nparam, Nintparam, paramG, weightsG, a, NrankAL, NrankAL)        
  call product_matrices (2*NmaxL, 2*NmaxC, 2*NmaxL, a, 2*NrankAL, 2*NrankAL,        &
       b, 2*NrankAL, 2*NrankAL)    
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,       &
       NmaxL, c)
  call product_matrix_vector (2*NmaxL, 2*NmaxL, a, 2*NrankAL, 2*NrankAL, c, c1)
  call extend_vector_positive (c1, cc, m, Mstart, Nrank, NmaxL, Nmaxmax)                                    
  call matrix_m_negativ (NmaxL, NmaxL, a, NrankAL, NrankAL)
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,      &
       NmaxL, c)
  call product_matrix_vector (2*NmaxL, 2*NmaxL, a, 2*NrankAL, 2*NrankAL, c, c1)
  call extend_vector_negative (c1, cc, m, Nrank, NmaxL, Nmaxmax)
  call DSCS (cc, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm,   &
      .false.,.true., h, v)
  call CQscat (cc, Mrank, Nrank, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,           &
       alfap, k, snorm, Cext, Qext)
  call delta_DSCS (Nteta, h, v, oldh, oldv, epsNrank, NthetaConv)
  print "(1x,'Nint = ',i5,',',1x,'Nrank = ',i3,',',1x,'Mrank = ',i3)", Nint,        &
          Nrank, m  
  do ipart = 1, Npart
    print "(1x,'Nrank for region ',i2,', Nrank = ',i3)", ipart, Nrankp1(ipart)  
  end do          
  print "(/,1x, a, /, 2x, a, 9x, a, 8x, a)",                                        &
 'normalized differential scattering cross section', 'theta', 'parallel',           &
 'perpendicular' 
  do i = 1, Nteta    
    teta = real(i - 1,O) * 180._O / real(Nteta - 1,O)                
    print "(1x,f6.2,5x,1pe13.4,5x,1pe13.4)", teta, h(i), v(i)   
  end do  
  print "(1x,'scattering efficiency = ',1pe13.4)",   Qscat
  print "(1x,'extinction efficiency = ',1pe13.4,/)", Qext 
  print "(1x, a, i2, a, 1f5.2, a, /)",                                               &
 'The solution converges in ', NthetaConv,' points with an relative error of ',      &
  100 * epsNrank, ' %;'              
  deallocate (a, b, c, c1, cc, h, v, oldh, oldv)
end subroutine convergence_NrankDSCOMP_m 
