subroutine TLAY
!------------------------------------------------------------------------------------
! 1. General Considerations                                                         !
! --------------------------                                                        !
! TLAY is a routine for computing the T matrix and the scattering characteristics   !
! of axisymmetric layered particles. The particle consists of Npart layers, each    !
! characterized by arbitrary but constant values of electric permittivity and       !
! magnetic permeability. Each layer is bounded by two surfaces. The enclosing       !
! (exterior) surface will be referred to as the layer surface, and obviously, the   !
! first layer surface is the largest surface.                                       ! 
!                                                                                   !
! The domain of analysis is divided into homogeneous layers. The electric and       !
! magnetic fields of a specific layer i can be approximated by using localized      !
! or distributed vector spherical wave functions, and Nrankp(i) will be referred    !
! to as the maximum expansion order of the layer i.                                 !
!                                                                                   !
! The transition matrix of a layered particle can be obtained by considering the    !
! null-field condition for the total electric field inside and outside each region. !
! Note that for axisymmetric layered particles, the scattering problem decouples    !
! over the azimuthal modes m and the T matrix can be computed separately for each m.!
!                                                                                   !
! Particle geometries currently supported include layered spheroids and cylinders.  !
! The following routines (from the file "GeomLib.f90") provide the required         !
! geometry parameters:                                                              !
! - the routine "interpolation_listLAY" computes the integration nodes and weights  !
!   for each layer generatrix,                                                      !
! - the routine "elem_geomLAY" provides the geometry parameters (position vector,   !
!   polar and azimuthal angles and normal unit vector in spherical coordinates)     !
!   at each integration node, and                                                   !
! - the routine "zDSLAY" generates the coordinates of the distributed sources.      !
! The generatrix of each layer is described with respect to a local coordinate      !
! system and the axial positions of the local coordinate systems are specified in   !
! the global coordinate system of the layered particle by the array OriginPart.     !
! The global coordinate system can be chosen as the local coordinate system of the  !
! first layer surface (host particle) by setting OriginPart(1) = 0.0. By convention,! 
! the COORDINATES OF THE DISTRIBUTED SOURCES corresponding to a specific layer are  !
! provided in the GLOBAL COORDINATE SYSTEM. The user can easily modify the above    !
! routines to generate particles with other geometries. Note that the list of       !
! parameters must be maintained and only geometries with an analytical description  !
! of the surface can be implemented.                                                !
!                                                                                   !
! 2. Convergence Test                                                               !
! --------------------                                                              !
! Convergence tests precede the scattering characteristics calculation. Three       !
! parameters control the T-matrix computation:                                      !
! - the number of integration points Nint,                                          !
! - the maximum expansion order Nrank, and                                          !
! - the maximum azimuthal order Mrank.                                              !
!                                                                                   !
! Before going any further, we will explain the significance of the parameters      !
! Nrank and Nint. Let Nrankp(i) be the maximum expansion order of the layer i,      !
! i = 1,2,...,Npart, where Npart is the number of layers. Nrank is defined as:      !
!                                                                                   !
!                             Nrank = Nrankp(1),                                    !
!                                                                                   !
! i.e., Nrank is the maximum expansion order corresponding to the host particle.    !
! The parameter Nrank determines the dimension of the T matrix and by convention,   !          
! Nrank will be referred to as the maximum expansion order of the layered particle. !
! Analogously, we define Nint as the number of integration points on the first      !
! generatrix, and Nint will be referred to as the global number of integration      !
! points for the layered particle. The number of integration points on the layer    !
! generatrix i is                                                                   !
!                                                                                   !
!                    Nint_layer(i) = l(i) * Nint / l(1),                            !
!                                                                                   !
! where l(i) is a characteristic length of the layer i. Each layer generatrix i     !
! can be divided into NparamPart(i) piecewise smooth curves and the number of       !
! integration points on the smooth curve j is                                       !
!                                                                                   !
!               Nint_layer_curve(i,j) = int { c(j) * Nint_layer(i) },               !
!                                                                                   !
! where c(j) is a multiplicative factor. The parameters c(j), do not appear as      !
! explicit parameters; they are hard coded in the routine "interpolation_listLAY"   !
! from the file "GeomLib.f90". For the supplied geometries,                         !
!                                                                                   !
!                         l(i) = max {a(i),b(i)},                                   !
!                                                                                   !
! where a(i) and b(i) are the semi-axes of the spheroidal surface i, or the         !
! half-length and the radius of the cylindrical surface i. The generatrix of a      !
! cylinder is divided into three straight lines. The first and third lines are      !
! perpendicular to the axis of symmetry, while the second line is parallel to the   !
! axis of symmetry. Denoting by a the half-length  of the cylinder and by b the     !
! cylinder radius, we have:                                                         !
!                                                                                   !
!                 Nint_layer_curve(cylinder,first line)                             !
!               = Nint_layer_curve(cylinder,third line)                             !
!               = int {[b/2/(a + b)] * Nint_layer(cylinder)}                        !
!                                                                                   !
! and                                                                               !
!                                                                                   !  
!               Nint_layer_curve(cylinder, second line)                             ! 
!      = Nint_layer(cylinder) - 2 * int {[b/2/(a + b)] * Nint_layer(cylinder) }.    !
!                                                                                   !
! For convergence tests, the incident wave is assumed to be a plane wave traveling  !
! along the Z-axis of the global coordinate system and the scattering               !
! characteristics  are computed in the azimuthal plane  phi = 0°. The convergence   !
! tests over Nint and Nrank are interactive, while the convergence test over Mrank  !
! is automatically performed.                                                       !  
!                                                                                   !
! For the integration and expansion order test, an axisymmetric orientation of the  !
! particle is considered, i.e., the incident wave propagates along the axis of      !
! symmetry of the particle. In this case, all azimuthal modes are zero except       !
! for m = - 1 and m = 1. For the convergence test over Nint, the scattering problem !
! is solved for Nint and Nint + dNint, while for the convergence test over Nrank,   !
! the scattering problem is solved for                                              !
!                                                                                   !
!                 Nrankp(i), i = 1,2,...,Npart,                                     !
!                                                                                   !
! and                                                                               !
!                                                                                   !
!                 Nrankp(i) - 1, i = 1,2,...,Npart.                                 !
!                                                                                   !  
! The normalized differential scattering cross section (DSCS) will be checked at    !
! 20° increments for convergence within epsX (epsNint or epsNrank) tolerance. If    !
! the calculated results converge within this tolerance at 80% of the scattering    !
! angles, then convergence is achieved.                                             !
!                                                                                   ! 
! After Nrank and Nint have been determined we pass to the azimuthal order test.    !
! The program automatically sets the particle to a more general orientation, i.e.,  !
! alpha = beta = 45°, and solves the scattering problem  for increasing m values    !
! until convergence of the angular scattering is achieved. The T matrix is stored   !
! for later use by other programs, and the values of Nrank and Mrank are printed to !
! the screen and to the T-matrix information file.  These values together with the  !
! T matrix serve as INPUT PARAMETERS for other programs.                            !
!                                                                                   !
! 3. Estimates of Nint and Nrank                                                    !
! -------------------------------                                                   !
! The above convergence tests require estimates for the number of integration       !
! points (Nint), and the maximum expansion orders for each layer (Nrankp(i),        !
! i = 1,2,...,Npart). These estimates are user-defined.                             !
!                                                                                   !
! A reasonable starting value of Nrankp is given by Wiscombe's truncation limit     !
! criterion [W. J. Wiscombe, Improved Mie scattering algorithms, Applied Optics,    !
! 19, 1505-1509, 1980], i.e.,                                                       !
!                                                                                   !                       
!                NrankpW = int(xp + 4.05 * xp**0.33 + 2),                           !
!                                                                                   !
! where xp is the size parameter of the layer, xp = k * lnormPart, k is the wave    !
! number and lnormPart is the radius of the smallest sphere circumscribing the      !
! layer.                                                                            !
!                                                                                   !
! The value of Nint depends on the type of discrete sources, the size parameter,    !
! and the optical properties of the layered particle. If                            !
!                                                                                   ! 
!                        Nrank = Nrankp(1),                                         !
!                                                                                   !
! we propose the following estimate                                                 !
!                                                                                   !
!                        Nint = Ndgs * Nrank,                                       !
!                                                                                   !
! where Ndgs = 15, 20, 25, .... For distributed sources, the value of Nint should   !
! be substantially larger than the value corresponding to localized sources.        !
!                                                                                   !
! 4. Strategies for Performing Convergence Tests                                    !
! -----------------------------------------------                                   !
! Two strategies for performing convergence tests over Nint and Nrank are           !
! recommended. The following technique is an extension of the approach used for     !
! homogeneous axisymmetric particles:                                               !
! 1. for each layer i, choose a value of Nrankp(i) close to the value predicted     !
!    by Wiscombe's truncation limit criterion;                                      !
! 2. for the set Nrankp(i), perform the expansion order test with                   !
!    Nint = Ndgs * Nrank, where e.g., Ndgs = 15;                                    !
! 3. if convergence is not achieved, set Nrankp(i) = Nrankp(i) + 1 for all layers   !
!    i, maintain the relation Nint = Ndgs * Nrank and go to Step 2;                 !
! 4. if convergence is achieved, perform several integration tests with the initial !
!    Nint-value, Nint = Ndgs * Nrank, until the DSCS converges.                     !
! At Step 3, it is not necessary to increase all Nrankp-values (for all layers).    !
! If the user has the feeling that some Nrankp-values are correct, these values     !
! need not to be changed.                                                           !  
!                                                                                   !
! For large and/or highly aspherical layered particles, T-matrix computations       !
! become poorly covergent or divergent. Frequently, a semi-convergent behaviour is  !
! attended: the errors decrease with increasing the maximum expansion order,        !
! attain a relative constant level and afterwards increase. The region of           !
! stability can be localized by using the following technique:                      !
! 1. for each layer i, choose a value of Nrankp(i) smaller than the value           !
!    predicted by Wiscombe's truncation limit criterion;                            !
! 2. for the set Nrankp(i), perform several integration tests with the initial      !
!    Nint-value, Nint = Ndgs * Nrank, until the DSCS converges;                     !
! 3. if convergence is not achieved, set Nrankp(i) = Nrankp(i) - 1 for all layers   !
!    i, or increase the tolerance epsNint, and go to Step 2;                        !
! 4. if convergence is achieved (for some Nint), perform several expansion order    !
!    tests for increasing Nrankp-values and monitor the errors of the               !
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
! logical variable DoConvTest to false. In this case, the values of Nint and        !
! Nrank must be specified in the input file.                                        !
!                                                                                   !
! Nint and Nrank have been determined for the azimuthal orders m = 1 and m = - 1.   ! 
! For higher azimuthal orders and sources distributed in the complex plane, the     !
! chosen values of Nint and Nrank may lead to false convergence. There are two      !
! strategies for improving the numerical stability.                                 !
! 1. The number of integration points can be increased at each azimuthal mode       !
! calculation with the integration step dNintMrank. The parameter dNintMrank         !
! should be found experimentally by repeating the azimuthal order test for          !
! different input values until the differential scattering cross section converges. !
! A conservative estimate is: dNintMrank = (0.1,...,0.2) * Nint.                    !
! 2. For sources generated automatically, an interactive convergence test over      !
! Nint and Nrank can be performed at each azimuthal mode calculation. The general   !
! strategy is to reduce the number of discrete sources (Nrankp(i), i = 1,2,...,     !
! Npart) and to increase the number of integration points (Nint) for increasing     !
! values of m. However, in the current version of the program, this test is not     !
! active (see the comments in the subroutine "convergence_MrankDSLAY").             !
!                                                                                   !
! The number of azimuthal modes is the same for all layers, and for localized       !
! sources, the relation Mrank <= Nrankp(i), for all i = 1,2,...,Npart, should be    !
! satisfied. However, the Nrankp(i), i = 1,2,...,Npart, are determined before Mrank !
! is computed. Therefore, for the azimuthal mode calculation m, with |m| > 1,       !
! we set Nrankp(i) = |m|, if Nrankp(i) < |m| for some i. In this case, the maximum  !
! expansion order of a layer is m-dependent and                                     !
!                                                                                   !
!                    n = |m|,|m|+1,...,Nrankp(i,|m|).                               !
!                                                                                   !
! 6. Input Parameters                                                               !
! ---------------------                                                             !
! The parameters specified in the input file "/INPUTFILES/InputLAY.dat" are listed  !
! below                                                                             !
!                                                                                   !
! - wavelength (real) - wavelength of the incident light in vacuo.                  !
!                                                                                   !
! - ind_refMed (real) - refractive index of the ambient medium (nonabsorbing        !
!   dielectric medium).                                                             !
!                                                                                   !
! - TypeGeom (integer) - parameter specifying the type of the particle geometry.    !
!                                                                                   !
! - Npart (integer) - number of layers.                                             !
!                                                                                   !
! - anorm (real) - characteristic length of the layered particle which is used      !
!   to normalize the differential scattering cross sections.                        !
!                                                                                   !
! - DoConvTest (logical) - if DoConvTest = t, the interactive convergence tests     !
!   over Nint and Nrank are invoked. Estimates of Nrank for all layers are          !
!   provided by Wiscombe's truncation limit criterion. If DoConvTest = f, the       !
!   values of Nint and Nrank must be supplied in the input file.                    !
!                                                                                   !
! - DS (logical) - if DS = t, distributed sources are used for T-matrix             !
!   calculation, otherwise localized sources are employed.                          !
!                                                                                   !
! - autGenDS (logical) - if autGenDS = t, the routine "zDSLAY" generates the        !
!   coordinates of the distributed sources. Otherwise, the coordinates of the       !
!   distributed sources must be provided by the user. If DS = f, the logical        ! 
!   variable autGenDS is ignored.                                                   ! 
!                                                                                   !
! - Nint (integer) - global number of integration points for the layered            !
!   particle. This parameter is used if the convergence tests are not performed     !
!   (DoConvTest = f).                                                               !
!                                                                                   !
! The next parameters (specified in a sequence of group statements) characterize     !
! each layer.                                                                       !
! - ind_refPartRel (complex) - relative refractive index of the actual layer        !
!   with respect to the ambient medium. The imaginary part of the relative          !
!   refractive index must be zero for nonabsorbing particles and positive for       !
!   absorbing particles.                                                            !
!                                                                                   !
! - NsurfPart (integer) - number of surface parameters of the actual layer.         !
!                                                                                   !
! - NparamPart (integer) - number of smooth curves forming the generatrix curve of  !
!   the actual layer.                                                               !
!                                                                                   !
! - OriginPart (real) - axial position of the local coordinate system of the        !
!   actual layer with respect to the global coordinate system of the layered        !
!   particle. This parameter can be positive or negative.                           !
!                                                                                   !
! - surfPart (real array: surfPart(1), surfPart(2),...,surfPart(NsurfPart)) -       !
!   surface parameters specifying the shape of the actual layer. The dimension of   !
!   the array is NsurfPD. The integer parameter NsurfPD is specified in the         !
!   routine  "Parameters.f90" and has the value NsurfPD = 10. If                    !
!   NsurfPart > NsurfPD, the execution is automatically terminated.                 !
!                                                                                   !
!   The permissive values of the surface parameters (for the supported geometries)  !
!   are summarized below.                                                           !
!                                                                                   !
!      Particle      TypeGeom   Nsurf   Nparam                surf                  ! 
!     concentric        1         2        1        surf(1) - length of the semi-   ! 
!     spheroids                                               axis along the        ! 
!                                                             symmetry axis         !  
!                                                   surf(2) - length of the second  ! 
!                                                             semi-axis             !     
!                                                                                   ! 
!     concentric        2         2        3        surf(1) - half-length of        ! 
!     cylinders                                               the cylinder          !  
!                                                   surf(2) - cylinder radius       ! 
!                                                                                   !
! - lnormPart (real) - characteristic length of the actual layer (usually the       !
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
! - EpsZReImPart (real) - input parameter of the routine "zDSLAY" which controls    !
!   the distribution of the discrete sources. This parameter is used if distributed !
!   sources are required and the coordinates of the distributed sources are         !
!   automatically generated (DS = t and autGenDS = t).                              !
!                                                                                   !
! - NrankPart (integer) - maximum expansion order for the actual layer. This        !
!   parameter is used if the convergence tests are not performed, or the coordinates!
!   of the distributed sources are user-defined. More specifically, NrankPart is    !
!   used if                                                                         !
!   - (DS = f, DoConvTest = f),                                                     !
!   - (DS = t, autGenDS = t, DoConvTest = f), or                                    ! 
!   - (DS = t, autGenDS = f, DoConvTest = t/f).                                     !
!                                                                                   !
! - zRePart, zImPart (real arrays: zXPart(1), zXPart(2),...,zXPart(Nrankp),         !
!   X = Re, Im) - coordinates of the distributed sources for the actual layer and   !
!   the expansion order NrankPart. These parameters are used if the coordinates of  !
!   the distributed sources are user-defined (DS = t and autGenDS = f). The         !
!   dimension of the arrays zRePart and zImPart is NrankPD and the inequality       !
!   NrankPart <= NrankPD must hold. The integer parameter NrankPD is specified in   !
!   the routine "Parameters.f90" and has the value NrankPD = 200. If                !
!   NrankPart > NrankPD, the execution is automatically terminated. Note that the   !
!   coordinates of the distributed sources are defined with respect to the GLOBAL   !
!   COORDINATE SYSTEM of the layered particle.                                      !
!                                                                                   !
! - zRePart1, zImPart1 (real arrays: zXPart1(1), zXPart1(2),...,zXPart1(Nrankp-1),  !
!   X = Re, Im) - coordinates of the distributed sources for the actual layer and   !
!   the expansion order NrankPart - 1. These parameters are used if the             !
!   coordinates of the distributed sources are user-defined and the expansion order !
!   test is performed (DS = t and autGenDS = f and TypeConvTest = 2). The dimension !
!   of the arrays zRePart1 and zImPart1 is NrankPD. As before, zRePart1 and zImPart1!
!   are defined with respect to the GLOBAL COORDINATE SYSTEM of the layered         !
!   particle.                                                                       !
!                                                                                   !
! NOTE: THE INPUT ARRAYS zRe, zIm AND zRe1, zIm1 MUST BE SPECIFIED IF (DS = t AND   !
! autGenDS = f), AND MUST BE SEPARATED BY A BLANK LINE. IF THE EXPANSION ORDER TEST !
! IS NOT PERFORMED (TypeConvTest /= 2), THE INPUT ARRAYS zRePart1 AND zImPart1 CAN  !
! BE OMITTED.                                                                       !
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
!   (DS = t) and THERE EXISTS AT LEST ONE LAYER with ComplexPlanePart = t). At      !
!   each azimuthal mode calculation, the number of integration points increases     !
!   with dNintMrank. For layers with sources distributed along the symmetry axis,   !
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
! NrankReg and SourceRegPosInp characterize a specific layer. THESE STATEMENTS      !
! MUST BE REPEATED FOR ALL Npart LAYERS.                                            !
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
!          < the code computes an estimate of Nrankp(i) for each layer i            !
!            by using Wiscombe's truncation criterion, and prompts for the          !
!            input value of Nrankp(i) >                                             !
!          < the code prompts for the estimated value of Nint for the layered       !
!            particle>                                                              !
!          < the code prompts for the type of convergence test:                     !
!            1 - over Nint, 2 - over Nrank or 3 - over Mrank >                      !
!       else if ( .not. do convergence test) then                                   !
!          < the input file provides the values of Nrankp(i) for each layer i >     !
!          < the input file provides the value of Nint for the layered particle>    !
!          type of convergence test = 3 (convergence test over Mrank)               !
!       end if                                                                      !
!                                                                                   !
!    else if ( use distributed sources ) then                                       !
!                                                                                   !
!       if ( use automatic sources generation ) then                                !
!          if ( do convergence test) then                                           !
!             < the code computes an estimate of Nrankp(i) for each layer i         !
!               by using Wiscombe's truncation criterion, and prompts for the       !
!               input value of Nrankp(i) >                                          !
!             < the code generates the coordinates of the distributed sources for   !
!               Nrank and all layers >                                              !
!             < the code prompts for the estimated value of Nint for the layered    !
!               particle>                                                           !
!             < the code prompts for the type of convergence test:                  !
!               1 - over Nint, 2 - over Nrank or 3 - over Mrank >                   !
!             if ( do maximum expansion order test ) then                           !
!                < the code generates the coordinates of the distributed sources    !
!                  for Nrank1 and all layers >                                      !
!             end if                                                                !
!          else if ( .not. do convergence test) then                                !
!             < the input file provides the values of Nrankp(i) for each layer i>   !
!             < the code generates the coordinates of the distributed sources for   !
!               Nrank and all layers >                                              !
!             < the input file provides the value of Nint for the layered particle> !
!             type of convergence test = 3 (convergence test over Mrank)            !
!       else if ( .not. use automatic sources generation ) then                     !
!          if ( do convergence test) then                                           !
!             < the input file provides the values of Nrankp(i) for each layer i >  !
!             < the input file provides the coordinates of the distributed sources  !
!               for Nrank and all layers >                                          !
!             < the code prompts for the estimated value of Nint for the layered    !
!               particle>                                                           !
!             < the code prompts for the type of convergence test:                  !
!               1 - over Nint, 2 - over Nrank or 3 - over Mrank >                   !
!             if ( do maximum expansion order test ) then                           !
!                < the input file provides the coordinates of the distributed       !
!                  sources for Nrank1 and all layers >                              ! 
!             end if                                                                ! 
!          else if ( .not. do convergence test) then                                !
!             < the input file provides the values of Nrankp(i) for each layer i >  !
!             < the input file provides the coordinates of the distributed sources  !
!               for Nrank and all layers >                                          !
!             < the input file provides the value of Nint for the layered particle> !
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
!       < the code computes the DSCS for Nrank and Nrank1 and                       !
!         write the results to the file "Output.dat" >                              !
!       < the execution is terminated >                                             !
!    end if                                                                         !
!    if ( do azimuthal order test, i.e., type of convergence test = 3 ) then        !
!       < the code computes the DSCS for increasing m values and                    !
!         write the results to the file "Output.dat"  >                             !
!       if ( use distributed sources .and. use automatic sources generation         !
!           .and. sources are distributed in the complex plane ) then               ! 
!          < an interactive convergence test over Nint and Nrank can be performed   !
!            at each azimuthal mode m >                                             !
!          Note: to make this test active comment out some lines in the subroutine  !
!          "convergence_MrankDSLAY".                                                ! 
!       end if                                                                      !
!       < the T matrix is stored in the file FileTmat and the T-matrix information  !
!         file is created >                                                         !
!       < Mrank and Nrank are printed to the screen and to the T-matrix             !
!         information file >                                                        !
!       < go to Step 3 >                                                            !
!    end if                                                                         !
! 3. < scattering characteristics calculation >                                     !
!------------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  use allocation, only: Nsurf, Nparam, Nrankp, Nrankp1, surf, zpart, zRe, zIm,      &
                        zRe1, zIm1, lnorm, ind_ref, ComplexPlane, EpsZReIm 
  implicit none 
  integer       :: TypeGeom, Npart, dNint, TypeConvTest, Nsurfmax, Nparammax,       &
                   Nrankpmax, Nrankpmax1, Nint, dNintMrank 
  real(O)       :: k, ind_refMed, wavelength, anorm, snorm, epsNint, epsNrank,      &
                   epsMrank                                           
  logical       :: DoConvTest, DS, autGenDS, PrnProgress 
  character(80) :: FileTmat  
! -----------------------------------------------------------------------------------
!                            Read the input file                                    ! 
! -----------------------------------------------------------------------------------       
  call readinputLAY ( wavelength, ind_refMed, TypeGeom, Npart, anorm,               &
       DoConvTest, DS, autGenDS, Nint, epsNint, epsNrank, epsMrank, dNint,          &
       dNintMrank, FileTmat, PrnProgress, k, snorm, Nsurfmax, Nparammax,            &
       Nrankpmax, Nrankpmax1, TypeConvTest )
! -----------------------------------------------------------------------------------
!                                      Main                                         !
! -----------------------------------------------------------------------------------   
  open (unit = iOutput, file = FileOutput, status = "replace") 
  call printinputLAY (TypeConvTest, TypeGeom, Nsurfmax, Nsurf, surf, Nparam, Npart, &
       lnorm, anorm, Nrankpmax, Nrankp, zpart, zRe, zIm, Nrankpmax1, Nrankp1, zRe1, &
       zIm1, ind_ref, dNintMrank, dNint, wavelength, ind_refMed, epsNint, epsNrank, &
       epsMrank, DS, autGenDS)  
  if (DoConvTest) then
    if (TypeConvTest == 1) then
      if (.not. DS) then
        call convergence_NintLAY (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,      &
             Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zpart, Nint, dNint,       &
             epsNint, PrnProgress)      
      else 
        call convergence_NintDSLAY (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,    &
             Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zRe, zIm, zpart, Nint,    &
             dNint, epsNint, PrnProgress)      
      end if
    else if (TypeConvTest == 2) then
      if (.not. DS) then
        call convergence_NrankLAY (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,     &
             Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zpart, Nint, epsNrank,    &
             PrnProgress)          
      else 
        call convergence_NrankDSLAY (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,   &
             Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zRe, zIm, Nrankpmax1,     &
             Nrankp1, zRe1, zIm1, zpart, Nint, epsNrank, PrnProgress)      
      end if
    else 
      if (.not. DS) then
        call convergence_MrankLAY (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,     &
             Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zpart, Nint, epsMrank,    &
             FileTmat, PrnProgress)          
      else 
        call convergence_MrankDSLAY (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,   &
             Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zRe, zIm, zpart, Nint,    &
             ComplexPlane, EpsZReIm, autGenDS, epsMrank, epsNrank, dNintMrank,      &
             FileTmat, PrnProgress)                              
      end if
    end if 
  else
    if (.not. DS) then
      call convergence_MrankLAY (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,       &
           Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zpart, Nint, epsMrank,      &
           FileTmat, PrnProgress)            
    else 
      call convergence_MrankDSLAY (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,     &
           Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zRe, zIm, zpart, Nint,      &
           ComplexPlane, EpsZReIm, autGenDS, epsMrank, epsNrank, dNintMrank,        &
           FileTmat, PrnProgress)                
    end if  
  end if             
  close (unit = iOutput)           
  deallocate (zRe, zIm, zRe1, zIm1)   
  deallocate (surf, Nsurf, Nparam, Nrankp, Nrankp1, ind_ref, zpart, lnorm) 
  deallocate (ComplexPlane, EpsZReIm)   
end subroutine TLAY
!***********************************************************************************
subroutine readinputLAY ( wavelength, ind_refMed, TypeGeom, Npart, anorm,           &
           DoConvTest, DS, autGenDS, Nint, epsNint, epsNrank, epsMrank, dNint,      &
           dNintMrank, FileTmat, PrnProgress, k, snorm, Nsurfmax, Nparammax,        &
           Nrankpmax, Nrankpmax1, TypeConvTest ) 
  use parameters
  use derived_parameters
  use allocation, only: Nsurf, Nparam, Nrankp, Nrankp1, surf, zpart, zRe, zIm,      &
                        zRe1, zIm1, lnorm, EpsZReIm, ind_ref, ComplexPlane  
  implicit none 
  integer       :: TypeGeom, Npart, NsurfPart, NparamPart, dNint, TypeConvTest,     &
                   i, ipart, isurf, NrankPart, NrankW, Nsurfmax, Nparammax,         &
                   Nrankpmax, Nrankpmax1, Nint, dNintMrank, ios
  real(O)       :: k, ind_refMed, wavelength, anorm, surfPart(NsurfPD), xpart,      &
                   snorm, epsNint, epsNrank, epsMrank, zRePart(NrankPD),            &
                   zImPart(NrankPD), zRePart1(NrankPD), zImPart1(NrankPD),          &
                   OriginPart, x, lnormPart, EpsZReImPart                        
  complex(O)    :: ind_refPartRel
  logical       :: DoConvTest, DS, autGenDS, InputDS, ComplexPlanePart,             &
                   PrnProgress, XFindPar 
  character(80) :: FileTmat, string  
! -----------------------------------------------------------------------------------
!                           Read the input file FileInputLAY                        ! 
! -----------------------------------------------------------------------------------
  call DrvParameters 
  open (unit = iInputLAY, file = FileInputLAY, status = "old", position = "rewind")   
  wavelength = 0.1_O * 2._O * Pi 
  ind_refMed = 1._O  
  string     = 'OptProp'
  if (XFindPar (iInputLAY, string)) then
    read (iInputLAY, *, iostat = ios) wavelength
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable wavelength;')"
      stop
    end if
    read (iInputLAY, *, iostat = ios) ind_refMed
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
  Npart = 2
  anorm = 1._O
  string   = 'GeomProp'
  if (XFindPar (iInputLAY, string)) then
    read (iInputLAY, *, iostat = ios) TypeGeom
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable TypeGeom;')"
      stop
    end if
    read (iInputLAY, *, iostat = ios) Npart
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Npart;')"
      stop
    end if
    read (iInputLAY, *, iostat = ios) anorm
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable anorm;')"
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
  if (XFindPar (iInputLAY, string)) then
    read (iInputLAY, *, iostat = ios) DoConvTest
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
  string   = 'Sources'
  if (XFindPar (iInputLAY, string)) then
    read (iInputLAY, *, iostat = ios) DS
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable DS;')"
      stop
    end if
    read (iInputLAY, *, iostat = ios) autGenDS
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
    if (XFindPar (iInputLAY, string)) then
      read (iInputLAY, *, iostat = ios) ind_refPartRel
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable ind_refPartRel;')"
        print "(  2x,'for the layer ',i3,';')", ipart
        stop
      end if               
    else
      print "(/,2x,'Group name OptRegProp not found;')"
      print "(  2x,'for the layer ',i3,';')", ipart
      stop  
    end if
    call check_ind_ref1 (ipart, ind_refPartRel)
    ind_ref(ipart) = ind_refPartRel
    NsurfPart  = 2
    NparamPart = 1
    OriginPart = 0._O
    do isurf = 1, NsurfPD
      surfPart(isurf) = 1._O
    end do
    lnormPart = 1._O
    string    = 'GeomRegProp'
    if (XFindPar (iInputLAY, string)) then
      read (iInputLAY, *, iostat = ios) NsurfPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable NsurfPart;')"
        print "(  2x,'for the layer ',i3,';')", ipart
        stop
      end if
      if (NsurfPart > NsurfPD) then
        print "(/,2x,'Input error: NsurfPart exceeds NsurfPD for the layer',i3)",   &
                ipart                                    
        stop
      end if      
      read (iInputLAY, *, iostat = ios) NparamPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable NparamPart;')"
        print "(  2x,'for the layer ',i3,';')", ipart
        stop
      end if
      read (iInputLAY, *, iostat = ios) OriginPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable OriginPart;')"
        print "(  2x,'for the layer ',i3,';')", ipart
        stop
      end if
      do isurf = 1, NsurfPart
        read (iInputLAY, *, iostat = ios) surfPart(isurf)
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable surfPart;')"
          print "(  2x,'for the layer ',i3,';')", ipart
          stop
        end if
      end do
      read (iInputLAY, *, iostat = ios) lnormPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable lnormPart;')"
        print "(  2x,'for the layer ',i3,';')", ipart
        stop
      end if
    else
      print "(/,2x,'Group name GeomRegProp not found;')"
      print "(  2x,'for the layer ',i3,';')", ipart
      stop  
    end if    
    Nsurf(ipart) = NsurfPart
    Nparam(ipart)= NparamPart
    zpart(ipart) = OriginPart
    lnorm(ipart) = lnormPart
    if (Nsurf(ipart)  > Nsurfmax)  Nsurfmax  = Nsurf(ipart)
    if (Nparam(ipart) > Nparammax) Nparammax = Nparam(ipart)
  end do
  rewind (unit = iInputLAY)
  allocate (surf(Npart, Nsurfmax))
  do ipart = 1, Npart  
    do isurf = 1, Nsurfmax
      surf(ipart,isurf) = 0._O
    end do   
    NsurfPart  = 2
    NparamPart = 1
    OriginPart = 0._O
    do isurf = 1, NsurfPD
      surfPart(isurf) = 1._O
    end do
    lnormPart = 1._O
    string    = 'GeomRegProp'
    if (XFindPar (iInputLAY, string)) then
      read (iInputLAY, *, iostat = ios) NsurfPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable NsurfPart;')"
        print "(  2x,'for the layer ',i3,';')", ipart
        stop
      end if
      read (iInputLAY, *, iostat = ios) NparamPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable NparamPart;')"
        print "(  2x,'for the layer ',i3,';')", ipart
        stop
      end if
      read (iInputLAY, *, iostat = ios) OriginPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable OriginPart;')"
        print "(  2x,'for the layer ',i3,';')", ipart
        stop
      end if
      do isurf = 1, NsurfPart
        read (iInputLAY, *, iostat = ios) surfPart(isurf)
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable surfPart;')"
          print "(  2x,'for the layer ',i3,';')", ipart
          stop
        end if
      end do
      read (iInputLAY, *, iostat = ios) lnormPart
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable lnormPart;')"
        print "(  2x,'for the layer ',i3,';')", ipart
        stop
      end if
    else
      print "(/,2x,'Group name GeomRegProp not found;')"
      print "(  2x,'for the layer ',i3,';')", ipart
      stop  
    end if
    do isurf = 1, NsurfPart
      surf(ipart,isurf) = surfPart(isurf)
    end do
  end do  
  call check_geomLAY (TypeGeom, Npart, Nsurf, Nparam)
!
  allocate (ComplexPlane(Npart), EpsZReIm(Npart))
  rewind (unit = iInputLAY)
  do ipart = 1, Npart
    ComplexPlane(ipart) = .false.
    EpsZReIm(ipart) = 0.95_O
    if (DS .and. autGenDS) then
      ComplexPlanePart = .false.
      EpsZReImPart = 0.95_O
      string     = 'SourceRegPosAut'
      if (XFindPar (iInputLAY, string)) then
        read (iInputLAY, *, iostat = ios) ComplexPlanePart
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable ComplexPlanePart;')"
          print "(  2x,'for the layer ',i3,';')", ipart
          stop
        end if
        read (iInputLAY, *, iostat = ios) EpsZReImPart
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable EpsZReImPart;')"
          print "(  2x,'for the layer ',i3,';')", ipart
          stop
        end if               
      else
        print "(/,2x,'Group name SourceRegPosAut not found;')"
        print "(  2x,'for the layer ',i3,';')", ipart
        stop  
      end if       
      ComplexPlane(ipart) = ComplexPlanePart
      EpsZReIm(ipart) = EpsZReImPart
    end if
  end do 
!  
  if (DoConvTest) then
    print "(/,2x,'Convergence Test for a Layered Particle')" 
    print "(  2x,'---------------------------------------')"         
  else
    print "(/,2x,'Convergence Test for a Layered Particle over Mrank')" 
    print "(  2x,'--------------------------------------------------')"    
  end if                          
  allocate (Nrankp(Npart), Nrankp1(Npart)) 
  rewind (unit = iInputLAY)
  Nrankpmax  = 0
  Nrankpmax1 = 0
  do ipart = 1, Npart
    if (ipart == 1) then
      x = k * lnorm(ipart)
    else
      x = k * lnorm(ipart) * abs(ind_ref(ipart)) ! or the real part of the ref. index
    end if
    NrankW = int(x + 4.05_O * x**0.33_O + 2._O)
    if (.not. DoConvTest .or. InputDS) then
      NrankPart = 17
      string    = 'NrankReg'
      if (XFindPar (iInputLAY, string)) then
        read (iInputLAY, *, iostat = ios) NrankPart
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable NrankPart;')"
          print "(  2x,'for the layer ',i3,';')", ipart
          stop
        end if
      else
        print "(/,2x,'Group name NrankReg not found;')"
        print "(  2x,'for the layer ',i3,';')", ipart
        stop  
      end if        
      Nrankp(ipart) = NrankPart
      if (ipart == 1) print "(/,2x,'Nrank input values:')"
      print "(2x,'the input value of Nrank for the layer ',i3,' is ', i4,',')",     &
              ipart, NrankPart
      print "(2x, a, i3, a)",                                                       &
     'while the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';' 
    else
      if (ipart == 1) print "(/,2x,'Nrank estimates:')"  
      print "(2x,'the estimated value of Nrank from Wiscombe''s criterion')" 
      print "(2x,'for layer ',i2,' is ',i3,';')", ipart, NrankW 
      print "(2x,'- enter the estimated value of Nrank for layer ',i2)", ipart      
      call read_integer (Nrankp(ipart))
    end if
    if (Nrankp(ipart) > Nrankpmax) Nrankpmax = Nrankp(ipart)
    Nrankp1(ipart) = Nrankp(ipart) - 1
    if (Nrankp1(ipart) > Nrankpmax1) Nrankpmax1 = Nrankp1(ipart)
  end do
!
  rewind (unit = iInputLAY)  
  if (DoConvTest) then
    print "(/,2x, a)",                                                              &
   '- enter the estimated values of Nint, where Nint = Ndgs * Nrankpmax,'        
    print "(  2x,'  Ndgs = 15,20,..., and Nrankpmax = ',i4,';')", Nrankpmax                 
    call read_integer (Nint) 
  else    
    Nint   = 100
    string = 'NintGlobal'
    if (XFindPar (iInputLAY, string)) then
      read (iInputLAY, *, iostat = ios) Nint
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Nint;')"
        stop
      end if          
    else
      print "(/,2x,'Group name NintGlobal not found;')"
      stop  
    end if 
    print "(/,2x,'Nint input value:')"
    print "(  2x, a, i4, a)",                                                       &
   'the input value of Nint is ', Nint, ', while Nint = Ndgs * Nrankpmax,'
    print "(  2x,'Ndgs = 15,20,..., and Nrankpmax = ',i4,';')", Nrankpmax                   
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
  rewind (unit = iInputLAY)
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
        call zDSLAY (TypeGeom, Npart, Nsurfmax, surf, Nrankpmax, Nrankp, zpart,     &
             ComplexPlane, EpsZReIm, zRe, zIm)
        if (TypeConvTest == 2) call zDSLAY (TypeGeom, Npart, Nsurfmax, surf,        &
                                    Nrankpmax1, Nrankp1, zpart, ComplexPlane,       &
                                    EpsZReIm, zRe1, zIm1)
      else
        do i = 1, NrankPD
          zRePart(i)  = 0._O
          zImPart(i)  = 0._O
          zRePart1(i) = 0._O
          zImPart1(i) = 0._O
        end do
        string = 'SourceRegPosInp'
        if (XFindPar (iInputLAY, string)) then
          do i = 1, Nrankp(ipart)
            read (iInputLAY, *, iostat = ios) zRePart(i), zImPart(i)
            if (ios /= 0) then
              print "(/,2x, a)",                                                    &
             'Error by reading the input variables zRePart and zImPart;'   	      
              print "(  2x,'for the layer ',i3,';')", ipart
              stop
            end if 
          end do
          if (TypeConvTest == 2) then	  
            read (iInputLAY, *)
            do i = 1, Nrankp(ipart) - 1
              read (iInputLAY, *, iostat = ios) zRePart1(i), zImPart1(i)
              if (ios /= 0) then
                print "(/,2x, a)",                                                  &
               'Error by reading the input variables zRePart1 and zImPart1;'		
                print "(  2x,'for the layer ',i3,';')", ipart
                stop
              end if 
            end do
          end if	    
        else
          print "(/,2x,'Group name SourceRegPosInp not found;')"
          print "(  2x,'for the layer ',i3,';')", ipart
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
  epsNint  = 5.e-2
  epsNrank = 5.e-2
  epsMrank = 5.e-2
  dNint = 4
  dNintMrank = 10
  string     = 'Errors'
  if (XFindPar (iInputLAY, string)) then
    read (iInputLAY, *, iostat = ios) epsNint
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsNint;')"
      stop
    end if
    read (iInputLAY, *, iostat = ios) epsNrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsNrank;')"
      stop
    end if
    read (iInputLAY, *, iostat = ios) epsMrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsMrank;')"
      stop
    end if 
    read (iInputLAY, *, iostat = ios) dNint
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable dNint;')"
      stop
    end if
    read (iInputLAY, *, iostat = ios) dNintMrank
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
  if (XFindPar (iInputLAY, string)) then
    read (iInputLAY, *, iostat = ios) FileTmat
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
  if (XFindPar (iInputLAY, string)) then
    read (iInputLAY, *, iostat = ios) PrnProgress
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable PrnProgress;')"
      stop
    end if             
  else
    print "(/,2x,'Group name PrintProgress not found;')"
    stop  
  end if  
  close (unit = iInputLAY)     
end subroutine readinputLAY    
!***********************************************************************************
subroutine printinputLAY (ic, TypeGeom, Nsurfmax, Nsurf, surf, Nparam, Npart,       &
           lnorm, anorm, Nrankpmax, Nrankp, zpart, zRe, zIm, Nrankpmax1, Nrankp1,   &
           zRe1, zIm1, ind_ref, dNintMrank, dNint, wavelength, ind_refMed, epsNint, &
           epsNrank, epsMrank, DS, autGenDS)
  use parameters
  implicit none  
  integer    :: ic, TypeGeom, Nsurfmax, Npart, Nrankpmax, Nrankpmax1, dNint,        &
                Nrankp(Npart), Nrankp1(Npart), Nparam(Npart), Nsurf(Npart),         &
                NsurfPart, Nrank, dNintMrank, i, j
  real(O)    :: wavelength, anorm, surf(Npart,Nsurfmax), zpart(Npart),              &
                lnorm(Npart), zRe(Npart,Nrankpmax), zIm(Npart,Nrankpmax),           &
                zRe1(Npart,Nrankpmax1), zIm1(Npart,Nrankpmax1),                     &
                epsNint, epsNrank, epsMrank, ind_refMed
  complex(O) :: ind_ref(Npart)  
  logical    :: DS, autGenDS  
!
  write (iOutput,"(/,2x,'Input Parameters of the Scattering Problem:',/)")
  write (iOutput,"(2x,'wavelength of the free space, wavelength = ',1pe13.4,';')")  &
         wavelength
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'refractive index of the ambient medium, ind_refMed = ', ind_refMed, ';'
  write (iOutput,"(2x, a)")                                                         &
 'relative refractive index of each layer with respect to the ambient medium:'
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
    write (iOutput,"(2x,'concentric spheroids;')")      
  else if (TypeGeom == 2) then
    write (iOutput,"(2x,'concentric cylinders;')")      
  end if
  write (iOutput,"(2x,'number of homogeneous layers, Npart = ',i2,';')") Npart
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'characteristic length of the layered particle, anorm = ', anorm, ';'
  do i = 1, Npart
    write (iOutput,*)
    write (iOutput,"(2x,'layer: ',i2)") i
    write (iOutput,"(2x,'number of surface parameters, Nsurf = ',i2,';')") Nsurf(i)
    write (iOutput,"(2x,'surface parameters:')")
    NsurfPart = Nsurf(i)
    do j = 1, NsurfPart     
      write (iOutput,"(2x,'surf(',i2,') = ',1pe10.3,',')") j, surf(i,j)     
    end do
    if (TypeGeom == 1) then
      write (iOutput,"(2x,'where surf(1) is the semi-axis along the symmetry axis')")
      write (iOutput,"(2x,'and   surf(2) is the second semi-axis;')")
    else if (TypeGeom == 2) then
      write (iOutput,"(2x,'where surf(1) is the half-length of the cylinder')")
      write (iOutput,"(2x,'and   surf(2) is the cylinder radius;')")
    end if       
    write (iOutput,"(2x,'number of integration surfaces, Nparam = ',i2,';')")       &
           Nparam(i)
    write (iOutput,"(2x, a, 1pe10.3, a)")                                           &
   'caracteristic length of the layer, lnorm = ', lnorm(i), ';'
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
      write (iOutput,"(2x,'layer: ',i2)") i
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
        write (iOutput,"(2x,'layer: ',i2)") i
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
  write (iOutput,"(2x, a, i3, a)")                                                  &
 'integration step for Mrank calculation, dNintMrank = ', dNintMrank, '.'
  write (iOutput,"(/)")                       
end subroutine printinputLAY
! **********************************************************************************
subroutine convergence_NintLAY (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,        &
           Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zpart, Nint, dNint,         &
           epsNint, PrnProgress)    
  use parameters 
  implicit none
  integer    :: TypeGeom, Nsurfmax, Nparammax, Npart, Nrankpmax, Nint, dNint,       &
                Nrankp(Npart), Nparam(Npart)
  real(O)    :: k, surf(Npart,Nsurfmax), zpart(Npart), snorm, epsNint                
  complex(O) :: ind_ref(Npart)
  logical    :: PrnProgress
!        
  integer    :: Mstart, Mrank, Nmax, Nmaxmax, Nmaxpmax, Nteta, i, ipart,            &
                m, iNint, NthetaConv
  real(O)    :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,        &
                Cext, Qext  
  integer,allocatable    :: Nintparam(:,:), Nmaxp(:)
  real(O),allocatable    :: paramG(:,:,:), weightsG(:,:,:), h(:), v(:), oldh(:),    &
                            oldv(:)
  complex(O),allocatable :: aa(:,:), bb(:,:), a(:,:), b(:,:), c(:,:), cv(:),        &
                            cv1(:), cc(:)
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
  call write_TypeConvHead (1)
  allocate (Nmaxp(Npart))
  Nmaxpmax = 0
  Nmax = 0
  do ipart = 1, Npart    
    Nmaxp(ipart) = Nrankp(ipart)    
    if (Nmaxp(ipart) > Nmaxpmax) Nmaxpmax = Nmaxp(ipart)
    if (ipart < Npart) then
      Nmax = Nmax + 2 * Nmaxp(ipart)
    else
      Nmax = Nmax + Nmaxp(ipart)
    endif
  end do
  allocate (aa(2*Nmax,2*Nmax), bb(2*Nmax,2*Nmax))            
  allocate (a(2*Nrankpmax,2*Nrankpmax), b(2*Nrankpmax,2*Nrankpmax),                 &
            c(2*Nrankpmax,2*Nrankpmax), cv(2*Nrankpmax),cv1(2*Nrankpmax))
  Nmaxmax = Nrankpmax + Mrank * (2 * Nrankpmax - Mrank + 1)    
  allocate (cc(2*Nmaxmax)) 
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do
  if (PrnProgress) call write_progress (.true., 1, 11)
  do iNint = 1, 2        
    allocate (paramG(Npart,Nparammax,Nint), weightsG(Npart,Nparammax,Nint),         &
              Nintparam(Npart,Nparammax))                 
    call interpolation_listLAY (TypeGeom, Npart, Nsurfmax, surf, Nint, Nparammax,   &
         Nintparam, paramG, weightsG)     
    call matrix_Q31_LAY (TypeGeom, k, ind_ref, Nsurfmax, surf, Npart, Nrankp,       &
         Nmaxpmax, Nmaxp, zpart, m, Nmax, Nint, Nparammax, Nparam, Nintparam,       &
         paramG, weightsG, aa, Nmax, Nmax)  
    if (PrnProgress) call write_progress (.false., 2+5*(iNint-1), 11)                           
    call inverse_matrix (aa, 2*Nmax, 2*Nmax, bb, 2*Nmax, 2*Nmax, 2*Nmax)
    if (PrnProgress) call write_progress (.false., 3+5*(iNint-1), 11)
    call matrix_Q1_LAY (TypeGeom, 1, k, ind_ref, Nsurfmax, surf, Npart, Nrankpmax,  &
         Nrankp, Nmaxpmax, Nmaxp, zpart, m, Nint, Nparammax, Nparam, Nintparam,     &
         paramG, weightsG, a, Nrankpmax, Nrankpmax)    
    call extract_matrix3 (1, Nmaxpmax, bb, Nmax, Nmax, b, Nrankpmax, Nrankpmax)  
    call product_matrices (2*Nmaxpmax, 2*Nmaxpmax, 2*Nmaxpmax,                      &
         a, 2*Nrankpmax, 2*Nrankpmax, b, 2*Nrankpmax, 2*Nrankpmax)     
    if (PrnProgress) call write_progress (.false., 4+5*(iNint-1), 11)
    call matrix_Q1_LAY (TypeGeom, 3, k, ind_ref, Nsurfmax, surf, Npart, Nrankpmax,  &
         Nrankp, Nmaxpmax, Nmaxp, zpart, m, Nint, Nparammax, Nparam, Nintparam,     &
         paramG, weightsG, c, Nrankpmax, Nrankpmax)            
    call extract_matrix3 (2, Nmaxpmax, bb, Nmax, Nmax, b, Nrankpmax, Nrankpmax)  
    call product_matrices (2*Nmaxpmax, 2*Nmaxpmax, 2*Nmaxpmax,                      &
         c, 2*Nrankpmax, 2*Nrankpmax, b, 2*Nrankpmax, 2*Nrankpmax)
    if (PrnProgress) call write_progress (.false., 5+5*(iNint-1), 11)
    call sum_matrices (2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax, 2*Nrankpmax,         &
         c, 2*Nrankpmax, 2*Nrankpmax)
    call incident_matrix_LAY (TypeGeom, k, Nsurfmax, surf, Npart, Nrankpmax,        &
         Nrankp, Nmaxpmax, Nmaxp, zpart, m, Nint, Nparammax, Nparam, Nintparam,     &
         paramG, weightsG, c, Nrankpmax, Nrankpmax)     
    call product_matrices (2*Nmaxpmax, 2*Nmaxpmax, 2*Nmaxpmax,                      &
         a, 2*Nrankpmax, 2*Nrankpmax, c, 2*Nrankpmax, 2*Nrankpmax)               
    if (PrnProgress) call write_progress (.false., 6+5*(iNint-1), 11)
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m,            &
         Nrankpmax, Nmaxpmax, cv)
    call product_matrix_vector (2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax,             &
         2*Nrankpmax, cv, cv1)      
    call extend_vector_positive (cv1, cc, m, Mstart, Nrankpmax, Nmaxpmax, Nmaxmax)                                   
    call matrix_m_negativ (Nmaxpmax, Nmaxpmax, a, Nrankpmax, Nrankpmax)
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m,           &
         Nrankpmax, Nmaxpmax, cv)
    call product_matrix_vector (2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax,             &
         2*Nrankpmax, cv, cv1)
    call extend_vector_negative (cv1, cc, m, Nrankpmax, Nmaxpmax, Nmaxmax)
    call DSCS (cc, Mrank, Nrankpmax, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k,    &
         snorm,.false.,.true., h, v)
    call CQscat (cc, Mrank, Nrankpmax, Nmaxmax, k, snorm, Cscat, Qscat)
    call CQext (cc, Mrank, Nrankpmax, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,     &
         alfap, k, snorm, Cext, Qext)
    call delta_DSCS (Nteta, h, v, oldh, oldv, epsNint, NthetaConv)
    call write_3ConvParamReg (Nint, m, Npart, Nrankp, .false.)
    call write_DSCS  (Nteta,.false., h, v)
    call write_Effic (Qscat, Qext)
    Nint = Nint + dNint
    deallocate (paramG, weightsG, Nintparam)
  end do
  call write_NintConvRes (NthetaConv, Nteta, epsNint)           
  deallocate (aa, bb, a, b, c, cv, cv1, cc, h, v, oldh, oldv, Nmaxp)
end subroutine convergence_NintLAY
! **********************************************************************************
subroutine convergence_NrankLAY (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,       &
           Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zpart, Nint, epsNrank,      &
           PrnProgress)    
  use parameters 
  implicit none
  integer    :: TypeGeom, Nsurfmax, Nparammax, Npart, Nrankpmax, Nint,              &
                Nrankp(Npart), Nparam(Npart)
  real(O)    :: k, surf(Npart,Nsurfmax), zpart(Npart), snorm, epsNrank                
  complex(O) :: ind_ref(Npart)
  logical    :: PrnProgress
!        
  integer    :: Mstart, Mrank, Nmax, Nmaxpmax, Nmaxmax, Nteta, i, m, ipart,         &
                NthetaConv
  real(O)    :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,        &
                Cext, Qext  
  integer,allocatable    :: Nintparam(:,:), Nmaxp(:)
  real(O),allocatable    :: paramG(:,:,:), weightsG(:,:,:), h(:), v(:), oldh(:),    &
                            oldv(:)
  complex(O),allocatable :: aa(:,:), bb(:,:), a(:,:), b(:,:), c(:,:), cv(:),        &
                            cv1(:), cc(:), aa1(:,:), b1(:,:), b3(:,:), a0(:,:)
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
  call write_TypeConvHead (2)
  allocate (Nmaxp(Npart))
  Nmaxpmax = 0
  Nmax = 0
  do ipart = 1, Npart    
    Nmaxp(ipart) = Nrankp(ipart)    
    if (Nmaxp(ipart) > Nmaxpmax) Nmaxpmax = Nmaxp(ipart)
    if (ipart < Npart) then
      Nmax = Nmax + 2 * Nmaxp(ipart)
    else
      Nmax = Nmax + Nmaxp(ipart)
    end if
  end do    
  allocate (aa1(2*Nmax,2*Nmax), b1(2*Nrankpmax,2*Nrankpmax),                        &
            b3(2*Nrankpmax,2*Nrankpmax), a0(2*Nrankpmax,2*Nrankpmax))
  allocate (aa(2*Nmax,2*Nmax), bb(2*Nmax,2*Nmax))            
  allocate (a(2*Nrankpmax,2*Nrankpmax), b(2*Nrankpmax,2*Nrankpmax),                 &
            c(2*Nrankpmax,2*Nrankpmax), cv(2*Nrankpmax), cv1(2*Nrankpmax))
  Nmaxmax = Nrankpmax + Mrank * (2 * Nrankpmax - Mrank + 1)    
  allocate (cc(2*Nmaxmax)) 
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))                        
  allocate (paramG(Npart,Nparammax,Nint), weightsG(Npart,Nparammax,Nint),           &
            Nintparam(Npart,Nparammax))           
  call interpolation_listLAY (TypeGeom, Npart, Nsurfmax, surf, Nint, Nparammax,     &
       Nintparam, paramG, weightsG)
  if (PrnProgress) call write_progress (.true., 1, 9)
! --- Nrank configuration ---  
  call matrix_Q31_LAY (TypeGeom, k, ind_ref, Nsurfmax, surf, Npart, Nrankp,         &
       Nmaxpmax, Nmaxp, zpart, m, Nmax, Nint, Nparammax, Nparam, Nintparam,         &
       paramG, weightsG, aa, Nmax, Nmax)
  if (PrnProgress) call write_progress (.false., 2, 9)
  call copy_matrix (2*Nmax, 2*Nmax, aa, 2*Nmax, 2*Nmax, aa1, 2*Nmax, 2*Nmax)
  call inverse_matrix (aa, 2*Nmax, 2*Nmax, bb, 2*Nmax, 2*Nmax, 2*Nmax)
  if (PrnProgress) call write_progress (.false., 3, 9)
  call matrix_Q1_LAY (TypeGeom, 1, k, ind_ref, Nsurfmax, surf, Npart, Nrankpmax,    &
       Nrankp, Nmaxpmax, Nmaxp, zpart, m, Nint, Nparammax, Nparam, Nintparam,       &
       paramG, weightsG, a, Nrankpmax, Nrankpmax)
  call copy_matrix (2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax, 2*Nrankpmax,            &
       b1, 2*Nrankpmax, 2*Nrankpmax)
  call extract_matrix3 (1, Nmaxpmax, bb, Nmax, Nmax, b, Nrankpmax, Nrankpmax)           
  call product_matrices (2*Nmaxpmax, 2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax,        &
       2*Nrankpmax, b, 2*Nrankpmax, 2*Nrankpmax)           
  if (PrnProgress) call write_progress (.false., 4, 9) 
  call matrix_Q1_LAY (TypeGeom, 3, k, ind_ref, Nsurfmax, surf, Npart, Nrankpmax,    &
       Nrankp, Nmaxpmax, Nmaxp, zpart, m, Nint, Nparammax, Nparam, Nintparam,       &
       paramG, weightsG, c, Nrankpmax, Nrankpmax)
  call copy_matrix (2*Nmaxpmax, 2*Nmaxpmax, c, 2*Nrankpmax, 2*Nrankpmax,            &
       b3, 2*Nrankpmax, 2*Nrankpmax)
  call extract_matrix3 (2, Nmaxpmax, bb, Nmax, Nmax, b, Nrankpmax, Nrankpmax)   
  call product_matrices (2*Nmaxpmax, 2*Nmaxpmax, 2*Nmaxpmax, c, 2*Nrankpmax,        &
       2*Nrankpmax, b, 2*Nrankpmax, 2*Nrankpmax)
  if (PrnProgress) call write_progress (.false., 5, 9)
  call sum_matrices (2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax, 2*Nrankpmax,           &
       c, 2*Nrankpmax, 2*Nrankpmax)
  call incident_matrix_LAY (TypeGeom, k, Nsurfmax, surf, Npart, Nrankpmax,          &
       Nrankp, Nmaxpmax, Nmaxp, zpart, m, Nint, Nparammax, Nparam, Nintparam,       &
       paramG, weightsG, c, Nrankpmax, Nrankpmax)
  call copy_matrix (2*Nmaxpmax, 2*Nmaxpmax, c, 2*Nrankpmax, 2*Nrankpmax,            &
       a0, 2*Nrankpmax, 2*Nrankpmax)
  call product_matrices (2*Nmaxpmax, 2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax,        &
       2*Nrankpmax, c, 2*Nrankpmax, 2*Nrankpmax)             
  if (PrnProgress) call write_progress (.false., 6, 9)
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m,              &
       Nrankpmax, Nmaxpmax, cv)
  call product_matrix_vector (2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax,               &
       2*Nrankpmax, cv, cv1)
  call extend_vector_positive (cv1, cc, m, Mstart, Nrankpmax, Nmaxpmax, Nmaxmax)                                 
  call matrix_m_negativ (Nmaxpmax, Nmaxpmax, a, Nrankpmax, Nrankpmax)
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m,             &
       Nrankpmax, Nmaxpmax, cv)
  call product_matrix_vector (2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax,               &
       2*Nrankpmax, cv, cv1)
  call extend_vector_negative (cv1, cc, m, Nrankpmax, Nmaxpmax, Nmaxmax)
  call DSCS (cc, Mrank, Nrankpmax, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k,      &
       snorm,.false.,.true., h, v)
  call CQscat (cc, Mrank, Nrankpmax, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, Nrankpmax, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,       &
       alfap, k, snorm, Cext, Qext)
  do i = 1, Nteta
    oldh(i) = h(i)
    oldv(i) = v(i)
  end do
  call write_3ConvParamReg (Nint, m, Npart, Nrankp,.false.)
  call write_DSCS (Nteta,.false., h, v)
  call write_Effic (Qscat, Qext)
! --- (Nrank - 1) configuration --- 
  call copy_matrix (2*Nmax, 2*Nmax, aa1, 2*Nmax, 2*Nmax, aa, 2*Nmax, 2*Nmax)
  call internal_matrix_LAY_Nrank_m (Npart, Nmaxp, aa, Nmax, Nmax)      
  if (PrnProgress) call write_progress (.false., 7, 9)
  call inverse_matrix (aa, 2*Nmax, 2*Nmax, bb, 2*Nmax, 2*Nmax, 2*Nmax) 
  if (PrnProgress) call write_progress (.false., 8, 9)
  call copy_matrix (2*Nmaxpmax, 2*Nmaxpmax, b1, 2*Nrankpmax, 2*Nrankpmax,           &
       a, 2*Nrankpmax, 2*Nrankpmax)
  call matrix_Nrank_m (Nmaxpmax, a, Nrankpmax, Nrankpmax)   
  call extract_matrix3 (1, Nmaxpmax, bb, Nmax, Nmax, b, Nrankpmax, Nrankpmax)         
  call product_matrices (2*Nmaxpmax, 2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax,        &
       2*Nrankpmax, b, 2*Nrankpmax, 2*Nrankpmax)
  call copy_matrix (2*Nmaxpmax, 2*Nmaxpmax, b3, 2*Nrankpmax, 2*Nrankpmax,           &
       c, 2*Nrankpmax, 2*Nrankpmax)
  call matrix_Nrank_m (Nmaxpmax, c, Nrankpmax, Nrankpmax)                  
  call extract_matrix3 (2, Nmaxpmax, bb, Nmax, Nmax, b, Nrankpmax, Nrankpmax)        
  call product_matrices (2*Nmaxpmax, 2*Nmaxpmax, 2*Nmaxpmax, c, 2*Nrankpmax,        &
       2*Nrankpmax, b, 2*Nrankpmax, 2*Nrankpmax)
  call sum_matrices (2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax, 2*Nrankpmax,           &
       c, 2*Nrankpmax, 2*Nrankpmax)
  call copy_matrix (2*Nmaxpmax, 2*Nmaxpmax, a0, 2*Nrankpmax, 2*Nrankpmax,           &
       c, 2*Nrankpmax, 2*Nrankpmax)
  call matrix_Nrank_m (Nmaxpmax, c, Nrankpmax, Nrankpmax)
  call product_matrices (2*Nmaxpmax, 2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax,        &
       2*Nrankpmax, c, 2*Nrankpmax, 2*Nrankpmax)             
  if (PrnProgress) call write_progress (.false., 9, 9)
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m,              &
       Nrankpmax, Nmaxpmax, cv)
  call product_matrix_vector (2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax,               &
       2*Nrankpmax, cv, cv1)
  call extend_vector_positive (cv1, cc, m, Mstart, Nrankpmax, Nmaxpmax, Nmaxmax)                                 
  call matrix_m_negativ (Nmaxpmax, Nmaxpmax, a, Nrankpmax, Nrankpmax)
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m,             &
       Nrankpmax, Nmaxpmax, cv)
  call product_matrix_vector (2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax,               &
       2*Nrankpmax, cv, cv1)
  call extend_vector_negative (cv1, cc, m, Nrankpmax, Nmaxpmax, Nmaxmax)
  call DSCS (cc, Mrank, Nrankpmax, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k,      &
       snorm,.false.,.true., h, v)
  call CQscat (cc, Mrank, Nrankpmax, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, Nrankpmax, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,       &
       alfap, k, snorm, Cext, Qext)
  call delta_DSCS (Nteta, h, v, oldh, oldv, epsNrank, NthetaConv)
  call write_3ConvParamReg (Nint, m, Npart, Nrankp, .true.)  
  call write_DSCS (Nteta,.false., h, v)
  call write_Effic (Qscat, Qext)
  call write_NrankConvRes (NthetaConv, Nteta, epsNrank)
  deallocate (aa1, b1, b3, a0)
  deallocate (aa, bb, a, b, c, cv, cv1, cc, h, v, oldh, oldv, Nmaxp, paramG,        &
              weightsG, Nintparam)
end subroutine convergence_NrankLAY
! **********************************************************************************
subroutine convergence_MrankLAY (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,       &
           Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zpart, Nint, epsMrank,      &
           FileTmat, PrnProgress)          
  use parameters 
  implicit none
  integer       :: TypeGeom, Nsurfmax, Nparammax, Npart, Nrankpmax, Nint,           &
                   Nrankp(Npart), Nparam(Npart)
  real(O)       :: k, surf(Npart,Nsurfmax), zpart(Npart), snorm, epsMrank                
  complex(O)    :: ind_ref(Npart)
  character(80) :: FileTmat
  logical       :: PrnProgress
!             
  integer       :: Mstart, Mrank, Nmaxpmax, Nrank, Nmax, Nmaxmax, Nteta,            &
                   i, m, ipart, NthetaConv
  real(O)       :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,     &
                   Cext, Qext
  integer,allocatable    :: Nintparam(:,:), Nmaxp(:)
  real(O),allocatable    :: paramG(:,:,:), weightsG(:,:,:), h(:), v(:), oldh(:),    &
                            oldv(:)
  complex(O),allocatable :: aa(:,:), bb(:,:), a(:,:), b(:,:), c(:,:), cv(:),        &
                            cv1(:), cc(:)
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
  Mrank  = Nrankpmax   
  open (unit = iTmat, file = FileTmat, status = 'replace')
  call write_HeadFileTmat (Nrankpmax, Nrankpmax)    
  call write_TypeConvHead (3)
  call write_2ConvParamReg (Nint, Npart, Nrankp)
  Nrank = 0
  do ipart = 1, Npart        
    if(ipart < Npart) then
      Nrank = Nrank + 2 * Nrankp(ipart)
    else
      Nrank = Nrank + Nrankp(ipart)
    endif
  end do
  allocate (aa(2*Nrank,2*Nrank), bb(2*Nrank,2*Nrank))            
  allocate (a(2*Nrankpmax,2*Nrankpmax), b(2*Nrankpmax,2*Nrankpmax),                 &
            c(2*Nrankpmax,2*Nrankpmax), cv(2*Nrankpmax),cv1(2*Nrankpmax))       
  Nmaxmax = Nrankpmax + Mrank * (2 * Nrankpmax - Mrank + 1)
  allocate (cc(2*Nmaxmax)) 
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do     
  allocate (paramG(Npart,Nparammax,Nint), weightsG(Npart,Nparammax,Nint),           &
            Nintparam(Npart,Nparammax))           
  call interpolation_listLAY (TypeGeom, Npart, Nsurfmax, surf, Nint, Nparammax,     &
       Nintparam, paramG, weightsG)     
  Mrank = - 1
  do m = Mstart, Nrankpmax
    call write_1ConvParam (m)
    Mrank = Mrank + 1    
    do ipart = 1, Npart     
      if (Nrankp(ipart) < m) Nrankp(ipart) = m      
    end do
    allocate (Nmaxp(Npart))
    Nmaxpmax = 0
    Nmax = 0
    do ipart = 1, Npart
      if (m == 0) then
        Nmaxp(ipart) = Nrankp(ipart)
      else
        Nmaxp(ipart) = Nrankp(ipart) - iabs(m) + 1
      end if
      if (Nmaxp(ipart) > Nmaxpmax) Nmaxpmax = Nmaxp(ipart)
      if (ipart < Npart) then
        Nmax = Nmax + 2 * Nmaxp(ipart)
      else
        Nmax = Nmax + Nmaxp(ipart)
      end if
    end do
    if (PrnProgress) call write_progress_m (.true., m, 1, 6)  
    call matrix_Q31_LAY (TypeGeom, k, ind_ref, Nsurfmax, surf, Npart, Nrankp,       &
         Nmaxpmax, Nmaxp, zpart, m, Nmax, Nint, Nparammax, Nparam, Nintparam,       &
         paramG, weightsG, aa, Nrank, Nrank)
    if (PrnProgress) call write_progress_m (.false., m, 2, 6)  
    call inverse_matrix (aa, 2*Nrank, 2*Nrank, bb, 2*Nrank, 2*Nrank, 2*Nmax)
    if (PrnProgress) call write_progress_m (.false., m, 3, 6)
    call matrix_Q1_LAY (TypeGeom, 1, k, ind_ref, Nsurfmax, surf, Npart, Nrankpmax,  &
         Nrankp, Nmaxpmax, Nmaxp, zpart, m, Nint, Nparammax, Nparam, Nintparam,     &
         paramG, weightsG, a, Nrankpmax, Nrankpmax)
    call extract_matrix3 (1, Nmaxpmax, bb, Nrank, Nrank, b, Nrankpmax, Nrankpmax)    
    call product_matrices (2*Nmaxpmax, 2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax,      &
         2*Nrankpmax, b, 2*Nrankpmax, 2*Nrankpmax)     
    if (PrnProgress) call write_progress_m (.false., m, 4, 6)
    call matrix_Q1_LAY (TypeGeom, 3, k, ind_ref, Nsurfmax, surf, Npart, Nrankpmax,  &
         Nrankp, Nmaxpmax, Nmaxp, zpart, m, Nint, Nparammax, Nparam, Nintparam,     &
         paramG, weightsG, c, Nrankpmax, Nrankpmax)        
    call extract_matrix3 (2, Nmaxpmax, bb, Nrank, Nrank, b, Nrankpmax, Nrankpmax)    
    call product_matrices (2*Nmaxpmax, 2*Nmaxpmax, 2*Nmaxpmax, c, 2*Nrankpmax,      &
         2*Nrankpmax, b, 2*Nrankpmax, 2*Nrankpmax)
    if (PrnProgress) call write_progress_m (.false., m, 5, 6)
    call sum_matrices (2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax, 2*Nrankpmax,         &
         c, 2*Nrankpmax, 2*Nrankpmax)
    call incident_matrix_LAY (TypeGeom, k, Nsurfmax, surf, Npart, Nrankpmax,        &
         Nrankp, Nmaxpmax, Nmaxp, zpart, m, Nint, Nparammax, Nparam, Nintparam,     &
         paramG, weightsG, c, Nrankpmax, Nrankpmax)
    call product_matrices (2*Nmaxpmax, 2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax,      &
         2*Nrankpmax, c, 2*Nrankpmax, 2*Nrankpmax)
    if (PrnProgress) call write_progress_m (.false., m, 6, 6) 
    call write_FileTmat (Nrankpmax, Nrankpmax, a)        
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m,            &
         Nrankpmax, Nmaxpmax, cv)
    call product_matrix_vector (2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax,             &
         2*Nrankpmax, cv, cv1)
    call extend_vector_positive (cv1, cc, m, Mstart, Nrankpmax, Nmaxpmax, Nmaxmax)                                 
    if (m /= 0) then
      call matrix_m_negativ (Nmaxpmax, Nmaxpmax, a, Nrankpmax, Nrankpmax)
      call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m,         &
           Nrankpmax, Nmaxpmax, cv)
      call product_matrix_vector (2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax,           &
           2*Nrankpmax, cv, cv1)
      call extend_vector_negative (cv1, cc, m, Nrankpmax, Nmaxpmax, Nmaxmax)
    end if      
    deallocate (Nmaxp)        
    call DSCS (cc, Mrank, Nrankpmax, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k,    &
         snorm,.false.,.true., h, v)    
    call delta_DSCS (Nteta, h, v, oldh, oldv, epsMrank, NthetaConv)    
    call write_DSCS (Nteta,.false., h, v)
    if (NthetaConv >= int(0.8d0*Nteta)) exit
  end do  
  close (unit = iTmat)  
  call CQscat (cc, Mrank, Nrankpmax, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, Nrankpmax, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,       &
       alfap, k, snorm, Cext, Qext)
  call write_Effic (Qscat, Qext)
  call write_MrankConvRes (NthetaConv, epsMrank)               
  if (NthetaConv >= int(0.8*Nteta)) then
    print "(/,2x,'Convergence criterion for Mrank is satisfied;')"                            
  else
    print "(/,2x,'Convergence criterion for Mrank is not satisfied;')"
  end if
  call write_InfoFileTmat (FileTmat, Mrank, Nrankpmax, .true., .false., .false.)
  call ScatCharact (k, FileTmat, Mrank, Nrankpmax, .true., .false., .false.)    
  print "(/,2x,'T matrix is stored in ',a50)", FileTmat
  print "(  2x,'The dimensions of the T matrix are given by:')"
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", Nrankpmax   
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank  
  deallocate (aa, bb, a, b, c, cv, cv1, cc, h, v, oldh, oldv, paramG, weightsG,     &
              Nintparam)
end subroutine convergence_MrankLAY    
! **********************************************************************************
subroutine convergence_NintDSLAY (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,      &
           Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zRe, zIm, zpart, Nint,      &
           dNint, epsNint, PrnProgress)    
  use parameters 
  implicit none
  integer    :: TypeGeom, Nsurfmax, Nparammax, Npart, Nrankpmax, Nint, dNint,       &
                Nrankp(Npart), Nparam(Npart)
  real(O)    :: k, surf(Npart,Nsurfmax), zpart(Npart), snorm, epsNint,              &
                zRe(Npart,Nrankpmax), zIm(Npart,Nrankpmax)                
  complex(O) :: ind_ref(Npart)
  logical    :: PrnProgress
!        
  integer    :: Mstart, Mrank, Nrank, Nmaxmax, Nmaxpmax, Nteta, i, ipart,           &
                m, iNint, NthetaConv
  real(O)    :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,        &
                Cext, Qext  
  integer,allocatable    :: Nintparam(:,:)
  real(O),allocatable    :: paramG(:,:,:), weightsG(:,:,:), h(:), v(:), oldh(:),    &
                            oldv(:)
  complex(O),allocatable :: aa(:,:), bb(:,:), a(:,:), b(:,:), c(:,:), cv(:),        &
                            cv1(:), cc(:)
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
  call write_TypeConvHead (1)
  Nmaxpmax = Nrankpmax
  Nrank = 0
  do ipart = 1,Npart - 1
    Nrank = Nrank + 2 * Nrankp(ipart)
  end do          
  Nrank = Nrank + Nrankp(Npart)    
  allocate (aa(2*Nrank,2*Nrank), bb(2*Nrank,2*Nrank))            
  allocate (a(2*Nrankpmax,2*Nrankpmax), b(2*Nrankpmax,2*Nrankpmax),                 &
            c(2*Nrankpmax,2*Nrankpmax), cv(2*Nrankpmax),cv1(2*Nrankpmax))
  Nmaxmax = Nrankpmax + Mrank * (2 * Nrankpmax - Mrank + 1)    
  allocate (cc(2*Nmaxmax)) 
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do
  if (PrnProgress) call write_progress (.true., 1, 11)
  do iNint = 1, 2        
    allocate (paramG(Npart,Nparammax,Nint), weightsG(Npart,Nparammax,Nint),         &
              Nintparam(Npart,Nparammax))                 
    call interpolation_listLAY (TypeGeom, Npart, Nsurfmax, surf, Nint, Nparammax,   &
         Nintparam, paramG, weightsG)     
    call matrix_Q31_DS_LAY (TypeGeom, k, ind_ref, Nsurfmax, surf, Npart,            &
         Nrankpmax, Nrankp, zRe, zIm, zpart, m, Nrank, Nint, Nparammax, Nparam,     &
         Nintparam, paramG, weightsG, aa, Nrank, Nrank)                                          
    if (PrnProgress) call write_progress (.false., 2+5*(iNint-1), 11)
    call inverse_matrix (aa, 2*Nrank, 2*Nrank, bb, 2*Nrank, 2*Nrank, 2*Nrank)     
    if (PrnProgress) call write_progress (.false., 3+5*(iNint-1), 11)
    call matrix_Q1_DS_LAY (TypeGeom, 1, k, ind_ref, Nsurfmax, surf, Npart,          &
         Nrankpmax, Nrankp, zRe, zIm, zpart, m, Nrankpmax, Nmaxpmax, Nint,          &
         Nparammax, Nparam, Nintparam, paramG, weightsG, a, Nrankpmax, Nrankpmax)    
    call extract_matrix3 (1, Nrankpmax, bb, Nrank, Nrank, b, Nrankpmax, Nrankpmax)        
    call product_matrices (2*Nmaxpmax, 2*Nrankpmax, 2*Nrankpmax, a, 2*Nrankpmax,    &
         2*Nrankpmax, b, 2*Nrankpmax, 2*Nrankpmax)     
    if (PrnProgress) call write_progress (.false., 4+5*(iNint-1), 11)
    call matrix_Q1_DS_LAY (TypeGeom, 3, k, ind_ref, Nsurfmax, surf, Npart,          &
         Nrankpmax, Nrankp, zRe, zIm, zpart, m, Nrankpmax, Nmaxpmax, Nint,          &
         Nparammax, Nparam, Nintparam, paramG, weightsG, c, Nrankpmax, Nrankpmax)             
    call extract_matrix3 (2, Nrankpmax, bb, Nrank, Nrank, b, Nrankpmax, Nrankpmax)           
    call product_matrices (2*Nmaxpmax, 2*Nrankpmax, 2*Nrankpmax, c, 2*Nrankpmax,    &
         2*Nrankpmax, b, 2*Nrankpmax, 2*Nrankpmax)    
    if (PrnProgress) call write_progress (.false., 5+5*(iNint-1), 11)
    call sum_matrices (2*Nmaxpmax, 2*Nrankpmax, a, 2*Nrankpmax, 2*Nrankpmax,        &
         c, 2*Nrankpmax, 2*Nrankpmax)
    call incident_matrix_DS_LAY (TypeGeom, k, Nsurfmax, surf, Npart, Nrankpmax,     &
         Nrankp, zRe, zIm, zpart, m, Nrankpmax, Nmaxpmax, Nint, Nparammax, Nparam,  &
         Nintparam, paramG, weightsG, c, Nrankpmax, Nrankpmax)     
    call product_matrices (2*Nmaxpmax, 2*Nrankpmax, 2*Nmaxpmax, a, 2*Nrankpmax,     &
         2*Nrankpmax, c, 2*Nrankpmax, 2*Nrankpmax)               
    if (PrnProgress) call write_progress (.false., 6+5*(iNint-1), 11)
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m,            &
         Nrankpmax, Nmaxpmax, cv)
    call product_matrix_vector (2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax,             &
         2*Nrankpmax, cv, cv1)
    call extend_vector_positive (cv1, cc, m, Mstart, Nrankpmax, Nmaxpmax, Nmaxmax)                                 
    call matrix_m_negativ (Nmaxpmax, Nmaxpmax, a, Nrankpmax, Nrankpmax)
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m,           &
         Nrankpmax, Nmaxpmax, cv)
    call product_matrix_vector (2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax,             &
         2*Nrankpmax, cv, cv1)
    call extend_vector_negative (cv1, cc, m, Nrankpmax, Nmaxpmax, Nmaxmax)
    call DSCS (cc, Mrank, Nrankpmax, Nmaxmax, Nteta, phiGS, alfa, beta, gama,k,     &
         snorm,.false.,.true., h, v)
    call CQscat (cc, Mrank, Nrankpmax, Nmaxmax, k, snorm, Cscat, Qscat)
    call CQext (cc, Mrank, Nrankpmax, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,     &
         alfap, k, snorm, Cext, Qext)
    call delta_DSCS (Nteta, h, v, oldh, oldv, epsNint, NthetaConv)    
    call write_3ConvParamReg (Nint, m, Npart, Nrankp,.false.)
    call write_DSCS (Nteta,.false., h, v)
    call write_Effic (Qscat, Qext)
    Nint = Nint + dNint
    deallocate (paramG, weightsG, Nintparam)
  end do
  call write_NintConvRes (NthetaConv, Nteta, epsNint)
  deallocate (aa, bb, a, b, c, cv, cv1, cc, h, v, oldh, oldv)
end subroutine convergence_NintDSLAY 
! **********************************************************************************   
subroutine convergence_NrankDSLAY (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,     &
           Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zRe, zIm, Nrankpmax1,       &
           Nrankp1, zRe1, zIm1, zpart, Nint, epsNrank, PrnProgress)        
  use parameters 
  implicit none
  integer    :: TypeGeom, Nsurfmax, Nparammax, Npart, Nrankpmax, Nrankpmax1, Nint,  &
                Nrankp(Npart), Nrankp1(Npart), Nparam(Npart)
  real(O)    :: k, surf(Npart,Nsurfmax), zpart(Npart), snorm, epsNrank,             &
                zRe(Npart,Nrankpmax), zIm(Npart,Nrankpmax),                         &
                zRe1(Npart,Nrankpmax1), zIm1(Npart,Nrankpmax1)
  complex(O) :: ind_ref(Npart)
  logical    :: PrnProgress
!        
  integer    :: Mstart, Mrank, Nteta, Nrank, Nrank1, Nmaxpmax, Nmaxmax, i,          &
                m, ipart, NthetaConv
  real(O)    :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,        &
                Cext, Qext  
  integer,allocatable    :: Nintparam(:,:)
  real(O),allocatable    :: paramG(:,:,:), weightsG(:,:,:), h(:), v(:), oldh(:),    &
                            oldv(:)
  complex(O),allocatable :: aa(:,:), bb(:,:), a(:,:), b(:,:), c(:,:), cv(:),        &
                            cv1(:), cc(:)
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
  call write_TypeConvHead (2)
  Nmaxpmax = Nrankpmax
  Nrank  = 0
  Nrank1 = 0
  do ipart = 1,Npart - 1
    Nrank  = Nrank  + 2 * Nrankp(ipart)
    Nrank1 = Nrank1 + 2 * Nrankp1(ipart)
  end do         
  Nrank  = Nrank  + Nrankp(Npart) 
  Nrank1 = Nrank1 + Nrankp1(Npart)   
  allocate (aa(2*Nrank,2*Nrank), bb(2*Nrank,2*Nrank))            
  allocate (a(2*Nrankpmax,2*Nrankpmax), b(2*Nrankpmax,2*Nrankpmax),                 &  
            c(2*Nrankpmax,2*Nrankpmax), cv(2*Nrankpmax), cv1(2*Nrankpmax))
  Nmaxmax = Nrankpmax + Mrank * (2 * Nrankpmax - Mrank + 1)    
  allocate (cc(2*Nmaxmax)) 
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))    
  allocate (paramG(Npart,Nparammax,Nint), weightsG(Npart,Nparammax,Nint),           &
            Nintparam(Npart,Nparammax))           
  call interpolation_listLAY (TypeGeom, Npart, Nsurfmax, surf, Nint, Nparammax,     &
       Nintparam, paramG, weightsG)
  if (PrnProgress) call write_progress (.true., 1, 11)          
! --- Nrank configuration ---  
  call matrix_Q31_DS_LAY (TypeGeom, k, ind_ref, Nsurfmax, surf, Npart, Nrankpmax,   &
       Nrankp, zRe, zIm, zpart, m, Nrank, Nint, Nparammax, Nparam, Nintparam,       &
       paramG, weightsG, aa, Nrank, Nrank)                           
  if (PrnProgress) call write_progress (.false., 2, 11)
  call inverse_matrix (aa, 2*Nrank, 2*Nrank, bb, 2*Nrank, 2*Nrank, 2*Nrank)
  if (PrnProgress) call write_progress (.false., 3, 11)
  call matrix_Q1_DS_LAY (TypeGeom, 1, k, ind_ref, Nsurfmax, surf, Npart, Nrankpmax, &
       Nrankp, zRe, zIm, zpart, m, Nrankpmax, Nmaxpmax, Nint, Nparammax, Nparam,    &
       Nintparam, paramG, weightsG, a, Nrankpmax, Nrankpmax)    
  call extract_matrix3 (1, Nrankpmax, bb, Nrank, Nrank, b, Nrankpmax, Nrankpmax)           
  call product_matrices (2*Nmaxpmax, 2*Nrankpmax, 2*Nrankpmax, a, 2*Nrankpmax,      &
       2*Nrankpmax, b, 2*Nrankpmax, 2*Nrankpmax)           
  if (PrnProgress) call write_progress (.false., 4, 11)
  call matrix_Q1_DS_LAY (TypeGeom, 3, k, ind_ref, Nsurfmax, surf, Npart, Nrankpmax, &
       Nrankp, zRe, zIm, zpart, m, Nrankpmax, Nmaxpmax, Nint, Nparammax, Nparam,    &
       Nintparam, paramG, weightsG, c, Nrankpmax, Nrankpmax)              
  call extract_matrix3 (2, Nrankpmax, bb, Nrank, Nrank, b, Nrankpmax, Nrankpmax)             
  call product_matrices (2*Nmaxpmax, 2*Nrankpmax, 2*Nrankpmax, c, 2*Nrankpmax,      &
       2*Nrankpmax, b, 2*Nrankpmax, 2*Nrankpmax)    
  if (PrnProgress) call write_progress (.false., 5, 11)
  call sum_matrices (2*Nmaxpmax, 2*Nrankpmax, a, 2*Nrankpmax, 2*Nrankpmax,          &
       c, 2*Nrankpmax, 2*Nrankpmax)
  call incident_matrix_DS_LAY (TypeGeom, k, Nsurfmax, surf, Npart, Nrankpmax,       &
       Nrankp, zRe, zIm, zpart, m, Nrankpmax, Nmaxpmax, Nint, Nparammax, Nparam,    &
       Nintparam, paramG, weightsG, c, Nrankpmax, Nrankpmax)   
  call product_matrices (2*Nmaxpmax, 2*Nrankpmax, 2*Nmaxpmax, a, 2*Nrankpmax,       &
       2*Nrankpmax, c, 2*Nrankpmax, 2*Nrankpmax)             
  if (PrnProgress) call write_progress (.false., 6, 11)
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrankpmax,   &
       Nmaxpmax, cv)
  call product_matrix_vector (2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax,               &
       2*Nrankpmax, cv, cv1)
  call extend_vector_positive (cv1, cc, m, Mstart, Nrankpmax, Nmaxpmax, Nmaxmax)                                 
  call matrix_m_negativ (Nmaxpmax, Nmaxpmax, a, Nrankpmax, Nrankpmax)
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m,             &
       Nrankpmax, Nmaxpmax, cv)
  call product_matrix_vector (2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax,               &
       2*Nrankpmax, cv, cv1)
  call extend_vector_negative (cv1, cc, m, Nrankpmax, Nmaxpmax, Nmaxmax)
  call DSCS (cc, Mrank, Nrankpmax, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k,      &
       snorm,.false.,.true., h, v)
  call CQscat (cc, Mrank, Nrankpmax, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, Nrankpmax, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,       &
       alfap, k, snorm, Cext, Qext)  
  do i = 1, Nteta
    oldh(i) = h(i)
    oldv(i) = v(i)
  end do
  call write_3ConvParamReg (Nint, m, Npart, Nrankp,.false.)
  call write_DSCS  (Nteta,.false., h, v)
  call write_Effic (Qscat, Qext)
! --- Nrank1 configuration ---  
  Nmaxpmax = Nrankpmax1
  Nmaxmax  = Nrankpmax1 + Mrank * (2 * Nrankpmax1 - Mrank + 1)       
  call matrix_Q31_DS_LAY (TypeGeom, k, ind_ref, Nsurfmax, surf, Npart, Nrankpmax1,  &
       Nrankp1, zRe1, zIm1, zpart, m, Nrank1, Nint, Nparammax, Nparam, Nintparam,   &
       paramG, weightsG, aa, Nrank, Nrank)                          
  if (PrnProgress) call write_progress (.false., 7, 11)
  call inverse_matrix (aa, 2*Nrank, 2*Nrank, bb, 2*Nrank, 2*Nrank, 2*Nrank1)
  if (PrnProgress) call write_progress (.false., 8, 11)
  call matrix_Q1_DS_LAY (TypeGeom, 1, k, ind_ref, Nsurfmax, surf, Npart,            &
       Nrankpmax1, Nrankp1, zRe1, zIm1, zpart, m, Nrankpmax1, Nmaxpmax, Nint,       &
       Nparammax, Nparam, Nintparam, paramG, weightsG, a, Nrankpmax, Nrankpmax)    
  call extract_matrix3 (1, Nrankpmax1, bb, Nrank, Nrank, b, Nrankpmax, Nrankpmax)       
  call product_matrices (2*Nmaxpmax, 2*Nrankpmax1, 2*Nrankpmax1, a, 2*Nrankpmax,    &
       2*Nrankpmax, b, 2*Nrankpmax, 2*Nrankpmax)           
  if (PrnProgress) call write_progress (.false., 9, 11)
  call matrix_Q1_DS_LAY (TypeGeom, 3, k, ind_ref, Nsurfmax, surf, Npart,            &
       Nrankpmax1, Nrankp1, zRe1, zIm1, zpart, m, Nrankpmax1, Nmaxpmax, Nint,       &
       Nparammax, Nparam, Nintparam, paramG, weightsG, c, Nrankpmax, Nrankpmax)           
  call extract_matrix3 (2, Nrankpmax1, bb, Nrank, Nrank, b, Nrankpmax, Nrankpmax)            
  call product_matrices (2*Nmaxpmax, 2*Nrankpmax1, 2*Nrankpmax1, c, 2*Nrankpmax,    &
       2*Nrankpmax, b, 2*Nrankpmax, 2*Nrankpmax)    
  if (PrnProgress) call write_progress (.false., 10, 11)
  call sum_matrices (2*Nmaxpmax, 2*Nrankpmax1, a, 2*Nrankpmax, 2*Nrankpmax,         &
       c, 2*Nrankpmax, 2*Nrankpmax)
  call incident_matrix_DS_LAY (TypeGeom, k, Nsurfmax, surf, Npart, Nrankpmax1,      &
       Nrankp1, zRe1, zIm1, zpart, m, Nrankpmax1, Nmaxpmax, Nint, Nparammax,        &
       Nparam, Nintparam, paramG, weightsG, c, Nrankpmax, Nrankpmax)   
  call product_matrices (2*Nmaxpmax, 2*Nrankpmax1, 2*Nmaxpmax, a, 2*Nrankpmax,      &
       2*Nrankpmax, c, 2*Nrankpmax, 2*Nrankpmax)             
  if (PrnProgress) call write_progress (.false., 11, 11)
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrankpmax1,  &
       Nmaxpmax, cv)
  call product_matrix_vector (2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax,               &
       2*Nrankpmax, cv, cv1)
  call extend_vector_positive (cv1, cc, m, Mstart, Nrankpmax1, Nmaxpmax, Nmaxmax)                                
  call matrix_m_negativ (Nmaxpmax, Nmaxpmax, a, Nrankpmax, Nrankpmax)
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m,             &
       Nrankpmax1, Nmaxpmax, cv)
  call product_matrix_vector (2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax,               &
       2*Nrankpmax, cv, cv1)
  call extend_vector_negative (cv1, cc, m, Nrankpmax1, Nmaxpmax, Nmaxmax)
  call DSCS (cc, Mrank, Nrankpmax1, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k,     &
       snorm,.false.,.true., h, v)
  call CQscat (cc, Mrank, Nrankpmax1, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, Nrankpmax1, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,      &
       alfap, k, snorm, Cext, Qext)  
  call delta_DSCS (Nteta, h, v, oldh, oldv, epsNrank, NthetaConv)  
  call write_3ConvParamReg (Nint, m, Npart, Nrankp1,.false.)
  call write_DSCS (Nteta,.false., h, v)
  call write_Effic (Qscat, Qext)
  call write_NrankConvRes (NthetaConv, Nteta, epsNrank)
  deallocate (aa, bb, a, b, c, cv, cv1, cc, h, v, oldh, oldv, paramG, weightsG,     &
              Nintparam)
end subroutine convergence_NrankDSLAY 
! **********************************************************************************
subroutine convergence_MrankDSLAY (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,     &
           Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zRe, zIm, zpart, Nint,      &
           ComplexPlane, EpsZReIm, autGenDS, epsMrank, epsNrank, dNintMrank,        &
           FileTmat, PrnProgress)    
  use parameters 
  implicit none
  integer       :: TypeGeom, Nsurfmax, Nparammax, Npart, Nrankpmax, Nint,           &
                   Nrankp(Npart), Nparam(Npart), dNintMrank
  real(O)       :: k, surf(Npart,Nsurfmax), zpart(Npart), snorm, epsMrank,          &
                   epsNrank, zRe(Npart,Nrankpmax), zIm(Npart,Nrankpmax),            &
                   EpsZReIm(Npart)
  complex(O)    :: ind_ref(Npart)
  character(80) :: FileTmat
  logical       :: autGenDS, PrnProgress, ComplexPlane(Npart)
!             
  integer       :: Mstart, Mrank, Nteta, Nrank, Nmaxmax, i, m, ipart, NthetaConv,   &
                   NrankAL, NrankpmaxAL, NrankG, NmaxG, NrankpmaxL, NrankpmaxL1
  real(O)       :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat,            &
                   Qscat, Cext, Qext
  logical       :: ComplexDS, ChangeDS, more, continuare, ConvTest
  integer,allocatable    :: Nintparam(:,:), Nrankp1(:)
  real(O),allocatable    :: paramG(:,:,:), weightsG(:,:,:), h(:), v(:), oldh(:),    &
                            oldv(:), zReL(:,:), zImL(:,:), zReL1(:,:), zImL1(:,:)
  complex(O),allocatable :: aa(:,:), bb(:,:), a(:,:), b(:,:), c(:,:), cv(:),        &
                            cv1(:), cc(:)
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
  Mrank  = Nrankpmax 
  Nrank  = 0
  do ipart = 1, Npart - 1
    Nrank = Nrank + 2 * Nrankp(ipart)
  end do     
  Nrank   = Nrank + Nrankp(Npart)      
  NrankAL = Nrank
  NrankpmaxAL = Nrankpmax  
  NrankG  = Nrankpmax 
  Nmaxmax = NrankG + Mrank * (2 * NrankG - Mrank + 1)
  open (unit = iTmat, file = FileTmat, status = 'replace')
  call write_HeadFileTmat (NrankpmaxAL, NrankpmaxAL) 
  call write_TypeConvHead (3)
  call write_2ConvParamReg (Nint, Npart, Nrankp)                  
  allocate (aa(2*NrankAL,2*NrankAL), bb(2*NrankAL,2*NrankAL))            
  allocate (a(2*NrankpmaxAL,2*NrankpmaxAL), b(2*NrankpmaxAL,2*NrankpmaxAL),         &
            c(2*NrankpmaxAL,2*NrankpmaxAL)) 
  allocate (cv(2*NrankG), cv1(2*NrankG))      
  allocate (cc(2*Nmaxmax))   
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do    
  allocate (paramG(Npart,Nparammax,Nint), weightsG(Npart,Nparammax,Nint),           &
            Nintparam(Npart,Nparammax))           
  call interpolation_listLAY (TypeGeom, Npart, Nsurfmax, surf, Nint, Nparammax,     &
       Nintparam, paramG, weightsG)
  ComplexDS = .false.
  do ipart = 1, Npart           
    do i = 1, Nrankp(ipart)                    
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
  do m = Mstart, NrankG
    call write_1ConvParam (m)
    Mrank = Mrank + 1
    if (m == 0) then
      NmaxG = NrankG
    else
      NmaxG = NrankG - m + 1
    end if
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
      if (.not. ChangeDS) then
        deallocate (paramG, weightsG, Nintparam)
        Nint = Nint + dNintMrank
        allocate (paramG(Npart,Nparammax,Nint), weightsG(Npart,Nparammax,Nint),     &
                  Nintparam(Npart,Nparammax))           
        call interpolation_listLAY (TypeGeom, Npart, Nsurfmax, surf, Nint,          &
             Nparammax, Nintparam, paramG, weightsG)     
      else
        print "(/,2x,'Azimuthal mode: m = ',i3)", m  
        more = .true.
        do while (more)
          print "(  2x, '- enter the number of integration points Nint')"
          call read_integer (Nint) 
          deallocate (paramG, weightsG, Nintparam)        
          allocate (paramG(Npart,Nparammax,Nint), weightsG(Npart,Nparammax,Nint),   &
                    Nintparam(Npart,Nparammax))           
          call interpolation_listLAY (TypeGeom, Npart, Nsurfmax, surf, Nint,        &
               Nparammax, Nintparam, paramG, weightsG)                 
          NrankpmaxL = 0  
          Nrank      = 0             
          do ipart = 1, Npart
            if (ComplexPlane(ipart)) then
              print "(2x,'- enter the estimated value of Nrank for region ',i2)",   &
	              ipart
              call read_integer (Nrankp(ipart))		
	    end if	 
            if (Nrankp(ipart)  > NrankpmaxL)  NrankpmaxL  = Nrankp(ipart)
            if (ipart < Npart) then
              Nrank = Nrank + 2 * Nrankp(ipart)
            else
              Nrank = Nrank + Nrankp(ipart)
            end if             	  
          end do 
          if (NrankpmaxL > Nrankpmax) then
            print "(/,2x,'Input error:')"
            print "(  2x, a)",                                                      &
           'the number of discrete sources for the first layer exceeds the number'	    
            print "(  2x, a)",                                                      &
           'of discrete sources corresponding to the initial configuration;'            
            stop   
	  end if                                     
          deallocate (zReL, zImL)	   	   
          allocate (zReL(Npart,NrankpmaxL), zImL(Npart,NrankpmaxL))               
          call zDSLAY (TypeGeom, Npart, Nsurfmax, surf, NrankpmaxL, Nrankp, zpart,  &
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
            call zDSLAY (TypeGeom, Npart, Nsurfmax, surf, NrankpmaxL1, Nrankp1,     &
                 zpart, ComplexPlane, EpsZReIm, zReL1, zImL1)  
            call convergence_NrankDSLAY_m (m, TypeGeom, k, ind_ref, snorm,          &
                 Nsurfmax, surf, Nparammax, Nparam, Npart, NrankpmaxL, Nrankp,      &
                 zReL, zImL, NrankpmaxL1, Nrankp1, zReL1, zImL1, zpart, NrankG,     &
                 NmaxG, Nint, Nintparam, paramG, weightsG, epsNrank)	    
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
      end if                                                                               
    end if
    if (Nrank > NrankAL) then
      NrankAL = Nrank
      deallocate (aa, bb)
      allocate (aa(2*NrankAL,2*NrankAL), bb(2*NrankAL,2*NrankAL))                                  
    end if    
!   --- end local convergence test ---     
    if (PrnProgress) call write_progress_m (.true., m, 1, 6)                    
    call matrix_Q31_DS_LAY (TypeGeom, k, ind_ref, Nsurfmax, surf, Npart,            &
         NrankpmaxL, Nrankp, zReL, zImL, zpart, m, Nrank, Nint, Nparammax,          &
         Nparam, Nintparam, paramG, weightsG, aa, NrankAL, NrankAL)         
    if (PrnProgress) call write_progress_m (.false., m, 2, 6)                                            
    call inverse_matrix (aa, 2*NrankAL, 2*NrankAL, bb, 2*NrankAL, 2*NrankAL, 2*Nrank)
    if (PrnProgress) call write_progress_m (.false., m, 3, 6)        
    call matrix_Q1_DS_LAY (TypeGeom, 1, k, ind_ref, Nsurfmax, surf, Npart,          &
         NrankpmaxL, Nrankp, zReL, zImL, zpart, m, NrankG, NmaxG, Nint, Nparammax,  &
         Nparam, Nintparam, paramG, weightsG, a, NrankpmaxAL, NrankpmaxAL)                               
    call extract_matrix3 (1, NrankpmaxL, bb, NrankAL, NrankAL, b,                   &    
         NrankpmaxAL, NrankpmaxAL)                
    call product_matrices (2*NmaxG, 2*NrankpmaxL, 2*NrankpmaxL, a, 2*NrankpmaxAL,   &
         2*NrankpmaxAL, b, 2*NrankpmaxAL, 2*NrankpmaxAL)                       
    if (PrnProgress) call write_progress_m (.false., m, 4, 6)
    call matrix_Q1_DS_LAY (TypeGeom, 3, k, ind_ref, Nsurfmax, surf, Npart,          &
         NrankpmaxL, Nrankp, zReL, zImL, zpart, m, NrankG, NmaxG, Nint, Nparammax,  &
         Nparam, Nintparam, paramG, weightsG, c, NrankpmaxAL, NrankpmaxAL)                      
    call extract_matrix3 (2, NrankpmaxL, bb, NrankAL, NrankAL, b,                   &
         NrankpmaxAL, NrankpmaxAL)               
    call product_matrices (2*NmaxG, 2*NrankpmaxL, 2*NrankpmaxL, c, 2*NrankpmaxAL,   &
         2*NrankpmaxAL, b, 2*NrankpmaxAL, 2*NrankpmaxAL)             
    if (PrnProgress) call write_progress_m (.false., m, 5, 6)    
    call sum_matrices (2*NmaxG, 2*NrankpmaxL, a, 2*NrankpmaxAL, 2*NrankpmaxAL,      &
         c, 2*NrankpmaxAL, 2*NrankpmaxAL)         
    call incident_matrix_DS_LAY (TypeGeom, k, Nsurfmax, surf, Npart, NrankpmaxL,    &
         Nrankp, zReL, zImL, zpart, m, NrankG, NmaxG, Nint, Nparammax, Nparam,      &
         Nintparam, paramG, weightsG, c, NrankpmaxAL, NrankpmaxAL)                                         
    call product_matrices (2*NmaxG, 2*NrankpmaxL, 2*NmaxG, a, 2*NrankpmaxAL,        &
         2*NrankpmaxAL, c, 2*NrankpmaxAL, 2*NrankpmaxAL)         
    if (PrnProgress) call write_progress_m (.false., m, 6, 6)
    call write_FileTmat (NrankpmaxAL, NrankpmaxAL, a)                    
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m,            &
         NrankG, NmaxG, cv)
    call product_matrix_vector (2*NmaxG, 2*NmaxG, a, 2*NrankpmaxAL, 2*NrankpmaxAL,  &
         cv, cv1)
    call extend_vector_positive (cv1, cc, m, Mstart, NrankG, NmaxG, Nmaxmax)                                 
    if (m /= 0) then
      call matrix_m_negativ (NmaxG, NmaxG, a, NrankpmaxAL, NrankpmaxAL)
      call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m,         &
           NrankG, NmaxG, cv)
      call product_matrix_vector (2*NmaxG, 2*NmaxG, a, 2*NrankpmaxAL,               &
           2*NrankpmaxAL, cv, cv1)
      call extend_vector_negative (cv1, cc, m, NrankG, NmaxG, Nmaxmax)
    end if      
    call DSCS (cc, Mrank, NrankG, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k,       &
         snorm,.false.,.true., h, v)                                        
    call delta_DSCS (Nteta, h, v, oldh, oldv, epsMrank, NthetaConv)    
    call write_DSCS (Nteta,.false., h, v)
    deallocate (zReL, zImL)     
    if (NthetaConv >= int(0.8*Nteta)) exit
  end do 
  close (unit = iTmat)  
  call CQscat (cc, Mrank, NrankG, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, NrankG, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,          &
       alfap, k, snorm, Cext, Qext)      
  call write_Effic (Qscat, Qext)             
  call write_MrankConvRes (NthetaConv, epsMrank)
  if (NthetaConv >= int(0.8*Nteta)) then
    print "(/,2x,'Convergence criterion for Mrank is satisfied;')"                                                              
  else
    print "(/,2x,'Convergence criterion for Mrank is not satisfied;')"
  end if 
  call write_InfoFileTmat (FileTmat, Mrank, NrankG, .true., .false., .false.)
  call ScatCharact (k, FileTmat, Mrank, NrankG, .true., .false., .false.) 
  print "(/,2x,'T matrix is stored in ',a50)", FileTmat
  print "(  2x,'The dimensions of the T matrix are given by:')"
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", NrankG
  print "(  2x,'- number of azimuthal modes, Mrank = ',i2,';')", Mrank            
  deallocate (aa, bb, a, b, c, cv, cv1, cc, h, v, oldh, oldv, paramG, weightsG,     &
              Nintparam)
end subroutine convergence_MrankDSLAY
! **********************************************************************************   
subroutine convergence_NrankDSLAY_m (m, TypeGeom, k, ind_ref, snorm, Nsurfmax,      &
           surf, Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zRe, zIm, Nrankpmax1, &
           Nrankp1, zRe1, zIm1, zpart, NrankG, NmaxG, Nint, Nintparam, paramG,      &
           weightsG, epsNrank)        
  use parameters 
  implicit none
  integer    :: m, TypeGeom, Nsurfmax, Nparammax, Npart, Nrankpmax, Nrankpmax1,     &
                Nint, Nrankp(Npart), Nrankp1(Npart), Nparam(Npart), NrankG, NmaxG,  &
                Nintparam(Npart,Nparammax)                
  real(O)    :: k, surf(Npart,Nsurfmax), zpart(Npart), snorm, epsNrank,             &
                zRe(Npart,Nrankpmax), zIm(Npart,Nrankpmax),                         &
                zRe1(Npart,Nrankpmax1), zIm1(Npart,Nrankpmax1),                     &
                paramG(Npart,Nparammax,Nint), weightsG(Npart,Nparammax,Nint) 
  complex(O) :: ind_ref(Npart)  
!        
  integer    :: Mstart, Mrank, Nteta, Nrank, Nrank1, Nmaxmax, i, ipart, NthetaConv, &
                NrankAL, NrankpmaxAL
  real(O)    :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,        &
                Cext, Qext, teta  
  real(O),allocatable    :: h(:), v(:), oldh(:), oldv(:)
  complex(O),allocatable :: aa(:,:), bb(:,:), a(:,:), b(:,:), c(:,:), cv(:),        &
                            cv1(:), cc(:)
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
  Nrank  = 0
  do ipart = 1,Npart - 1
    Nrank  = Nrank  + 2 * Nrankp(ipart)    
  end do         
  Nrank  = Nrank  + Nrankp(Npart)     
  NrankpmaxAL = max(Nrankpmax,NrankG)
  NrankAL     = Nrank  
  allocate (aa(2*NrankAL,2*NrankAL), bb(2*NrankAL,2*NrankAL))            
  allocate (a(2*NrankpmaxAL,2*NrankpmaxAL), b(2*NrankpmaxAL,2*NrankpmaxAL),         &  
            c(2*NrankpmaxAL,2*NrankpmaxAL)) 
  allocate (cv(2*NrankG), cv1(2*NrankG))
  Nmaxmax = NrankG + Mrank * (2 * NrankG - Mrank + 1)    
  allocate (cc(2*Nmaxmax)) 
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))                  
! --- Nrank configuration ---  
  call matrix_Q31_DS_LAY (TypeGeom, k, ind_ref, Nsurfmax, surf, Npart, Nrankpmax,   &
       Nrankp, zRe, zIm, zpart, m, Nrank, Nint, Nparammax, Nparam, Nintparam,       &
       paramG, weightsG, aa, NrankAL, NrankAL)                                    
  call inverse_matrix (aa, 2*NrankAL, 2*NrankAL, bb, 2*NrankAL, 2*NrankAL, 2*Nrank)   
  call matrix_Q1_DS_LAY (TypeGeom, 1, k, ind_ref, Nsurfmax, surf, Npart,            &
       Nrankpmax, Nrankp, zRe, zIm, zpart, m, NrankG, NmaxG, Nint, Nparammax,       &
       Nparam, Nintparam, paramG, weightsG, a, NrankpmaxAL, NrankpmaxAL)           
  call extract_matrix3 (1, Nrankpmax, bb, NrankAL, NrankAL, b,                      &
       NrankpmaxAL, NrankpmaxAL)             
  call product_matrices (2*NmaxG, 2*Nrankpmax, 2*Nrankpmax, a, 2*NrankpmaxAL,       &
       2*NrankpmaxAL, b, 2*NrankpmaxAL, 2*NrankpmaxAL)                   
  call matrix_Q1_DS_LAY (TypeGeom, 3, k, ind_ref, Nsurfmax, surf, Npart,            &
       Nrankpmax, Nrankp, zRe, zIm, zpart, m, NrankG, NmaxG, Nint, Nparammax,       &
       Nparam, Nintparam, paramG, weightsG, c, NrankpmaxAL, NrankpmaxAL)                     
  call extract_matrix3 (2, Nrankpmax, bb, NrankAL, NrankAL, b,                      &
       NrankpmaxAL, NrankpmaxAL)               
  call product_matrices (2*NmaxG, 2*Nrankpmax, 2*Nrankpmax, c, 2*NrankpmaxAL,       &
       2*NrankpmaxAL, b, 2*NrankpmaxAL, 2*NrankpmaxAL)            
  call sum_matrices (2*NmaxG, 2*Nrankpmax, a, 2*NrankpmaxAL, 2*NrankpmaxAL,         &
       c, 2*NrankpmaxAL, 2*NrankpmaxAL)       
  call incident_matrix_DS_LAY (TypeGeom, k, Nsurfmax, surf, Npart, Nrankpmax,       &
       Nrankp, zRe, zIm, zpart, m, NrankG, NmaxG, Nint, Nparammax, Nparam,          &
       Nintparam, paramG, weightsG, c, NrankpmaxAL, NrankpmaxAL)          
  call product_matrices (2*NmaxG, 2*Nrankpmax, 2*NmaxG, a, 2*NrankpmaxAL,           &
       2*NrankpmaxAL, c, 2*NrankpmaxAL, 2*NrankpmaxAL)                    
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, NrankG,      &
       NmaxG, cv)       
  call product_matrix_vector (2*NmaxG, 2*NmaxG, a, 2*NrankpmaxAL,                   &
       2*NrankpmaxAL, cv, cv1)       
  call extend_vector_positive (cv1, cc, m, Mstart, NrankG, NmaxG, Nmaxmax)
  call matrix_m_negativ (NmaxG, NmaxG, a, NrankpmaxAL, NrankpmaxAL)  
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m,             &
       NrankG, NmaxG, cv)       
  call product_matrix_vector (2*NmaxG, 2*NmaxG, a, 2*NrankpmaxAL,                   &
       2*NrankpmaxAL, cv, cv1)       
  call extend_vector_negative (cv1, cc, m, NrankG, NmaxG, Nmaxmax)  
  call DSCS (cc, Mrank, NrankG, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k,         &
       snorm,.false.,.true., h, v)
  call CQscat (cc, Mrank, NrankG, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, NrankG, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,          &
       alfap, k, snorm, Cext, Qext)  
  do i = 1, Nteta
    oldh(i) = h(i)
    oldv(i) = v(i)
  end do    
  print "(1x,'Convergence Test over Nrank')"
  print "(1x,'---------------------------')"  
  print "(1x,'Nint = ',i5,',',1x,'Nrank = ',i3,',',1x,'Mrank = ',i3)", Nint,        &
          NrankG, m  
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
  Nrank1 = 0
  do ipart = 1, Npart - 1   
    Nrank1 = Nrank1 + 2 * Nrankp1(ipart)
  end do         
  Nrank1 = Nrank1 + Nrankp1(Npart)  
  call matrix_Q31_DS_LAY (TypeGeom, k, ind_ref, Nsurfmax, surf, Npart, Nrankpmax1,  &
       Nrankp1, zRe1, zIm1, zpart, m, Nrank1, Nint, Nparammax, Nparam, Nintparam,   &
       paramG, weightsG, aa, NrankAL, NrankAL)                                    
  call inverse_matrix (aa, 2*NrankAL, 2*NrankAL, bb, 2*NrankAL, 2*NrankAL, 2*Nrank1)   
  call matrix_Q1_DS_LAY (TypeGeom, 1, k, ind_ref, Nsurfmax, surf, Npart,            &
       Nrankpmax1, Nrankp1, zRe1, zIm1, zpart, m, NrankG, NmaxG, Nint, Nparammax,   &
       Nparam, Nintparam, paramG, weightsG, a, NrankpmaxAL, NrankpmaxAL)           
  call extract_matrix3 (1, Nrankpmax1, bb, NrankAL, NrankAL, b,                     &
       NrankpmaxAL, NrankpmaxAL)             
  call product_matrices (2*NmaxG, 2*Nrankpmax1, 2*Nrankpmax1, a, 2*NrankpmaxAL,     &
       2*NrankpmaxAL, b, 2*NrankpmaxAL, 2*NrankpmaxAL)                   
  call matrix_Q1_DS_LAY (TypeGeom, 3, k, ind_ref, Nsurfmax, surf, Npart,            &
       Nrankpmax1, Nrankp1, zRe1, zIm1, zpart, m, NrankG, NmaxG, Nint, Nparammax,   &
       Nparam, Nintparam, paramG, weightsG, c, NrankpmaxAL, NrankpmaxAL)                     
  call extract_matrix3 (2, Nrankpmax1, bb, NrankAL, NrankAL, b,                     &
       NrankpmaxAL, NrankpmaxAL)               
  call product_matrices (2*NmaxG, 2*Nrankpmax1, 2*Nrankpmax1, c, 2*NrankpmaxAL,     &
       2*NrankpmaxAL, b, 2*NrankpmaxAL, 2*NrankpmaxAL)            
  call sum_matrices (2*NmaxG, 2*Nrankpmax1, a, 2*NrankpmaxAL, 2*NrankpmaxAL,        &
       c, 2*NrankpmaxAL, 2*NrankpmaxAL)       
  call incident_matrix_DS_LAY (TypeGeom, k, Nsurfmax, surf, Npart, Nrankpmax1,      &
       Nrankp1, zRe1, zIm1, zpart, m, NrankG, NmaxG, Nint, Nparammax, Nparam,       &
       Nintparam, paramG, weightsG, c, NrankpmaxAL, NrankpmaxAL)          
  call product_matrices (2*NmaxG, 2*Nrankpmax1, 2*NmaxG, a, 2*NrankpmaxAL,          &
       2*NrankpmaxAL, c, 2*NrankpmaxAL, 2*NrankpmaxAL)                    
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, NrankG,      &
       NmaxG, cv)       
  call product_matrix_vector (2*NmaxG, 2*NmaxG, a, 2*NrankpmaxAL,                   &
       2*NrankpmaxAL, cv, cv1)       
  call extend_vector_positive (cv1, cc, m, Mstart, NrankG, NmaxG, Nmaxmax)  
  call matrix_m_negativ (NmaxG, NmaxG, a, NrankpmaxAL, NrankpmaxAL)  
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m,             &
       NrankG, NmaxG, cv)       
  call product_matrix_vector (2*NmaxG, 2*NmaxG, a, 2*NrankpmaxAL,                   &
       2*NrankpmaxAL, cv, cv1)       
  call extend_vector_negative (cv1, cc, m, NrankG, NmaxG, Nmaxmax)  
  call DSCS (cc, Mrank, NrankG, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k,         &
       snorm,.false.,.true., h, v)
  call CQscat (cc, Mrank, NrankG, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, NrankG, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,          &
       alfap, k, snorm, Cext, Qext) 
  call delta_DSCS (Nteta, h, v, oldh, oldv, epsNrank, NthetaConv)   
  print "(1x,'Nint = ',i5,',',1x,'Nrank = ',i3,',',1x,'Mrank = ',i3)", Nint,        &
          NrankG, m  
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
  print "(1x, a, i2, a, 1f5.2, a, /)",                                              &
 'The solution converges in ', NthetaConv, ' points with an relative error of ',    &
  100 * epsNrank, ' %;'                  
  deallocate (aa, bb, a, b, c, cv, cv1, cc, h, v, oldh, oldv)
end subroutine convergence_NrankDSLAY_m 
