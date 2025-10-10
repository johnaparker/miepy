subroutine TAXSYM
!------------------------------------------------------------------------------------
! 1. General Considerations                                                         !
! --------------------------                                                        !
! TAXSYM is a routine for computing the T matrix and the scattering characteristics !
! of homogeneous, dielectric (isotropic, chiral) and perfectly conducting,          !
! axisymmetric particles.                                                           !
!                                                                                   !
! The transition matrix of a homogeneous particle can be written as:                !
!                                                                                   !
!                        T = - Q11 * Q31^(-1),                                      !
!                                                                                   ! 
! where the expressions of the matrices Q31 and Q11 follow from the extinction      !
! theorem and Huygens principle, respectively. The elements of the matrices Q31     !
! and Q11 are expressed as integrals over the particle surface. For axisymmetric    !
! particles, the T-matrix computation decouples over the azimuthal modes and the    !
! surface integrals reduce to line integrals over the generatrix. The integrals     !
! are evaluated by means of the Gauss-Legendre quadrature method, and note that for !
! particles with a plane of symmetry perpendicular to the axis of rotation (mirror  !
! symmetric particles), the integration can be performed along the half-generatrix  !
! curve. The numerical stability and accuracy of T-matrix computations for an       !
! nonsmooth generatrix can be enhanced by dividing the generatrix into piecewise    !
! smooth curves and applying separate quadrature formulas to each curve.            !
!                                                                                   !
! Localized or distributed sources can be used for T-matrix calculation and note    !
! that distributed sources are preferable for highly aspherical particles. For      !
! highly elongated particles, the sources are distributed along the axis of         !
! symmetry of the particle, while for highly flattened particles, the sources       !
! are distributed in the complex plane. 			                    !
!                                                                                   !
! Particle geometries currently supported include: spheroids, cylinders and         !
! rounded oblate cylinders. The following routines (from the file "GeomLib.f90")    !
! provide the required geometry parameters:                                         !
! - the routine "interpolation_listAXSYM" computes the integration nodes and        !
!   weights for an axisymmmetric surface,                                           !
! - the routine "elem_geomAXSYM" provides the geometry parameters (position vector, !
!   zenith and azimuthal angles, and normal unit vector in spherical coordinates)   !
!   at each integration node, and                                                   !
! - the routine "zDSAXSYM" generates the coordinates of the distributed sources.    !
! The user can easily modify these routines to generate particles with other        !
! geometries. Note that the list of parameters must be maintained and only          !
! geometries with an analytical description of the surface can be implemented.      !
! Alternatively, the option FileGeom  = .true., causes the code to read the         !
! particle geometry information from the file FileFEM. The file FileFEM is read     !
! by the routine "read_FileFEMAxsym" from the file "InputOutput.f90". The user can  !
! customize the routine "read_FileFEMAxsym" as needed to conform to the manner in   !
! which the particle geometry is stored in the file FileFEM. However, as supplied,  !
! the routine "read_FileFEMAxsym" expects the file FileFEM to have the following    !
! structure:                                                                        !
! - one line providing the number of smooth curves Nlines; the format for           !
!   reading Nlines is i7, i.e.,                                                     !
!                   read (iFEM, "(i7)", iostat = ios)  Nlines;                      ! 
! - a set of Nlines-sequences of data specifying the geometry parameters on each     !
!   smooth curve; each sequence includes:                                           !
!   - one line containing the number of curve elements (vertices) NVvr; the         !
!     format for reading NVvr is i7, i.e.,                                          !
!                   read (iFEM, "(i7)", iostat = ios) NVvr;                         !                                         
!   - NVvr lines containing the current index ielem, the Cartesian coordinates      !   
!     of the curve element center: x = r(1) and z = r(2), the Cartesian components  !
!     of the unit normal vector at the curve element center: nx = n(1) and          !
!     nz = n(2), and the area of the surface element; the format for reading the    !
!     data is:                                                                      !
!                   read (iFEM, "(i7,2x,5(e15.7,2x))", iostat = ios)                !
!                         ielem, r(1), r(2), n(1), n(2), area.                      !
! As provided, the file FileFEM is contained in the directory "GEOMFILES".          !
! If the particle geometry is read from file, the coordinates of the distributed    !
! sources (if used) must be supplied by the user.                                   ! 
!                                                                                   !
! 2. Convergence Test                                                               !
! -------------------                                                               !
! Convergence tests precede the scattering characteristics calculation. Three       !
! parameters control the T-matrix computation:                                      !
! - the number of integration points for computing integrals over the generatrix,   !
!   Nint,                                                                           !
! - the maximum expansion order Nrank, and                                          !
! - the maximum azimuthal order Mrank.                                              ! 
! Nint is the global number of integration points on the generatrix curve. As       ! 
! mentioned before, the generatrix of an nonsmooth surface is divided into Nparam   !
! piecewise smooth curves and the number of integration points on the curve j is    !
!                   Nint_curve(j) = int { c(j) * Nint },                            !
! where c(j) is a multiplicative factor and j = 1, 2,.., Nparam. The parameters     !
! c(j), j = 1, 2,..., Nparam, do not appear as explicit parameters; they are hard   !
! coded in the routine "interpolation_listAXSYM" from the file "GeomLib.f90".       !
! For example, the generatrix of a cylinder is divided into three straight lines.   !
! The first and third lines are perpendicular to the axis of symmetry, while the    ! 
! second line is parallel to the axis of symmetry. Denoting by a the                !
! half-length of the cylinder and by b the cylinder radius, we have:                !
!                                                                                   !
!     Nint_curve(first line) = Nint_curve(third line) = int {[b/2/(a+b)] * Nint}    !
!                                                                                   !
! and                                                                               !
!                                                                                   ! 
!     Nint_curve(second line) = Nint - 2 * int {[b/2/(a+b)] * Nint }.               !
!                                                                                   !
! For convergence tests, the incident wave is assumed to be a vector plane wave     !
! traveling along the Z-axis of the global coordinate system and the scattering     !
! characteristics are computed in the azimuthal plane phi = 0°. The convergence     !
! tests over Nint and Nrank are interactive, while the convergence test over Mrank  !
! is automatically performed.                                                       ! 
!                                                                                   !
! The convergence tests discussed in the literature are based on the analysis of    !
! the differential scattering cross section or extinction and scattering cross      !
! sections. The convergence test based on differential scattering cross section     !
! analysis is more conservative and is implemented in our code.                     !
!                                                                                   !
! For the integration and expansion order test an axisymmetric orientation of the   !
! particle is considered, i.e., the incident wave propagates along the axis of      !
! symmetry of the particle. In this case, all azimuthal modes are zero except       !
! for m = - 1 and m = 1. For the convergence test over Nint the scattering problem  !
! is solved for Nint and Nint + dNint, while for the convergence test over Nrank    !
! the scattering problem is solved for Nrank and Nrank - 1. The normalized          !
! differential scattering cross section (DSCS) will be checked at 20° increments    !
! for convergence within epsX (epsNint or epsNrank) tolerance. If the calculated    !
! results converge within this tolerance at 80% of the scattering angles, then      !
! convergence is achieved. At each calculation, the extinction and scattering       !
! cross sections Cext and Cscat are also computed. Although the convergence of Cext !
! and Cscat does not guarantee that the DSCS converges, the divergence of Cext and  !
! Cscat implies the divergence of the DSCS. The results of the convergence          !
! analysis are written to the file "/OUTPUTFILES/Output.dat". The above             !
! convergence tests (based on the analysis of the differential scattering cross     !
! section) have been proposed by Barber and Hill [P. W. Barber and S. C. Hill,      !
! Light scattering by particles: Computational Methods, World Scientific,           !
! Singapore, 1990].                                                                 ! 
!                                                                                   ! 
! After Nrank and Nint have been determined we pass to the azimuthal order test.    !
! The program automatically sets the particle to a more general orientation, i.e.,  !
! alpha = beta = 45°, and solves the scattering problem  for increasing m values    !
! until convergence of the angular scattering is achieved. The T matrix is stored   !
! for later use by other programs, and the values of Nrank and Mrank are printed to !
! the screen and to the T-matrix information file (see "Description.txt"). These    !
! values together with the T matrix serve as INPUT PARAMETERS for other programs.   ! 
!                                                                                   !
! 3. Estimates of Nint and Nrank                                                    !
! -------------------------------                                                   !
! The above convergence tests require estimates of Nint and Nrank. These            !
! estimates must be supplied by the user, and for this purpose, the automatic       !
! convergence procedure proposed by Mishchenko [M. I. Mishchenko, Light scattering  !
! by a size-shape distributions of randomly oriented axially symmetric particles    !
! of a size comparable to the wavelength, Applied Optics, 32, 4652-4666, 1993]      !
! or Wiscombe's truncation limit criterion [W. J. Wiscombe, Improved Mie            !
! scattering algorithms, Applied Optics, 19, 1505-1509, 1980] can be used.          !
!                                                                                   !
! The procedure developed by Mishchenko finds reliable a priori estimates of        !
! both Nint and Nrank by using only the zeroth order block matrices of the T        !
! matrix. The estimates are computed by checking the convergence of the quantities: !
!                                                                                   !
!        Ce = - (2Pi / k**2) Re {SUM (2n + 1) * [T11_0n0n + T22_0n0n] }             !
!                                                                                   !
! and                                                                               !
!                                                                                   !
!        Cs =   (2Pi / k**2) SUM (2n + 1) * [|T11_0n0n|**2 + |T22_0n0n|**2]         !
!                                                                                   ! 
! where T11_0n0n' and T22_0n0n' are the zeroth order block matrices of the T        !
! matrix. Note that for spherical particles, Ce and Cs are the extinction and       !
! scattering cross sections. An important numerical parameter of the convergence    !
! procedure is the multiplicity factor Ndgs,                                        !
!                                                                                   !
!                          Nint = Ndgs * Nrank,                                     !
!                                                                                   !
! which controls the number of integration points and must be optimized for each    !
! particle shape. The main steps of the convergence procedure are                   !
! 1. choose Ndgs;                                                                   !
! 2. compute a starting value of Nrank, as for instance                             !
!                                                                                   !
!                 NrankW = int(x + 4.05 * x**0.33),                                 !
!                                                                                   !
!    where x is the size parameter, x = k * Rcirc, and Rcirc is the radius of the   !
!    smallest circumscribing sphere;                                                !
! 3. solve the scattering problem for increasing Nrank-values (in unit steps) and   !
!    Nint = Ndgs * Nrank until the convergence criterion                            !
!                                                                                   !
!          max { [Cs(Nrank) - Cs(Nrank-1) ] / Cs(Nrank),                            !
!                [Ce(Nrank)  - Ce(Nrank-1)] / Ce(Nrank) } < delta                   !
!                                                                                   !
!    is satisfied, where delta is the desired accuracy;                             !
! 4. after Nrank has been determined, increase Nint in dNint steps until Cs(Nrank)  !
!    and Ce(Nrank) converge within delta.                                           !
!                                                                                   !
! If the convergence procedure is finished WE CAN PASS DIRECTLY TO THE AZIMUTHAL    !
! ORDER TEST. However, it is recommended to check at least the value of Nrank by    !
! using the convergence procedures based on the differential scattering cross       !
! section analysis. The convergence procedure proposed by Mishchenko can not be     !
! used for chiral particles, if the coordinates of the discrete sources are         !
! user-defined and if the particle geometry is read from file.                      !
!                                                                                   !
! A reasonable starting value of Nrank is given by Wiscombe's truncation limit      !
! criterion, i.e.,                                                                  !
!                                                                                   !    
!                NrankW = int(x + 4.05 * x**0.33 + 2),                              !
!                                                                                   !
! and note that the convergence procedure proposed by Mishchenko also uses NrankW   !
! as initial Nrank-value.                                                           !
!                                                                                   !
! The value of Nint depends on the type of discrete sources, the size parameter,    !
! the particle shape and the relative refractive index. For localized sources,      !
! size parameters smaller than 15.0, and particle aspect ratios ranging between     !
! 0.33 and 3.0, the recommended values are:                                         !
!                for spheroids:       Ndgs = (3...5), and                           !
!                for cylinders:       Ndgs = (4...7).                               !
! These values should be substantially increased as the size parameter increases,   !
! the particles become less aspherical or distributed sources are used.             !
! For spheroids with a size parameter smaller than 30.0, the following routine      !
! can be used to compute a conservative estimate of Nint:                           !
!                                                                                   !
! a = surf(1), length of the semi-axis along the axis of symmetry                   !
! b = surf(2), knight of the second semi-axis                                       !
! e = a / b                                                                         !
! xmax  = k * max(a,b)                                                              !
! Nrank = int(xmax + 4.05 * xmax**0.33 + 2.0)                                       !
! if (xmax <= 30.0) then                                                            !
!     if (e == 1.0) then                                                            !
!         Nint = 4 * Nrank                                                          !
!     else if (e < 1.0) then                                                        !
!         kb = k * b                                                                !
!         if (kb < 5.0) then                                                        !
!             Nint = 100                                                            !
!         else                                                                      !
!             Nint = 20 * int(kb)                                                   !
!         end if                                                                    !
!         if (e < 0.5) then                                                         !
!             dNr  = 25.0 + 10.0 * kb                                               !
!             Nint = Nint + int(dNr * (1.0 / e - 2.0))                              !
!         end if                                                                    ! 
!     else                                                                          !
!         ka = k * a                                                                !
!         if (ka < 10.0) then                                                       !
!             Nint = 40                                                             !
!         else                                                                      !
!             Nint = 4 * int(ka)                                                    !
!         end if                                                                    !
!         if (e > 2) then                                                           !
!             if (ka < 10.0) then                                                   !
!                 dNr = 30.0                                                        !
!             else                                                                  !
!                 dNr = 30.0 + 3.5 * (ka - 10.0)                                    !
!             end if                                                                !
!             Nint = Nint + int(dNr * (e - 2.0))                                    !
!         end if                                                                    !
!     end if                                                                        !
! end if                                                                            !
!                                                                                   !
! 4. Strategies for Performing Convergence Tests                                    !
! ----------------------------------------------                                    !
! Two strategies for performing convergence tests over Nint and Nrank are           !
! recommended. The following technique has been proposed by Barber and Hill:        !
! 1. choose a value of Nrank close to the value predicted by Wiscombe's truncation  !
!    limit criterion;                                                               !
! 2. for Nrank, perform the expansion order test with Nint = Ndgs * Nrank (in the   !
!    original code of Barber and Hill, Ndgs = 2);                                   !
! 3. if convergence is not achieved, set Nrank = Nrank + 1, maintain the relation   !
!    Nint = Ndgs * Nrank and go to Step 2;                                          !
! 4. if convergence is achieved, perform several integration tests with the         !
!    initial Nint-value, Nint = Ndgs * Nrank, until the DSCS converges.             !
!                                                                                   !
! For large and/or highly aspherical particles, T-matrix computations become        !
! poorly convergent or even divergent. Frequently, a semi-convergent behaviour is   !
! attended: the errors decrease with increasing the maximum expansion order,        !
! attain a relative constant level and afterwards increase. The region of           !
! stability can be localized by using the following technique:                      !
! 1. choose a value of Nrank smaller (or significantly smaller) than the value      !
!    predicted by Wiscombe's truncation limit criterion;                            !
! 2. for Nrank, perform several integration tests with the initial Nint-value,      !
!    Nint = Ndgs * Nrank, until the DSCS converges;                                 !
! 3. if convergence can not be achieved, set Nrank = Nrank - 1 or increase the      !
!    tolerance epsNint, and go to Step 2;                                           !
! 4. if convergence is achieved (for some Nint), perform several expansion order    !
!    tests for increasing Nrank-values and monitor the errors of the differential   !
!    scattering cross sections.                                                     !
!                                                                                   !
! For particles which are too extreme in size and/or aspect ratio, the errors of    !
! the extinction and scattering cross sections (instead of the differential         !
! scattering cross section) can be analyzed.                                        ! 
!                                                                                   !
! 5. Additional Comments                                                            !
! ----------------------                                                            !
! For the expansion order test and distributed sources, two configurations of       !
! sources (with Nrank and Nrank - 1 discrete sources) are considered. For both      !
! configurations, the coordinates of the distributed sources can be generated       !
! automatically or can be specified in the input file.                              !
!                                                                                   !
! The convergence tests over Nint and Nrank can be switched off by setting the      !
! logical variable DoConvTest to false. In this case, the values of Nint and        !
! Nrank must be specified in the input file.                                        !
!                                                                                   !
! The number of integration points has been determined for the azimuthal orders     !
! m = 1 and m = - 1. For higher azimuthal orders and sources distributed in the     !
! complex plane, the chosen value of Nint may lead to false convergence. Therefore, !
! the number of integration points should be increased at each azimuthal mode       !
! calculation. The integation step is dNintMrank. The parameter dNintMrank should   !
! be found experimentally by repeating the azimuthal order test for different       !
! input values until the differential scattering cross section converges.           !
! A conservative estimate is: dNintMrank = (0.1,...,0.2) * Nint.                    !
!                                                                                   !
! 6. Input Parameters                                                               !
! --------------------                                                              !
! The parameters specified in the input file "/INPUTFILES/InputAXSYM.dat" are       !
! listed below.                                                                     !
!                                                                                   !
! - wavelength (real) - wavelength of the incident light in vacuo.                  !
!                                                                                   !
! - ind_refMed (real) - refractive index of the ambient medium (nonabsorbing        !
!   dielectric medium).                                                             !
!                                                                                   !
! - ind_refRel (complex) - relative refractive index of the particle with respect   !
!   to the ambient medium. The imaginary part of the relative refractive index      !
!   must be zero for nonabsorbing particles and positive for absorbing particles.   !
!                                                                                   !
! - perfectcond (logical) - if perfectcond = t, the particle is perfectly           !
!   conducting.                                                                     !
!                                                                                   !
! - chiral (logical) - if chiral = t, the particle is optical active (chiral).      !
!                                                                                   !
!   The permissive values of the logical variables perfectcond and chiral are       !
!   given below.                                                                    !
!                                                                                   !
!       particle              perfectcond            chiral                         !
!      dielectric                  f                   f                            !
!   perfectly conducting           t                   f                            !
!        chiral                    f                   t                            !
!                                                                                   !
! - kb (real) - parameter of chirality. This parameter is used for chiral           !
!   particles and is ignored for dielectric and perfectly conducting particles.     !
!   The wave numbers for the left- and right-handed circularly polarized waves      !
!   are given by:                                                                   !
!               kil = ki / (1.0 - kb)  and   kir = ki / (1.0 + kb)                  !
!   where ki is the wave number in the chiral medium.                               !
!                                                                                   !
! - FileGeom (logical) - if FileGeom = t, the particle geometry is supplied by the  !
!   input file FileFEM,                                                             ! 
!                                                                                   !
! - FileFEM (character(80)) - name of the file containing the particle geometry.    !
!                                                                                   !
! - TypeGeom (integer) - parameter specifying the type of the particle geometry.    !
!                                                                                   !
! - Nsurf (integer) - number of surface parameters.                                 !
!                                                                                   !
! - surf (real array: surf(1), surf(2),...,surf(Nsurf)) - surface parameters        !
!   specifying the shape of the particle. The dimension of the array surf is        !
!   NsurfPD. The integer parameter NsurfPD is specified in the routine              !
!   "Parameters.f90" and has the value NsurfPD = 10. If Nsurf > NsurfPD, the        !
!   execution is automatically terminated.                                          !
!                                                                                   !
! - Nparam (integer) - number of smooth curves forming the generatrix curve.        !
!                                                                                   !
!   The permissive values of the surface parameters (for the supported geometries)  !
!   are summarized below.                                                           !
!                                                                                   !
!   Particle      TypeGeom   Nsurf   Nparam                surf                     !
!   spheroid         1         2       1         surf(1) - length of the semi-      !
!                                                          axis along the           !
!                                                          axis of symmetry         !
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
! - anorm (real) - characteristic length of the particle which is used to           !
!   normalize the differential scattering cross sections.                           !
!                                                                                   !
! - Rcirc (real) - characteristic length of the particle (usually the radius of     !
!   the smallest circumscribing sphere) which is used to compute an estimate of the !
!   maximum expansion order by using Wiscombe's truncation limit criterion (the     !
!   size parameter is x = k * Rcirc, where k is the wave number in the ambient      !
!   medium). Alternatively, Rcirc can be chosen as the equal-volume sphere radius   !
!   or the surface-equivalent-sphere radius. This parameter is used if the          !
!   interactive convergence tests over Nint and Nrank are performed                 !
!  (DoConvTest = t). Otherwise, Rcirc is ignored.                                   !
!                                                                                   !
! - miror (logical) - if miror = t, the particle is mirror symmetric (the plane of  !
!   symmetry or the plane of reflection is perpendicular to the axis of rotation).  !
!                                                                                   ! 
! NOTE: FOR CHIRAL PARTICLES AND DISTRIBUTED SOURCES SET miror = f.                 ! 
!                                                                                   !
! - DoConvTest (logical) - if DoConvTest = t, the interactive convergence tests     !
!   over Nint and Nrank are invoked. Estimates of Nint and Nrank can be obtained    !
!   by using Mishchenko's convergence procedure. Alternatively, an estimate of      !
!   Nrank is given by Wiscombe's truncation limit criterion.                        !
!                                                                                   !
! - MishConvTest (logical) - if MishConvTest = t, estimates of Nint and Nrank are   !
!   computed with the convergence procedure proposed by Mishchenko. The code        !
!   prompts for the integer parameter Ndgs which controls the initial value of the  !
!   number of integration points. Some orientative values of Ndgs are printed and   !
!   the user can enter his option. This convergence test can be performed if:       !
!   - localized sources are used, the interactive convergence tests are performed,  !
!     the particle is not optical active, and the particle geometry is not read from!
!     file (DS = f and DoConvTest = t and chiral = f and FileGeom = f), or          !
!   - distributed sources are used, the coordinates of the distributed sources are  !
!     not user-defined, the interactive convergence tests are performed, the        !
!     particle is not optical active, and the particle geometry is not read from    ! 
!     file (DS = t and autGenDS = t and DoConvTest = t and chiral = f and           !
!     FileGeom = f).                                                                !
!                                                                                   !
! NOTE: IF THE PARTICLE IS OPTICAL ACTIVE (chiral = t) OR THE PARTICLE GEOMETRY IS  !
! SUPPLIED BY THE FILE FileFEM (FileGeom = t), THE CODE SETS MishConvTest = f.      !
!                                                                                   !
! - DS (logical) - if DS = t, distributed sources are used for T-matrix             !
!   calculation, otherwise localized sources are employed.                          !
!                                                                                   !
! - autGenDS (logical) - if autGenDS = t, the routine "zDSAXSYM" generates the      !
!   coordinates of the distributed sources. Otherwise, the coordinates of the       !
!   distributed sources must be provided by the user. If DS = f, the logical        !
!   variable autGenDS is ignored. If the particle geometry is read from file        !  
!   (FileGeom = t), the coordinates of the distributed sources are user-defined,    !
!   and the code sets autgenDS = f.                                                 !
!                                                                                   !
! - ComplexPlane (logical) - if ComplexPlane = t, the distributed sources are       !
!   situated in the complex plane. This parameter is used if distributed sources    !
!   are required and the coordinates of the distributed sources are automatically   !
!   generated (DS = t and autGenDS = t).                                            !
!                                                                                   !
! - epsZReIm (real) - input parameter of the routine "zDSAXSYM" which controls      !
!   the distribution of the discrete sources. This parameter is used if distributed !
!   sources are required and the coordinates of the distributed sources are         !
!   automatically generated (DS = t and autGenDS = t).                              !
!                                                                                   !
! - Nint (integer) - number of integration points in computing integrals over       !
!   the generatrix curve. This parameter is used (for localized and distributed     !
!   sources) if the convergence tests are not performed and the particle geometry   !
!   is not supplied by the file FileFEM. More specifically, Nint is used if         !
!   - (DS = f and DoConvTest = f and FileGeom = f) or,                              !
!   - (DS = t and autGenDS = t and DoConvTest = f) or,                              !
!   - (DS = t and autGenDS = f and DoConvTest = f and FileGeom = f).                !
!                                                                                   !
! - Nrank (integer) - maximum expansion order. This parameter is used (for localized!
!   and distributed sources) if the convergence tests are not performed, or the     !
!   coordinates of the distributed sources are user-defined. More specifically,     !
!   Nrank is used if                                                                !
!   - (DS = f and DoConvTest = f), or                                               !
!   - (DS = t and autGenDS = t and DoConvTest = f), or                              !
!   - (DS = t and autGenDS = f and DoConvTest = t), or                              !
!   - (DS = t and autGenDS = f and DoConvTest = f).                                 !
!                                                                                   !
! NOTE: IN ORDER TO SIMPLIFY THE ROUTINE FOR READING THE DATA, THE VARIABLES Nint   !
! AND Nrank MUST BE PROVIDED (IN THE INPUT FILE "InptAXSYM.dat") IF                 !
! ((DoConvTest = f) OR (DS = t AND autgenDS = f)). BECAUSE THE SIGNIFICANCE OF      !
! THESE PARAMETERS IS AS ABOVE, THE MAIN PROGRAM WILL IGNORE THE IRRELEVANT DATA    !
! INFORMATIONS.                                                                     !
!                                                                                   !
! - zRe, zIm (real arrays: zX(1), zX(2),...,zX(Nrank), X = Re, Im) - coordinates    !
!   of the distributed sources for the expansion order Nrank. These parameters      !
!   are used if the coordinates of the distributed sources are user-defined         !
!   (DS = t and autGenDS = f). The dimension of the arrays zRe and zIm is NrankPD   !
!   and the inequality Nrank <= NrankPD must hold. The integer parameter NrankPD is !
!   specified in the routine "Parameters.f90" and has the value NrankPD = 200. If   !
!   Nrank > NrankPD, the execution is automatically terminated.                     !
!                                                                                   !
! - zRe1, zIm1 (real arrays: zX1(1), zX1(2),...,zX1(Nrank-1), X = Re, Im) -         !
!   coordinates of the distributed sources for the expansion order Nrank - 1.       !
!   These parameters are used if the coordinates of the distributed sources are     !
!   user-defined and the expansion order test is performed (DS = t and autGenDS = f !
!   and TypeConvTest = 2). As before, the dimension of the arrays zRe1 and zIm1     !
!   is NrankPD.                                                                     !
!                                                                                   !
! NOTE: IN ORDER TO SIMPLIFY THE ROUTINE FOR READING THE DATA, THE REAL ARRAYS      !
! zRe, zIm AND zRe1, zIm1 MUST BE SPECIFIED (IN THE INPUT FILE "InptAXSYM.dat") IF  !
! (DS = t AND autGenDS = f). FURTHER NOTE THAT THE INPUT ARRAYS MUST BE SEPARATED   !
! BY A BLANK LINE, AND THAT THE INPUT ARRAYS zRe1 AND zIm1 MUST BE SPECIFIED AS     !
! VECTORS WITH Nrank - 1 COMPONENTS.                                                !
!                                                                                   !
! - epsNint (real) - error tolerance for the integration test.                      !
!                                                                                   !
! - epsNrank (real) - error tolerance for the expansion order test.                 !
!                                                                                   !
! - epsMrank (real) - error tolerance for the azimuthal order test.                 !
!                                                                                   !
! - dNint (integer) - number of division points for the integration test and        !
!   Mishchenko's convergence procedure. Note that the scattering problem is solved  !
!   for Nint and Nint + dNint.                                                      ! 
!                                                                                   !
! - dNintMrank (integer) - number of division points for azimuthal mode             !
!   calculation. This parameter is used if the distributed sources are situated     !
!   in the complex plane (DS = t and ComplexPlane = t). At each azimuthal mode      !
!   calculation, the number of integration points increases with dNintMrank.        !
!                                                                                   !
! - FileTmat (character(80)) - name of the file to which the T matrix is written.   !
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
! 1.  < read input data >                                                           !
! 2.  < T-matrix calculation >                                                      !
!                                                                                   !
!    if ( read geometry from file ) then                                            !
!       < do not use Mishchenko's convergence procedure >                           !
!       < do not use automatic sources generation >                                 !
!    end if                                                                         !
!    if ( chiral particle ) then                                                    !
!       < do not use Mishchenko's convergence procedure >                           !
!    end if                                                                         !
!    if ( use localized sources ) then                                              !
!                                                                                   !
!       if ( do convergence test) then                                              !
!          if ( use Mishchenko's convergence procedure ) then                       !
!             < the code prompts for the integer parameter Ndgs and computes        !
!               the estimated values of Nint and Nrank by using Mishchenko's        !
!               convergence procedure >                                             !
!          else                                                                     !
!             < the code computes an estimate of Nrank by using                     !
!               Wiscombe's truncation criterion >                                   !
!          end if                                                                   !  
!          if ( .not. read geometry from file) then                                 !
!            < the code prompts for the estimated values of Nint and Nrank >        !
!            < the code prompts for the type of convergence test:                   !
!              1 - over Nint, 2 - over Nrank or 3 - over Mrank >                    !
!          else                                                                     !
!            < the code prompts for the estimated value of Nrank >                  !
!            < the code prompts for the type of convergence test:                   !
!              2 - over Nrank or 3 - over Mrank >                                   !           
!          end if                                                                   !
!       else if ( .not. do convergence test) then                                   !
!          if ( .not. read geometry from file) then                                 !
!             < the input file provides the values of Nint and Nrank >              !
!          else                                                                     ! 
!             < the input file provides the value of Nrank >                        !
!          end if                                                                   !
!          type of convergence test = 3 (convergence test over Mrank)               !
!       end if                                                                      !
!                                                                                   !
!    else if ( use distributed sources ) then                                       !
!                                                                                   !
!       if ( use automatic sources generation ) then                                !
!          if ( do convergence test) then                                           !
!             if ( use Mishchenko's convergence procedure ) then                    !
!                < the code prompts for the integer parameter Ndgs and computes     !
!                  the estimated values of Nint and Nrank by using Mishchenko's     !
!                  convergence procedure >                                          ! 
!             else                                                                  !
!                < the code computes an estimate of Nrank by using                  !
!                  Wiscombe's truncation criterion >                                !
!             end if                                                                !
!             < the code prompts for the estimated values of Nint and Nrank >       !
!             < the code generates the coordinates of the distributed sources       !
!               zRe and zIm corresponding to Nrank >                                !
!             < the code prompts for the type of convergence test:                  !
!               1 - over Nint, 2 - over Nrank or 3 - over Mrank >                   !
!             if ( do maximum expansion order test ) then                           !
!                < the code generate the coordinates of the distributed sources     !
!                  zRe1 and zIm1 corresponding to Nrank - 1 >                       !
!             end if                                                                !
!          else if ( .not. do convergence test) then                                !
!             < the input file provides the values of Nint and Nrank >              !
!             < the code generate the coordinates of the distributed sources        !
!               zRe and zIm corresponding to Nrank >                                !
!             type of convergence test = 3 (convergence test over Mrank)            ! 
!          end if                                                                   !
!       else if ( .not. use automatic sources generation ) then                     !
!          if ( do convergence test) then                                           !
!             < the input file provides the value of Nrank >                        !
!             < the input file provides the coordinates of the distributed sources  !
!               zRe and zIm corresponding to Nrank >                                !
!             if ( .not. read geometry from file ) then                             !
!                < the code prompts for the estimated value of Nint >               !
!                < the code prompts for the type of convergence test:               !
!                  1 - over Nint, 2 - over Nrank or 3 - over Mrank >                !
!             else                                                                  !
!                < the code prompts for the type of convergence test:               !
!                  2 - over Nrank or 3 - over Mrank >                               !
!             end if                                                                !
!             if ( do expansion order test ) then                                   !
!                < the input file provides the coordinates of the distributed       !
!                  sources zRe1 and zIm1 corresponding to Nrank - 1 >               !
!             end if                                                                !
!          else if ( .not. do convergence test) then                                !
!             if ( .not. read geometry from file ) then                             !
!                < the input file provides the values of Nint and Nrank >           !
!             else                                                                  !
!                < the input file provides the value of Nrank >                     !
!             end if                                                                ! 
!             < the input file provides the coordinates of the distributed sources  !
!               zRe and zIm corresponding to Nrank >                                !
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
!       < the code computes the DSCS for Nrank and Nrank - 1 and                    !
!         write the results to the file "Output.dat" >                              !
!       < the execution is terminated >                                             !
!    end if                                                                         !
!    if ( do azimuthal order test, i.e., type of convergence test = 3 ) then        !
!       < the code computes the DSCS for increasing m values and                    !
!         write the results to the file "Output.dat" >                              !
!       < the T matrix is stored in the file FileTmat and the T-matrix information  !
!         file is created >                                                         !
!       < Mrank and Nrank are printed to the screen and to the T-matrix             !
!         information file >                                                        !
!       < go to Step 3 >                                                            !
!    end if                                                                         !
! 3. < scattering characteristics calculation >                                     !
!------------------------------------------------------------------------------------
  use parameters  
  implicit none 
  integer       :: TypeGeom, Nsurf, Nparam, TypeConvTest, dNint, Nrank, Nint,       &
                   Nrank1, dNintMrank, Ndgs, NrankMax, NrankW, Nface                     
  real(O)       :: k, ind_refMed, wavelength, anorm, surf(NsurfPD), snorm,          &
                   kb, epsNint, epsNrank, epsMrank, zRe(NrankPD), zIm(NrankPD),     &
                   zRe1(NrankPD), zIm1(NrankPD), Rcirc, x, delta, Cscat1, Cext1,    &
                   epsZReIm, rp(2,NfacePD), np(2,NfacePD), area(NfacePD)                          
  complex(O)    :: ind_refRel
  logical       :: FileGeom, miror, perfectcond, DoConvTest, DS, chiral, autGenDS,  &
                   ComplexPlane, MishConvTest, PrnProgress
  character(80) :: FileTmat, FileFEM
! -----------------------------------------------------------------------------------
!                               Read the input file                                 ! 
! -----------------------------------------------------------------------------------       
  call readinputAXSYM ( wavelength, ind_refMed, ind_refRel, perfectcond,            &
       chiral, kb, FileGeom, TypeGeom, FileFEM, Nsurf, surf, Nparam, anorm,         &
       Rcirc, miror, DoConvTest, MishConvTest, DS, autGenDS, ComplexPlane,          &
       epsZReIm, Nint, Nrank, zRe, zIm, zRe1, zIm1, epsNint, epsNrank, epsMrank,    &
       dNint, dNintMrank, FileTmat, PrnProgress, k, snorm, Nface, rp, np, area )                   
! -----------------------------------------------------------------------------------
!         Select the type of convergence test and the values of Nint and Nrank      !
! -----------------------------------------------------------------------------------
  x = k * Rcirc
  NrankW = int(x + 4.05_O * x**0.33_O + 2._O)
  if (.not. DS) then
    if (DoConvTest) then                      
      if (MishConvTest) then                         
        print "(/,2x, a)",                                                          &
       'Estimates of Nint and Nrank Using Mishchenko''s Convergence Procedure' 
        print "(  2x, a)",                                                          &
       '---------------------------------------------------------------------'                                         
        call estimateNrankMishchenko (TypeGeom, k, ind_refRel, Nsurf, surf,         &
             zRe, zIm, Nparam, miror, perfectcond, DS, ComplexPlane,                &
             epsZReIm, x, delta, Ndgs, Nint, Nrank, NrankMax, Cscat1, Cext1)
        call estimateNintMishchenko (TypeGeom, k, ind_refRel, Nsurf, surf,          &
             zRe, zIm, Nparam, miror, perfectcond, DS, x, delta, Ndgs, Nint,        &
             dNint, Nrank, NrankMax, Cscat1, Cext1)
        print "(/,2x,'Convergence Test for an Axisymmetric Particle')"
        print "(  2x,'---------------------------------------------')"	       
        print "(/,2x,'- enter the estimated values of Nint and Nrank;')"            
        call read_integer2 (Nint, Nrank)                                
      else 
        print "(/,2x,'Convergence Test for an Axisymmetric Particle')"
        print "(  2x,'---------------------------------------------')"             
        print "(/,2x,'Nrank estimate:')"                                                            
        print "(  2x, a, i3, a)",                                                   &  
       'the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';'
        if (.not. FileGeom) then
          print "(/,2x, a)",                                                        &
         '- enter the estimated values of Nint and Nrank, where Nint = Ndgs * Nrank'
          print "(  2x,'  and Ndgs = 3,4,...;')"
          call read_integer2 (Nint, Nrank)                                                                      
        else
          print "(/,2x,'- enter the estimated value of Nrank;')"              
          call read_integer (Nrank)     
        end if  
      end if
      if (.not. FileGeom) then        
        print "(/,2x, a)",                                                          &
       '- enter the type of convergence test: 1 - Nint, 2 - Nrank, 3 - Mrank;' 
        call read_integerbound (TypeConvTest, 1, 3)
      else
        print "(/,2x,'- enter the type of convergence test: 2 - Nrank, 3 - Mrank;')"
        call read_integerbound (TypeConvTest, 2, 3)
      end if
    else
      print "(/,2x,'Convergence Test for an Axisymmetric Particle over Mrank')" 
      print "(  2x,'--------------------------------------------------------')"
      print "(/,2x,'Input values:')"
      if (.not. FileGeom) then
        print "(  2x, a, i4, a, i4, a)",                                            & 
       'the input values of Nint and Nrank are ', Nint, ' and ', Nrank,             &
       ', respectively,'  
        print "(  2x, a, i3, a)",                                                   &
       'while the estimated value of Nrank from Wiscombe''s criterion is ', NrankW,';'
      else
        print "(  2x,'the input value of Nrank is ', i4,', while')", Nrank
        print "(  2x, a, i3, a)",                                                   &
       'the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';' 
      end if                                                                          
      TypeConvTest = 3      
    end if
    Nrank1 = Nrank - 1 ! redundant            
  else
    if (autGenDS) then
      if (DoConvTest) then                        
        if (MishConvTest) then           
          print "(/,2x, a)",                                                        &
         'Estimates of Nint and Nrank Using Mishchenko''s Convergence Procedure'
          print "(  2x, a)",                                                        &
         '---------------------------------------------------------------------'	  	 	 	                    
          call estimateNrankMishchenko (TypeGeom, k, ind_refRel, Nsurf, surf,       &
               zRe, zIm, Nparam, miror, perfectcond, DS, ComplexPlane,              &
               epsZReIm, x, delta, Ndgs, Nint, Nrank, NrankMax, Cscat1, Cext1)
          call estimateNintMishchenko (TypeGeom, k, ind_refRel, Nsurf, surf,        &
               zRe, zIm, Nparam, miror, perfectcond, DS, x, delta, Ndgs, Nint,      &
               dNint, Nrank, NrankMax, Cscat1, Cext1)   
          print "(/,2x,'Convergence Test for an Axisymmetric Particle')"
          print "(  2x,'---------------------------------------------')"	       
          print "(/,2x,'- enter the estimated values of Nint and Nrank;')"           
          call read_integer2 (Nint, Nrank)                               
        else
          print "(/,2x,'Convergence Test for an Axisymmetric Particle')"
          print "(  2x,'---------------------------------------------')"	
          print "(/,2x,'Nrank estimate:')"                                                                
          print "(  2x, a, i3, a)",                                                 &
         'the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';'
          print "(/,2x, a)",                                                        &
         '- enter the estimated values of Nint and Nrank, where Nint = Ndgs * Nrank'
          print "(  2x,'  and Ndgs = 5,6,...;')"
          call read_integer2 (Nint, Nrank)                                                                                 
        end if      
        call check_MaxNrank (Nrank)             
        call zDSAXSYM (TypeGeom, Nsurf, surf, Nrank, ComplexPlane, epsZReIm, zRe, zIm)     
        print "(/,2x, a)",                                                          &
       '- enter the type of convergence test: 1 - Nint, 2 - Nrank, 3 - Mrank;'
        call read_integerbound (TypeConvTest, 1, 3)
        if (TypeConvTest == 2) then
          Nrank1 = Nrank - 1
          call zDSAXSYM (TypeGeom, Nsurf, surf, Nrank1, ComplexPlane, epsZReIm,     &
               zRe1, zIm1)         
        end if
      else
        print "(/,2x,'Convergence Test for an Axisymmetric Particle over Mrank')" 
        print "(  2x,'--------------------------------------------------------')"
        print "(/,2x,'Input values:')"
        print "(  2x, a, i4, a, i4, a)",                                            &
       'the input values of Nint and Nrank are ', Nint, ' and ', Nrank,             &
       ', respectively,' 
        print "(  2x, a, i3, a)",                                                   &
       'while the estimated value of Nrank from Wiscombe''s criterion is ', NrankW,';'         
        call check_MaxNrank (Nrank)
        call zDSAXSYM (TypeGeom, Nsurf, surf, Nrank, ComplexPlane, epsZReIm, zRe, zIm) 
        TypeConvTest = 3    
      end if
      Nrank1 = Nrank - 1 ! redundant  
    else 
      if (DoConvTest) then      
        print "(/,2x,'Convergence Test for an Axisymmetric Particle')"
        print "(  2x,'---------------------------------------------')"
        print "(/,2x,'Input values:')"
        print "(  2x,'the input value of Nrank is ', i4,', while')", Nrank
        print "(  2x, a, i3, a)",                                                   &
       'the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';'
        if (.not. FileGeom) then                                            
          print "(/,2x, a)",                                                        &
         '- enter the estimated value of Nint, where Nint = Ndgs * Nrank'
          print "(  2x,'  and Ndgs = 5,6,...;')"
          call read_integer (Nint) 
          print "(/,2x, a)",                                                        &
         '- enter the type of convergence test: 1 - Nint, 2 - Nrank, 3 - Mrank;'
          call read_integerbound (TypeConvTest, 1, 3)
        else
          print "(/,2x,'- enter the type of convergence test: 2 - Nrank, 3 - Mrank;')"
          call read_integerbound (TypeConvTest, 2, 3)
        end if
        if (TypeConvTest == 2) Nrank1 = Nrank - 1
      else
        print "(/,2x,'Convergence Test for an Axisymmetric Particle over Mrank')" 
        print "(  2x,'--------------------------------------------------------')" 
        print "(/,2x,'Input values:')" 
        if (.not. FileGeom) then                  
          print "(  2x, a, i4, a, i4, a)",                                          &
         'the input values of Nint and Nrank are ', Nint, ' and ', Nrank,           &
         ', respectively,' 
          print "(  2x, a, i3, a)",                                                 &
         'while the estimated value of Nrank from Wiscombe''s criterion is ',       &
          NrankW,';' 
        else
          print "(  2x,'the input value of Nrank is ', i4,', while')", Nrank
          print "(  2x, a, i3, a)",                                                 &
         'the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';' 
        end if
        TypeConvTest = 3              
      end if
      Nrank1 = Nrank - 1 ! redundant          
    end if
  end if                                                                            
! -----------------------------------------------------------------------------------
!                               Convergence test                                    !
! -----------------------------------------------------------------------------------
  open (unit = iOutput, file = FileOutput, status = "replace") 
  call printinputAXSYM (TypeConvTest, FileGeom, TypeGeom, FileFEM, Nsurf, Nparam,   &
       Nrank, Nrank1, dNint, dNintMrank, ind_refMed, wavelength, anorm, Rcirc,      &
       surf, kb, epsNint, epsNrank, epsMrank, zRe, zIm, zRe1, zIm1, ind_refRel,     &
       miror, perfectcond, DS, chiral, autGenDS)                     
  if (DoConvTest) then              
    if (TypeConvTest == 1) then
      if (.not. DS) then                  
        call convergence_NintAXSYM (FileGeom, TypeGeom, k, ind_refRel, snorm,       &
             Nsurf, surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint,       &
             dNint, miror, perfectcond, DS, chiral, kb, epsNint, PrnProgress)
      else 
        call convergence_NintDSAXSYM (FileGeom, TypeGeom, k, ind_refRel, snorm,     &
             Nsurf, surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint,       &
             dNint, miror, perfectcond, DS, chiral, kb, epsNint, PrnProgress)
      end if
    else if (TypeConvTest == 2) then
      if (.not. DS) then        
        call convergence_NrankAXSYM (FileGeom, TypeGeom, k, ind_refRel, snorm,      &
             Nsurf, surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint,       &
             miror, perfectcond, DS, chiral, kb, epsNrank, PrnProgress)
      else      
        call convergence_NrankDSAXSYM (FileGeom, TypeGeom, k, ind_refRel, snorm,    &
             Nsurf, surf, rp, np, area, Nface, zRe, zIm, zRe1, zIm1, Nparam, Nrank, &
             Nrank1, Nint, miror, perfectcond, DS, chiral, kb, epsNrank,            &
             PrnProgress)
      end if
    else 
      if (.not. DS) then    
        call convergence_MrankAXSYM (FileGeom, TypeGeom, k, ind_refRel, snorm,      &
             Nsurf, surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint,       &
             miror, perfectcond, DS, chiral, kb, epsMrank, FileTmat, PrnProgress)
      else 
        call convergence_MrankDSAXSYM (FileGeom, TypeGeom, k, ind_refRel, snorm,    &
             Nsurf, surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint,       &
             miror, perfectcond, DS, chiral, kb, epsMrank, dNintMrank, FileTmat,    &
             PrnProgress)
      end if
    end if  
  else      
    if (.not. DS) then    
      call convergence_MrankAXSYM (FileGeom, TypeGeom, k, ind_refRel, snorm, Nsurf, &
           surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint, miror,         &
           perfectcond, DS, chiral, kb, epsMrank, FileTmat, PrnProgress)
    else          
      call convergence_MrankDSAXSYM (FileGeom, TypeGeom, k, ind_refRel, snorm,      &
           Nsurf, surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint, miror,  &
           perfectcond, DS, chiral, kb, epsMrank, dNintMrank, FileTmat, PrnProgress)
    end if  
  end if 
  close (unit = iOutput)
end subroutine TAXSYM
!************************************************************************************
subroutine readinputAXSYM ( wavelength, ind_refMed, ind_refRel, perfectcond,        &
           chiral, kb, FileGeom, TypeGeom, FileFEM, Nsurf, surf, Nparam, anorm,     &
           Rcirc, miror, DoConvTest, MishConvTest, DS, autGenDS, ComplexPlane,      &
           epsZReIm, Nint, Nrank, zRe, zIm, zRe1, zIm1, epsNint, epsNrank,          &
           epsMrank, dNint, dNintMrank, FileTmat, PrnProgress, k, snorm, Nface,     &
           rp, np, area )
  use parameters
  use derived_parameters
  implicit none 
  integer       :: TypeGeom, Nsurf, Nparam, dNint, Nrank, Nint, i, dNintMrank, ios, &
                   Nface, j                     
  real(O)       :: k, ind_refMed, wavelength, anorm, surf(NsurfPD), xpart, snorm,   &
                   kb, epsNint, epsNrank, epsMrank, zRe(NrankPD), zIm(NrankPD),     &
                   zRe1(NrankPD), zIm1(NrankPD), Rcirc, epsZReIm, rp(2,NfacePD),    &
                   np(2,NfacePD), area(NfacePD), dp                          
  complex(O)    :: ind_refRel
  logical       :: FileGeom, miror, perfectcond, DoConvTest, DS, chiral, autGenDS,  &
                   ComplexPlane, InputDS, MishConvTest, PrnProgress, XFindPar
  character(80) :: FileTmat, FileFEM, string       
! -----------------------------------------------------------------------------------
!                        Read the input file FileInputAXSYM                         ! 
! -----------------------------------------------------------------------------------    
  call DrvParameters  
  open (unit = iInputAXSYM, file = FileInputAXSYM, status = "old",                  &
        position = "rewind")   
  wavelength = 0.1_O * 2._O * Pi  
  ind_refMed = 1._O                                                                 
  ind_refRel = (1.5_O,0._O) 
  string     = 'OptProp'    
  if (XFindPar (iInputAXSYM, string)) then
    read (iInputAXSYM, *, iostat = ios) wavelength
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable wavelength;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) ind_refMed
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refMed;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) ind_refRel
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refRel;')"
      stop
    end if      
  else
    print "(/,2x,'Group name OptProp not found;')"
    stop  
  end if
  call check_ind_ref (ind_refRel)
  k = 2._O * Pi * ind_refMed / wavelength 
!
  perfectcond = .false.  
  chiral = .false.
  kb     = 0._O  
  string = 'MatProp'
  if (XFindPar (iInputAXSYM, string)) then
    read (iInputAXSYM, *, iostat = ios) perfectcond
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable perfectcond;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) chiral
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable chiral;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) kb
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable kb;')"
      stop
    end if      
  else
    print "(/,2x,'Group name MatProp not found;')"
    stop  
  end if      
  call check_MatPropAXSYM (perfectcond, chiral, kb)   
  if (chiral) call check_chirality (kb)
!
  FileGeom = .false.
  FileFEM  = ' '  
  TypeGeom = 1  
  Nsurf    = 2
  do i = 1, NsurfPD
    surf(i)= 1._O
  end do
  Nparam = 1
  anorm  = 1._O
  Rcirc  = 1._O 
  miror  = .false.
  string = 'GeomProp'
  if (XFindPar (iInputAXSYM, string)) then
    read (iInputAXSYM, *, iostat = ios) FileGeom
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable FileGeom;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) FileFEM
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable FileFEM;')"
      stop
    end if       
    read (iInputAXSYM, *, iostat = ios) TypeGeom
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable TypeGeom;')"
      stop
    end if       
    read (iInputAXSYM, *, iostat = ios) Nsurf
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nsurf;')"
      stop
    end if
    if (Nsurf > NsurfPD) then
      print "(/,2x,'Input error: Nsurf exceeds NsurfPD;')"                                    
      stop
    end if            
    do i = 1, Nsurf
      read (iInputAXSYM, *, iostat = ios) surf(i)
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable surf;')"
        stop
      end if
    end do    
    read (iInputAXSYM, *, iostat = ios) Nparam
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nparam;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) anorm
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable anorm;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) Rcirc
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Rcirc;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) miror
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable miror;')"
      stop
    end if                                                
  else
    print "(/,2x,'Group name GeomProp not found;')"
    stop  
  end if                          
  call check_geomAXSYM (TypeGeom, Nsurf, Nparam)
  call check_geomAXSYMOblate (TypeGeom, Nsurf, surf)  
  call check_anorm (anorm)
  xpart = k * anorm
  snorm = Pi * xpart * xpart 
  Nface = 1
  do i = 1, NfacePD
    do j = 1, 2
      rp(j,i) = 0._O
      np(j,i) = 0._O
    end do
    area(i) = 0._O
  end do   
  if (FileGeom) then    
    call read_FileFEMAxsym (FileFEM, Nface, rp, np, area)     
    Rcirc = 0._O    
    do i = 1, Nface
      dp = sqrt(rp(1,i)**2 + rp(2,i)**2)
      if (dp > Rcirc) Rcirc = dp
    end do    
  end if          
!
  DoConvTest   = .true.
  MishConvTest = .true.
  string       = 'ConvTest'
  if (XFindPar (iInputAXSYM, string)) then
    read (iInputAXSYM, *, iostat = ios) DoConvTest
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable DoConvTest;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) MishConvTest
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable MishConvTest;')"
      stop
    end if         
  else
    print "(/,2x,'Group name ConvTest not found;')"
    stop  
  end if   
  if (chiral)   MishConvTest = .false.     
  if (FileGeom) MishConvTest = .false.  
!
  DS = .false.
  autGenDS = .true.
  string   = 'Sources'
  if (XFindPar (iInputAXSYM, string)) then
    read (iInputAXSYM, *, iostat = ios) DS
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable DS;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) autGenDS
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable autGenDS;')"
      stop
    end if         
  else
    print "(/,2x,'Group name Sources not found;')"
    stop  
  end if          
  call check_inputAXSYM (miror, chiral, DS)
  if (FileGeom) autGenDS = .false.  
  InputDS = .false.
  if (DS .and. .not. autGenDS) InputDS = .true. 
!
  ComplexPlane = .false.
  epsZReIm = 0.95_O 
  if (DS .and. autGenDS) then
    string = 'SourcePosAut'
    if (XFindPar (iInputAXSYM, string)) then
      read (iInputAXSYM, *, iostat = ios) ComplexPlane
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable ComplexPlane;')"
        stop
      end if
      read (iInputAXSYM, *, iostat = ios) epsZReIm
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable epsZReIm;')"
        stop
      end if         
    else
      print "(/,2x,'Group name SourcePosAut not found;')"
      stop  
    end if              
  end if
!
  Nint   = 100
  Nrank  =  17 
  if (.not. DoConvTest .or. InputDS) then       
    string = 'NintNrank'
    if (XFindPar (iInputAXSYM, string)) then
      read (iInputAXSYM, *, iostat = ios) Nint
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Nint;')"
        stop
      end if
      read (iInputAXSYM, *, iostat = ios) Nrank
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Nrank;')"
        stop
      end if         
    else
      print "(/,2x,'Group name NintNrank not found;')"
      stop  
    end if                              
  end if    
!
  do i = 1, NrankPD
    zRe(i)  = 0._O
    zIm(i)  = 0._O
    zRe1(i) = 0._O
    zIm1(i) = 0._O
  end do
  if (InputDS) then   
    call check_MaxNrank (Nrank) 
    string = 'SourcePosInp'
    if (XFindPar (iInputAXSYM, string)) then
      do i = 1, Nrank
        read (iInputAXSYM, *, iostat = ios) zRe(i), zIm(i)
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variables zRe and zIm;')"
          stop
        end if
      end do      
      read (iInputAXSYM, *)
      do i = 1, Nrank - 1
        read (iInputAXSYM, *, iostat = ios) zRe1(i), zIm1(i)
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variables zRe1 and zIm1;')"
          stop
        end if
      end do
    else
      print "(/,2x,'Group name SourcePosInp not found;')"
      stop  
    end if                       
  end if         
!
  epsNint  = 5.e-2_O
  epsNrank = 5.e-2_O
  epsMrank = 5.e-2_O
  dNint    = 4
  dNintMrank = 10
  string   = 'Errors'
  if (XFindPar (iInputAXSYM, string)) then
    read (iInputAXSYM, *, iostat = ios) epsNint
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsNint;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) epsNrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsNrank;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) epsMrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsMrank;')"
      stop
    end if 
    read (iInputAXSYM, *, iostat = ios) dNint
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable dNint;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) dNintMrank
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
  if (XFindPar (iInputAXSYM, string)) then
    read (iInputAXSYM, *, iostat = ios) FileTmat
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
  if (XFindPar (iInputAXSYM, string)) then
    read (iInputAXSYM, *, iostat = ios) PrnProgress
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable PrnProgress;')"
      stop
    end if             
  else
    print "(/,2x,'Group name PrintProgress not found;')"
    stop  
  end if       
  close (unit = iInputAXSYM) 
end subroutine readinputAXSYM   
!************************************************************************************
subroutine printinputAXSYM (ic, FileGeom, TypeGeom, FileFEM, Nsurf, Nparam, Nrank,  &
           Nrank1, dNint, dNintMrank,ind_refMed, wavelength, anorm, Rcirc, surf,    &
           kb, epsNint, epsNrank, epsMrank, zRe, zIm, zRe1, zIm1, ind_refRel,       &
           miror, perfectcond, DS, chiral, autGenDS)
  use parameters
  implicit none
  integer       :: ic, TypeGeom, Nsurf, Nparam, Nrank, Nrank1, dNint, dNintMrank,   &
                   i, LenString                     
  real(O)       :: ind_refMed, wavelength, anorm, Rcirc, surf(Nsurf), kb, epsNint,  &
                   epsNrank, epsMrank, zRe(Nrank), zIm(Nrank), zRe1(Nrank1),        &
                   zIm1(Nrank1)
  complex(O)    :: ind_refRel
  character(80) :: FileFEM, FileFEMWrite
  logical       :: FileGeom, miror, perfectcond, DS, chiral, autGenDS  
!
  write (iOutput,"(/,2x,'Input Parameters of the Scattering Problem:',/)")
  write (iOutput,"(2x,'wavelength of the free space, wavelength = ',1pe13.4,';')")  &
         wavelength
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'refractive index of the ambient medium, ind_refMed = ', ind_refMed, ';'
  write (iOutput,"(2x, a, 1pe10.3, ',', 1pe10.3, a)")                               &
 'relative refractive index of the particle, ind_refRel = (', ind_refRel, ');' 
  write (iOutput,*)
  if (FileGeom) then
    FileFEMWrite = FileFEM(14:LenString(FileFEM))
    write (iOutput,"(2x, a, a30)")                                                  &
   'name of the file containing the particle geometry, FileFEM = ', FileFEMWrite
  else
    write (iOutput,"(2x,'index of the particle geometry, TypeGeom = ',i2,';')")     &
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
      write (iOutput,"(2x,'where surf(1) is the semi-axis along the symmetry axis,')")
      write (iOutput,"(2x,'and   surf(2) is the second semi-axis;')")
    else if (TypeGeom == 2 .or. TypeGeom == 3) then
      write (iOutput,"(2x,'where surf(1) is the half-length of the cylinder')") 
      write (iOutput,"(2x,'and   surf(2) is the cylinder radius;')")
    end if               
    write (iOutput,"(2x,'number of integration surfaces, Nparam = ',i2,';')") Nparam
  end if  
  write (iOutput,"(2x,a, 1pe10.3, a)")                                              &
 'characteristic length of the particle, anorm = ', anorm, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'radius of the circumscribing sphere, Rcirc = ', Rcirc, ';'
  if (miror) write (iOutput,"(2x,'mirror symmetric particle;')")  
  write (iOutput,*)
  if (perfectcond) then
    write (iOutput,"(2x,'perfectly conducting particle;')")
  else if (chiral) then
    write (iOutput,"(2x,'chiral particle;')")
    write (iOutput,"(2x,'characteristic of chirality, kb = ',1pe10.3,';')") kb  
  else 
    write (iOutput,"(2x,'dielectric particle;')")
  end if 
  write (iOutput,*)
  if (.not.DS) then
    write (iOutput,"(2x,'computation with localized sources;')")
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
    do i = 1, Nrank
      if (i /= Nrank) then
        write (iOutput,"(2x,'zRe(',i3,') = ',1pe10.3,', zIm(',i3,') = ',1pe10.3,',')")&
               i, zRe(i), i, zIm(i)
      else
        write (iOutput,"(2x,'zRe(',i3,') = ',1pe10.3,', zIm(',i3,') = ',1pe10.3,';')")&
               i, zRe(i), i, zIm(i)
      end if
    end do
    if (ic == 2) then
      write (iOutput,*)
      write (iOutput,"(2x,'second configuration:')")
      do i = 1,Nrank1
        if (i /= Nrank1) then
          write (iOutput,"(2x,'zRe(',i3,') = ',1pe10.3,', zIm(',i3,') = ',1pe10.3,',')")&
                 i, zRe1(i), i, zIm1(i)
        else
          write (iOutput,"(2x,'zRe(',i3,') = ',1pe10.3,', zIm(',i3,') = ',1pe10.3,';')")&
                 i, zRe1(i), i, zIm1(i)
        end if
      end do
    end if
  end if
  write (iOutput,*)
  write (iOutput,"(2x,'integration tolerance, epsNint = ',1pe10.3,';')") epsNint
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum expansion order tolerance, epsNrank = ', epsNrank, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum azimuthal order tolerance, epsMrank = ', epsMrank, ';'          
  write (iOutput,"(2x,'integration step, dNint = ',i3,'.')") dNint
  if (DS) then
    write (iOutput,"(2x, a, i3, a)")                                                &
   'integration step for Mrank calculation, dNintMrank = ',     dNintMrank, '.'
  end if
  write (iOutput,"(/)")                       
end subroutine printinputAXSYM
!***********************************************************************************
subroutine convergence_NintAXSYM (FileGeom, TypeGeom, k, ind_ref, snorm, Nsurf,     &
           surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint, dNint, miror,  &
           perfectcond, DS, chiral, kb, epsNint, PrnProgress)
  use parameters                                                                           
  implicit none 
  integer    :: TypeGeom, Nsurf, Nface, Nparam, Nrank, Nint, dNint
  real(O)    :: k, surf(Nsurf), snorm, epsNint, kb, zRe(Nrank), zIm(Nrank),         &
                rp(2,NfacePD), np(2,NfacePD), area(NfacePD)
  complex(O) :: ind_ref
  logical    :: FileGeom, miror, perfectcond, DS, chiral, PrnProgress
!      
  integer    :: Nteta, Mstart, Mrank, Nmax, Nmaxmax, i, m, iNint, NthetaConv,       &
                iprog, Nprog
  real(O)    :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,        &
                Cext, Qext
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG(:,:), weightsG(:,:), h(:), v(:), oldh(:), oldv(:)
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
  Nmax   = Nrank        
  m = 1            
  call write_TypeConvHead (1)
  allocate (a(2*Nrank,2*Nrank), b(2*Nrank,2*Nrank), c(2*Nrank), c1(2*Nrank))
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  allocate (cc(2*Nmaxmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))  
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do
  if (.not. chiral) then
    Nprog = 7
    iprog = 3
  else
    Nprog = 13
    iprog =  6
  end if
  if (PrnProgress) call write_progress (.true., 1, Nprog)
  do iNint = 1, 2         
    allocate (paramG(Nparam,Nint), weightsG(Nparam,Nint), Nintparam(Nparam))
    call interpolation_listAXSYM (TypeGeom, Nsurf, surf, Nint, Nparam,              &
         Nintparam, paramG, weightsG, miror)       
    call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,    &
         weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
    if (PrnProgress) call write_progress (.false., 2+iprog*(iNint-1), Nprog)
    call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,    &
         weightsG, miror, perfectcond, DS, chiral, kb, b, Nrank, Nrank)
    if (PrnProgress) call write_progress (.false., 3+iprog*(iNint-1), Nprog)
    call LU_SYSTEM (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nmax)
    if (PrnProgress) call write_progress (.false., 4+iprog*(iNint-1), Nprog)
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,     &
         Nmax, c)
    call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, c, c1)
    call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)
    if (.not. chiral) then
      call matrix_m_negativ (Nmax, Nmax, b, Nrank, Nrank)
    else 
      call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np,   &
           area, Nface, zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG, &
           weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
      if (PrnProgress) call write_progress (.false., 5+iprog*(iNint-1), Nprog)
      call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np,   &
           area, Nface, zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG, &
           weightsG, miror, perfectcond, DS, chiral, kb, b, Nrank, Nrank)
      if (PrnProgress) call write_progress (.false., 6+iprog*(iNint-1), Nprog)
      call LU_SYSTEM (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nmax)
      if (PrnProgress) call write_progress (.false., 7+iprog*(iNint-1), Nprog)
    end if
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,    &
         Nmax, c)
    call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, c, c1)
    call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)
    call DSCS (cc, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm, &
        .false.,.true., h, v)    
    call delta_DSCS (Nteta, h, v, oldh, oldv, epsNint, NthetaConv)
    call CQscat (cc, Mrank, Nrank, Nmaxmax, k, snorm, Cscat, Qscat)
    call CQext (cc, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,         &
         alfap, k, snorm, Cext, Qext)
    call write_3ConvParam (Nint, Nrank, m)
    call write_DSCS(Nteta,.false., h, v)
    call write_Effic (Qscat, Qext)  
    Nint = Nint + dNint
    deallocate (paramG, weightsG, Nintparam)
  end do  
  call write_NintConvRes (NthetaConv, Nteta, epsNint)                 
  deallocate (a, b, c, c1, cc, h, v, oldh, oldv)
end subroutine convergence_NintAXSYM 
!***********************************************************************************
subroutine convergence_NrankAXSYM (FileGeom, TypeGeom, k, ind_ref, snorm, Nsurf,    &
           surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint, miror,         &
           perfectcond, DS, chiral, kb, epsNrank, PrnProgress)                  
  use parameters
  implicit none  
  integer    :: TypeGeom, Nsurf, Nface, Nparam, Nrank, Nint
  real(O)    :: k, surf(Nsurf), snorm, epsNrank, kb, zRe(Nrank), zIm(Nrank),        &
                rp(2,NfacePD), np(2,NfacePD), area(NfacePD)
  complex(O) :: ind_ref
  logical    :: FileGeom, miror, perfectcond, DS, chiral, PrnProgress
!      
  integer    :: Nteta, Mstart, Mrank, Nmax, Nmaxmax, i, m, NthetaConv, iprog,       &
                Nprog, j   
  real(O)    :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,        &
                Cext, Qext
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG(:,:), weightsG(:,:), h(:), v(:), oldh(:), oldv(:)
  complex(O),allocatable :: a(:,:), b(:,:), c(:), c1(:), cc(:), ap(:,:), bp(:,:),   &
                            am(:,:), bm(:,:)
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
  Nmax   = Nrank    
  m = 1 
  call write_TypeConvHead (2)
  allocate (a(2*Nrank,2*Nrank), b(2*Nrank,2*Nrank), c(2*Nrank), c1(2*Nrank))     
  allocate (ap(2*Nrank,2*Nrank), bp(2*Nrank,2*Nrank))
  allocate (am(2*Nrank,2*Nrank), bm(2*Nrank,2*Nrank))  
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  allocate (cc(2*Nmaxmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))  
  allocate (paramG(Nparam,Nint), weightsG(Nparam,Nint), Nintparam(Nparam))
  if (.not. FileGeom) then
    call interpolation_listAXSYM (TypeGeom, Nsurf, surf, Nint, Nparam,              &
         Nintparam, paramG, weightsG, miror)
  else
    do i = 1, Nparam
      do j = 1, Nint
        paramG(i,j)   = 0._O 
        weightsG(i,j) = 0._O
      end do
      Nintparam(i) = 1
    end do
  end if
  if (.not. chiral) then
    Nprog = 5
  else
    Nprog = 8
  end if
  if (PrnProgress) call write_progress (.true., 1, Nprog)             
! --- Nrank configuration ---                    
  call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np, area, &
       Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG, weightsG,  &
       miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
  if (PrnProgress) call write_progress (.false., 2, Nprog)
  call copy_matrix (2*Nmax, 2*Nmax, a, 2*Nrank, 2*Nrank, ap, 2*Nrank, 2*Nrank)
  call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np, area, &
       Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG, weightsG,  &
       miror, perfectcond, DS, chiral, kb, b, Nrank, Nrank)
  if (PrnProgress) call write_progress (.false., 3, Nprog)
  call copy_matrix (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, bp, 2*Nrank, 2*Nrank)
  call LU_SYSTEM (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nmax)
  if (PrnProgress) call write_progress (.false., 4, Nprog)
  iprog = 4
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,       &
       Nmax, c)
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, c, c1)
  call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)
  if (.not. chiral) then
    call matrix_m_negativ (Nmax, Nmax, b, Nrank, Nrank)
  else 
    call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,   &
         weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
    if (PrnProgress) call write_progress (.false., 4, Nprog)
    call copy_matrix (2*Nmax, 2*Nmax, a, 2*Nrank, 2*Nrank, am, 2*Nrank, 2*Nrank)        
    call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,   &
         weightsG, miror, perfectcond, DS, chiral, kb, b, Nrank, Nrank)
    if (PrnProgress) call write_progress (.false., 5, Nprog)
    call copy_matrix (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, bm, 2*Nrank, 2*Nrank)    
    call LU_SYSTEM (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nmax)
    if (PrnProgress) call write_progress (.false., 6, Nprog)
    iprog = 6
  end if
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,      &
       Nmax, c)
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, c, c1)
  call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)
  call DSCS (cc, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm,   &
      .false.,.true., h, v) 
  call CQscat (cc, Mrank, Nrank, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,           &
       alfap, k, snorm, Cext, Qext)
  do i = 1, Nteta
    oldh(i) = h(i)
    oldv(i) = v(i)
  end do  
  call write_3ConvParam (Nint, Nrank, m)
  call write_DSCS (Nteta,.false., h, v)
  call write_Effic (Qscat, Qext)
! --- (Nrank - 1) configuration ---
  call copy_matrix (2*Nmax, 2*Nmax, ap, 2*Nrank, 2*Nrank, a, 2*Nrank, 2*Nrank)
  call matrix_Nrank_m_left (Nmax, a, Nrank, Nrank)
  call copy_matrix (2*Nmax, 2*Nmax, bp, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank)
  call matrix_Nrank_m (Nmax, b, Nrank, Nrank)
  call LU_SYSTEM (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nmax)
  if (PrnProgress) call write_progress (.false., iprog+1, Nprog)
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,       &
       Nmax, c)
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, c, c1)
  call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)
  if (.not. chiral) then
    call matrix_m_negativ (Nmax, Nmax, b, Nrank, Nrank)
  else 
    call copy_matrix (2*Nmax, 2*Nmax, am, 2*Nrank, 2*Nrank, a, 2*Nrank, 2*Nrank)
    call matrix_Nrank_m_left (Nmax, a, Nrank, Nrank)
    call copy_matrix (2*Nmax, 2*Nmax, bm, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank)        
    call matrix_Nrank_m (Nmax, b, Nrank, Nrank)
    call LU_SYSTEM (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nmax)
    if (PrnProgress) call write_progress (.false., iprog+2, Nprog)
  end if
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,      &
       Nmax, c)
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, c, c1)
  call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)
  call DSCS (cc, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm,   &
      .false.,.true., h, v)
  call CQscat (cc, Mrank, Nrank, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,           &
       alfap, k, snorm, Cext, Qext)
  call delta_DSCS (Nteta, h, v, oldh, oldv, epsNrank, NthetaConv)  
  call write_3ConvParam (Nint, Nrank - 1, m)
  call write_DSCS (Nteta,.false., h, v)
  call write_Effic (Qscat, Qext)
  call write_NrankConvRes (NthetaConv, Nteta, epsNrank) 
  deallocate (ap, bp, am, bm)
  deallocate (a, b, c, c1, cc, h, v, oldh, oldv, paramG, weightsG, Nintparam)
end subroutine convergence_NrankAXSYM
!***********************************************************************************
subroutine convergence_MrankAXSYM (FileGeom, TypeGeom, k, ind_ref, snorm, Nsurf,    &
           surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint, miror,         &
           perfectcond, DS, chiral, kb, epsMrank, FileTmat, PrnProgress)                
  use parameters
  implicit none   
  integer       :: TypeGeom, Nsurf, Nface, Nparam, Nrank, Nint
  real(O)       :: k, surf(Nsurf), snorm, epsMrank, kb, zRe(Nrank), zIm(Nrank),     &
                   rp(2,NfacePD), np(2,NfacePD), area(NfacePD)
  complex(O)    :: ind_ref
  logical       :: FileGeom, miror, perfectcond, DS, chiral, PrnProgress
  character(80) :: FileTmat
!      
  integer       :: Nteta, Mstart, Mrank, Nmax, Nmaxmax, i, m, NthetaConv, Nprog, j
  real(O)       :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,     &
                   Cext, Qext  
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG(:,:), weightsG(:,:), h(:), v(:), oldh(:), oldv(:)
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
  open (unit = iTmat, file = FileTmat, status = 'replace')  
  call write_HeadFileTmat (Nrank, Nrank) 
  call write_TypeConvHead (3)
  call write_2ConvParamAxsym (Nint, Nrank)
  allocate (a(2*Nrank,2*Nrank), b(2*Nrank,2*Nrank), c(2*Nrank), c1(2*Nrank))
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  allocate (cc(2*Nmaxmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do  
  allocate (paramG(Nparam,Nint), weightsG(Nparam,Nint), Nintparam(Nparam))
  if (.not. FileGeom) then
    call interpolation_listAXSYM (TypeGeom, Nsurf, surf, Nint, Nparam,              &
         Nintparam, paramG, weightsG, miror)
  else
    do i = 1, Nparam
      do j = 1, Nint
        paramG(i,j)   = 0._O 
        weightsG(i,j) = 0._O
      end do
      Nintparam(i) = 1
    end do
  end if        
  if (.not. chiral) then
    Nprog = 4
  else
    Nprog = 7
  end if
  Mrank = - 1
  do m = Mstart, Nrank    
    call write_1ConvParam (m)
    Mrank = Mrank + 1
    if (m == 0) then
      Nmax = Nrank
    else
      Nmax = Nrank - m + 1
    end if  
    if (PrnProgress) call write_progress_m (.true., m, 1, Nprog)    
    call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,    &
         weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank) 
    if (PrnProgress) call write_progress_m (.false., m, 2, Nprog)       
    call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,    &
         weightsG, miror, perfectcond, DS, chiral, kb, b, Nrank, Nrank) 
    if (PrnProgress) call write_progress_m (.false., m, 3, Nprog)                           
    call LU_SYSTEM (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nmax)
    if (PrnProgress) call write_progress_m (.false., m, 4, Nprog)                   
    call write_FileTmat (Nrank, Nrank, b)
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,     &
         Nmax, c)
    call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, c, c1)
    call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)       
    if (m /= 0) then
      if (.not. chiral) then
        call matrix_m_negativ (Nmax, Nmax, b, Nrank, Nrank)     
      else 
        call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np, &
             area, Nface, zRe, zIm, -m, Nrank,Nmax, Nint, Nparam, Nintparam,        &
             paramG, weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
        if (PrnProgress) call write_progress_m (.false., m, 5, Nprog)                     
        call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np, &
             area, Nface, zRe, zIm, -m, Nrank,Nmax, Nint, Nparam, Nintparam,        &
             paramG, weightsG, miror, perfectcond, DS, chiral, kb, b, Nrank, Nrank)    
        if (PrnProgress) call write_progress_m (.false., m, 6, Nprog)
        call LU_SYSTEM (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nmax)       
        if (PrnProgress) call write_progress_m (.false., m, 7, Nprog)
        call write_FileTmat (Nrank, Nrank, b)   
      end if
      call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,  &
           Nmax, c) 
      call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, c, c1)
      call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)     
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
  call write_InfoFileTmat (FileTmat, Mrank, Nrank, .true., .false., chiral)
  call ScatCharact (k, FileTmat, Mrank, Nrank, .true., .false., chiral)
  print "(/,2x,'T matrix is stored in ',a50)", FileTmat
  print "(  2x,'The dimensions of the T matrix are given by:')"      
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", Nrank
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank                     
  deallocate (a, b, c, c1, cc, h, v, oldh, oldv, paramG, weightsG, Nintparam)
end subroutine convergence_MrankAXSYM
!***********************************************************************************
subroutine convergence_NintDSAXSYM (FileGeom, TypeGeom, k, ind_ref, snorm, Nsurf,   &
           surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint, dNint, miror,  &
           perfectcond, DS, chiral, kb, epsNint, PrnProgress)
  use parameters                                                                           
  implicit none 
  integer    :: TypeGeom, Nsurf, Nface, Nparam, Nrank, Nint, dNint
  real(O)    :: k, surf(Nsurf), snorm, epsNint, kb, zRe(Nrank), zIm(Nrank),         &
                rp(2,NfacePD), np(2,NfacePD), area(NfacePD)
  complex(O) :: ind_ref
  logical    :: FileGeom, miror, perfectcond, DS, chiral, PrnProgress
!      
  integer    :: Nteta, Mstart, Mrank, Nmax, Nmaxmax, i, m, iNint, NthetaConv,       &
                iprog, Nprog 
  real(O)    :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,        &
                Cext, Qext
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG(:,:), weightsG(:,:), h(:), v(:), oldh(:), oldv(:)
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
  Nmax   = Nrank        
  m = 1           
  call write_TypeConvHead (1)
  allocate (a(2*Nrank,2*Nrank), b(2*Nrank,2*Nrank), c(2*Nrank), c1(2*Nrank))
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  allocate (cc(2*Nmaxmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))  
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do
  if (.not. chiral) then
    Nprog = 9
    iprog = 4
  else
    Nprog = 17
    iprog =  8
  end if
  if (PrnProgress) call write_progress (.true., 1, Nprog)
  do iNint = 1, 2         
    allocate(paramG(Nparam,Nint), weightsG(Nparam,Nint), Nintparam(Nparam))
    call interpolation_listAXSYM (TypeGeom, Nsurf, surf, Nint, Nparam,              &
         Nintparam, paramG, weightsG, miror)       
    call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,    &
         weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
    if (PrnProgress) call write_progress (.false., 2+iprog*(iNint-1), Nprog)
    call incident_matrix_m (FileGeom, TypeGeom, k, Nsurf, surf, rp, np, area, Nface,& 
         zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG, weightsG, b,    &
         Nrank, Nrank)
    if (PrnProgress) call write_progress (.false., 3+iprog*(iNint-1), Nprog)
    call LU_SYSTEM_DIRECT (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nrank, 2*Nmax)
    if (PrnProgress) call write_progress (.false., 4+iprog*(iNint-1), Nprog)
    call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,    &
         weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
    call product_matrices (2*Nmax, 2*Nrank, 2*Nmax, a, 2*Nrank, 2*Nrank, b, 2*Nrank,&
         2*Nrank)    
    if (PrnProgress) call write_progress (.false., 5+iprog*(iNint-1), Nprog)
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,     &
         Nmax, c)
    call product_matrix_vector (2*Nmax, 2*Nmax, a, 2*Nrank, 2*Nrank, c, c1)
    call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)
    if (.not. chiral) then
      call matrix_m_negativ (Nmax, Nmax, a, Nrank, Nrank)
    else 
      call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np,   &
           area, Nface, zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG, &
           weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
      if (PrnProgress) call write_progress (.false., 6+iprog*(iNint-1), Nprog)
      call incident_matrix_m (FileGeom, TypeGeom, k, Nsurf, surf, rp, np, area,     &
           Nface, zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,       &
           weightsG, b, Nrank, Nrank)
      if (PrnProgress) call write_progress (.false., 7+iprog*(iNint-1), Nprog)
      call LU_SYSTEM_DIRECT (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nrank,     &
           2*Nmax)
      if (PrnProgress) call write_progress (.false., 8+iprog*(iNint-1), Nprog)
      call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np,   &
           area, Nface, zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG, &
           weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)    
      call product_matrices (2*Nmax, 2*Nrank, 2*Nmax, a, 2*Nrank, 2*Nrank,          &
           b, 2*Nrank, 2*Nrank)
      if (PrnProgress) call write_progress (.false., 9+iprog*(iNint-1), Nprog)
    end if
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,    &
         Nmax, c)
    call product_matrix_vector (2*Nmax, 2*Nmax, a, 2*Nrank, 2*Nrank, c, c1)
    call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)
    call DSCS (cc, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm, &
        .false.,.true., h, v)
    call CQscat (cc, Mrank, Nrank, Nmaxmax, k, snorm, Cscat, Qscat)
    call CQext (cc, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,         &
         alfap, k, snorm, Cext, Qext)    
    call delta_DSCS (Nteta, h, v, oldh, oldv, epsNint, NthetaConv)
    call write_3ConvParam (Nint, Nrank, m)
    call write_DSCS (Nteta,.false., h, v)
    call write_Effic (Qscat, Qext)
    Nint = Nint + dNint
    deallocate (paramG, weightsG, Nintparam)
  end do
  call write_NintConvRes (NthetaConv, Nteta, epsNint)                 
  deallocate (a, b, c, c1, cc, h, v, oldh, oldv)
end subroutine convergence_NintDSAXSYM
!***********************************************************************************
subroutine convergence_NrankDSAXSYM (FileGeom, TypeGeom, k, ind_ref, snorm, Nsurf,  &
           surf, rp, np, area, Nface, zRe, zIm, zRe1, zIm1, Nparam, Nrank, Nrank1,  &
           Nint, miror, perfectcond, DS, chiral, kb, epsNrank, PrnProgress)              
  use parameters
  implicit none  
  integer    :: TypeGeom, Nsurf, Nface, Nparam, Nrank, Nrank1, Nint
  real(O)    :: k, surf(Nsurf), snorm, epsNrank, kb, zRe(Nrank), zIm(Nrank),        &
                zRe1(Nrank1), zIm1(Nrank1), rp(2,NfacePD), np(2,NfacePD),           &
                area(NfacePD)
  complex(O) :: ind_ref
  logical    :: FileGeom, miror, perfectcond, DS, chiral, PrnProgress
!      
  integer    :: Nteta, Mstart, Mrank, Nmax, Nmaxmax, i, m, NthetaConv,              &
                iprog, Nprog, j
  real(O)    :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,        &
                Cext, Qext
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG(:,:), weightsG(:,:), h(:), v(:), oldh(:), oldv(:)
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
  Nmax   = Nrank  
  m = 1
  call write_TypeConvHead (2)
  allocate (a(2*Nrank,2*Nrank), b(2*Nrank,2*Nrank), c(2*Nrank), c1(2*Nrank))
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  allocate (cc(2*Nmaxmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))  
  allocate (paramG(Nparam,Nint), weightsG(Nparam,Nint), Nintparam(Nparam))
  if (.not. FileGeom) then
    call interpolation_listAXSYM (TypeGeom, Nsurf, surf, Nint, Nparam,              &
         Nintparam, paramG, weightsG, miror)
  else
    do i = 1, Nparam
      do j = 1, Nint
        paramG(i,j)   = 0._O 
        weightsG(i,j) = 0._O
      end do
      Nintparam(i) = 1
    end do
  end if         
  if (.not. chiral) then
    Nprog = 9
  else
    Nprog = 17
  end if
  if (PrnProgress) call write_progress (.true., 1, Nprog)
! --- Nrank configuration ---                    
  call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np, area, &
       Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG, weightsG,  &
       miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
  if (PrnProgress) call write_progress (.false., 2, Nprog)
  call incident_matrix_m (FileGeom, TypeGeom, k, Nsurf, surf, rp, np, area, Nface,  &
       zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG, weightsG, b,      &
       Nrank, Nrank)                              
  if (PrnProgress) call write_progress (.false., 3, Nprog)
  call LU_SYSTEM_DIRECT (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nrank, 2*Nmax)
  if (PrnProgress) call write_progress (.false., 4, Nprog)
  call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np, area, &
       Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG, weightsG,  &
       miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
  call product_matrices (2*Nmax, 2*Nrank, 2*Nmax, a, 2*Nrank, 2*Nrank, b, 2*Nrank,  &
       2*Nrank)    
  if (PrnProgress) call write_progress (.false., 5, Nprog)
  iprog = 5
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank, Nmax, c)
  call product_matrix_vector (2*Nmax, 2*Nmax, a, 2*Nrank, 2*Nrank, c, c1)
  call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)
  if (.not. chiral) then
    call matrix_m_negativ (Nmax, Nmax, a, Nrank, Nrank)
  else 
    call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,   &
         weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
    if (PrnProgress) call write_progress (.false., 6, Nprog)
    call incident_matrix_m (FileGeom, TypeGeom, k, Nsurf, surf, rp, np, area, Nface,&
         zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG, weightsG, b,   &
         Nrank, Nrank)
    if (PrnProgress) call write_progress (.false., 7, Nprog)
    call LU_SYSTEM_DIRECT (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nrank, 2*Nmax)
    if (PrnProgress) call write_progress (.false., 8, Nprog)
    call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,   &
         weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
    call product_matrices (2*Nmax, 2*Nrank, 2*Nmax, a, 2*Nrank, 2*Nrank, b, 2*Nrank,&
         2*Nrank) 
    if (PrnProgress) call write_progress (.false., 9, Nprog)
        iprog = 9
  end if
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,      &
       Nmax, c)
  call product_matrix_vector (2*Nmax, 2*Nmax, a, 2*Nrank, 2*Nrank, c, c1)
  call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)  
  call DSCS (cc, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm,   &
      .false.,.true., h, v)
  call CQscat (cc, Mrank, Nrank, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,           &
       alfap, k, snorm, Cext, Qext)
  do i = 1, Nteta
    oldh(i) = h(i)
    oldv(i) = v(i)
  end do  
  call write_3ConvParam (Nint, Nrank, m)
  call write_DSCS (Nteta,.false., h, v)
  call write_Effic (Qscat, Qext)
! --- Nrank1 configuration ---
  Nmax   = Nrank1
  Nmaxmax= Nrank1 + Mrank*(2*Nrank1 - Mrank + 1)  
  call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np, area, &
       Nface, zRe1, zIm1, m, Nrank1, Nmax, Nint, Nparam, Nintparam, paramG,         &
       weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
  if (PrnProgress) call write_progress (.false., iprog+1, Nprog)
  call incident_matrix_m (FileGeom, TypeGeom, k, Nsurf, surf, rp, np, area, Nface,  &
       zRe1, zIm1, m, Nrank1, Nmax, Nint, Nparam, Nintparam, paramG, weightsG, b,   &
       Nrank, Nrank)
  if (PrnProgress) call write_progress (.false., iprog+2, Nprog)
  call LU_SYSTEM_DIRECT (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nrank1, 2*Nmax)
  if (PrnProgress) call write_progress (.false., iprog+3, Nprog)
  call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np, area, &
       Nface, zRe1, zIm1, m, Nrank1, Nmax, Nint, Nparam, Nintparam, paramG,         &
       weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
  call product_matrices (2*Nmax, 2*Nrank1, 2*Nmax, a, 2*Nrank, 2*Nrank, b, 2*Nrank, &
       2*Nrank)    
  if (PrnProgress) call write_progress (.false., iprog+4, Nprog)
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank1,      &
       Nmax, c)
  call product_matrix_vector (2*Nmax, 2*Nmax, a, 2*Nrank, 2*Nrank, c, c1)
  call extend_vector_positive (c1, cc, m, Mstart, Nrank1, Nmax, Nmaxmax)  
  if (.not. chiral) then 
    call matrix_m_negativ (Nmax, Nmax, a, Nrank, Nrank)
  else 
    call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe1, zIm1, -m, Nrank1, Nmax, Nint, Nparam, Nintparam,        &
         paramG, weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)    
    if (PrnProgress) call write_progress (.false., iprog+5, Nprog)
    call incident_matrix_m (FileGeom, TypeGeom, k, Nsurf, surf, rp, np, area, Nface,& 
         zRe1, zIm1, -m, Nrank1, Nmax, Nint, Nparam, Nintparam, paramG, weightsG,   &
         b, Nrank, Nrank)    
    if (PrnProgress) call write_progress (.false., iprog+6, Nprog)
    call LU_SYSTEM_DIRECT (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nrank1,      &
         2*Nmax)    
    if (PrnProgress) call write_progress (.false., iprog+7, Nprog)
    call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe1, zIm1, -m, Nrank1, Nmax, Nint, Nparam, Nintparam,        &
         paramG, weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)    
    call product_matrices (2*Nmax, 2*Nrank1, 2*Nmax, a, 2*Nrank, 2*Nrank,           &
         b, 2*Nrank, 2*Nrank)       
    if (PrnProgress) call write_progress (.false., iprog+8, Nprog)
  end if
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank1,     &
       Nmax, c)
  call product_matrix_vector (2*Nmax, 2*Nmax, a, 2*Nrank, 2*Nrank, c, c1)
  call extend_vector_negative (c1, cc, m, Nrank1, Nmax, Nmaxmax)  
  call DSCS (cc, Mrank, Nrank1, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm,  &
      .false.,.true., h, v)
  call CQscat (cc, Mrank, Nrank1, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, Nrank1, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,          &
       alfap, k, snorm, Cext, Qext)
  call delta_DSCS (Nteta, h, v, oldh, oldv, epsNrank, NthetaConv)  
  call write_3ConvParam (Nint, Nrank1, m)
  call write_DSCS (Nteta,.false., h, v)
  call write_Effic (Qscat, Qext)
  call write_NrankConvRes (NthetaConv, Nteta, epsNrank)
  deallocate (a, b, c, c1, cc, h, v, oldh, oldv, paramG, weightsG, Nintparam)
end subroutine convergence_NrankDSAXSYM 
!***********************************************************************************
subroutine convergence_MrankDSAXSYM (FileGeom, TypeGeom, k, ind_ref, snorm, Nsurf,  &
           surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint, miror,         &
           perfectcond, DS, chiral, kb, epsMrank, dNintMrank, FileTmat, PrnProgress)            
  use parameters
  implicit none   
  integer       :: TypeGeom, Nsurf, Nface, Nparam, Nrank, Nint, dNintMrank
  real(O)       :: k, surf(Nsurf), snorm, epsMrank, kb, zRe(Nrank), zIm(Nrank),     &
                   rp(2,NfacePD), np(2,NfacePD), area(NfacePD)
  complex(O)    :: ind_ref
  logical       :: FileGeom, miror, perfectcond, DS, chiral, PrnProgress
  character(80) :: FileTmat
!      
  integer       :: Nteta, Mstart, Mrank, Nmax, Nmaxmax, i, m, NthetaConv, Nprog, j
  real(O)       :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,     &
                   Cext, Qext  
  logical       :: ComplexPlane
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG(:,:), weightsG(:,:), h(:), v(:), oldh(:), oldv(:)
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
  open (unit = iTmat, file = FileTmat, status = 'replace') 
  call write_HeadFileTmat (Nrank, Nrank) 
  call write_TypeConvHead (3)
  call write_2ConvParamAxsym (Nint, Nrank)
  allocate (a(2*Nrank,2*Nrank), b(2*Nrank,2*Nrank), c(2*Nrank), c1(2*Nrank))
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  allocate (cc(2*Nmaxmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do  
  allocate (paramG(Nparam,Nint), weightsG(Nparam,Nint), Nintparam(Nparam))
  if (.not. FileGeom) then
    call interpolation_listAXSYM (TypeGeom, Nsurf, surf, Nint, Nparam,              &
         Nintparam, paramG, weightsG, miror)
  else
    do i = 1, Nparam
      do j = 1, Nint
        paramG(i,j)   = 0._O 
        weightsG(i,j) = 0._O
      end do
      Nintparam(i) = 1
    end do
  end if       
  ComplexPlane = .false.        
  do i = 1, Nrank                  
    if (zIm(i) /= 0._O) ComplexPlane = .true.
  end do
  if (.not. chiral) then
    Nprog = 5
  else
    Nprog = 9
  end if                        
  Mrank = - 1
  do m = Mstart, Nrank          
    call write_1ConvParam (m)
    Mrank = Mrank + 1
    if (m == 0) then
      Nmax = Nrank
    else
      Nmax = Nrank - m + 1
    end if  
    if (ComplexPlane) then
      deallocate (paramG, weightsG, Nintparam)
      Nint = Nint + dNintMrank
      allocate (paramG(Nparam,Nint), weightsG(Nparam,Nint), Nintparam(Nparam))
      call interpolation_listAXSYM (TypeGeom, Nsurf, surf, Nint, Nparam,            &
           Nintparam, paramG, weightsG, miror)          
    end if 
    if (PrnProgress) call write_progress_m (.true., m, 1, Nprog)  
    call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,    &
         weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
    if (PrnProgress) call write_progress_m (.false., m, 2, Nprog) 
    call incident_matrix_m (FileGeom, TypeGeom, k, Nsurf, surf, rp, np, area, Nface,&
         zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG, weightsG, b,    &
         Nrank, Nrank)
    if (PrnProgress) call write_progress_m (.false., m, 3, Nprog)
    call LU_SYSTEM_DIRECT (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nrank, 2*Nmax)
    if (PrnProgress) call write_progress_m (.false., m, 4, Nprog)
    call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,    &
         weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
    call product_matrices (2*Nmax, 2*Nrank, 2*Nmax, a, 2*Nrank, 2*Nrank, b, 2*Nrank,&
         2*Nrank)   
    if (PrnProgress) call write_progress_m (.false., m, 5, Nprog)
    call write_FileTmat (Nrank, Nrank, a)
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,     &
         Nmax, c)           
    call product_matrix_vector (2*Nmax, 2*Nmax, a, 2*Nrank, 2*Nrank, c, c1)
    call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)       
    if (m /= 0) then
      if (.not. chiral) then
        call matrix_m_negativ (Nmax, Nmax, a, Nrank, Nrank)
      else 
        call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np, &
             area, Nface, zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam,       &
             paramG, weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
        if (PrnProgress) call write_progress_m (.false., m, 6, Nprog)
        call incident_matrix_m (FileGeom, TypeGeom, k, Nsurf, surf, rp, np, area,   &
             Nface, zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,     &
             weightsG, b, Nrank, Nrank)
        if (PrnProgress) call write_progress_m (.false., m, 7, Nprog)
        call LU_SYSTEM_DIRECT (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nrank,   &
             2*Nmax)
        if (PrnProgress) call write_progress_m (.false., m, 8, Nprog)
        call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np, &
             area, Nface, zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam,       &
             paramG, weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
        call product_matrices (2*Nmax, 2*Nrank, 2*Nmax, a, 2*Nrank, 2*Nrank,        &
             b, 2*Nrank, 2*Nrank)               
        if (PrnProgress) call write_progress_m (.false., m, 9, Nprog)
        call write_FileTmat (Nrank, Nrank, a)   
      end if
      call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,  &
           Nmax, c)
      call product_matrix_vector (2*Nmax, 2*Nmax, a, 2*Nrank, 2*Nrank, c, c1)
      call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)
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
  call write_InfoFileTmat (FileTmat, Mrank, Nrank, .true., .false., chiral)     
  call ScatCharact (k, FileTmat, Mrank, Nrank, .true., .false., chiral) 
  print "(/,2x,'T matrix is stored in ',a50)", FileTmat
  print "(  2x,'The dimensions of the T matrix are given by:')"
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", Nrank 
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank             
  deallocate (a, b, c, c1, cc, h, v, oldh, oldv, paramG, weightsG, Nintparam)
end subroutine convergence_MrankDSAXSYM
