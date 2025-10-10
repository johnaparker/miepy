subroutine EFMED
!-----------------------------------------------------------------------------------
! 1. General Considerations                                                        !
! --------------------------                                                       !
! EFMED is a code for computing the effective wave number of a medium with         !
! randomly distributed spheroidal particles. The particles are identical and are   !
! illuminated by a vector plane wave traveling along the Z-axis of the global      !
! coordinate system (normal incidence). The T-matrix method together with the      !
! quasi-crystalline approximation is used to generate a set of homogeneous         !
! equations (the generalized Lorentz-Lorenz law):                                  ! 
!                                                                                  !
!                         [I - n * T * (J1 + J2)] * s = 0,                         !
!                                                                                  !
! where I is the identity matrix, T is the transition matrix, n is the number      !
! of particles per unit volume (n = N / V), J1 and J2 are the averaged translation ! 
! matrices and s is the unknown vector of the scattered field coefficients.        !
! Note that the elements of the averaged translation matrices are expressed as     !
! integrals over the zenith angle beta. For a nontrivial solution of the system of ! 
! homogeneous equations, the determinant must vanishes giving an equation for      ! 
! the effective wave number. The code computes the effective wave numbers for      !
! size parameters ranging from xmin to xmax in dx steps. For a spheroidal          ! 
! particles with semi-axes a and b, the size parameter is defined as               !
!                                                                                  !
!                             x = k * max(a,b),                                    !
!                                                                                  !
! where k is the wave number in the ambient medium. The root of the equation       !
!                                                                                  !
!                   det {I - n * T * (J1 + J2)} = 0,                               !
!                                                                                  !
! is searched in the complex plane using Mueller's method. Mueller's method        ! 
! requires an initial guess and two starting values for the effective wave         !
! number. For x = xmin, the initial guess is computed by using the first and       !
! second order approximations of the long-wavelength T-matrix method. The          !
! second and third starting values of the wave number are given by:                !
!                                                                                  !
!                           K2 = K1 + delta,                                       !
!                                                                                  ! 
! and                                                                              !
!                                                                                  ! 
!                           K3 = K2 + delta,                                       !
!                                                                                  !
! respectively, where K1 is the initial guess and delta is a parameter specified   !
! in the namelist input. The process of calculation is iterative, i.e., the wave   !
! number corresponding to the size parameter x serves as initial guess for         !
! computing the wave number corresponding to the size parameter x + dx.            !
!                                                                                  !
! For axisymmetric particles with the axis of symmetry directed along the Z-axis,  !
! the transition matrix is diagonal with respect to azimuthal index. Because the   !
! non-vanishing incident field coefficients correspond to the azimuthal modes      !
! m = 1 and m = -1, only these modes are required for computing the effective      !
! wave number.                                                                     !   
!                                                                                  !
! The elements of the averaged translation matrix J2 are expressed in terms of     !
! the quantities:                                                                  !
!                                                                                  !
!                   |inf                                                           !
!             G_n = |  [g(2Rx) - 1] * h_n(2kRx) * j_n(2Krx) x**2 dx,               !                  
!                   |0                                                             !
!                                                                                  !
! where g is the pair distribution function, R is the radius of a sphere           !
! enclosing the particle (for spherical particles, R is the particle radius),      !
! k is the wave number in the ambient medium and h_n and j_n are the spherical     !
! Hankel and Bessel functions. For the Percus-Yevick pair distribution function,   !
! we have:                                                                         !
!                                                                                  !
!                             Pi  |inf          sin(ux)                            !
!              g(2Rx) - 1 = ------|  H(u/2R) * --------- u**2 du                   !  
!                           2R**2 |0              ux                               !
!                                                                                  !
! where H is the Fourier transform of the total correlation function. Note that    !
! Percus-Yevick pair distribution function is the solution of the Ornstein-Zernike ! 
! equation for the total and direct correlation functions of two particles h(r)    !
! and d(r) under the assumption that h(r) - d(r) = y(r) - 1. For the case of       !
! hard-sphere potential, the Ornstein-Zernicke equation with Percus-Yevick         !
! approximation has a closed-form solution H. The function g(2Rx) - 1 is computed  !
! by sintransform                                                                  !
!                                                                                  !
!                              Pi  |L                                              !
!              g(2Rx) - 1 = -------|  H(u/2R) * u * sin(ux) du,                    !
!                           2xR**2 |0                                              !
!                                                                                  !
! while the G_n are given by                                                       !
!                                                                                  !  
!                   |Xmax                                                          !
!             G_n = |  [g(2Rx) - 1] * h_n(2kRx) * j_n(2Krx) x**2 dx.               !
!                   |0                                                             !
!                                                                                  !
! The number of discrete points N for computing the pair distribution function     ! 
! by sintransform is                                                               !
!                                                                                  !
!                         N = 2**pN,                                               !
!                                                                                  !
! while the upper limits L and Xmax are                                            !
!                                                                                  !
!                         L ~ 2**(pN - p)                                          !
!                                                                                  !
! and                                                                              !
!                                                                                  !
!                         Xmax ~ 2**p,                                             !
!                                                                                  !
! respectively. The integers p and pN are supplied in the namelist input, and the  !
! coefficients G_n are computed by using the trapez integration method.            !
!                                                                                  !
! For spherical particles with small size parameters, the effective wave numbers   !
! can be also computed by using the long-wavelength quasi-crystalline              !
! approximations with and without coherent potential.                              !
!                                                                                  !   
! The input data are supplied by the file "/INPUTFILES/InputEFMED.dat", while      !
! the output data are written to the file "/OUTPUTFILES/Output.dat".               !
!                                                                                  !
! 2. Input Parameters                                                              !
! --------------------                                                             !
! The parameters specified in the input file "/INPUTFILES/InputEFMED.dat" are      !  
! listed below.                                                                    !
!                                                                                  !
! - ind_refMed (real) - refractive index of the ambient medium.                    !
!                                                                                  !
! - ind_refRel (complex) - relative refractive index of the spheroidal particles   !
!   with respect to the ambient medium.                                            !
!                                                                                  !
! - xmax (real) - maximum size parameter. The size parameter is defined as         !
!   x = k * max(a,b), where k is the wave number in the ambient medium, a is       !
!   the semi-axis along the axis of symmetry and b is the second half-axis.        !
!                                                                                  !
! - xmin (real) - minimum size parameter.                                          !
!                                                                                  !
! - dx (real) - size parameter step. This parameter MUST BE POSITIVE.              !
!                                                                                  !
! - e (real) - eccentricity of the spheroidal particle. e is defined as            !
!   e = a / b, where as before, a is the semi-axis along the axis of symmetry      !
!   and b is the second half-axis. For spherical particles, set e = 1.0.           !
!                                                                                  !
! - c (real) - volume concentration, c = N * Vp / V, where N is the number of      !
!   particles, Vp is the volume of the particle and V is the volume occupied by    !      
!   the particles.                                                                 !
!                                                                                  !
! - Nint (integer) - number of integration points on the particle surface for      !
!   computing the T matrix.                                                        !
!                                                                                  !
! - RandomOrientation (logical) - if RandomOrientation = t, the spheroidal         !
!   particle is randomly oriented, otherwise the axis of symmetry is along the     ! 
!   Z-axis (the incident direction).                                               !
!                                                                                  !
! - ComputeQCA (logical) - if ComputeQCA = t, the long-wavelength                  !  
!   quasicrystalline approximations with and without coherent potential are used   !
!   to compute the effective wave number. These methods can be used for spherical  !
!   particles with small size parameter.                                           !
!                                                                                  !
! - Nbeta (integer) - number of integration points for computing the elements of   !
!   the averaged translation matrices J1 and J2.                                   !
!                                                                                  !
! - p, pN (integer variables) - determine the numbers of discretization points     !
!   for computing the pair distribution function by sintransform,                  !
!                                                                                  !
!                                |L                                                !
!                         f(x) = |  F(u)*sin(u*x) du,                              !
!                                |0                                                !
!                                                                                  !
!   The number of discrete points is N = 2**pN, L ~ Pi * 2**(pN - p), and the      !
!   upper limit for the pair-distribution function integration is Xmax ~ 2**p.     !
!   For low volume concentration we recommend to use p = 3, while for high         !
!   volume concentration we recommend to use p = 4. Here is a table giving some    !
!   approximate values of the parameters.                                          !
!                                                                                  !
!                         N         L          Xmax                                !
!                                                                                  !
!     p = 3, pN =  8     256       100           8                                 !
!     p = 3, pN =  9     512       200           8                                 ! 
!     p = 3, pN = 10    1024       400           8                                 !
!     p = 3, pN = 11    2048       800           8                                 ! 
!     p = 3, pN = 12    4096      1600           8                                 !
!     p = 3, pN = 13    8192      3200           8                                 ! 
!                                                                                  !
!     p = 4, pN =  8     256        50          16                                 !
!     p = 4, pN =  9     512       100          16                                 ! 
!     p = 4, pN = 10    1024       200          16                                 !
!     p = 4, pN = 11    2048       400          16                                 ! 
!     p = 4, pN = 12    4096       800          16                                 !
!     p = 4, pN = 13    8192      1600          16                                 !
!                                                                                  !
!     p = 5, pN =  8     256        25          32                                 !
!     p = 5, pN =  9     512        50          32                                 ! 
!     p = 5, pN = 10    1024       100          32                                 !
!     p = 5, pN = 11    2048       200          32                                 ! 
!     p = 5, pN = 12    4096       400          32                                 !
!     p = 5, pN = 13    8192       800          32                                 !  
!                                                                                  !
!  - epsX, epsY (real variables) - tolerance specifying the accuracy of root     !
!    computation.  For Mueller's method we impose a conservative criterion, i.e.,  !
!    the convergence is achieved if the iterates converge (epsX cratered) and     !
!    the residual is smaller than a given tolerance (epsY criterion).             !
!                                                                                  !
! - Niter (integer) - maximum number of iterations for Mueller's method.           !
!                                                                                  !
! WARNING: OUR FEELING IS THAT THERE IS AN ERROR IN THE ROUTINE FOR COMPUTING THE  !
! EFFECTIVE WAVE NUMBER OF SPHEROIDAL PARTICLES. IN THIS REGARD, WE RECOMMEND TO   !
! USE THE CODE ONLY FOR SPHERICAL PARTICLES.                                       ! 
!-----------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  implicit none 
  integer    :: Nint, Nbeta, Niter, p, pN, ios, N
  real(O)    :: xmax, xmin, dx, e, c, wavelength, epsX, epsY, ind_refMed,           &
                lambdaMed, wavenumberMed, du, L, lmax 
  complex(O) :: ind_refRel
  logical    :: RandomOrientation, ComputeQCA
! -----------------------------------------------------------------------------------
!                               Read the Input File                                 ! 
! ----------------------------------------------------------------------------------- 
  call readinputEFMED ( wavelength, ind_refMed, ind_refRel, xmax, xmin, dx, e,      &
       c, Nint, Nbeta, p, pN, epsX, epsY, Niter, RandomOrientation, ComputeQCA )       
  lambdaMed     = wavelength / ind_refMed
  wavenumberMed = 2 * Pi / lambdaMed        
! ----------------------------------------------------------------------------------
!                         Write Informations to Output File                        ! 
! ----------------------------------------------------------------------------------               
  open (unit = iOutput, file = FileOutput, status = "replace") 
  write (iOutput,"(/,2x,'Input Parameters of the Scattering Problem:',/)")
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'refractive index of the ambient medium, ind_refMed = ', ind_refMed, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'wave number in the ambient medium, wave_numberMed  = ', wavenumberMed, ';'
  write (iOutput,"(2x, a, 1pe10.3, ',', 1pe10.3, a)")                               &
 'relative refractive index of the particle, ind_refRel = (', ind_refRel, ');'  
  write (iOutput,"(2x,'volume concentration,   c    = ',1pe10.3,';')") c
  write (iOutput,"(2x,'maximum size parameter, xmax = ',1pe10.3,';')") xmax
  write (iOutput,"(2x,'minimum size parameter, xmin = ',1pe10.3,';')") xmin
  write (iOutput,"(2x,'size parameter step,    dx   = ',1pe10.3,';')") dx
  write (iOutput,"(2x,'spheroid eccentricity,  e    = ',1pe10.3,';')") e                       
  write (iOutput,"(2x, a, i3,';')")                                                 &
 'number of integration points on the particle surface, Nint = ', Nint
  if (RandomOrientation) then
    write (iOutput,"(2x,'randomly oriented spheroidal particles;')")
  else
    write (iOutput,"(2x,'particles in fixed orientation;')")
  end if
  N    = 2**pN
  du   = Pi / 2**p
  L    = real(N - 1,O) * du  
  lmax = real(N - 1,O) * Pi / real(N,O) / du
  write (iOutput,"(2x, a, i2, a, i2, a)")                                           &
 'parameters controlling the discretization process, p = ', p, ', pN = ', pN, ';'
  write (iOutput,"(2x,'number of discrete points, N = ',i6,';')") N
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'upper limit for sintransform integration, L = ', L, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'upper limit for pair-distribution function integration, Xmax = ', lmax, ';'
  write (iOutput,"(2x, a, i6, a)")                                                  &
 'maximum number of iterations for Mueller''s method, Niter = ', Niter, ';'
  write (iOutput,"(2x, a, 1pe10.3, a, 1pe10.3, a)")                                 &
 'tolerences controlling Mueller''s iterations, epsX = ', epsX, ', epsY = ', epsY, ';'
!-----------------------------------------------------------------------------------
!                                   Main                                           !
!----------------------------------------------------------------------------------- 
  print "(/,2x,'Effective Wavenumber of a Medium with Spheroidal Particles',/)"
  call IterativeSpheroidMullerSolver ( wavenumberMED, ind_refRel, xmax, xmin, dx,   &
       e, c, Nint, Nbeta, p, pN, Niter, epsX, epsY, RandomOrientation, ComputeQCA )
  close (unit = iOutput)  
end subroutine EFMED
!***********************************************************************************
subroutine readinputEFMED ( wavelength, ind_refMed, ind_refRel, xmax, xmin, dx, e,  &
           c, Nint, Nbeta, p, pN, epsX, epsY, Niter, RandomOrientation, ComputeQCA ) 
  use parameters
  use derived_parameters
  implicit none 
  integer       :: Nint, Nbeta, Niter, p, pN, ios
  real(O)       :: wavelength, xmax, xmin, dx, e, c, epsX, epsY, ind_refMed
  complex(O)    :: ind_refRel
  character(80) :: string
  logical       :: RandomOrientation, ComputeQCA, XFindPar    
!
  call DrvParameters  
  open (unit = iInputEFMED, file = FileInputEFMED, status = "old", position = "rewind")
  wavelength = 0.2_O * Pi        
  ind_refMed = 1._O                                                         
  ind_refRel = (1.668_O,0.510_O)
  string     = 'OptProp'
  if (XFindPar (iInputEFMED, string)) then   
    read (iInputEFMED, *, iostat = ios) wavelength
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable wavelength;')"
      stop
    end if   
    read (iInputEFMED, *, iostat = ios) ind_refMed
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refMed;')"
      stop
    end if
    read (iInputEFMED, *, iostat = ios) ind_refRel
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refRel;')"
      stop
    end if      
  else
    print "(/,2x,'Group name OptProp not found;')"
    stop  
  end if
!
  xmax = 0.5_O
  xmin = 0.1_O  
  dx   = 0.1_O  
  e    = 1.0_O
  c    = 0.1_O 
  Nint = 60
  RandomOrientation = .true.
  ComputeQCA        = .true.      
  string = 'TmatPart'
  if (XFindPar (iInputEFMED, string)) then    
    read (iInputEFMED, *, iostat = ios) xmax
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable xmax;')"
      stop
    end if
    read (iInputEFMED, *, iostat = ios) xmin    
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable xmin;')"
      stop
    end if    
    read (iInputEFMED, *, iostat = ios) dx    
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable dx;')"
      stop
    end if 
    read (iInputEFMED, *, iostat = ios) e    
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable e;')"
      stop
    end if                         
    read (iInputEFMED, *, iostat = ios) c
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable c;')"
      stop
    end if 
    read (iInputEFMED, *, iostat = ios) Nint
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nint;')"
      stop
    end if 
    read (iInputEFMED, *, iostat = ios) RandomOrientation    
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable RandomOrientation;')"
      stop
    end if 
    read (iInputEFMED, *, iostat = ios) ComputeQCA
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ComputeQCA;')"
      stop
    end if                    
  else
    print "(/,2x,'Group name TmatPart not found;')"
    stop  
  end if
  if (e == 1._O) then
    RandomOrientation = .false.
  else
    ComputeQCA = .false. 
  end if
!
  Nbeta = 180
  p  = 5
  pN = 10            
  epsX   = 1.e-6_O
  epsY   = 1.e-6_O 
  Niter  = 1000
  string = 'AlgProp'
  if (XFindPar (iInputEFMED, string)) then    
    read (iInputEFMED, *, iostat = ios) Nbeta
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nbeta;')"
      stop
    end if
    read (iInputEFMED, *, iostat = ios) p
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable p;')"
      stop
    end if 
    read (iInputEFMED, *, iostat = ios) pN
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable pN;')"
      stop
    end if    
    read (iInputEFMED, *, iostat = ios) epsX
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsX;')"
      stop
    end if
    read (iInputEFMED, *, iostat = ios) epsY
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsY;')"
      stop
    end if
    read (iInputEFMED, *, iostat = ios) Niter
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Niter;')"
      stop
    end if
  else
    print "(/,2x,'Group name AlgProp not found;')"
    stop  
  end if
  close (unit = iInputEFMED)
end subroutine readinputEFMED
!***********************************************************************************
subroutine IterativeSpheroidMullerSolver ( wavenumberMED, ind_refRel, xmax, xmin,   &
           dx, e, c, Nint, Nbeta, p, pN, Niter, epsX, epsY, RandomOrientation,      &
           ComputeQCA )
  use parameters
  implicit none 
  integer    :: Nint, Nbeta, p, pN, Niter
  real(O)    :: xmax, xmin, dx, e, c, epsX, epsY, wavenumberMed
  complex(O) :: ind_refRel            
  logical    :: RandomOrientation, ComputeQCA           
!
  integer    :: Nx, ix, Nrank
  real(O)    :: Nxreal, x, SemiAxisA, SemiAxisB, a, x3
  complex(O) :: kX, wavenumberEFMED, epsRel, tamp, fct
  logical    :: sphere, Conv
!
  sphere = .false.
  if ( e == 1._O ) sphere = .true.
  if ( xmin > xmax ) then
    Nx = 1
  else
    Nxreal = ( xmax - xmin ) / dx
    Nx     = int(Nxreal)
    if (Nx + 0.1_O > Nxreal) then
      Nx = Nx + 1
    else
      Nx = Nx + 2
    end if
  end if
  kX = cmplx(wavenumberMED,0.0,O)
  print "(2x,'Iterative calculation over the size parameter:')"
! ----------------------------------------------------------------------------------
!                              T-Matrix Method                                     !
! ----------------------------------------------------------------------------------
  write (iOutput,"(/)")
  write (iOutput,"(2x,'Results',/)")
  write (iOutput,"(2x,'T-Matrix Method')")
  write (iOutput,"(2x,'size parameter          effective wave number',/)") 
  do ix = 1, Nx
    if (ix < Nx) then
      x = xmin + (ix - 1) * dx
    else
      x = xmax
    end if
    if (e >= 1._O) then
      SemiAxisA = x / wavenumberMED
      SemiAxisB = SemiAxisA / e
    else
      SemiAxisB = x / wavenumberMED
      SemiAxisA = SemiAxisB * e
    end if
    Nrank = int(x + 4.05_O * x**0.33_O + 2._O)
    if (Nrank < 4) Nrank = Nrank + 4
    a = max (SemiAxisA, SemiAxisB) 
    print "(2x,'iteration      = ', i3, ';')", ix
    print "(2x,'size parameter =', 1pe13.4, ';')", x 
    if ( ix == 1) then
      call InitialGuessSpheroid (wavenumberMED, ind_refRel, a, SemiAxisA,           &
           SemiAxisB, sphere, RandomOrientation, c, Nrank, Nint, epsX, epsY,        &
           Niter, Conv, kX)
      if (Conv) then        
        print "(2x,'initial guess  = (',1pe13.4,',',1pe13.4,'   );')", kX
      else
        print "(/,2x,'convergence is not achieved at the initial guess;')" 
        print "(  2x, a)",                                                          &
       'the tolerances epsX/epsY, or the iteration number Niter are too low;'	             
        stop
      end if		         
    end if
    wavenumberEFMED = kX
    call MullerSpheroid (wavenumberEFMED, wavenumberMED, ind_refRel, a, SemiAxisA,  &
         SemiAxisB, sphere, RandomOrientation, c, Nrank, Nint, Nbeta, p, pN, Niter, &
         epsX, epsY, Conv, kX)	  				               		                                                        
    if (Conv) then	  	  
      print "(2x,'initial  effective wave number = (',1pe13.4,',',1pe13.4,'   );')",  &
	      wavenumberEFMED
      print "(2x,'computed effective wave number = (',1pe13.4,',',1pe13.4,'   );',/)",&
              kX
      write (iOutput,"(2x, 1pe13.4, '      (',1pe13.4,',',1pe13.4,'   )')") x, kX
     !write (iOutput,"(2x, 1pe13.4, 2x, 1pe13.4, 2x, 1pe13.4)") x, wavenumberMED/real(kX), &
     !	     2.d0*aimag(kX)/real(kX)
    else
      print "(/,2x,'convergence is not achieved;')" 
      print "(  2x, a)",                                                            &
     'the tolerances epsX/epsY, or the iteration number Niter are too low;'	             
      stop
    end if
  end do
  if (sphere .and. ComputeQCA) then
    write (iOutput,"(/,2x,'QCA and QCA-CP')")
    write (iOutput,"(2x, a,/)")                                                     &
   'size parameter         effective wave number QCA          effective wave number QCA-CP'
    epsRel = ind_refRel**2   
    tamp   = (epsRel - 1._O) / (epsRel + 2._O)
    do ix = 1, Nx
      if (ix < Nx) then
        x = xmin + (ix - 1) * dx
      else
        x = xmax
      end if	  	  
      a = x / wavenumberMED
!     ------------------------------------------------------------------------------
!                 Long-Wavelength Quasicrystalline Approximation                   !
!     ------------------------------------------------------------------------------   
      x3  = x**3
      fct = 3._O * c * tamp  / (1._O - c * tamp)
      fct = fct * (1._O + im * 2._O * x3 * tamp * (1._O - c)**4 /                   &
          ((1._O - c * tamp) * (1._O + 2._O * c)**2 * 3._O))
      wavenumberEFMED = wavenumberMed * sqrt(1._O + fct)            
!     ------------------------------------------------------------------------------
!         Long-Wavelength Quasicrystalline Approximation with Coherent Potential   !
!     ------------------------------------------------------------------------------       
      call MullerQCACP (wavenumberEFMED, wavenumberMED, ind_refRel, a, c, Niter,    &
           epsX, epsY, Conv, kX)      
      if (Conv) then	  	  	  
        write (iOutput,"(2x,1pe13.4,a,1pe13.4,',',1pe13.4,a,1pe13.4,',',1pe13.4,a)")&
               x,'      (', wavenumberEFMED, '   )     (', kX, '   )'   
      else
        print "(/,2x,'convergence is not achieved for QCA-CP;')" 
        print "(  2x, a)",                                                          &
       'the tolerances epsX/epsY, or the iteration number Niter are too low;'	             
        stop
      end if
    end do
  end if
end subroutine IterativeSpheroidMullerSolver  
!***********************************************************************************
subroutine MullerSpheroid (wavenumberEFMED, wavenumberMED, ind_refRel, a, SemiAxisA,&
           SemiAxisB, sphere, RandomOrientation, c, Nrank, Nint, Nbeta, p, pN,      &
           Niter, epsX, epsY, Conv, x)
  use parameters
  implicit none 
  integer    :: Nrank, Nint, Nbeta, p, pN, Niter
  real(O)    :: a, SemiAxisA, SemiAxisB, c, epsX, epsY, wavenumberMed
  complex(O) :: wavenumberEFMED, ind_refRel, x
  logical    :: sphere, RandomOrientation, Conv 
!  
  integer    :: iter
  real(O)    :: deltaX
  complex(O) :: delta, x0, x1, x2, P0, P1, P2, s, t, u, den, den1, d, q         
! 
  delta = 0.01_O * wavenumberEFMED 
  x0 = wavenumberEFMED
  x1 = x0 + delta
  x2 = x1 + delta
  call DeterminantEFMEDSpheroid (x0, wavenumberMED, ind_refRel, a, SemiAxisA,       &
       SemiAxisB, sphere, RandomOrientation, c, Nrank, Nint, Nbeta, p, pN, P0)
  call DeterminantEFMEDSpheroid (x1, wavenumberMED, ind_refRel, a, SemiAxisA,       &
       SemiAxisB, sphere, RandomOrientation, c, Nrank, Nint, Nbeta, p, pN, P1)
  call DeterminantEFMEDSpheroid (x2, wavenumberMED, ind_refRel, a, SemiAxisA,       &
       SemiAxisB, sphere, RandomOrientation, c, Nrank, Nint, Nbeta, p, pN, P2)        
  iter = 0
  Conv = .false.        
  print "(2x,'iteration summary of Mueller''s method:')"
  print "(2x,'determinants at starting points:')"
  print "(2x,'                  x0              x1              x2')"
  print "(2x,'determinant ', 1pe12.3, 4x, 1pe12.3, 4x, 1pe12.3)", abs(P0), abs(P1), abs(P2)
  print "(2x,'iteration',1x,'|x(n) - x(n-1)|',6x,'|det(x(n))|')"  
  do while ((iter < Niter) .and. (.not.Conv))
    iter = iter + 1
    q = (x2 - x1) / (x1 - x0)                   
    s = q * P2 - q * (1._O + q) * P1 + q**2 * P0
    t = (2._O * q + 1._O) * P2 - (1._O + q)**2 * P1 + q**2 * P0
    u = (1._O + q) * P2
    d = sqrt(t**2 - 4._O * s * u)
    den  = t + d
    den1 = t - d
    if (abs(den) < abs(den1)) den = den1 
    x = x2 - (x2 - x1) * 2._O * u / den    
    deltaX = abs(x - x2)
    if (abs(x - x1) < abs(x - x0)) then
      s  = x1
      x1 = x0
      x0 = s
      t  = P1
      P1 = P0
      P0 = t
    else if (abs(x - x1) > abs(x - x2)) then
      s  = x2  
      x2 = x1  
      x1 = s
      t  = P2  
      P2 = P1  
      P1 = t
    end if
    x2 = x    
    call DeterminantEFMEDSpheroid (x2, wavenumberMED, ind_refRel, a, SemiAxisA,     &
         SemiAxisB, sphere, RandomOrientation, c, Nrank, Nint, Nbeta, p, pN, P2)          
    print "(5x,i3,4x,1pe12.3,7x,1pe12.3)", iter, deltaX, abs(P2)    
    if ((deltaX < epsX) .and. (abs(P2) < epsY)) Conv = .true.           
  end do  
  x = x2
  if (imag(x) < 0._O) x = wavenumberEFMED  
end subroutine MullerSpheroid 
!***********************************************************************************
subroutine InitialGuessSpheroid (k, ind_refRel, a, SemiAxisA, SemiAxisB, sphere,    &
           RandomOrientation, c, Nrank, Nint, epsX, epsY, Niter, Conv, kX) 
  use parameters
  implicit none
  integer    :: Nrank, Nint, Niter
  real(O)    :: a, SemiAxisA, SemiAxisB, c, k, epsX, epsY
  complex(O) :: kX, ind_refRel
  logical    :: sphere, RandomOrientation, Conv
!  
  integer    :: iter
  real(O)    :: deltaX  
  complex(O) :: delta, T11, T22, T12, T21, x, alpha, x0, x1, x2, P0, P1, P2, s, v,  &
                u, den, den1, d, q 
  complex(O),allocatable :: T(:,:), tv(:), tv11(:), tv12(:), tv21(:), tv22(:)
! --- first order approximation ---        
  if (.not. sphere) then
    if (.not. RandomOrientation) then
      allocate (T(2*Nrank,2*Nrank)) 
      call TmatrixSpheroid (k, ind_refRel, SemiAxisA, SemiAxisB, Nrank, Nint, T)
      T11 = T(1,1)
      T22 = T(Nrank+1,Nrank+1)  
      T12 = T(1,Nrank+1)
      T21 = T(Nrank+1,1)
      deallocate (T)
    else
      allocate (tv11(Nrank), tv12(Nrank), tv21(Nrank), tv22(Nrank))      
      call AverageTmatrixSpheroid (k, ind_refRel, SemiAxisA, SemiAxisB, Nrank,      &
           Nint, tv11, tv12, tv21, tv22)
      T11 = tv11(1)
      T22 = tv22(1)
      T12 = tv12(1)
      T21 = tv21(1)
      deallocate (tv11, tv12, tv21, tv22)
    end if
  else
    allocate (tv(2*Nrank)) 
    call coefficients_fg_m (k, SemiAxisA, ind_refRel, 1, Nrank, Nrank, tv)
    T11 = tv(1)
    T22 = tv(Nrank+1)
    T12 = 0._O
    T21 = 0._O
    deallocate (tv)
  end if
  x     = cmplx(k * a,0.0,O)
  alpha = x * x * x * ( (1._O - c)**4 / (1._O + 2._O * c)**2 - 1._O / T22 )
  kX    = k * sqrt ( (2._O * alpha + 6._O * im * c) / (2._O * alpha - 3._O * im * c) ) 
  if ( imag(kX) < 0._O ) kX = - kX
! --- second order approximation --- 
  delta = 0.01_O * kX
  x0 = kX
  x1 = x0 + delta
  x2 = x1 + delta
  call DeterminantLWLSpheroid (x0, k, a, c, T11, T22, T12, T21, P0)
  call DeterminantLWLSpheroid (x1, k, a, c, T11, T22, T12, T21, P1)
  call DeterminantLWLSpheroid (x2, k, a, c, T11, T22, T12, T21, P2)
  iter = 0
  Conv = .false.    
  do while ((iter < Niter) .and. (.not. Conv))
    iter = iter + 1
    q = (x2 - x1) / (x1 - x0)                   
    s = q * P2 - q * (1._O + q) * P1 + q**2 * P0
    v = (2._O * q + 1._O) * P2 - (1._O + q)**2 * P1 + q**2 * P0
    u = (1._O + q) * P2
    d = sqrt(v**2 - 4._O * s * u)
    den  = v + d
    den1 = v - d
    if (abs(den) < abs(den1)) den = den1 
    x = x2 - (x2 - x1) * 2._O * u / den    
    deltaX = abs(x - x2)
    if (abs(x - x1) < abs(x - x0)) then
      s  = x1
      x1 = x0
      x0 = s
      v  = P1
      P1 = P0
      P0 = v
    else if (abs(x - x1) > abs(x - x2)) then
      s  = x2  
      x2 = x1  
      x1 = s
      v  = P2  
      P2 = P1  
      P1 = v
    end if
    x2 = x    
    call DeterminantLWLSpheroid (x2, k, a, c, T11, T22, T12, T21, P2)    
    if ((deltaX < epsX) .and. (abs(P2) < epsY)) Conv = .true.           
  end do  
  if (imag(x2) >= 0._O) kX  = x2
end subroutine InitialGuessSpheroid
!***********************************************************************************
subroutine DeterminantEFMEDSpheroid (kX, k, ind_refRel, a, SemiAxisA, SemiAxisB,    &
           sphere, RandomOrientation, c, Nrank, Nint, Nbeta, p, pN, Pv) 
  use parameters
  implicit none
  integer    :: Nrank, Nint, Nbeta, p, pN
  real(O)    :: a, SemiAxisA, SemiAxisB,c, k
  complex(O) :: kX, ind_refRel, Pv
  logical    :: sphere, RandomOrientation
!  
  integer    :: i, j
  complex(O) :: fact, det
  complex(O),allocatable :: AB(:,:), T(:,:), tv(:), tv11(:), tv12(:), tv21(:),      &
                            tv22(:)
! 
  allocate (AB(2*Nrank,2*Nrank), T(2*Nrank,2*Nrank))
  call matrix_translationAB_mn_m1n1_EFMED (k, kX, a, c, Nrank, Nbeta, p, pN,        &
       AB, Nrank, Nrank)
  if ( .not. sphere) then
    if (.not. RandomOrientation) then   	           
      call TmatrixSpheroid (k, ind_refRel, SemiAxisA, SemiAxisB, Nrank, Nint, T)  
      call product_matrices (2*Nrank, 2*Nrank, 2*Nrank, T, 2*Nrank, 2*Nrank,        &
           AB, 2*Nrank, 2*Nrank)  
    else 		 
      allocate (tv11(Nrank), tv12(Nrank), tv21(Nrank), tv22(Nrank))
      call AverageTmatrixSpheroid (k, ind_refRel, SemiAxisA, SemiAxisB, Nrank,      &
           Nint, tv11, tv12, tv21, tv22)
      call product_AverageTmatrix (Nrank, tv11, tv12, tv21, tv22, Nrank,            &
           AB, 2*Nrank, 2*Nrank, T, 2*Nrank, 2*Nrank)	
      deallocate (tv11, tv12, tv21, tv22)		     		 
    end if		   
  else
    allocate (tv(2*Nrank))
    call coefficients_fg_m (k, SemiAxisA, ind_refRel, 1, Nrank, Nrank, tv)                  
    call product_matrix_vector1 (2*Nrank, 2*Nrank, tv, AB, 2*Nrank, 2*Nrank,        &
         T, 2*Nrank, 2*Nrank)
    deallocate (tv)		    
  end if
  fact = 24._O * c / ((k * a)**2 - (kX * a)**2)
  do i = 1, 2*Nrank
    do j = 1, 2*Nrank
      if (i == j) then
        T(i,j) = one - fact * T(i,j)
      else
        T(i,j) = - fact * T(i,j)
      end if
    end do
  end do
  call determinant (T, 2*Nrank, 2*Nrank, 2*Nrank, det)
  Pv = det   
  deallocate (AB, T)
end subroutine DeterminantEFMEDSpheroid
!***********************************************************************************
subroutine DeterminantLWLSpheroid (kX, k, r, c, T11, T22, T12, T21, P)
  use parameters
  implicit none
  real(O)    :: r, c, k
  complex(O) :: kX, T11, T22, T12, T21, P
!
  complex(O) :: x, xx, x2, ci, JH0, JH1, JH2, A, B, a11, a12, a21, a22
!    
  x   = cmplx(k * r,0.0,O) 
  xx  = kX * r
  x2  = x**2 - xx**2
  ci  = 3._O * im * c
  JH0 = ci / x / x2 + (1._O - c)**4 / (1._O + 2._O * c)**2
  JH1 = ci * xx / x / x / x2
  JH2 = ci * xx**2 / x / x / x / x2
  A   = JH0 + 0.5_O * JH2
  B   = 1.5_O * JH1
  a11 = one - ( T11 * A + T12 * B )
  a12 = - ( T11 * B + T12 * A )
  a21 = - ( T21 * A + T22 * B )
  a22 = one - ( T21 * B + T22 * A )
  P   = a11 * a22 - a12 * a21
end subroutine DeterminantLWLSpheroid
!***********************************************************************************
subroutine TmatrixSpheroid (k, ind_ref, SemiAxisA, SemiAxisB, Nrank, Nint, T)                
  use parameters
  implicit none   
  integer       :: Nrank, Nint
  real(O)       :: k, SemiAxisA, SemiAxisB
  complex(O)    :: ind_ref, T(2*Nrank,2*Nrank)  
!      
  integer       :: Nmax, TypeGeom, Nsurf, Nparam, Nface, i, j, m
  real(O)       :: rp(2,NfacePD), np(2,NfacePD), area(NfacePD), kb
  logical       :: FileGeom, miror, perfectcond, DS, chiral
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: surf(:), zRe(:), zIm(:), paramG(:,:), pondereG(:,:)
  complex(O),allocatable :: a(:,:)
!       
  m      = 1
  Nmax   = Nrank
  DS     = .false.          
  chiral = .false.       
  kb     =  0.1  
  miror  = .false. 
  perfectcond = .false.          
  FileGeom = .false.
  TypeGeom = 1                                                
  Nsurf    = 2  
  allocate (surf(Nsurf))               
  surf(1) = SemiAxisA              
  surf(2) = SemiAxisB               
  Nparam  = 1                               
  Nface   = 1
  do i = 1, NfacePD
    do j = 1, 2
      rp(j,i) = 0._O
      np(j,i) = 0._O
    end do
    area(i) = 0._O
  end do 			    
  allocate (zRe(Nrank), zIm(Nrank)) 
  do i = 1, Nrank
    zRe(i)  = 0._O
    zIm(i)  = 0._O   
  end do              
  allocate (a(2*Nrank,2*Nrank))  
  allocate (paramG(Nparam,Nint), pondereG(Nparam,Nint), Nintparam(Nparam))
  call interpolation_listAXSYM (TypeGeom, Nsurf, surf, Nint, Nparam,                &
       Nintparam, paramG, pondereG, miror)	   	   	         
  call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np,       &
       area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,      &
       pondereG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)      
  call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np,       &
       area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,      &
       pondereG, miror, perfectcond, DS, chiral, kb, T, Nrank, Nrank)                           
  call LU_SYSTEM (a, 2*Nrank, 2*Nrank, T, 2*Nrank, 2*Nrank, 2*Nmax)                                   
  deallocate (a, surf, zRe, zIm, paramG, pondereG, Nintparam)
end subroutine TMatrixSpheroid
!***********************************************************************************
subroutine AverageTmatrixSpheroid (k, ind_ref, SemiAxisA, SemiAxisB, Nrank, Nint,   &
           tv11, tv12, tv21, tv22)                
  use parameters
  implicit none   
  integer       :: Nrank, Nint
  real(O)       :: k, SemiAxisA, SemiAxisB
  complex(O)    :: ind_ref, tv11(Nrank), tv12(Nrank), tv21(Nrank), tv22(Nrank)   
!      
  integer       :: Nmax, TypeGeom, Nsurf, Nparam, Nface, i, j, m, n, kn
  real(O)       :: rp(2,NfacePD), np(2,NfacePD), area(NfacePD), kb, fact
  logical       :: FileGeom, miror, perfectcond, DS, chiral
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: surf(:), zRe(:), zIm(:), paramG(:,:), pondereG(:,:)
  complex(O),allocatable :: a(:,:), T(:,:)
!       
  DS     = .false.          
  chiral = .false.       
  kb     =  0.1  
  miror  = .false. 
  perfectcond = .false.          
  FileGeom = .false.
  TypeGeom = 1                                                
  Nsurf    = 2  
  allocate (surf(Nsurf))               
  surf(1) = SemiAxisA              
  surf(2) = SemiAxisB               
  Nparam  = 1                               
  Nface   = 1
  do i = 1, NfacePD
    do j = 1, 2
      rp(j,i) = 0._O
      np(j,i) = 0._O
    end do
    area(i) = 0._O
  end do 			    
  allocate (zRe(Nrank), zIm(Nrank)) 
  do i = 1, Nrank
    zRe(i)  = 0._O
    zIm(i)  = 0._O   
  end do              
  allocate (a(2*Nrank,2*Nrank), T(2*Nrank,2*Nrank))  
  allocate (paramG(Nparam,Nint), pondereG(Nparam,Nint), Nintparam(Nparam))
  call interpolation_listAXSYM (TypeGeom, Nsurf, surf, Nint, Nparam,                &
       Nintparam, paramG, pondereG, miror)	   	   	                                           
  do n = 1, Nrank
    tv11(n) = zero
    tv12(n) = zero
    tv21(n) = zero
    tv22(n) = zero
  end do
  do m = 0, Nrank        
    if (m == 0) then
      Nmax = Nrank
    else
      Nmax = Nrank - m + 1
    end if  
    call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,    &
         pondereG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank) 
    call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,    &
         pondereG, miror, perfectcond, DS, chiral, kb, T, Nrank, Nrank)                       
    call LU_SYSTEM (a, 2*Nrank, 2*Nrank, T, 2*Nrank, 2*Nrank, 2*Nmax)
    do kn = 1, Nmax
      if (m == 0) then
        n    = kn
        fact = 1._O
      else
        n    = m + kn - 1
        fact = 2._O
      end if
      tv11(n) = tv11(n) + fact * T(kn,kn)
      tv12(n) = tv12(n) + fact * T(kn,kn+Nmax)
      tv21(n) = tv21(n) + fact * T(kn+Nmax,kn)
      tv22(n) = tv22(n) + fact * T(kn+Nmax,kn+Nmax)
    end do  
  end do
  do n = 1, Nrank
    fact    = 1._O / (2._O * n + 1._O)
    tv11(n) = fact * tv11(n)
    tv12(n) = fact * tv12(n)
    tv21(n) = fact * tv21(n)
    tv22(n) = fact * tv22(n)			
  end do  
  deallocate (a, T, surf, zRe, zIm, paramG, pondereG, Nintparam)
end subroutine AverageTMatrixSpheroid
!***********************************************************************************
subroutine matrix_translationAB_mn_m1n1_EFMED (wavenumber, wavenumberX, a, c, Nrank,&
           Nbeta, pow, powN, AB, NP, MP)
  use parameters
  implicit none
  integer    :: Nrank, Nbeta, pow, powN, NP, MP
  real(O)    :: a, c, wavenumber 
  complex(O) :: wavenumberX, AB(2*NP,2*MP)
!  
  integer    :: Nsup, n, n1, n2, pint
  real(O)    :: abeta, bbeta, beta, f, nn1
  complex(O) :: ka, kaX, sum, Fn2, fact, nm, hder, jder, fpp, ftt, fpt, ftp 
  real(O),allocatable    :: bt(:), wt(:), Pn(:), dPn(:), pin(:), taun(:),           &
                            Pnm(:), dPnm(:), pinm(:), taunm(:)
  complex(O),allocatable :: j(:), jd(:), h(:), hd(:), Fn2C(:)
!
  Nsup = 4 * Nrank ! we can use a lower limit for summation, i.e., Nsup = 2*Nrank;
  do n1 = 1, 2*Nrank
    do n = 1, 2*Nrank
      AB(n1,n) = zero
    end do  
  end do 
  allocate (j(0:Nsup+1), jd(0:Nsup+1), h(0:Nsup+1), hd(0:Nsup+1))
  allocate (bt(Nbeta), wt(Nbeta))
  allocate (Pn(0:Nsup), dPn(0:Nsup), pin(0:Nsup), taun(0:Nsup))
  allocate (Pnm(0:Nrank), dPnm(0:Nrank), pinm(0:Nrank), taunm(0:Nrank)) 
  ka  = cmplx(2._O * wavenumber * a,0.0,O)
  kaX = 2._O * wavenumberX * a     
  call besel_h (ka, Nsup + 1, h, hd)
  call besel_j (kaX, Nsup + 1, j, jd)  
  allocate (Fn2C(0:Nsup))
  call Fn2Cor (wavenumber, wavenumberX, a, c, pow, powN, Nsup, Fn2C)
  abeta = 0._O
  bbeta = Pi
  call Gauss_Legendre (abeta, bbeta, Nbeta, wt, bt)   
  do pint = 1, Nbeta
    beta = bt(pint)
    call Leg_normalized (beta, 0, Nsup, Pn, dPn, pin, taun)
    sum = zero
    do n2 = 0, Nsup
      hder = n2 * h(n2) / ka - h(n2+1)
      jder = n2 * j(n2) / kaX - j(n2+1)
      Fn2  = 0.5_O * (ka * hder * j(n2) - kaX * jder * h(n2))         
      Fn2  = Fn2 + Fn2C(n2)	  
      f = (- 1._O)**n2 * sqrt(real(2 * n2 + 1,O) / 2._O)
      sum = sum + f * Fn2 * Pn(n2)
    end do      
    fact = sum * wt(pint) * sin(beta)
    call Leg_normalized (beta, 1, Nrank, Pnm, dPnm, pinm, taunm)          
    do n = 1, Nrank      
      do n1 = 1, Nrank        
        nn1 = real(n * (n + 1) * n1 * (n1 + 1),O)
        nm  = im**(n1 - n) * fact / sqrt(nn1)		
        fpp = pinm(n)  * pinm(n1)  * nm
        ftt = taunm(n) * taunm(n1) * nm
        fpt = pinm(n)  * taunm(n1) * nm
        ftp = taunm(n) * pinm(n1)  * nm 
        AB(n,n1)       = AB(n,n1) + (fpp + ftt)        
		AB(n,n1+Nrank) = AB(n,n1+Nrank) - (fpt + ftp)	
      end do  
    end do 	 	                                           
  end do
  do n = 1, Nrank
    do n1 = 1, Nrank
      AB(n+Nrank,n1)       = AB(n,n1+Nrank)
      AB(n+Nrank,n1+Nrank) = AB(n,n1)
     end do  
  end do                           
  deallocate (j, jd, h, hd, bt, wt, Pn, dPn, pin, taun, Pnm, dPnm, pinm,            &
              taunm, Fn2C)
end subroutine matrix_translationAB_mn_m1n1_EFMED
!***********************************************************************************
subroutine MullerQCACP (wavenumberEFMED, wavenumberMED, ind_refRel, a, c, Niter,    &
           epsX, epsY, Conv, kX)
  use parameters
  implicit none 
  integer    :: Niter
  real(O)    :: a, c, epsX, epsY, wavenumberMed
  complex(O) :: wavenumberEFMED, ind_refRel, kX 
  logical    :: Conv                       
!  
  integer    :: iter
  real(O)    :: deltaX
  complex(O) :: delta, x0, x1, x2, x, P0, P1, P2, s, t, u, den, den1, d, q         
!  
  delta = 0.01_O * wavenumberEFMED
  x0 = wavenumberEFMED
  x1 = x0 + delta
  x2 = x1 + delta  
  call FuncQCACP (x0, wavenumberMED, ind_refRel, a, c, P0) 
  call FuncQCACP (x1, wavenumberMED, ind_refRel, a, c, P1)
  call FuncQCACP (x2, wavenumberMED, ind_refRel, a, c, P2)                       
  iter = 0
  Conv = .false.    
  do while ((iter < Niter) .and. (.not.Conv))
    iter = iter + 1
    q = (x2 - x1) / (x1 - x0)                   
    s = q * P2 - q * (1._O + q) * P1 + q**2 * P0
    t = (2._O * q + 1._O) * P2 - (1._O + q)**2 * P1 + q**2 * P0
    u = (1._O + q) * P2
    d = sqrt(t**2 - 4._O * s * u)
    den  = t + d
    den1 = t - d
    if (abs(den) < abs(den1)) den = den1 
    x = x2 - (x2 - x1) * 2._O * u / den    
    deltaX = abs(x - x2)
    if (abs(x - x1) < abs(x - x0)) then
      s  = x1
      x1 = x0
      x0 = s
      t  = P1
      P1 = P0
      P0 = t
    else if (abs(x - x1) > abs(x - x2)) then
      s  = x2  
      x2 = x1  
      x1 = s
      t  = P2  
      P2 = P1  
      P1 = t
    end if
    x2 = x        
    call FuncQCACP (x2, wavenumberMED, ind_refRel, a, c, P2)    
    if ((deltaX < epsX) .and. (abs(P2) < epsY)) Conv = .true.           
  end do
  if (imag(x2) >= 0._O) then    
    kX = x2 
  else
    kX = wavenumberEFMED
  end if
end subroutine MullerQCACP 
!***********************************************************************************
subroutine FuncQCACP (kX, k, ind_refRel, a, c, PQCA) 
  use parameters
  implicit none
  real(O)    :: a, c, k
  complex(O) :: kX, ind_refRel, PQCA
!
  complex(O) :: ks, k2, ka3, dk, fct
!
  ks  = k * ind_refRel
  k2  = ks * ks - k * k
  ka3 = k2 * kX * a**3
  dk  = k2 / (3._O * kX * kX)
  fct = c * k2 / (1._O + dk * (1._O - c))
  fct = fct * (1._O + im * 2._O * ka3 * (1._O - c)**4 /                             &
       (9._O * (1._O + dk * (1._O - c)) * (1._O + 2._O * c)**2))
  PQCA = fct + k**2 - kX**2
end subroutine FuncQCACP
!***********************************************************************************
subroutine Fn2Cor (wavenumber, wavenumberX, a, c, p, pN, Nsup, Fn2C)
  use parameters
  implicit none
  integer    :: p, pN, Nsup
  real(O)    :: a, c, wavenumber
  complex(O) :: wavenumberX, Fn2C(0:Nsup) 
!
  integer    :: N, kmin, k, n2
  real(O)    :: dx, xl, gPY, w
  complex(O) :: sum, ka, kaX
  real(O),allocatable    :: x(:), G(:)
  complex(O),allocatable :: j(:), jd(:), h(:), hd(:)
!
  N = 2**pN
  allocate (x(N), G(N))
  allocate (j(0:Nsup), jd(0:Nsup), h(0:Nsup), hd(0:Nsup))
  call PYVct (c, p, N, x, G, kmin)  
  dx = 1._O / 2**(pN - p)  
  do n2 = 0, Nsup
    sum = zero
    do k = kmin, N
      if (k == kmin .or. k == N) then
        w = 0.5_O * dx
      else
        w = dx
      end if
      xl  = x(k)                  
      gPY = G(k)
      ka  = cmplx(2._O * wavenumber * a * xl,0.0,O)
      kaX = 2._O * wavenumberX * a * xl     
      call besel_h (ka, Nsup, h, hd)
      call besel_j (kaX, Nsup, j, jd) 
      sum = sum + (gPY - 1._O) * h(n2) * j(n2) * xl * xl * w      
    end do            
    Fn2C(n2) = 2._O * ((wavenumber * a)**2 - (wavenumberX * a)**2) * sum
  end do
  deallocate (x, G, j, jd, h, hd)
end
!***********************************************************************************   
subroutine PYVct (c, p, N, x, G, kmin)
  use parameters
  implicit none
  integer :: p, N, kmin
  real(O) :: c, x(N), G(N)
!
  integer :: j, k
  real(O) :: Nr, N1r, L, du, uj, dx, f, max 
!
  Nr  = real(N,O)
  N1r = Nr - 1._O
  L   = Pi * N1r / 2**p  
  du  = L / N1r 
  do j = 1, N    
    uj   = real(j - 1,O) * du
    G(j) = f(c,uj) * du
  end do 
  call SINFT (G, N)
  dx = Pi / Nr / du
  do k = 1, N
    x(k) = real(k - 1,O) * dx
    if (x(k) >= 1._O) then
      G(k) = 1._O + G(k) / x(k)
    else
      G(k) = 0._O
    end if
  end do
! --- compute kmin for which the maximum of G is attained ---
  max  = - 1.e+10_O
  kmin =   0
  do k = 1, N
    if (G(k) > max) then
      max  = G(k)
      kmin = k
    end if
  end do                          
end subroutine PYVct   
!***********************************************************************************
function f(c,u) result(func)
  use parameters
  real(O),intent(in) :: c, u
  real(O) :: func
  real(O) :: a, b, d, g
!
  a =   (1._O + 2._O * c)**2 / (1._O - c)**4
  b = - 6._O * c * (1._O + 0.5_O * c)**2 / (1._O - c)**4
  d =   0.5_O * c * a  
  g =   (a + b + d) * u**4 * cos(u) - (a + 2 * b + 4 * d) * u**3 * sin(u) -         &
        2._O * (b + 6 * d) * u**2 * cos(u) + 2._O * b * u**2 +                      &
        24._O * d * u * sin(u) + 24._O * d * (cos(u) - 1._O)   
  g =   (u**6 / (u**6 - 24._O * c * g) - 1._O) / 24._O / c           
  func = 2._O * g * u / Pi
end function f
! **********************************************************************************
subroutine product_AverageTmatrix (Nrank, tv11, tv12, tv21, tv22, ntp, a, nap, map, &
           b, nbp, mbp)
!-----------------------------------------------------------------------------------
!     The routine computes the matrix product:                                     ! 
!         | t11   t12 |   | a11   a12 |   | t11*a11 + t12*a21   t11*a12 + t12*a22 |!
!     b = |           | * |           | = |                                       |!
!         | t21   t22 |   | a21   a22 |   | t21*a11 + t22*a21   t21*a12 + t22*a22 |!
!     where t11, t12, t21, t22 are diagonal matrices.                              !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: Nrank, ntp, nap, map, nbp, mbp, j, i
  complex(O) :: tv11(ntp), tv12(ntp), tv21(ntp), tv22(ntp), a(nap,map), b(nbp,mbp)
!
  do i = 1, Nrank
    do j = 1, Nrank
      b(i,j)             = tv11(i) * a(i,j)       + tv12(i) * a(i+Nrank,j)
	  b(i,j+Nrank)       = tv11(i) * a(i,j+Nrank) + tv12(i) * a(i+Nrank,j+Nrank)
      b(i+Nrank,j)       = tv21(i) * a(i,j)       + tv22(i) * a(i+Nrank,j)
	  b(i+Nrank,j+Nrank) = tv21(i) * a(i,j+Nrank) + tv22(i) * a(i+Nrank,j+Nrank)
    end do
  end do      
end subroutine product_AverageTmatrix
! ***********************************************************************************
! *             FAST SIN TRANSFORM ROUTINE FROM NUMERICAL RECIPES                   *
! ***********************************************************************************
subroutine SINFT (Y, N)
  use parameters
  implicit none
  integer :: N, M, j
  real(O) :: Y(N), wr, wi, wpr, wpi, wtemp, theta, Y1, Y2, sum
!
  theta = Pi / N
  wr  =   1._O
  wi  =   0._O
  wpr = - 2._O * sin(0.5_O * theta)**2
  wpi =   sin(theta)
  M   =   N / 2
  Y(1) = 0._O  
  do j = 1, M
    wtemp = wr
    wr = wr * wpr - wi * wpi + wr
    wi = wi * wpr + wtemp * wpi + wi
    Y1 = wi * (Y(j+1) + Y(N-j+1))
    Y2 = 0.5_O * (Y(j+1) - Y(N-j+1))
    Y(j+1)   = Y1 + Y2
    Y(N-j+1) = Y1 - Y2
  end do
  call REALFT (Y, M, 1)
  sum  = 0._O
  Y(1) = 0.5_O * Y(1)
  Y(2) = 0._O
  do j = 1, N - 1, 2
    sum    = sum + Y(j)
    Y(j)   = Y(j+1)
    Y(j+1) = sum
  end do
end subroutine SINFT
!*********************************************************************************
subroutine REALFT (Y, N, isign)
  use parameters
  implicit none
  integer :: N, isign, n2p3, i, i1, i2, i3, i4
  real(O) :: Y(2*N), wr, wi, wpr, wpi, wtemp, c1, c2, h1r, h1i, h2r, h2i, theta
!    
  theta = Pi / N
  c1    = 0.5_O
  if (isign == 1) then
    c2 = - 0.5_O
    call FOUR (Y, N, +1)
  else
    c2    =   0.5_O
    theta = - theta
  end if
  wpr  = - 2._O * sin(0.5_O * theta)**2
  wpi  =   sin(theta)
  wr   =   1._O + wpr
  wi   =   wpi
  n2p3 =   2 * N + 3
  do i = 2, N/2 + 1
    i1  = 2 * i - 1
    i2  = i1 + 1
    i3  = n2p3 - i2
    i4  = i3 + 1
    h1r =   c1 * (Y(i1) + Y(i3))
    h1i =   c1 * (Y(i2) - Y(i4))
    h2r = - c2 * (Y(i2) + Y(i4))
    h2i =   c2 * (Y(i1) - Y(i3))
    Y(i1) =   h1r + wr * h2r - wi * h2i
    Y(i2) =   h1i + wr * h2i + wi * h2r
    Y(i3) =   h1r - wr * h2r + wi * h2i
    Y(i4) = - h1i + wr * h2i + wi * h2r
    wtemp = wr
    wr = wr * wpr - wi * wpi + wr
    wi = wi * wpr + wtemp * wpi + wi
  end do     
  if (isign == 1) then   
    h1r  = Y(1)
    Y(1) = h1r + Y(2)
    Y(2) = h1r - Y(2)       
  else
    h1r  = Y(1)
    Y(1) = c1 * (h1r + Y(2))
    Y(2) = c1 * (h1r - Y(2))                
    call FOUR (Y, N, -1)
  end if
end subroutine REALFT  
!*********************************************************************************
subroutine FOUR (Y, NN, isign)
  use parameters
  implicit none
  integer :: NN, isign, n, j, i, m, mmax, istep
  real(O) :: Y(2*NN), tempr, tempi, wr, wi, wpr, wpi, wtemp, theta   
!      
  n = 2 * nn              
  j = 1
  do i = 1, n, 2
    if (j > i) then
      tempr  = Y(j)
      tempi  = Y(j+1)
      Y(j)   = Y(i)
      Y(j+1) = Y(i+1)
      Y(i)   = tempr
      Y(i+1) = tempi
    end if
    m = n / 2  
    do while (m >= 2 .and. j > m)
      j = j - m
      m = m / 2   
    end do
    j = j + m
  end do            
  mmax = 2
  do while (n > mmax)
    istep =   2 * mmax
    theta =   2._O * Pi / (isign * mmax)
    wpr   = - 2._O * sin(0.5_O * theta)**2
    wpi   =   sin(theta)
    wr = 1._O
    wi = 0._O
    do m = 1, mmax, 2
      do i = m, n, istep
        j = i + mmax
        tempr  = wr * Y(j) - wi * Y(j+1)
        tempi  = wr * Y(j+1) + wi * Y(j)
        Y(j)   = Y(i) - tempr
        Y(j+1) = Y(i+1) - tempi
        Y(i)   = Y(i) + tempr
        Y(i+1) = Y(i+1) + tempi
      end do
      wtemp = wr
      wr = wr * wpr - wi * wpi + wr
      wi = wi * wpr + wtemp * wpi + wi
    end do
    mmax = istep
  end do                                                
end subroutine FOUR 
