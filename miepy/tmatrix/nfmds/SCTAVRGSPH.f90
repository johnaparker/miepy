subroutine SCTAVRGSPH
!------------------------------------------------------------------------------------ 
! 1. General Considerations                                                         !
! --------------------------                                                        !     
! The routine computes the scattering characteristics of an ensemble of             !
! polydisperse, homogeneous spherical particles. The external excitation is a       !
! vector plane wave propagating along the Z-axis of the global coordinate system    !
! and the scattering plane is the XZ-plane. The code computes the following average !
! quantities:                                                                       ! 
! - the extinction and scattering cross sections,                                   !
! - the asymmetry parameter,                                                        ! 
! - the scattering matrix at a set of NthetaAvrg scattering angles,                 !
! - the differential scattering cross sections (DSCS) at a set of NthetaGSAvrg      !
!   scattering angles and for a specified incident polarization state               !
!   (alphapAvrg).                                                                   !
! The scattering matrix is computed at scattering angles which are uniformly        !
! spaced in the interval (thetaminAvrg, thetamaxAvrg), while the differential       !
! scattering cross sections are computed at scattering angles which are uniformly   !
! spaced in the interval (0, 180°).                                                 !
!                                                                                   !
! The basic optical characteristics of a homogeneous spherical particles can be     !
! directly expressed in terms of the Lorenz-Mie coefficients Fn and Gn, and         !
! simple formulas are available for the average quantities. For example, the        !
! average scattering cross section is given by                                      !
!                                                                                   ! 
!                       2*Pi                                                        !
!              Cscat = ------ SUM { (2n + 1) * [<|Fn|**2> + <|Gn|**2>] }            !
!                       k**2                                                        !
!                                                                                   !
! and the the scattering efficiency by                                              !
!                                                                                   ! 
!                         Cscat                                                     ! 
!              Qscat = -------------,                                               ! 
!                       Pi*anorm**2                                                 !
!                                                                                   !
! where anorm is the normalization length and k is the wave number. The average     !
! extinction cross section can be computed as                                       !
!                                                                                   !           
!                         2*Pi                                                      !
!               Cext = - ------ RE { SUM { (2n + 1) * [<Fn> + <Gn>] }}              !
!                         k**2                                                      !
!                                                                                   !
! and the the extinction efficiency as                                              !
!                                                                                   !
!                          Cext                                                     !
!               Qext = -------------.                                               !
!                       Pi*anorm**2                                                 !
!                                                                                   !
! The average asymmetry parameter is given by                                       ! 
!                                                                                   !
!                          2*Pi                                                     !
!                  g = ------------ RE{ SUM { a(n) * [<Fn * conjg(Fn+1)             ! 
!                       Cscat*k**2                                                  !
!                      + Gn * conjg(Gn+1)>] + b(n) * < Fn * conjg(Gn) > } },        !
!                                                                                   !
! where                                                                             !
!                       n * (n + 2)                                                 !
!               a(n) = -------------,                                               !
!                          n + 1                                                    !
! and                                                                               !
!                        2 * n + 1                                                  !
!               b(n) = -------------.                                               ! 
!                       n * (n + 1)                                                 !
!                                                                                   !
! The following analytical size distribution are considered:                        !
! - the modified gamma distribution                                                 !
!                                                                                   !    
!         n(a) ~ exp { alpha * log(a) - alpha * (a / ac)**gamma / gamma },          !
!                                                                                   !
! - the log normal distribution                                                     !
!                                                                                   !
!         n(a) ~ exp { - log(a) - 0.5_O * (log(a/ag) / log(sg))**2 },               !
!                                                                                   ! 
! - the gamma distribution                                                          !
!                                                                                   !
!         n(a) ~ exp {(1._O - 3._O * gamma) * log(a) / gamma - a / alpha / gamma }, !
!                                                                                   !
! - the power law distribution                                                      !
!                                                                                   !
!         n(a)  ~   exp { - alpha * log(a) },                                       !
!                                                                                   !
! - the modified bimodal log normal distribution                                    !
!                                                                                   !
!         n(a)  ~ exp { - 4._O * log(a) - 0.5_O * (log(a/ag)  / log(sg))**2 }       !
!                 + gamma * exp { - 4._O * log(a) - 0.5_O * (log(a/ag1)             !
!                 / log(sg1))**2},                                                  !
!                                                                                   !
! where a is the particle radii and alpha, gamma, ag, sg, ag1 and sg1 are the       !
! parameters of the size distributions. The constant for each size distribution     !
! is chosen such that the size distribution satisfies the normalization             !
! condition. For analytical size distributions, the particle radii varies           !
! between zero and infinity. In the code, we consider truncated size                !
! distributions, i.e., the particle radii varies in the interval (amin, amax).      !
! The numerical integration of scattering characteristics over a size               !
! distribution is performed with Simpson's rule and the number of integration       !
! points Nint must be an odd number.                                                !
!                                                                                   !
! The code also computes the following characteristics of the size distribution:    !
! - the effective radius                                                            !
!                                                                                   !
!         a_eff = < a * Pi * a**2 > / area_avrg,                                    !
!                                                                                   !
! - the effective variance                                                          !
!                                                                                   ! 
!         v_eff = < (a - aeff)**2 * Pi * a**2 > / (area_avrg * aeff**2>),           !
!                                                                                   !
! - the average area of the geometric projection                                    !
!                                                                                   !
!         area_avrg  = < Pi * a**2 >,                                               !
!                                                                                   !
! - the average radius                                                              !
!                                                                                   !
!         a_avrg = < a >,                                                           !
!                                                                                   ! 
! - the average volume                                                              !
!                                                                                   !
!         v_avrg = < 4/3 * Pi * a**3 >,                                             !
!                                                                                   ! 
! - the volume-weighted average radius                                              !
!                                                                                   !
!         avw_avrg = < a * 4/3 * Pi * a**3 > / v_avrg.                              ! 
!                                                                                   !
! The average quantities are expressed in terms of three average quantities         !
! <S_pq(e_r,e_z) S_p1q1(e_r,e_z)*>, where S_pq are the elements of the              ! 
! amplitude matrix, e_r is the unit vector along the scattering direction           !
! (theta,0) and e_z is the unit vector along the Z axis. The average quantities     !
! are computed at NthetaGSAvrg scattering angles, which are uniformly spaced in     !  
! the interval (0,180°). The scattering matrix is computed at the same sample       !
! angles and polynomial interpolation is used to evaluate the scattering matrix     ! 
! at any zenith angle theta (in the range (thetaminAvrg, thetamaxAvrg)). The        !
! same scattering angles are used for computing the DSCS and the average            !
! quantities <S_pq S_p1q1*>.                                                        !
!                                                                                   !
! 2. Input Parameters                                                               !
! --------------------                                                              !
! The input parameters specified in the file "/INPUTFILES/InputSCTAVRGSPH.dat"      !
! are listed below.                                                                 !
!                                                                                   !
! - wavelength (real) - wavelenght of the ambient medium.                           !
!                                                                                   !
! - ind_ref (complex) - refractive index of the spherical particles.                !
!                                                                                   !
! - anormAvrg (real) - characteristic length of the ensemble which is used          !
!   to normalize the average differential scattering cross sections.                !
!                                                                                   ! 
! - TypeDist (integer) - specifies the type of the size distribution. The           !
!   permissive values are:                                                          !
!   - 1 for the modified gamma distribution,                                        !
!   - 2 for the log normal distribution,                                            !
!   - 3 for the gamma distribution,                                                 !
!   - 4 for the power law distribution,                                             !   
!   - 5 for the modified bimodal log normal distribution.                           !
!                                                                                   ! 
! - amin, amax (real variables) - specify the minimal and maximal values of the     !
!   particle radii for the truncated size distributions. These values are the       !  
!   integration limits for particle size averaging.                                 !
!                                                                                   !
! - Npar (integer) - number of parameters characterizing the size distributions.    !
!   The permissive values are:                                                      !
!   - 3 for the modified gamma distribution,                                        !
!   - 2 for the log normal distribution,                                            !
!   - 2 for the gamma distribution,                                                 !
!   - 1 for the power law distribution,                                             !   
!   - 5 for the modified bimodal log normal distribution.                           !
!                                                                                   ! 
! - par (real array) - parameters of the size distributions. The significance of    !
!   the parameters is as follows:                                                   !
!   - par(1) = alpha, par(2) = gamma and par(3) = ac for the modified gamma         !
!     distribution,                                                                 !
!   - par(1) = ag and par(2) = sg for the log normal distribution,                  !
!   - par(1) = alpha and par(2) = gamma for the gamma distribution,                 !
!   - par(1) = alpha for the power law distribution,                                !  
!   - par(1) = ag, par(2) = sg, par(3) = ag1, par(4) = sg1 and par(5) = gamma       !
!     for the modified bimodal log normal distribution.                             !
!                                                                                   !
! - Nint (integer) - number of integration points for the numerical integration     !
!   of the scattering characteritics over a size distribution. Because the          !
!   numerical integration is performed with Simpson's rule, Nint should be an       !
!   ODD number.                                                                     !
!                                                                                   !
! - alphapAvrg (real) - polarization angle of the incident wave. In general, if     ! 
!  (e_r,e_theta,e_phi) are the spherical unit vectors of the incident               !
!   direction (thetaGI,phiGI), then alphapAvrg is the angle between the electric    ! 
!   field vector and the unit vector e_theta. In the present case, the incident     !
!   direction is (0,0) and alphapAvrg is the angle between the electric field       !
!   vector and the X-axis of the global coordinate system. This parameter is        !
!   used for computing the average differential scattering cross sections.          !
!                                                                                   !
! - ComputeDSCSAvrg (logical) - if ComputeDSCSAvrg = t, the average differential    ! 
!   scattering cross sections are computed in the azimuthal plane phi = 0°.         !
!                                                                                   !
! - NthetaGSAvrg (integer) - number of zenith angles at which the average           !
!   differential scattering cross section are computed. The scattering angles       !
!   are uniformly spaced in the interval (0, 180°).                                 !
!                                                                                   !
! - normalizedAvrg (logical) - if normalizedAvrg = t, the average differential      ! 
!   scattering cross section are normalized by Pi * anormAvrg**2, where             !
!   anormAvrg is the characteristic lenght of the ensemble.                         !
!                                                                                   !
! - FileDSCSAvrg (character(80)) - name of the file to which the average            !
!   differential scattering cross sections are written.                             !
!                                                                                   !
! - ComputeScatParAvrg (logical) - if ComputeScatParAvrg = t, the scattering        !
!   characteristics are computed at specified scattering directions.                !
!   The following average quantities are computed: the extinction and               !
!   scattering cross sections, the asymmetry parameter and the scattering           !
!   matrix at NthetaAvrg scattering angles.                                         !
!                                                                                   !
! - NthetaAvrg (integer) - number of zenith angle values in the interval            !
!  (thetaminAvrg, thetamaxAvrg) at which the scattering matrix is computed.         !
!   The zenith angles are uniformly spaced.                                         !    
!                                                                                   !
! - thetaminAvrg, thetamaxAvrg (real variables) - specify the minimal and           !
!   maximal values of the zenith angle theta at which the scattering matrix         !
!   is computed.                                                                    !
!                                                                                   !
! - FileSCATAvrg (character(80)) - name of the file to which the scattering         !
!   characteristics are written.                                                    !
!                                                                                   ! 
! - PrnProgress (logical) - if PrnProgress = t, the progress of calculation         ! 
!   is printed.                                                                     !
!                                                                                   !  
! - WriteInputInfo (logical) - if WriteInputInfo = t, the input parameters of       !
!   the scattering problem are written to the output files FileDSCSAvrg and         !
!   FileSCATAvrg.                                                                   !
!------------------------------------------------------------------------------------
  use parameters
  implicit none
  integer       :: TypeDist, Nrank, NthetaGSAvrg, NthetaAvrg, Nint, Npar, i                   
  real(O)       :: wavelength, k, anormAvrg, snormAvrg, Cext, Cscat, Qext, Qscat,   &
                   AsymPar, thetaGS, Z(4,4), alphapAvrg, amin, amax, par(NparPD),   &
                   Cnst, thetaminAvrg, thetamaxAvrg, dtheta, AvrgArea, EffRadius,   &
                   AvrgRadius, AvrgVolume, VWAvrgRadius, EffVar
  complex(O)    :: ind_ref
  character(80) :: FileDSCSAvrg, FileScatAvrg 
  logical       :: ComputeDSCSAvrg, normalizedAvrg, ComputeScatParAvrg,             &
                   PrnProgress, WriteInputInfo
  real(O),allocatable    :: ZE(:,:)
  complex(O),allocatable :: SS(:,:)
! -----------------------------------------------------------------------------------
!                               Read the input file                                 !
! -----------------------------------------------------------------------------------
  call readinputSCTAVRGSPH ( wavelength, ind_ref, anormAvrg, TypeDist,              &
       amin, amax, Npar, par, Nint, alphapAvrg, ComputeDSCSAvrg, NthetaGSAvrg,      &
       normalizedAvrg, FileDSCSAvrg, ComputeScatParAvrg, NthetaAvrg,                &
       thetaminAvrg, thetamaxAvrg, FileScatAvrg, PrnProgress, WriteInputInfo,       &
       k, snormAvrg, Nrank, Cnst, AvrgArea, EffRadius, AvrgRadius, AvrgVolume,      &
       VWAvrgRadius, EffVar ) 
! -----------------------------------------------------------------------------------
!                                       Main                                        !
! -----------------------------------------------------------------------------------       
  print "(/,2x,'Scattering Characteristics of an Ensemble of Spheres')"
  print "(  2x,'----------------------------------------------------',/)"     
  allocate (SS(3,NthetaGSAvrg))        
  if (PrnProgress)                                                                  &
      print "(/,2x,'progress of main calculation for the average matrix <SS*>:')"    
  call SijSCpqSPH (k, ind_ref, amin, amax, TypeDist, Npar, par, Cnst, Nrank,        &
       Nint, NthetaGSAvrg, SS, PrnProgress)                   
  call AvCscatCextSPH (k, snormAvrg, ind_ref, amin, amax, TypeDist, Npar, par,      &
       Cnst, Nrank, Nint, Cext, Cscat, Qext, Qscat, AsymPar)
  if (ComputeDSCSAvrg) then     
    open (unit = iDSCS, file = FileDSCSAvrg, status = "replace")                         
    if (WriteInputInfo) call inputDSCS_SCATAvrgSPH (.true., amin, amax, TypeDist,   &
                             Npar, par, Nrank, Nint, alphapAvrg, wavelength,        &
                             anormAvrg, normalizedAvrg)
    call DiffScatCrossSectAvrgSPH (anormAvrg, alphapAvrg, NthetaGSAvrg, SS,         &
         normalizedAvrg, Cext, Cscat, Qext, Qscat)                                                                             
    close (unit = iDSCS)                       
  end if        
  if (ComputeScatParAvrg) then
    if (NthetaAvrg /= 1) then
      dtheta = (thetamaxAvrg - thetaminAvrg) / (NthetaAvrg - 1)
    else
      dtheta = 0._O
    end if
    open (unit = iSCAT, file = FileScatAvrg, status = "replace")
    if (WriteInputInfo) call inputDSCS_SCATAvrgSPH (.false., amin, amax, TypeDist,  &
                             Npar, par, Nrank, Nint, alphapAvrg, wavelength,        &
                             anormAvrg, normalizedAvrg)     
    allocate (ZE(4,NthetaGSAvrg))
    call ExtendedScatteringMatrixSPH (NthetaGSAvrg, SS, ZE)     
    do i = 1, NthetaAvrg
      thetaGS = thetaminAvrg + (i - 1) * dtheta                        
      call ScatteringMatrixSPH (thetaGS, NthetaGSAvrg, ZE, Z)
      call write_ScatMatrixAvrg (thetaGS, Z, i, Cscat, Cext, Qscat, Qext, AsymPar,  &
           AvrgArea, EffRadius, AvrgRadius, AvrgVolume, VWAvrgRadius, EffVar)
    end do                   
    close (unit = iSCAT)                                         
    deallocate (ZE)                                      
  end if  
  deallocate (SS)            
end subroutine SCTAVRGSPH
!***********************************************************************************
subroutine readinputSCTAVRGSPH ( wavelength, ind_ref, anormAvrg, TypeDist,          &
           amin, amax, Npar, par, Nint, alphapAvrg, ComputeDSCSAvrg, NthetaGSAvrg,  &
           normalizedAvrg, FileDSCSAvrg, ComputeScatParAvrg, NthetaAvrg,            &
           thetaminAvrg, thetamaxAvrg, FileScatAvrg, PrnProgress, WriteInputInfo,   &
           k, snormAvrg, Nrank, Cnst, AvrgArea, EffRadius, AvrgRadius, AvrgVolume,  &
           VWAvrgRadius, EffVar )
  use parameters
  use derived_parameters
  implicit none
  integer       :: TypeDist, Nrank, NthetaGSAvrg, NthetaAvrg, Nint, Npar, i, ios                   
  real(O)       :: wavelength, k, anormAvrg, xnorm, snormAvrg, xmax, alphapAvrg,    &
                   amin, amax, par(NparPD), Cnst, thetaminAvrg, thetamaxAvrg, grd,  &
                   AvrgArea, EffRadius, AvrgRadius, AvrgVolume, VWAvrgRadius, EffVar
  complex(O)    :: ind_ref
  character(80) :: FileDSCSAvrg, FileScatAvrg, string 
  logical       :: ComputeDSCSAvrg, normalizedAvrg, ComputeScatParAvrg,             &
                   PrnProgress, WriteInputInfo, XFindPar
! -----------------------------------------------------------------------------------
!                        Read the input file FileInputSCTAVRGSPH                    !
! -----------------------------------------------------------------------------------
  call DrvParameters   
  open (unit = iInputSCTAVRGSPH, file = FileInputSCTAVRGSPH, status = "old",        &
        position = "rewind")   
  wavelength = 0.1_O * 2._O * Pi  
  string     = 'OptProp'
  if (XFindPar (iInputSCTAVRGSPH, string)) then
    read (iInputSCTAVRGSPH, *, iostat = ios) wavelength
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable wavelength;')"
      stop
    end if            
  else
    print "(/,2x,'Group name OptProp not found;')"
    stop  
  end if  
  k = 2._O * Pi / wavelength
!
  ind_ref = (1.5_O,0._O)
  string  = 'Tmat'
  if (XFindPar (iInputSCTAVRGSPH, string)) then
    read (iInputSCTAVRGSPH, *, iostat = ios) ind_ref
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_ref;')"
      stop
    end if            
  else
    print "(/,2x,'Group name Tmat not found;')"
    stop  
  end if  
  call check_ind_ref (ind_ref) 
!
  anormAvrg = 1._O
  string    = 'GeomProp'
  if (XFindPar (iInputSCTAVRGSPH, string)) then
    read (iInputSCTAVRGSPH, *, iostat = ios) anormAvrg
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable anormAvrg;')"
      stop
    end if            
  else
    print "(/,2x,'Group name GeomProp not found;')"
    stop  
  end if  
  call check_anorm (anormAvrg)             
  xnorm     = k * anormAvrg
  snormAvrg = Pi * xnorm * xnorm  
!
  TypeDist = 1
  amin = 1._O
  amax = 20._O
  Npar = 3
  do i = 1, NparPD
    par(i) = 1._O
  end do
  Nint   = 181
  string = 'SizeDst'
  if (XFindPar (iInputSCTAVRGSPH, string)) then
    read (iInputSCTAVRGSPH, *, iostat = ios) TypeDist
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable TypeDist;')"
      stop
    end if
    read (iInputSCTAVRGSPH, *, iostat = ios) amin
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable amin;')"
      stop
    end if
    read (iInputSCTAVRGSPH, *, iostat = ios) amax
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable amax;')"
      stop
    end if
    read (iInputSCTAVRGSPH, *, iostat = ios) Npar
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Npar;')"
      stop
    end if
    do i = 1, Npar
      read (iInputSCTAVRGSPH, *, iostat = ios) par(i)
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable par;')"
        stop
      end if
    end do
    read (iInputSCTAVRGSPH, *, iostat = ios) Nint
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nint;')"
      stop
    end if
  else
    print "(/,2x,'Group name SizeDst not found;')"
    stop  
  end if 
  call check_inputSCTAVRGSPH (amin, amax, TypeDist, Npar)
  call check_NSimpson (Nint)
  xmax  = k * amax
  Nrank = int(xmax + 4.05_O * xmax**0.33_O + 2._O)                  
  call GeometryPars (amin, amax, TypeDist, Npar, par, Nint, Cnst, AvrgArea,         &
       EffRadius, AvrgRadius, AvrgVolume, VWAvrgRadius, EffVar)
!
  alphapAvrg = 45._O
  string     = 'IncWaveAvrg'
  if (XFindPar (iInputSCTAVRGSPH, string)) then
    read (iInputSCTAVRGSPH, *, iostat = ios) alphapAvrg
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable alphapAvrg;')"
      stop
    end if            
  else
    print "(/,2x,'Group name IncWaveAvrg not found;')"
    stop  
  end if  
  call check_polarization_angle (alphapAvrg)
  grd = Pi / 180._O
  alphapAvrg = alphapAvrg * grd
!
  ComputeDSCSAvrg = .true. 
  NthetaGSAvrg    = 91 
  normalizedAvrg  = .true. 
  FileDSCSAvrg    = '../OUTPUTFILES/dscsavrg_sph.dat'
  string          = 'DSCSAvrg'
  if (XFindPar (iInputSCTAVRGSPH, string)) then
    read (iInputSCTAVRGSPH, *, iostat = ios) ComputeDSCSAvrg
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ComputeDSCSAvrg;')"
      stop
    end if 
    read (iInputSCTAVRGSPH, *, iostat = ios) NthetaGSAvrg
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable NthetaGSAvrg;')"
      stop
    end if
    read (iInputSCTAVRGSPH, *, iostat = ios) normalizedAvrg
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable normalizedAvrg;')"
      stop
    end if
    read (iInputSCTAVRGSPH, *, iostat = ios) FileDSCSAvrg
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable FileDSCSAvrg;')"
      stop
    end if
  else
    print "(/,2x,'Group name DSCSAvrg not found;')"
    stop  
  end if 
!
  ComputeScatParAvrg = .true.
  NthetaAvrg    = 1
  thetaminAvrg  = 72._O 
  thetamaxAvrg  = 72._O    
  FileScatAvrg  = '../OUTPUTFILES/scat_avrg_sph.dat'
  string        = 'ScatParAvrg'
  if (XFindPar (iInputSCTAVRGSPH, string)) then
    read (iInputSCTAVRGSPH, *, iostat = ios) ComputeScatParAvrg
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ComputeScatParAvrg;')"
      stop
    end if 
    read (iInputSCTAVRGSPH, *, iostat = ios) NthetaAvrg 
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable NthetaAvrg;')"
      stop
    end if
    read (iInputSCTAVRGSPH, *, iostat = ios) thetaminAvrg
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable thetaminAvrg;')"
      stop
    end if
    read (iInputSCTAVRGSPH, *, iostat = ios) thetamaxAvrg
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable thetamaxAvrg;')"
      stop
    end if
    read (iInputSCTAVRGSPH, *, iostat = ios) FileScatAvrg
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable FileScatAvrg;')"
      stop
    end if
  else
    print "(/,2x,'Group name ScatParAvrg not found;')"
    stop  
  end if         
  call check_tetaminmax (thetaminAvrg, thetamaxAvrg, NthetaAvrg)
  thetaminAvrg = thetaminAvrg * grd
  thetamaxAvrg = thetamaxAvrg * grd
!
  PrnProgress    = .true.
  WriteInputInfo = .true.
  string         = 'PrintInfo'
  if (XFindPar (iInputSCTAVRGSPH, string)) then
    read (iInputSCTAVRGSPH, *, iostat = ios) PrnProgress
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable PrnProgress;')"
      stop
    end if 
    read (iInputSCTAVRGSPH, *, iostat = ios) WriteInputInfo
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable WriteInputInfo;')"
      stop
    end if
  else
    print "(/,2x,'Group name PrintInfo not found;')"
    stop  
  end if        
  close (unit = iInputSCTAVRGSPH)
end subroutine readinputSCTAVRGSPH
