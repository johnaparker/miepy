! ***********************************************************************************
! *          SCATTERING CHARACTERISTICS: INCOMPLETE UNIFORM DISTRIBUTION            *
! *    ---------------------------------------------------------------------------  *
! *    Partial list of subroutines:                                                 *
! *     ScatCharact,             DSCS_SCAT,               DSCS,                     *          
! *     F_azimuthal,             F_azimuthal_tetaGS,      CQscat,                   *
! *     CQext,                   CQextGen,                CQextPolGen,              *
! *     CQextPolPlane(UNUSED)    PhaseMatrix,             ExtinctionMatrix,         *  
! *     extend_vector_positive,  extend_vector_negative,  delta_DSCS,               *  
! *     NodesAlphaBetaGamma,     IndexElements,           DiffScatCrossSectPARTSUB, *
! *     DSCSPARTSUB,             EMFPARTSUB,              delta_DSCSPARTSUB,        *
! *     readinputSCT,            readinputPARTSUB1                                  *
! ***********************************************************************************
subroutine ScatCharact (k, FileTmat, Mrank, Nrank, axsym, sphere, chiral)
!------------------------------------------------------------------------------------
! The subroutine computes the scattering characteristics of a particle using the    !
! T matrix stored in file FileTmat. The input parameters are:                       !
! - k (real) - wave number of the ambient medium.                                   !
! - FileTmat (character(80)) - name of the file containing the T matrix.            !
! - Mrank, Nrank (integer variables) - Mrank is the number of azimuthal modes       !
!   (maximum degree), while Nrank is the maximum expansion order.                   !
! - axsym (logical) - if axsym = t, the scatterer is a rotationally symmetric       !
!   particle (axisymmetric particle).                                               !
! - sphere (logical) - if sphere = t, the scatterer is a spherical particle.        !
! - chiral (logical) - if chiral = t, the scatterer is an optical active particle   !
!  (chiral particle).                                                               !
!------------------------------------------------------------------------------------
  use parameters
  use allocation, only: IndI, IndJ, NameElem
  implicit none
  integer        :: Nrank, Mrank, Nalpha, Nbeta, Ngamma, NthetaGS, NthetaRND, Nphi, &
                    Ntheta(NphiMax), Nelem, MatrixElem(16), NthetaAsym, NphiAsym,   &
                    Nrankeff, Mrankeff, ifail, itheta, Nfail, ialpha, ibeta,        &
                    igamma, i, j, kelem, iphi, NthetaAL, NphiAL
  real(O)        :: k, wavelength, alphamin, alphamax, betamin, betamax, gammamin,  &
                    gammamax, anorm, snorm, thetaGI, phiGI, alphapGauss, x0, y0,    &
                    z0, w0, deltaOrder, phiGS, thetaminRND, thetamaxRND,            &
                    phi(NphiMax), thetamin(Nphimax), thetamax(NphiMax), theta,      &
                    dtheta, AsymPar, AsymParV, Cext, Cscat, Qext, Qscat, Z(4,4),    &
                    CextAv, CscatAv, QextAv, QscatAv, alpha, beta, gamma, fact,     &
                    norm, CscatPAv, CscatSAv, CextPAv, CextSAv, QscatPAv, QscatSAv, &
                    QextPAv, QextSAv, gP_xAv, gP_yAv, gP_zAv, gS_xAv, gS_yAv,       &
                    gS_zAv, gP_x, gP_y, gP_z, gS_x, gS_y, gS_z, gP_t, gP_p, gP_r,   &
                    gS_t, gS_p, gS_r, KeAv(4,4), Ke(4,4), CscatP, CscatS, CextP,    &
                    CextS, CscatV, CextV, QscatV, QextV, aphi, bphi, atheta, btheta
  complex(O)     :: epol_beta, epol_alpha		    
  character(5)   :: TypeExcit  
  character(80)  :: FileTmat, FileDSCS, FileSCAT
  character(256) :: FailMessage(10)
  logical        :: MirorSym, axsym, sphere, chiral, RandomOrientation, DoNumAvrg,  &
                    UseSimpson, ReducedOrder, ExtThetaDom, normalized, ComputeDSCS, &
                    ComputeScatPar, FailHovTest, ComputeAsymPar, PrnProgress,       &
                    WriteInputInfo
  real(O),allocatable      :: ZE(:,:,:), walpha(:), xalpha(:), wbeta(:), xbeta(:),  &
                              wgamma(:), xgamma(:), hAv(:), vAv(:), h(:), v(:),     &
                              ZTPAv(:,:,:), ZTP(:,:,:), wphi(:), xphi(:), wtheta(:),&
                              xtheta(:)
  complex(O),allocatable   :: SS(:,:)
! -----------------------------------------------------------------------------------
!                                Read the input file                                !
! -----------------------------------------------------------------------------------           
  call readinputSCT ( .false., wavelength, k, FileTmat, Nrank, Mrank,               &
       axsym, sphere, chiral, RandomOrientation, MirorSym, DoNumAvrg,               &
       UseSimpson, ReducedOrder, deltaOrder, alphamin, alphamax, Nalpha,            &
       betamin, betamax, Nbeta, gammamin, gammamax, Ngamma, anorm, TypeExcit,       &
       thetaGI, phiGI, epol_beta, epol_alpha, alphapGauss, x0, y0, z0, w0,          &
       ComputeDSCS, phiGS, NthetaGS, ExtThetaDom, normalized, FileDSCS,             &
       ComputeScatPar, NthetaRND, thetaminRND, thetamaxRND, Nphi, phi, Ntheta,      &
       thetamin, thetamax, FileScat, Nelem, MatrixElem, ComputeAsymPar,             &
       NthetaAsym, NphiAsym, PrnProgress, WriteInputInfo, snorm )
  if (.not. ComputeDSCS .and. .not. ComputeScatPar) return     
! -----------------------------------------------------------------------------------
!                                       Main                                        !
! -----------------------------------------------------------------------------------              
  if (RandomOrientation) then          
    print "(/,2x,'Scattering Characteristics of a Randomly Oriented Particle')"  
    print "(  2x,'----------------------------------------------------------')"     
    allocate (SS(10,NthetaGS))    
    if (PrnProgress)                                                                &
    print "(/,2x,'progress of main calculation for the average matrix <SS*>:')"
    if (.not. DoNumAvrg) then                   
      call AvCscatCextEff (ReducedOrder, k, snorm, FileTmat, axsym, chiral,         &
           deltaOrder, Nrank, Mrank, Nrankeff, Mrankeff, Cext, Cscat,               &
           Qext, Qscat)                                                                       
      print "(2x,'effective number of azimuthal modes, Mrankeff = ',i3,';')",       &
              Mrankeff
      print "(2x,'effective maximum expansion order,   Nrankeff = ',i3,';')",       &
              Nrankeff   								                                        
      call SijSCpqEff (MirorSym, k, axsym, chiral, FileTmat, Nrank, Mrank,          &
           Nrankeff, Mrankeff, NthetaGS, SS, PrnProgress)                       
    else                                                               
      call SijSCpqInt (MirorSym, k, snorm, FileTmat, axsym, chiral, Nrank, Mrank,   &
           Nalpha, Nbeta, Ngamma, NthetaGS, UseSimpson, Cext, Cscat, Qext, Qscat,   &
           SS, PrnProgress)
    end if                                                                  
    CscatV = 0._O
    CextV  = 0._O
    QscatV = 0._O
    QextV  = 0._O    
    if (ComputeScatPar) then
      if (NthetaRND /= 1) then
        dtheta = (thetamaxRND - thetaminRND) / (NthetaRND - 1)
      else
        dtheta = 0._O
      end if  
      open (unit = iSCAT, file = FileScat, status = "replace")
      if (WriteInputInfo) call inputDSCS_SCATAvrg (.false., wavelength, axsym,      &
                               chiral, FileTmat, Nrank, Mrank, anorm, normalized,   &                              
                               DoNumAvrg, Nalpha, Nbeta, Ngamma, NthetaGS,          &
                               epol_beta, epol_alpha) 
      call RandomExtinctionMatrix (k, FileTmat, MirorSym, axsym, chiral, Nrank,     &
           Mrank, KeAv)     	   			                                            
      allocate (ZE(4,4,NthetaGS))
      call ExtendedScatteringMatrix (MirorSym, NthetaGS, SS, ZE)
      if (.not.MirorSym) call AvCscatCextV (k, snorm, NthetaGS, ZE, KeAv, CscatV,   &
                              CextV, QscatV, QextV)            
      call AsymmetryParameter (MirorSym, NthetaGS, ZE, Cscat, AsymPar, AsymParV)
      FailHovTest = .false.
      Nfail = 0
      do ifail = 1, 10
        FailMessage(ifail) = ' '
      end do                                
      do itheta = 1, NthetaRND
        theta = thetaminRND + (itheta - 1) * dtheta                          
        call ScatteringMatrix (MirorSym, theta, NthetaGS, ZE, Z)                      
        if (axsym) then
          FailHovTest = .false. 
          call HovTest (Z, FailHovTest, Nfail, FailMessage)
        end if          
        call write_ScatMatrix (MirorSym, axsym, FailHovTest, Nfail, FailMessage,    &
             KeAv, theta, Z, itheta, Nelem, IndI, IndJ, NameElem, Cscat, Cext,      &
             Qscat, Qext, CscatV, CextV, QscatV, QextV, AsymPar, AsymParV)
      end do        
      close (unit = iSCAT)    
      deallocate (ZE)                                            
    end if
    if (ComputeDSCS) then
      open (unit = iDSCS, file = FileDSCS, status = "replace")       			      
      if (WriteInputInfo) call inputDSCS_SCATAvrg (.true., wavelength, axsym,       &
                               chiral, FileTmat, Nrank, Mrank, anorm, normalized,   &                              
                               DoNumAvrg, Nalpha, Nbeta, Ngamma, NthetaGS,          &
                               epol_beta, epol_alpha)
      if (.not.MirorSym)  then      
        call RandomExtinctionMatrix (k, FileTmat, MirorSym, axsym, chiral, Nrank,   &
             Mrank, KeAv)     	   			                                            
        allocate (ZE(4,4,NthetaGS))
        call ExtendedScatteringMatrix (MirorSym, NthetaGS, SS, ZE)
        call AvCscatCextV (k, snorm, NthetaGS, ZE, KeAv, CscatV, CextV, QscatV,     &
	     QextV)                  
      end if  			       
      call DiffScatCrossSectAvrg (MirorSym, SS, NthetaGS, anorm, epol_beta,         &
           epol_alpha, phiGS, normalized, Cscat, Cext, Qscat, Qext, CscatV, CextV,  &
           QscatV, QextV)
      close (unit = iDSCS) 
      if (.not.MirorSym) deallocate (ZE)                          
    end if                       
    deallocate (SS)   
  else      
    print "(/,2x,'Scattering Characteristics of a Nonrandomly Oriented Particle')" 
    print "(  2x,'-------------------------------------------------------------')" 
    allocate (walpha(Nalpha), xalpha(Nalpha), wbeta(Nbeta), xbeta(Nbeta),           &
              wgamma(Ngamma), xgamma(Ngamma))                                   
    NphiAL   = Nphi
    NthetaAL = 0
    do iphi = 1, Nphi
      if (Ntheta(iphi) > NthetaAL) NthetaAL = Ntheta(iphi)
    end do
    allocate (h(NthetaGS), v(NthetaGS), ZTP(NphiAL,NthetaAL,Nelem))    
    allocate (xphi(NphiAsym), wphi(NphiAsym), xtheta(NthetaAsym), wtheta(NthetaAsym))                
    call NodesAlphaBetaGamma (alphamin, alphamax, Nalpha, betamin, betamax, Nbeta,  &
         gammamin, gammamax, Ngamma, UseSimpson, xalpha, walpha, xbeta, wbeta,      &
         xgamma, wgamma)	    	                                         
    if (ComputeDSCS)    open (unit = iDSCS, file = FileDSCS, status = "replace")       
    if (ComputeScatPar) open (unit = iSCAT, file = FileScat, status = "replace")                                            
    if (ComputeDSCS .and. WriteInputInfo)                                           &                      
      call inputDSCS_SCAT (.true., axsym, sphere, chiral, Mrank, Nrank, phiGS,      &
           thetaGI, phiGI, alphamin, alphamax, Nalpha, betamin, betamax, Nbeta,     &
           gammamin, gammamax, Ngamma, epol_beta, epol_alpha, alphapGauss, x0,      &
           y0, z0, w0, TypeExcit, wavelength, anorm, normalized, FileTmat)
    if (ComputeScatPar .and. WriteInputInfo)                                        &   
      call inputDSCS_SCAT (.false., axsym, sphere, chiral, Mrank, Nrank, phiGS,     &
           thetaGI, phiGI, alphamin, alphamax, Nalpha, betamin, betamax, Nbeta,     &
           gammamin, gammamax, Ngamma, epol_beta, epol_alpha, alphapGauss, x0,      &
           y0, z0, w0, TypeExcit, wavelength, anorm, normalized, FileTmat)                  
    if (ComputeDSCS) then
      allocate (hAv(NthetaGS), vAv(NthetaGS))
      CextAv  = 0._O
      CscatAv = 0._O
      do itheta = 1, NthetaGS
        hAv(itheta) = 0._O
        vAv(itheta) = 0._O
      end do    
    end if
    if (ComputeScatPar) then
      allocate (ZTPAv(NphiAL,NthetaAL,Nelem))  
      CextPAv  = 0._O
      CextSAv  = 0._O         
      CscatPAv = 0._O
      CscatSAv = 0._O                  
      do iphi = 1, NphiAL 
        do itheta = 1, NthetaAL               
          do kelem = 1, Nelem        
            ZTPAv(iphi,itheta,kelem) = 0._O
          end do
        end do
      end do
      do i = 1, 4
        do j = 1, 4
          KeAv(i,j) = 0._O
        end do
      end do
      if (ComputeAsymPar) then                    
        gP_xAv = 0._O
        gP_yAv = 0._O
        gP_zAv = 0._O      
        gS_xAv = 0._O
        gS_yAv = 0._O
        gS_zAv = 0._O  
!        
        aphi = 0._O
        bphi = 2._O * Pi
        call Simpson (aphi, bphi, NphiAsym, xphi, wphi)   
        if (.not. UseSimpson) then
          atheta = 0._O
          btheta = Pi   
          call Gauss_Legendre (atheta, btheta, NthetaAsym, wtheta, xtheta)   
          do itheta = 1, NthetaAsym
            theta = xtheta(itheta)
            wtheta(itheta) = wtheta(itheta) * sin(theta)      
          end do        
        else
          atheta = - 1._O
          btheta =   1._O
          call Simpson (atheta, btheta, NthetaAsym, xtheta, wtheta)     
        end if     	                  
      end if 	
    end if
    if (PrnProgress .and. Nalpha > 1)                                               &
        print "(/,2x,'progress of main calculation for the averaging procedure:')"   
    do ialpha = 1, Nalpha
      alpha = xalpha(ialpha)
      if (PrnProgress .and. Nalpha > 1) call write_progress (.false., ialpha, Nalpha)
      do ibeta = 1, Nbeta          
        if (.not. UseSimpson) then 
          beta = xbeta(ibeta)
        else
          beta = acos(xbeta(ibeta))
        end if
        do igamma = 1, Ngamma 
          gamma = xgamma(igamma)           
          fact  = walpha(ialpha) * wbeta(ibeta) * wgamma(igamma)          	       	       	       
          call DSCS_SCAT (ComputeDSCS, ComputeScatPar, ComputeAsymPar, axsym,       &
               sphere, chiral, Mrank, Nrank, NthetaGS, phiGS, Nphi, phi, Ntheta,    &
               thetamin, thetamax, NthetaAsym, NphiAsym, UseSimpson, xphi, wphi,    &
               xtheta, wtheta, thetaGI, phiGI, alpha, beta, gamma, epol_beta,       &
               epol_alpha, alphapGauss, x0, y0, z0, w0, TypeExcit, k, snorm,        &
               FileTmat, ExtThetaDom, normalized, Nelem, IndI, IndJ, h, v, ZTP,     &
               NphiAL, NthetaAL, Ke, Cscat, Cext, CscatP, CscatS, CextP, CextS,     &
               gP_x, gP_y, gP_z, gS_x, gS_y, gS_z)  	       	        
          if (ComputeDSCS) then                                        
            CextAv  = CextAv  + Cext  * fact
            CscatAv = CscatAv + Cscat * fact
            do itheta = 1, NthetaGS
              hAv(itheta) = hAv(itheta) + h(itheta) * fact
              vAv(itheta) = vAv(itheta) + v(itheta) * fact
            end do
          end if   
          if (ComputeScatPar) then
            CextPAv  = CextPAv  + CextP  * fact
            CextSAv  = CextSAv  + CextS  * fact               	    
            CscatPAv = CscatPAv + CscatP * fact
            CscatSAv = CscatSAv + CscatS * fact             	        
            do iphi = 1, Nphi 
              do itheta = 1, Ntheta(iphi)                           
                do kelem = 1, Nelem        
                  ZTPAv(iphi,itheta,kelem) = ZTPAv(iphi,itheta,kelem) +             &
                                             ZTP(iphi,itheta,kelem) * fact
                end do
              end do
            end do
            if (TypeExcit(1:5) == 'PLANE') then  
              do i = 1, 4
                do j = 1, 4
                  KeAv(i,j) = KeAv(i,j) + Ke(i,j) * fact
                end do
              end do           
            end if
	    if (ComputeAsymPar) then
              gP_xAv = gP_xAv + gP_x * fact
              gP_yAv = gP_yAv + gP_y * fact
              gP_zAv = gP_zAv + gP_z * fact	                
              gS_xAv = gS_xAv + gS_x * fact
              gS_yAv = gS_yAv + gS_y * fact
              gS_zAv = gS_zAv + gS_z * fact	                	    	    	    
            end if 	      
          end if
        end do
      end do
    end do      
    norm = Pi * anorm * anorm
    norm = 1._O / norm                 
    if (ComputeDSCS) then
      QextAv  = CextAv  * norm  
      QscatAv = CscatAv * norm      
      call write_DSCS1 (NthetaGS, normalized, ExtThetaDom, hAv, vAv, CscatAv,       &
           CextAv, QscatAv, QextAv)                                            
      close (unit = iDSCS)
      deallocate (hAv, vAv)
    end if
    if (ComputeScatPar) then
      QextPAv  = CextPAv  * norm 
      QextSAv  = CextSAv  * norm                   
      QscatPAv = CscatPAv * norm 
      QscatSAv = CscatSAv * norm 
      gP_t = 0._O
      gP_p = 0._O
      gP_r = 0._O  
      gS_t = 0._O 
      gS_p = 0._O
      gS_r = 0._O      
      if (ComputeAsymPar) then                                           
        call T_cartesian_global_local (gP_xAv, gP_yAv, gP_zAv, phiGI, thetaGI,      &
             0._O, gP_t, gP_p, gP_r)
        call T_cartesian_global_local (gS_xAv, gS_yAv, gS_zAv, phiGI, thetaGI,      &
             0._O, gS_t, gS_p, gS_r)      
      end if	     
      call write_SCAT (Nelem, NameElem, Nphi, phi, Ntheta, thetamin, thetamax,      &
           ZTPAv, NphiAL, NthetaAL, KeAv, CscatPAv, CscatSAv, CextPAv, CextSAv,     &
           QscatPAv, QscatSAv, QextPAv, QextSAv, ComputeAsymPar, gP_t, gP_p, gP_r,  &
           gS_t, gS_p, gS_r, TypeExcit)
      close (unit = iSCAT)             
      deallocate (ZTPAv)      
    end if
    deallocate (walpha, xalpha, wbeta, xbeta, wgamma, xgamma)
    deallocate (h, v, ZTP)
    deallocate (xphi, wphi, xtheta, wtheta)    
  end if
  deallocate (IndI, IndJ, NameElem)  
end subroutine ScatCharact
! **********************************************************************************
subroutine readinputSCT ( CompleteFile, wavelength, k, FileTmat, Nrank, Mrank,      &
           axsym, sphere, chiral, RandomOrientation, MirorSym, DoNumAvrg,           &
	   UseSimpson, ReducedOrder, deltaOrder, alphamin, alphamax, Nalpha,        &
	   betamin, betamax, Nbeta, gammamin, gammamax, Ngamma, anorm, TypeExcit,   &
           thetaGI, phiGI, epol_beta, epol_alpha, alphapGauss, x0, y0, z0, w0,      &
           ComputeDSCS, phiGS, NthetaGS, ExtThetaDom, normalized, FileDSCS,         &
           ComputeScatPar, NthetaRND, thetaminRND, thetamaxRND, Nphi, phi, Ntheta,  &
           thetamin, thetamax, FileScat, Nelem, MatrixElem, ComputeAsymPar,         &
           NthetaAsym, NphiAsym, PrnProgress, WriteInputInfo, snorm )
  use parameters
  use allocation, only: IndI, IndJ, NameElem
  implicit none
  integer        :: Nrank, Mrank, Nalpha, Nbeta, Ngamma, NthetaGS, NthetaRND, Nphi, &
                    Ntheta(NphiMax), Nelem, MatrixElem(16), NthetaAsym, NphiAsym,   &
                    i, ios, iphi
  real(O)        :: k, wavelength, alphamin, alphamax, betamin, betamax, gammamin,  &
                    gammamax, anorm, xpart, snorm, thetaGI, phiGI, alphapGauss,     &
                    x0, y0, z0, w0, deltaOrder, phiGS, thetaminRND, thetamaxRND,    &
                    phi(NphiMax), thetamin(Nphimax), thetamax(NphiMax), grd   
  complex(O)     :: epol_beta, epol_alpha		    
                    
  character(5)   :: TypeExcit  
  character(80)  :: FileDSCS, FileSCAT, FileTmat, string 
  logical        :: CompleteFile, axsym, sphere, chiral, RandomOrientation,         &
                    MirorSym, DoNumAvrg, UseSimpson, ReducedOrder, ExtThetaDom,     &
                    normalized, ComputeDSCS, ComputeScatPar, ComputeAsymPar,        &
                    PrnProgress, WriteInputInfo, XFindPar  
! -----------------------------------------------------------------------------------
!                           Read the input file FileInputSCT                        !
! -----------------------------------------------------------------------------------      
  open (unit = iInputSCT, file = FileInputSCT, status = "old", position = "rewind") 
!
  if (CompleteFile) then
    call DrvParameters           
!
    wavelength = 0.1_O * 2._O * Pi 
    string     = 'OptProp'
    if (XFindPar (iInputSCT, string)) then
      read (iInputSCT, *, iostat = ios) wavelength
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
    FileTmat = '../TMATFILES/T.dat'
    Nrank  = 7
    Mrank  = 4  
    axsym  = .true.
    sphere = .false.
    chiral = .false.
    string = 'Tmat'
    if (XFindPar (iInputSCT, string)) then
      read (iInputSCT, *, iostat = ios) FileTmat
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable FileTmat;')"
        stop
      end if
      read (iInputSCT, *, iostat = ios) Nrank
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Nrank;')"
        stop
      end if
      read (iInputSCT, *, iostat = ios) Mrank
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Mrank;')"
        stop
      end if
      read (iInputSCT, *, iostat = ios) axsym
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable axsym;')"
        stop
      end if
      read (iInputSCT, *, iostat = ios) sphere
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable sphere;')"
        stop
      end if
      read (iInputSCT, *, iostat = ios) chiral
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable chiral;')"
        stop
      end if
    else
      print "(/,2x,'Group name Tmat not found;')"
      stop  
    end if     
  else   
    wavelength = 2._O * Pi / k 
  end if
  grd = Pi / 180._O
! --- initialization of input variables ---                             
  RandomOrientation = .true.
  MirorSym     = .true.
  DoNumAvrg    = .false.
  UseSimpson   = .true.
  ReducedOrder = .true.
  deltaOrder   = 1.e-4_O
  alphamin     = 0._O
  alphamax     = 0._O
  Nalpha       = 1
  betamin      = 0._O
  betamax      = 0._O
  Nbeta        = 1  
  gammamin     = 0._O
  gammamax     = 0._O
  Ngamma       = 1  
  anorm        = 1._O 
  TypeExcit    = 'PLANE'
  thetaGI      =  0._O
  phiGI        =  0._O
  epol_beta    = one 
  epol_alpha   = zero
  alphapGauss  =  0._O
  x0 =  1._O
  y0 =  1._O
  z0 =  1._O
  w0 = 10._O  
  ComputeDSCS  = .true. 
  phiGS        =  0._O 
  NthetaGS     = 91 
  ExtThetaDom  = .true. 
  normalized   = .true.   
  FileDSCS     = '../OUTPUTFILES/dscs.dat'
  ComputeScatPar = .true.
  NthetaRND    = 1
  thetaminRND  = 0._O 
  thetamaxRND  = 0._O 
  Nphi = 1
  do iphi = 1, NphiMax
    phi(iphi)      = 0._O
    Ntheta(iphi)   = 1
    thetamin(iphi) = 0._O
    thetamax(iphi) = 0._O  
  end do      
  FileScat = '../OUTPUTFILES/scat.dat'
  Nelem         =  6
  MatrixElem(1) = 11
  MatrixElem(2) = 11
  MatrixElem(3) = 11
  MatrixElem(4) = 11
  MatrixElem(5) = 11
  MatrixElem(6) = 11
  do i = 7, 16
    MatrixElem(i) = 11 
  end do
  ComputeAsymPar = .false. 
  NthetaAsym = 61
  NphiAsym   = 31  
  PrnProgress    = .true.
  WriteInputInfo = .true. 
  allocate (IndI(Nelem), IndJ(Nelem), NameElem(Nelem))
  do i = 1, Nelem
    IndI(i) = 1
    IndJ(i) = 1   
    NameElem(i) = '11'  
  end do  
! ..................................................................................   
!                                start reading data                                !
! ..................................................................................
  string = 'TypePDF'
  if (XFindPar (iInputSCT, string)) then
    read (iInputSCT, *, iostat = ios) RandomOrientation
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable RandomOrientation;')"
      stop
    end if
  else
    print "(/,2x,'Group name TypePDF not found;')"
    stop  
  end if 
  if (sphere) then
    if (.not. axsym) axsym = .true.
    if (Randomorientation) RandomOrientation = .false.     
    Mrank = Nrank  ! redundant      
  end if               
! ..................................................................................  
!                         read data for a complete PDF                             !  
! ..................................................................................
  if (RandomOrientation) then
    string = 'CompletePDF - TypeMedium'
    if (XFindPar (iInputSCT, string)) then
      read (iInputSCT, *, iostat = ios) MirorSym
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable MirorSym;')"
        stop
      end if
    else
      print "(/,2x,'Group name CompletePDF - TypeMedium not found;')"
      stop  
    end if
    if (axsym .and. .not. MirorSym) then
      MirorSym = .true.  
      print "(/,2x,'Warning: the variable MirorSym has been setted to .true.')"         	
    end if            
!
    string = 'CompletePDF - AvrgMtrSS'
    if (XFindPar (iInputSCT, string)) then
      read (iInputSCT, *, iostat = ios) DoNumAvrg
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable DoNumAvrg;')"
        stop
      end if
      read (iInputSCT, *, iostat = ios) NthetaGS
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable NthetaGS;')"
        stop
      end if                                                                                                  
    else
      print "(/,2x,'Group name CompletePDF - AvrgMtrSS not found;')"
      stop  
    end if   
    call check_NthetaGS (NthetaGS)        
!
    if (.not. DoNumAvrg) then
      string = 'CompletePDF - AnalytComputAvrgMtrSS'
      if (XFindPar (iInputSCT, string)) then      
        read (iInputSCT, *, iostat = ios) ReducedOrder
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable ReducedOrder;')"
          stop
        end if
        if (ReducedOrder) then
          read (iInputSCT, *, iostat = ios) deltaOrder
          if (ios /= 0) then
            print "(/,2x,'Error by reading the input variable deltaOrder;')"
            stop
          end if                                                                                                
        end if	
      else
        print "(/,2x,'Group name CompletePDF - AnalytComputAvrgMtrSS not found;')"
        stop  
      end if      
    end if        
!
    if (DoNumAvrg) then
      string = 'CompletePDF - NumComputAvrgMtrSS'
      if (XFindPar (iInputSCT, string)) then      
        read (iInputSCT, *, iostat = ios) UseSimpson
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable UseSimpson;')"
          stop
        end if             
        read (iInputSCT, *, iostat = ios) Nalpha
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable Nalpha;')"
          stop
        end if
        read (iInputSCT, *, iostat = ios) Nbeta
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable Nbeta;')"
          stop
        end if
        read (iInputSCT, *, iostat = ios) Ngamma
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable Ngamma;')"
          stop
        end if 
        call check_NalphaNbetaNgamma (1, Nalpha)
        call check_NalphaNbetaNgamma (2, Nbeta)
        if (.not. axsym) call check_NalphaNbetaNgamma (3, Ngamma)
	call check_NSimpson (Nalpha)
        if (UseSimpson)  call check_NSimpson (Nbeta)
        if (.not. axsym) call check_NSimpson (Ngamma)      	
      else
        print "(/,2x,'Group name CompletePDF - NumComputAvrgMtrSS not found;')"
        stop  
      end if                        
    end if        
!
    string = 'CompletePDF - NormConst'
    if (XFindPar (iInputSCT, string)) then
      read (iInputSCT, *, iostat = ios) anorm
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable anorm;')"
        stop
      end if
    else
      print "(/,2x,'Group name CompletePDF - NormConst not found;')"
      stop  
    end if
!
    string  = 'CompletePDF - DSCS'
    if (XFindPar (iInputSCT, string)) then
      read (iInputSCT, *, iostat = ios) ComputeDSCS 
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable ComputeDSCS;')"
        stop
      end if
      if (ComputeDSCS) then      
        read (iInputSCT, *, iostat = ios) epol_beta
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable epol_beta;')"
          stop
        end if
        read (iInputSCT, *, iostat = ios) epol_alpha
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable epol_alpha;')"
          stop
        end if                  
        read (iInputSCT, *, iostat = ios) phiGS
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable phiGS;')"
          stop
        end if            
        read (iInputSCT, *, iostat = ios) normalized
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable normalized;')"
          stop
        end if
        read (iInputSCT, *, iostat = ios) FileDSCS
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable FileDSCS;')"
          stop
        end if
        call check_azimuthal_plane (phiGS)      
        phiGS = phiGS  * grd	
      end if	
    else
      print "(/,2x,'Group name CompletePDF - DSCS not found;')"
      stop  
    end if          
!
    string = 'CompletePDF - ScatPars'
    if (XFindPar (iInputSCT, string)) then
      read (iInputSCT, *, iostat = ios) ComputeScatPar 
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable ComputeScatPar;')"
        stop
      end if
      if (ComputeScatPar) then            
        read (iInputSCT, *, iostat = ios) NthetaRND
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable NthetaRND;')"
          stop
        end if
        read (iInputSCT, *, iostat = ios) thetaminRND
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable thetaminRND;')"
          stop
        end if
        read (iInputSCT, *, iostat = ios) thetamaxRND
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable thetamaxRND;')"
          stop
        end if                                         
        read (iInputSCT, *, iostat = ios) FileScat
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable FileScat;')"
          stop
        end if
        read (iInputSCT, *, iostat = ios) Nelem
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable Nelem;')"
          stop
        end if                    
        do i = 1, Nelem
          read (iInputSCT, *, iostat = ios) MatrixElem(i)                        
          if (ios /= 0) then
            print "(/,2x,'Error by reading the input variable MatrixElem;')"
            stop
          end if
        end do
        call check_tetaminmax (thetaminRND, thetamaxRND, NthetaRND)
        thetaminRND = thetaminRND * grd
        thetamaxRND = thetamaxRND * grd 						
        call check_MatrixElem (Nelem, MatrixElem)
	deallocate (indI, indJ, NameElem)
	allocate (IndI(Nelem), IndJ(Nelem), NameElem(Nelem))     
        call IndexElements (Nelem, MatrixElem, IndI, IndJ, NameElem)	
      end if	
    else
      print "(/,2x,'Group name CompletePDF - ScatPars not found;')"
      stop  
    end if                    
! ..................................................................................  
!                       read data for an incomplete PDF                            !  
! ..................................................................................    
  else
    string = 'IncompletePDF - TypeIntegration'
    if (XFindPar (iInputSCT, string)) then      
      read (iInputSCT, *, iostat = ios) UseSimpson
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable UseSimpson;')"
        stop
      end if      
    else
      print "(/,2x,'Group name IncompletePDF - TypeIntegration not found;')"
      stop  
    end if  
!      
    string  = 'IncompletePDF - OrientationAngles'
    if (XFindPar (iInputSCT, string)) then
      read (iInputSCT, *, iostat = ios) alphamin
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable alphamin;')"
        stop
      end if
      read (iInputSCT, *, iostat = ios) alphamax
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable alphamax;')"
        stop
      end if
      read (iInputSCT, *, iostat = ios) Nalpha
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Nalpha;')"
        stop
      end if
      read (iInputSCT, *, iostat = ios) betamin
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable betamin;')"
        stop
      end if
      read (iInputSCT, *, iostat = ios) betamax
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable betamax;')"
        stop
      end if
      read (iInputSCT, *, iostat = ios) Nbeta
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Nbeta;')"
        stop
      end if
      read (iInputSCT, *, iostat = ios) gammamin
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable gammamin;')"
        stop
      end if
      read (iInputSCT, *, iostat = ios) gammamax
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable gammamax;')"
        stop
      end if
      read (iInputSCT, *, iostat = ios) Ngamma
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Ngamma;')"
        stop
      end if
      if (sphere) then    
        alphamin = 0._O
        alphamax = 0._O
        Nalpha   = 1
        betamin  = 0._O
        betamax  = 0._O
        Nbeta    = 1  
        gammamin = 0._O
        gammamax = 0._O
        Ngamma   = 1  
      end if
      call check_alfaminmax (alphamin, alphamax, Nalpha)
      call check_betaminmax (betamin, betamax, Nbeta)     
      if (.not. axsym) then
        call check_gamaminmax (gammamin, gammamax, Ngamma)              
      else
        gammamin = 0._O
        gammamax = 0._O
        Ngamma   = 1 
      end if
      alphamin = alphamin * grd
      alphamax = alphamax * grd
      betamin  = betamin  * grd
      betamax  = betamax  * grd
      gammamin = gammamin * grd
      gammamax = gammamax * grd      
      call check_NSimpson (Nalpha)
      if (UseSimpson)  call check_NSimpson (Nbeta)
      if (.not. axsym) call check_NSimpson (Ngamma)                
    else
      print "(/,2x,'Group name IncompletePDF - OrientationAngles not found;')"
      stop  
    end if  
!              
    string = 'IncompletePDF - NormConst'
    if (XFindPar (iInputSCT, string)) then
      read (iInputSCT, *, iostat = ios) anorm
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable anorm;')"
        stop
      end if
    else
      print "(/,2x,'Group name IncompletePDF - NormConst not found;')"
      stop  
    end if   
!              
    string = 'IncompletePDF - IncWave'
    if (XFindPar (iInputSCT, string)) then
      read (iInputSCT, *, iostat = ios) TypeExcit
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable TypeExcit;')"
        stop
      end if
      read (iInputSCT, *, iostat = ios) thetaGI
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable thetaGI;')"
        stop
      end if
      read (iInputSCT, *, iostat = ios) phiGI
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable phiGI;')"
        stop
      end if
      if (TypeExcit(1:5) == 'GAUSS') then       
        read (iInputSCT, *, iostat = ios) x0
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable x0;')"
          stop
        end if
        read (iInputSCT, *, iostat = ios) y0
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable y0;')"
          stop
        end if
        read (iInputSCT, *, iostat = ios) z0
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable z0;')"
          stop
        end if
        read (iInputSCT, *, iostat = ios) w0
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable w0;')"
          stop
        end if 
      end if	
      call check_TypeExcit (TypeExcit) 
      call check_incident_direction (thetaGI, phiGI)  
      thetaGI = thetaGI * grd
      phiGI   = phiGI   * grd                      
    else
      print "(/,2x,'Group name IncompletePDF - IncWave not found;')"
      stop  
    end if                           
!
    string = 'IncompletePDF - DSCS'
    if (XFindPar (iInputSCT, string)) then
      read (iInputSCT, *, iostat = ios) ComputeDSCS 
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable ComputeDSCS;')"
        stop
      end if                                                            
      if (ComputeDSCS) then   	
        read (iInputSCT, *, iostat = ios) epol_beta
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable epol_beta;')"
          stop
        end if
        read (iInputSCT, *, iostat = ios) epol_alpha
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable epol_alpha;')"
          stop
        end if                    
        read (iInputSCT, *, iostat = ios) alphapGauss
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable alphapGauss;')"
          stop
        end if        	        
	read (iInputSCT, *, iostat = ios) phiGS
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable phiGS;')"
          stop
        end if
        read (iInputSCT, *, iostat = ios) NthetaGS
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable NthetaGS;')"
          stop
        end if
        read (iInputSCT, *, iostat = ios) ExtThetaDom
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable ExtThetaDom;')"
          stop
        end if
        read (iInputSCT, *, iostat = ios) normalized
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable normalized;')"
          stop
        end if
        read (iInputSCT, *, iostat = ios) FileDSCS
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable FileDSCS;')"
          stop
        end if									
	call check_azimuthal_plane (phiGS)
	phiGS = phiGS  * grd        		
        call check_polarization_angle (alphapGauss)        
        alphapGauss = alphapGauss * grd			
      end if                                                	
    else
      print "(/,2x,'Group name IncompletePDF - DSCS not found;')"
      stop  
    end if   
!    
    string = 'IncompletePDF - ScatPars'
    if (XFindPar (iInputSCT, string)) then
      read (iInputSCT, *, iostat = ios) ComputeScatPar 
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable ComputeScatPar;')"
        stop
      end if            
      if (ComputeScatPar) then              
        read (iInputSCT, *, iostat = ios) Nphi
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable Nphi;')"
          stop
        end if
        do iphi = 1, Nphi
          read (iInputSCT, *, iostat = ios) phi(iphi)
          if (ios /= 0) then
            print "(/,2x,'Error by reading the input variable phi;')"
            stop
          end if      
          read (iInputSCT, *, iostat = ios) Ntheta(iphi)
          if (ios /= 0) then
            print "(/,2x,'Error by reading the input variable Ntheta;')"
            stop
          end if      
          read (iInputSCT, *, iostat = ios) thetamin(iphi)
          if (ios /= 0) then
            print "(/,2x,'Error by reading the input variable thetamin;')"
            stop
          end if        
          read (iInputSCT, *, iostat = ios) thetamax(iphi)
          if (ios /= 0) then
            print "(/,2x,'Error by reading the input variable thetamax;')"
            stop
          end if                
        end do                    
        read (iInputSCT, *, iostat = ios) FileScat
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable FileScat;')"
          stop
        end if
        read (iInputSCT, *, iostat = ios) Nelem
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable Nelem;')"
          stop
        end if                    
        do i = 1, Nelem
          read (iInputSCT, *, iostat = ios) MatrixElem(i)                        
          if (ios /= 0) then
            print "(/,2x,'Error by reading the input variable MatrixElem;')"
            stop
          end if
        end do        	
	call check_Nphi (Nphi)
        call check_teta_phiminmax (.true., Nphi, phi, Ntheta, thetamin, thetamax)
        do iphi = 1, Nphi    
          phi(iphi)      = phi(iphi)      * grd
          thetamin(iphi) = thetamin(iphi) * grd
          thetamax(iphi) = thetamax(iphi) * grd      
        end do
        call check_MatrixElem (Nelem, MatrixElem)  
	deallocate (indI, indJ, NameElem)
	allocate (IndI(Nelem), IndJ(Nelem), NameElem(Nelem))    
        call IndexElements (Nelem, MatrixElem, IndI, IndJ, NameElem)	
      end if      	
    else
      print "(/,2x,'Group name IncompletePDF - ScatPars not found;')"
      stop  
    end if
!           
    string = 'IncompletePDF - MeanDirPropagatScatWave'
    if (XFindPar (iInputSCT, string)) then
      read (iInputSCT, *, iostat = ios) ComputeAsymPar
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable ComputeAsymPar;')"
        stop
      end if
      if (ComputeAsymPar .and. .not. ComputeScatPar) then
        ComputeAsymPar = .false.         
        print                                                                       &
       "(/,2x,'Warning: the variable ComputeAsymPar has been setted to .false.')"	
      end if
      if (ComputeAsymPar) then
        read (iInputSCT, *, iostat = ios) NthetaAsym
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable NthetaAsym;')"
          stop
        end if
        read (iInputSCT, *, iostat = ios) NphiAsym
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable NphiAsym;')"
          stop
        end if
        call check_NSimpson (NphiAsym)
        if (UseSimpson)  call check_NSimpson (NthetaAsym) 	
      end if   	
    else
      print "(/,2x,'Group name IncompletePDF - MeanDirPropagatScatWave not found;')"
      stop  
    end if                      
  end if
  call check_anorm (anorm)
  xpart = k * anorm
  snorm = Pi * xpart * xpart 
!
  string = 'PrintInfo'
  if (XFindPar (iInputSCT, string)) then
    read (iInputSCT, *, iostat = ios) PrnProgress
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable PrnProgress;')"
      stop
    end if
    read (iInputSCT, *, iostat = ios) WriteInputInfo
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable WriteInputInfo;')"
      stop
    end if
  else
    print "(/,2x,'Group name PrintInfo not found;')"
    stop  
  end if 
  close (unit = iInputSCT) 
end subroutine readinputSCT  
! **********************************************************************************
subroutine DSCS_SCAT (ComputeDSCS, ComputeScatPar, ComputeAsymPar, axsym, sphere,   &
           chiral, Mrank, Nrank, NthetaGS, phiAZIMUT, Nphi, phi, Ntheta, thetamin,  &
           thetamax, NthetaAsym, NphiAsym, UseSimpson, xphi, wphi, xtheta, wtheta,  &
           thetaGI, phiGI, alpha, beta, gamma, epol_beta, epol_alpha, alphapGauss,  &
           x0, y0, z0, w0, TypeExcit, wavenumber, snorm, FileTmat, ExtThetaDom,     &
           normalized, Nelem, IndI, IndJ, h, v, ZTP, NphiAL, NthetaAL, Ke, Cscat,   &
           Cext, CscatX, CscatY, CextX, CextY, gX_x, gX_y, gX_z, gY_x, gY_y, gY_z )                       
! -----------------------------------------------------------------------------------
! The routine computes:                                                             !
! - the optical cross sections for incident parallel and perpendicular linear       !
!   polarizations,                                                                  !
! - the mean direction of propagation of the scattered radiation for incident       !
!   parallel and perpendicular linear polarizations,                                !
! - the extinction matrix for a plane wave incidence, and                           !   
! - the phase matrix at Nphi scattering planes.                                     !
! Note that the optical cross sections and the mean direction of propagation of the !
! scattered field are computed for linearly polarized waves (vector plane waves and ! 
! Gaussian beams propagating in the direction (thetaGI,phiGI)). The calculations    !
! are performed for incident parallel and perpendicular polarizations, that is for  !
! alphap = 0° and alphap = 90°. The code also determines the differential           !
! scattering cross sections at a set of NthetaGS scattering angles in the azimuthal !
! plane phiAZIMUT. The calculations are performed for elliptically polarized plane  !
! waves and linearly polarized Gaussian beams. For a plane wave incidence, the      !
! differential scattering cross sections correspond to the complex polarization     !
! unit vector                                                                       !      
!                 epol = epol_beta * e_beta + epol_alpha * e_alpha.                 !
! In the case of a Gaussian beam illumination, the calculations are performed for a !
! linearly polarized wave and the state of polarization of the Gaussian beam is     !
! described by the polarization angle alphapGauss. Note that all scattering         !
! characteristics are computed for a fixed particle orientation characterized by    !
! the Euler rotation angles alpha, beta and gamma.                                  !
! -----------------------------------------------------------------------------------
  use parameters
  implicit none      
  integer       :: Mrank, Nrank, NthetaGS, Ntheta(NphiMax), Nphi, Nelem,            &
                   IndI(Nelem), IndJ(Nelem), NthetaAsym, NphiAsym, NphiAL, NthetaAL
  real(O)       :: phiAZIMUT, thetamin(NphiMax), thetamax(NphiMax), phi(NphiMax),   &
                   thetaGI, phiGI, alpha, beta, gamma, alphapGauss, x0, y0, z0, w0, &
                   wavenumber, snorm, Cscat, Cext, CscatX, CscatY, CextX, CextY,    &
                   Ke(4,4), h(NthetaGS), v(NthetaGS), ZTP(NphiAL,NthetaAL,Nelem),   &
                   gX_x, gX_y, gX_z, gY_x, gY_y, gY_z, xphi(NphiAsym),              &
                   wphi(NphiAsym), xtheta(NthetaAsym), wtheta(NthetaAsym)		                   
  complex(O)    :: epol_beta, epol_alpha		    
  character(5)  :: TypeExcit
  character(80) :: FileTmat
  logical       :: ComputeDSCS, ComputeScatPar, ComputeAsymPar, sphere, axsym,      &
                   chiral, UseSimpson, ExtThetaDom, normalized
!     
  integer       :: Mstart, m, Nmax, Nmaxmax, ntg, mtg, ntl, mtl, itheta, iphi,      &
                   i, j, kelem  
  real(O)       :: alphapX, alphapY, QscatX, QscatY, QextX, QextY, thetaGSL,        &
                   phiGSL, thetaGS, phiGS, dtheta, Z(4,4), fact, Qscat, Qext, norm   
  complex(O)    :: FthetaX, FphiX, FthetaY, FphiY, S(2,2), epol1_beta, epol1_alpha
  complex(O),allocatable :: T(:,:), tvg(:), tv(:), cX(:), cY(:), c1X(:), c1Y(:),    &
                            ccX(:), ccY(:), cc(:)
!
  alphapX  = 0._O
  alphapY  = Pi / 2._O 
  Nmaxmax  = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  allocate (ccX(2*Nmaxmax), ccY(2*Nmaxmax)) 
  open (unit = iTmat, file = FileTmat, status = 'unknown')
  if (.not. axsym) then
    call read_HeadFileTmat (ntg, mtg)
    call check_dimensionMat (ntg, mtg, Nmaxmax)           
    allocate (T(2*ntg,2*mtg), cX(2*Nmaxmax), cY(2*Nmaxmax))    
    call read_FileTmat (ntg, mtg, T) 
    if (TypeExcit(1:5) == 'PLANE') then                   
      call PWcoefficients_ab (thetaGI, phiGI, alpha, beta, gamma, alphapX, Mrank,   &
           Nrank, Nmaxmax, cX)
      call PWcoefficients_ab (thetaGI, phiGI, alpha, beta, gamma, alphapY, Mrank,   &
           Nrank, Nmaxmax, cY)             	   	   	   	   
    else if (TypeExcit(1:5) == 'GAUSS') then                                             
      call GBcoefficients_ab (wavenumber, x0, y0, z0, w0, thetaGI, phiGI,           &
           alpha, beta, gamma, alphapX, Mrank, Nrank, Nmaxmax, cX) 
      call GBcoefficients_ab (wavenumber, x0, y0, z0, w0, thetaGI, phiGI,           &
           alpha, beta, gamma, alphapY, Mrank, Nrank, Nmaxmax, cY)                         
    end if                                                                              
    call product_matrix_vector (2*Nmaxmax, 2*Nmaxmax, T, 2*ntg, 2*mtg, cX, ccX)
    call product_matrix_vector (2*Nmaxmax, 2*Nmaxmax, T, 2*ntg, 2*mtg, cY, ccY)         
    deallocate (T, cX, cY)
  else 
    gamma = 0._O ! redundant
    if (.not. sphere) then      
      call read_HeadFileTmat (ntl, mtl)
      call check_dimensionMat (ntl, mtl, Nrank)             
      allocate (T(2*ntl, 2*mtl))
    else 
      call read_HeadFileTmatVct (ntl)
      call check_dimensionVct (ntl, Nrank)
      allocate (tvg(2*ntl), tv(2*ntl))
      call read_FileTmatVct (ntl, tvg)       
    end if
    allocate (cX(2*Nrank), cY(2*Nrank), c1X(2*Nrank), c1Y(2*Nrank))                                                          
    Mstart = 0
    do m = Mstart, Mrank
      if (m == 0) then
        Nmax = Nrank
      else
        Nmax = Nrank - m + 1
      end if
      if (.not. sphere) then        
        call read_FileTmat (ntl, mtl, T)
      else
        call form_Tvector (Nrank, Nmax, m, tvg, tv)
      end if
      if (TypeExcit(1:5) == 'PLANE') then             
        call PWcoefficients_ab_m (thetaGI, phiGI, alpha, beta, gamma, alphapX, m,   &
             Nrank, Nmax, cX)
        call PWcoefficients_ab_m (thetaGI, phiGI, alpha, beta, gamma, alphapY, m,   &
             Nrank, Nmax, cY)        	     	     	     	     
      else if (TypeExcit(1:5) == 'GAUSS') then   
        call GBcoefficients_ab_m (wavenumber, x0, y0, z0, w0, thetaGI, phiGI,       &
             alpha, beta, gamma, alphapX, m, Nmax, cX) 
        call GBcoefficients_ab_m (wavenumber, x0, y0, z0, w0, thetaGI, phiGI,       &
             alpha, beta, gamma, alphapY, m, Nmax, cY)                                                                                 
      end if                                       
      if (.not. sphere) then        
        call product_matrix_vector (2*Nmax, 2*Nmax, T, 2*ntl, 2*mtl, cX, c1X)
        call product_matrix_vector (2*Nmax, 2*Nmax, T, 2*ntl, 2*mtl, cY, c1Y)		  
      else 
        call product_vector_vector (2*Nmax, tv, cX, c1X)
        call product_vector_vector (2*Nmax, tv, cY, c1Y)		  
      end if                     
      call extend_vector_positive (c1X, ccX, m, Mstart, Nrank, Nmax, Nmaxmax)
      call extend_vector_positive (c1Y, ccY, m, Mstart, Nrank, Nmax, Nmaxmax)      
      if (m /= 0) then
        if (.not. sphere) then
          if (.not. chiral) then
            call matrix_m_negativ (Nmax, Nmax, T, ntl, mtl)
          else              
            call read_FileTmat (ntl, mtl, T)        
          end if    
        end if  
        if (TypeExcit(1:5) == 'PLANE') then      
          call PWcoefficients_ab_m (thetaGI, phiGI, alpha, beta, gamma, alphapX,    &
               - m, Nrank, Nmax, cX)
          call PWcoefficients_ab_m (thetaGI, phiGI, alpha, beta, gamma, alphapY,    &
               - m, Nrank, Nmax, cY)          	       	       	       	       	       	               
        else if (TypeExcit(1:5) == 'GAUSS') then         
          call GBcoefficients_ab_m (wavenumber, x0, y0, z0, w0, thetaGI, phiGI,     &
               alpha, beta, gamma, alphapX, - m, Nmax, cX) 
          call GBcoefficients_ab_m (wavenumber, x0, y0, z0, w0, thetaGI, phiGI,     &
               alpha, beta, gamma, alphapY, - m, Nmax, cY)
        end if 
        if (.not. sphere) then        
          call product_matrix_vector (2*Nmax, 2*Nmax, T, 2*ntl, 2*mtl, cX, c1X)
          call product_matrix_vector (2*Nmax, 2*Nmax, T, 2*ntl, 2*mtl, cY, c1Y)           	        
        else 
          call product_vector_vector (2*Nmax, tv, cX, c1X)
          call product_vector_vector (2*Nmax, tv, cY, c1Y)          	  
        end if                      
        call extend_vector_negative (c1X, ccX, m, Nrank, Nmax, Nmaxmax)
        call extend_vector_negative (c1Y, ccY, m, Nrank, Nmax, Nmaxmax)	
      end if 
    end do 
    if (.not. sphere) then
      deallocate (T)
    else 
      deallocate (tvg, tv)      
    end if                                   
    deallocate (cX, cY, c1X, c1Y)
  end if
  close (unit = iTmat)   
! --- DSCS and Cross Sections for plane waves and Gaussian beams ---
  Cscat = 0._O
  Cext  = 0._O
  do itheta = 1, NthetaGS
    h(itheta) = 0._O
    v(itheta) = 0._O
  end do
  if (ComputeDSCS) then        
    norm        =   sqrt ( abs(epol_beta)**2 + abs(epol_alpha)**2 )  
    epol1_beta  =   epol_beta  / norm
    epol1_alpha =   epol_alpha / norm  
    allocate (cc(2*Nmaxmax))
    if (TypeExcit(1:5) == 'PLANE') then                                                                        
      do i = 1, 2*Nmaxmax
        cc(i) = ccX(i) * epol1_beta + ccY(i) * epol1_alpha        
      end do    
    else if (TypeExcit(1:5) == 'GAUSS') then
      do i = 1, 2*Nmaxmax
        cc(i) = ccX(i) * cos(alphapGauss) + ccY(i) * sin(alphapGauss)        
      end do  
    end if             
    call DSCS (cc, Mrank, Nrank, Nmaxmax, NthetaGS, phiAZIMUT, alpha, beta, gamma,  &
         wavenumber, snorm, ExtThetaDom, normalized, h, v)
    call CQscat (cc, Mrank, Nrank, Nmaxmax, wavenumber, snorm, Cscat, Qscat)
    call CQextPolGen (cc, Mrank, Nrank, Nmaxmax, thetaGI, phiGI, alpha, beta,       &
         gamma, epol1_beta, epol1_alpha, alphapGauss, x0, y0, z0, w0, TypeExcit,    &
         wavenumber, snorm, Cext, Qext)            	      
    deallocate (cc)	 
  end if
! --- Scattering Characteristics for linearly polarized waves ---
  do i = 1, 4
    do j = 1, 4
      Ke(i,j) = 0._O
    end do
  end do
  do iphi = 1, NphiAL     
    do itheta = 1, NthetaAL               
      do kelem = 1, Nelem                           
        ZTP(iphi,itheta,kelem) = 0._O
      end do
    end do
  end do
  CextX  = 0._O
  CextY  = 0._O
  CscatX = 0._O
  CscatY = 0._O                    
  gX_x   = 0._O
  gX_y   = 0._O
  gX_z   = 0._O
  gY_x   = 0._O
  gY_y   = 0._O
  gY_z   = 0._O  
  if (ComputeScatPar) then 
!   --- extinction matrix --- 
    if (TypeExcit(1:5) == 'PLANE') then       
      thetaGSL = thetaGI
      phiGSL   = phiGI
      call F_azimuthal_tetaGS (ccX, Mrank, Nrank, Nmaxmax, thetaGSL, phiGSL, alpha, &
           beta, gamma, wavenumber, FthetaX, FphiX)
      call F_azimuthal_tetaGS (ccY, Mrank, Nrank, Nmaxmax, thetaGSL, phiGSL, alpha, &
           beta, gamma, wavenumber, FthetaY, FphiY)
      S(1,1) = FthetaX
      S(2,1) = FphiX
      S(1,2) = FthetaY
      S(2,2) = FphiY 
      call ExtinctionMatrix (wavenumber, S, Ke)        
    end if       
!   --- phase matrix ---
    do iphi = 1, Nphi
      phiGS = phi(iphi)
      if (Ntheta(iphi) /= 1) then
        dtheta = (thetamax(iphi) - thetamin(iphi)) / (Ntheta(iphi) - 1)
      else
        dtheta = 0._O
      end if        
      do itheta = 1, Ntheta(iphi)
        thetaGS = thetamin(iphi) + (itheta - 1) * dtheta      
        call F_azimuthal_tetaGS (ccX, Mrank, Nrank, Nmaxmax, thetaGS, phiGS,        &
             alpha, beta, gamma, wavenumber, FthetaX, FphiX)
        call F_azimuthal_tetaGS (ccY, Mrank, Nrank, Nmaxmax, thetaGS, phiGS,        &
             alpha, beta, gamma, wavenumber, FthetaY, FphiY)               
        S(1,1) = FthetaX
        S(2,1) = FphiX
        S(1,2) = FthetaY
        S(2,2) = FphiY     
        call PhaseMatrix (S, Z)   
        do kelem = 1, Nelem
          i = IndI(kelem)
          j = IndJ(kelem)       
          ZTP(iphi,itheta,kelem) = Z(i,j)
        end do
      end do
    end do                    
!   --- scattering and extinction cross sections for linearly polarized waves ---              
    call CQscat (ccX, Mrank, Nrank, Nmaxmax, wavenumber, snorm, CscatX, QscatX)
    call CQscat (ccY, Mrank, Nrank, Nmaxmax, wavenumber, snorm, CscatY, QscatY)  
    call CQextGen (ccX, Mrank, Nrank, Nmaxmax, thetaGI, phiGI, alpha, beta, gamma,  &
         alphapX, x0, y0, z0, w0, TypeExcit, wavenumber, snorm, CextX, QextX) 
    call CQextGen (ccY, Mrank, Nrank, Nmaxmax, thetaGI, phiGI, alpha, beta, gamma,  &
         alphapY, x0, y0, z0, w0, TypeExcit, wavenumber, snorm, CextY, QextY)     
!   --- mean direction of propagation of the scattered radiation ---    
    if (ComputeAsymPar) then        
      do itheta = 1, NthetaAsym
        if (.not. UseSimpson) then 
          thetaGS = xtheta(itheta)
        else
          thetaGS = acos(xtheta(itheta))
        end if
        do iphi = 1, NphiAsym
          phiGS = xphi(iphi)
          fact  = wtheta(itheta) * wphi(iphi)
          call F_azimuthal_tetaGS (ccX, Mrank, Nrank, Nmaxmax, thetaGS, phiGS,      &
               alpha, beta, gamma, wavenumber, FthetaX, FphiX) 
          call F_azimuthal_tetaGS (ccY, Mrank, Nrank, Nmaxmax, thetaGS, phiGS,      &
               alpha, beta, gamma, wavenumber, FthetaY, FphiY)
          gX_x = gX_x + (abs(FthetaX)**2 + abs(FphiX)**2) * sin(thetaGS) *          &
                         cos(phiGS) * fact 
          gX_y = gX_y + (abs(FthetaX)**2 + abs(FphiX)**2) * sin(thetaGS) *          &
                         sin(phiGS) * fact
          gX_z = gX_z + (abs(FthetaX)**2 + abs(FphiX)**2) * cos(thetaGS) * fact    
          gY_x = gY_x + (abs(FthetaY)**2 + abs(FphiY)**2) * sin(thetaGS) *          &
                         cos(phiGS) * fact 
          gY_y = gY_y + (abs(FthetaY)**2 + abs(FphiY)**2) * sin(thetaGS) *          &
                         sin(phiGS) * fact
          gY_z = gY_z + (abs(FthetaY)**2 + abs(FphiY)**2) * cos(thetaGS) * fact
        end do
      end do 
      gX_x = gX_x / CscatX
      gX_y = gX_y / CscatX
      gX_z = gX_z / CscatX
      gY_x = gY_x / CscatY
      gY_y = gY_y / CscatY
      gY_z = gY_z / CscatY                        	      
    end if
  end if
  deallocate (ccX, ccY)       
end subroutine DSCS_SCAT
! **********************************************************************************
subroutine DSCS (c, Mrank, Nrank, Nmax, Ntheta, phiAZIMUT, alpha, beta, gamma,      &
           wavenumber, snorm, ExtThetaDom, normalized, h, v)
!------------------------------------------------------------------------------------       
! Let                                                                               !
!           ES = [exp(jkR) / (kR)] * ESinf = [exp(jkR) / R] * FS,                   !
! be the scattered field in the direction (theta,phi), where FS is the far-field    !
! pattern, ESinf = k * FS and k is the wave number. If we consider the decomposition!
!           ESinf = ESinf_theta * e_theta + ESinf_phi * e_phi,                      !
! the routine computes the differential scattering cross sections for parallel and  !
! perpendicular polarizations                                                       !
!           h = |ESinf_theta|**2 / (k**2) = |FS_theta|**2  and                      !
!           v = |ESinf_phi|**2   / (k**2) = |FS_phi|**2,                            !
! and the normalized differential scattering cross sections for parallel and        !
! perpendicular polarizations                                                       !
!           hn = |ES_theta|**2 / [Pi * (k*a)**2] = |FS_theta|**2 / (Pi * a**2) and  !
!           vn = |ES_phi|**2   / [Pi * (k*a)**2] = |FS_phi|**2   / (Pi * a**2),     !
! where a is the characteristic dimension of the scatterer. Thus,                   !
!           hn = h / (Pi * a**2)   and  vn = v / (Pi * a**2).                       !
! Note that the normalization constant snorm is snorm = Pi * (k*a)**2. The          !
! quantities are computed at a set of Ntheta scattering angles in the azimuthal     !
! plane phiAZIMUT.                                                                  !                           
!------------------------------------------------------------------------------------
  use parameters
  implicit none
  integer      :: Mrank, Nrank, Nmax, Ntheta
  real(O)      :: phiAZIMUT, alpha, beta, gamma, wavenumber, snorm,                 &
                  h(Ntheta), v(Ntheta)
  complex(O)   :: c(2*Nmax)
  logical      :: ExtThetaDom, normalized
!
  integer      :: i, m, k, N0, l
  real(O)      :: fact, thetaGS, phiGS, thetaLS, phiLS, cosu, sinu
  complex(O)   :: sum(3)
  complex(O),allocatable :: Minf(:,:), Ninf(:,:)
!
  if (normalized) then
    fact = snorm
  else 
    fact = wavenumber * wavenumber
  end if           
  allocate (Minf(3,Nmax), Ninf(3,Nmax))
  do i = 1, Ntheta
    if (ExtThetaDom) then      
      thetaGS = real(i - 1,O) * 2._O * Pi / real(Ntheta - 1,O)
      phiGS   = phiAZIMUT
      if (thetaGS > Pi) then
        thetaGS = 2._O * Pi - thetaGS
        phiGS   = phiAZIMUT + Pi
      end if
    else 
      thetaGS = real(i - 1,O) * Pi / real(Ntheta - 1,O)
      phiGS   = phiAZIMUT    
    end if  
    call T_spherical_global_local (thetaGS, phiGS, alpha, beta, gamma, thetaLS, phiLS)         
    call MN_infinit_complete (thetaLS, phiLS, Mrank, Nrank, Nmax, .true., Minf, Ninf)
    sum(1) = zero
    sum(2) = zero
    sum(3) = zero
    do m = 0, Mrank
      if (m == 0) then
        do k = 1, Nrank
          sum(1) = sum(1) + (Minf(1,k) * c(k) + Ninf(1,k) * c(Nmax+k))
          sum(2) = sum(2) + (Minf(2,k) * c(k) + Ninf(2,k) * c(Nmax+k))
          sum(3) = sum(3) + (Minf(3,k) * c(k) + Ninf(3,k) * c(Nmax+k))
        end do
      else                    
        N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
        do l = 1, 2
          do k = 1, Nrank - m + 1
            sum(1) = sum(1) + (Minf(1,N0+k) * c(N0+k) + Ninf(1,N0+k) * c(Nmax+N0+k))
            sum(2) = sum(2) + (Minf(2,N0+k) * c(N0+k) + Ninf(2,N0+k) * c(Nmax+N0+k))
            sum(3) = sum(3) + (Minf(3,N0+k) * c(N0+k) + Ninf(3,N0+k) * c(Nmax+N0+k))
          end do
          N0 = N0 + Nrank - m + 1
        end do                        
      end if
    end do
    call angle_unitvct_ItL_ItG (thetaGS, phiGS, thetaLS, phiLS, alpha, beta, gamma, &
         cosu, sinu)        
    h(i) = (abs(cosu * sum(2) - sinu * sum(3)))**2 / fact
    v(i) = (abs(sinu * sum(2) + cosu * sum(3)))**2 / fact
  end do
  deallocate (Minf, Ninf)      
end subroutine DSCS
! **********************************************************************************
subroutine F_azimuthal (c, Mrank, Nrank, Nmax, Ntheta, phiAZIMUT, alpha, beta,      &
           gamma, wavenumber, ExtThetaDom, Ftheta, Fphi)
!------------------------------------------------------------------------------------ 
! Let                                                                               !
!           ES = [exp(jkR) / (kR)] * ESinf = [exp(jkR) / R] * FS,                   !
! be the scattered field in the direction (theta,phi), where FS is the far-field    !
! pattern, ESinf = k * FS and k is the wave number. If we consider the decomposition!
!           ESinf = ESinf_theta * e_theta + ESinf_phi * e_phi,                      !
! the routine computes                                                              ! 
!           FS_theta = ESinf_theta / k, and                                         !
!           FS_phi   = ESinf_phi   / k                                              !
! at Ntheta scattering angles in the azimuthal plane phiAZIMUT.                     !
!------------------------------------------------------------------------------------                  
  use parameters
  implicit none
  integer      :: Mrank, Nrank, Nmax, Ntheta
  real(O)      :: phiAZIMUT, alpha, beta, gamma, wavenumber             
  complex(O)   :: c(2*Nmax), Ftheta(Ntheta), Fphi(Ntheta)
  logical      :: ExtThetaDom
!
  integer      :: i, m, k, N0, l
  real(O)      :: thetaGS, phiGS, thetaLS, phiLS, cosu, sinu
  complex(O)   :: sum(3)
  complex(O),allocatable :: Minf(:,:), Ninf(:,:)
!
  allocate (Minf(3,Nmax), Ninf(3,Nmax))
  do i = 1, Ntheta
    if (ExtThetaDom) then         
      thetaGS = real(i - 1,O) * 2._O * Pi / real(Ntheta - 1,O)
      phiGS  = phiAZIMUT
      if (thetaGS > Pi) then
        thetaGS = 2._O * Pi - thetaGS
        phiGS  = phiAZIMUT + Pi
      end if
    else 
      thetaGS = real(i - 1,O) * Pi / real(Ntheta - 1,O)
      phiGS  = phiAZIMUT    
    end if  
    call T_spherical_global_local (thetaGS, phiGS, alpha, beta, gamma, thetaLS, phiLS)            
    call MN_infinit_complete (thetaLS, phiLS, Mrank, Nrank, Nmax, .true., Minf, Ninf)
    sum(1) = zero
    sum(2) = zero
    sum(3) = zero
    do m = 0, Mrank
      if (m == 0) then
        do k = 1, Nrank
          sum(1) = sum(1) + (Minf(1,k) * c(k) + Ninf(1,k) * c(Nmax+k))
          sum(2) = sum(2) + (Minf(2,k) * c(k) + Ninf(2,k) * c(Nmax+k))
          sum(3) = sum(3) + (Minf(3,k) * c(k) + Ninf(3,k) * c(Nmax+k))
        end do
      else                    
        N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
        do l = 1, 2
          do k = 1, Nrank - m + 1
            sum(1) = sum(1) + (Minf(1,N0+k) * c(N0+k) + Ninf(1,N0+k) * c(Nmax+N0+k))
            sum(2) = sum(2) + (Minf(2,N0+k) * c(N0+k) + Ninf(2,N0+k) * c(Nmax+N0+k))
            sum(3) = sum(3) + (Minf(3,N0+k) * c(N0+k) + Ninf(3,N0+k) * c(Nmax+N0+k))
          end do
          N0 = N0 + Nrank - m + 1
        end do                    
      end if
    end do
    call angle_unitvct_ItL_ItG (thetaGS, phiGS, thetaLS, phiLS, alpha, beta, gamma, &
         cosu, sinu)
    Ftheta(i) = (cosu * sum(2) - sinu * sum(3)) / wavenumber
    Fphi(i)  = (sinu * sum(2) + cosu * sum(3)) / wavenumber                              
  end do
  deallocate (Minf, Ninf)      
end subroutine F_azimuthal      
! **********************************************************************************
subroutine F_azimuthal_tetaGS (c, Mrank, Nrank, Nmax, thetaGS, phiGS, alpha, beta,  &
           gamma, wavenumber, Ftheta, Fphi)
!------------------------------------------------------------------------------------
! Let                                                                               !
!           ES = [exp(jkR) / (kR)] * ESinf = [exp(jkR) / R] * FS,                   !
! be the scattered field in the direction (theta,phi), where FS is the far-field    !
! pattern, ESinf = k * FS and k is the wave number. If we consider the decomposition!
!           ESinf = ESinf_theta * e_theta + ESinf_phi * e_phi,                      !
! the routine computes                                                              !  
!           FS_theta = ESinf_theta / k, and                                         !
!           FS_phi   = ESinf_phi   / k                                              !
! at the scattering direction (thetaGS, phiGS).                                     !     
!------------------------------------------------------------------------------------  
  use parameters
  implicit none
  integer      :: Mrank, Nrank, Nmax
  real(O)      :: thetaGS, phiGS, alpha, beta, gamma, wavenumber
  complex(O)   :: c(2*Nmax), Ftheta, Fphi
!
  integer      :: m, k, N0, l
  real(O)      :: thetaLS, phiLS, cosu, sinu
  complex(O)   :: sum(3)
  complex(O),allocatable :: Minf(:,:), Ninf(:,:)
!
  allocate (Minf(3,Nmax), Ninf(3,Nmax))            
  call T_spherical_global_local (thetaGS, phiGS, alpha, beta, gamma, thetaLS, phiLS)            
  call MN_infinit_complete (thetaLS, phiLS, Mrank, Nrank, Nmax, .true., Minf, Ninf)           
  sum(1) = zero
  sum(2) = zero
  sum(3) = zero
  do m = 0, Mrank
    if (m == 0) then
      do k = 1, Nrank
        sum(1) = sum(1) + (Minf(1,k) * c(k) + Ninf(1,k) * c(Nmax+k))
        sum(2) = sum(2) + (Minf(2,k) * c(k) + Ninf(2,k) * c(Nmax+k))
        sum(3) = sum(3) + (Minf(3,k) * c(k) + Ninf(3,k) * c(Nmax+k))
      end do
    else        
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      do l = 1, 2
        do k = 1, Nrank - m + 1
          sum(1) = sum(1) + (Minf(1,N0+k) * c(N0+k) + Ninf(1,N0+k) * c(Nmax+N0+k))
          sum(2) = sum(2) + (Minf(2,N0+k) * c(N0+k) + Ninf(2,N0+k) * c(Nmax+N0+k))
          sum(3) = sum(3) + (Minf(3,N0+k) * c(N0+k) + Ninf(3,N0+k) * c(Nmax+N0+k))
        end do
        N0 = N0 + Nrank - m + 1
      end do               
    end if
  end do
  call angle_unitvct_ItL_ItG (thetaGS, phiGS, thetaLS, phiLS, alpha, beta, gamma,   &
       cosu, sinu)   
  Ftheta = (cosu * sum(2) - sinu * sum(3)) / wavenumber
  Fphi  = (sinu * sum(2) + cosu * sum(3)) / wavenumber            
  deallocate (Minf, Ninf)      
end subroutine F_azimuthal_tetaGS       
! **********************************************************************************
subroutine CQscat (c, Mrank, Nrank, Nmax, wavenumber, snorm, Cscat, Qscat)
!------------------------------------------------------------------------------------
! The routine computes the scattering cross section                                 !
!        C_scat = (Pi / k**2) * SUM { |Fmn|**2 + |Gmn|**2 }                         !
! end the scattering efficiency                                                     !
!        Q_scat = C_scat / (Pi * a**2)                                              !
!               = [1 / (ka)**2] * SUM { |Fmn|**2 + |Gmn|**2 }                       ! 
!               = [Pi / Pi * (ka)**2] * SUM { |Fmn|**2 + |Gmn|**2 }                 !
! where Fmn and Gmn are the expansion coefficients of the scattered field, k is the !
! wave number and a is the characteristic dimension of the particle. Note that the  !
! normalization constant snorm is snorm = Pi * (k*a)**2 and that the scattered      !
! field coefficients (which are stored in the vector c) depend on the direction and !
! polarization of the incident wave.                                                ! 
!------------------------------------------------------------------------------------                 
  use parameters
  implicit none        
  integer    :: Mrank, Nrank, Nmax
  real(O)    :: wavenumber, snorm, Cscat, Qscat
  complex(O) :: c(2*Nmax)
!
  integer    :: k, m, N0, l
  real(O)    :: sum, fct, wavenumber2
!
  sum = 0._O
  do m = 0, Mrank
    if (m == 0) then
      do k = 1, Nrank
        sum = sum + abs(c(k))**2 + abs(c(k+Nmax))**2              
      end do
    else        
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      do l = 1, 2
        do k = 1, Nrank - m + 1
          sum = sum + abs(c(k+N0))**2 + abs(c(k+N0+Nmax))**2               
        end do
        N0 = N0 + Nrank - m + 1
      end do         
    end if
  end do
  fct = sum * Pi
  wavenumber2 = wavenumber * wavenumber
  Cscat = fct / wavenumber2
  Qscat = fct / snorm      
end subroutine CQscat
! **********************************************************************************
subroutine CQext (c, Mrank, Nrank, Nmax, thetaGI, phiGI, alpha, beta, gamma,        &  
           alphap, wavenumber, snorm, Cext, Qext)
!------------------------------------------------------------------------------------ 
! The routine computes the extinction cross section (according to optical theorem)  !
!        C_ext = (4 * Pi / k**2) * IM{ ES * e0 }                                    !
! and the extinction efficiency                                                     ! 
!        Q_ext = C_ext / (Pi * a**2)                                                !
!              = [4 / (ka)**2] * IM{ ES * e0 }                                      !                                
!              = [4Pi / Pi * (ka)**2] * IM{ ES * e0 }                               !
! where ES is the scattered field in the direction (thetaGI,phiGI), e0 is the unit  !
! vector along the incident polarization direction, i.e.,                           !
!        e0 = cos(alphap) * e_theta0 + sin(alphap) * e_phi0,                        !
! alphap, is the polarization angle, k is the wave number and a is the              !
! characteristic length of the particle. The incident wave is a linearly polarized  !
! plane wave traveling in the direction (thetaGI,phiGI). Note that the normalization! 
! constant snorm is snorm = Pi * (k * a)**2 and that the scattered field            !
! coefficients depend on the direction and polarization of the incident wave.       !
! The incident wave is a linearly polarized vector plane wave.                      !
!------------------------------------------------------------------------------------
  use parameters
  implicit none
  integer     :: Mrank, Nrank, Nmax
  real(O)     :: thetaGI, phiGI, alpha, beta, gamma, alphap, wavenumber, snorm,     &
                 Cext, Qext
  complex(O)  :: c(2*Nmax)
!
  integer     :: m, k, N0, l
  real(O)     :: thetaGS, phiGS, thetaLS, phiLS, cosu, sinu, Img, wavenumber2
  complex(O)  :: sum(3)
  complex(O),allocatable :: Minf(:,:), Ninf(:,:)
!
  allocate (Minf(3,Nmax), Ninf(3,Nmax))
  thetaGS = thetaGI
  phiGS  = phiGI
  call T_spherical_global_local (thetaGS, phiGS, alpha, beta, gamma, thetaLS, phiLS)          
  call MN_infinit_complete (thetaLS, phiLS, Mrank, Nrank, Nmax, .true., Minf, Ninf)            
  sum(1) = zero
  sum(2) = zero
  sum(3) = zero
  do m = 0, Mrank
    if (m == 0) then
      do k = 1, Nrank
        sum(1) = sum(1) + (Minf(1,k) * c(k) + Ninf(1,k) * c(Nmax+k))
        sum(2) = sum(2) + (Minf(2,k) * c(k) + Ninf(2,k) * c(Nmax+k))
        sum(3) = sum(3) + (Minf(3,k) * c(k) + Ninf(3,k) * c(Nmax+k))
      end do
    else        
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      do l = 1, 2
        do k = 1, Nrank - m + 1
          sum(1) = sum(1) + (Minf(1,N0+k) * c(N0+k) + Ninf(1,N0+k) * c(Nmax+N0+k))
          sum(2) = sum(2) + (Minf(2,N0+k) * c(N0+k) + Ninf(2,N0+k) * c(Nmax+N0+k))
          sum(3) = sum(3) + (Minf(3,N0+k) * c(N0+k) + Ninf(3,N0+k) * c(Nmax+N0+k))
        end do
        N0 = N0 + Nrank - m + 1
      end do          
    end if
  end do
  call angle_unitvct_ItL_ItG (thetaGS, phiGS, thetaLS, phiLS, alpha, beta, gamma,   &
       cosu, sinu)      
  Img  = aimag((cosu * sum(2) - sinu * sum(3)) * cos(alphap)) 
  Img  = Img + aimag((sinu * sum(2) + cosu * sum(3)) * sin(alphap))      
  Img  = 4._O * Img * Pi
  wavenumber2 = wavenumber * wavenumber
  Cext = Img / wavenumber2
  Qext = Img / snorm
  deallocate (Minf, Ninf)      
end subroutine CQext
! **********************************************************************************
subroutine CQextGen (c, Mrank, Nrank, Nmax, thetaGI, phiGI, alpha, beta, gamma,     &
           alphap, x0, y0, z0, w0, TypeExcit, wavenumber, snorm, Cext, Qext)
!------------------------------------------------------------------------------------ 
! The routine computes the extinction cross section                                 !
!       C_ext = - Pi * SUM { RE{ Amn * conjg(Fmn) + conjg(Bmn) * Gmn } } / k**2     !
! and the extinction efficiency                                                     !
!       Q_ext = C_ext / (Pi * a**2)                                                 !
! where Fmn and Gmn are the expansion coefficients of the scattered field, Amn and  !
! Bmn are the expansion coefficients of the incident wave, k is the wave number and !
! a is the characteristic length of the particle. Note that the normalization       !
! constant snorm is snorm = Pi*(k*a)**2. The incident wave is a linearly polarized  !
! vector plane wave.                                                                !
!------------------------------------------------------------------------------------
  use parameters
  implicit none
  integer      :: Mrank, Nrank, Nmax
  real(O)      :: thetaGI, phiGI, alpha, beta, gamma, alphap, x0, y0, z0, w0,       &
                  wavenumber, snorm, Cext, Qext 
  character(5) :: TypeExcit                  
  complex(O)   :: c(2*Nmax)
!
  integer      :: m, k, N0, l
  real(O)      :: fct, wavenumber2
  complex(O)   :: sum
  complex(O),allocatable :: cc(:)
!
  allocate (cc(2*Nmax))
  if (TypeExcit(1:5) == 'PLANE') then
    call PWcoefficients_ab (thetaGI, phiGI, alpha, beta, gamma, alphap, Mrank,      &
         Nrank, Nmax, cc)
  else if (TypeExcit(1:5) == 'GAUSS') then
    call GBcoefficients_ab (wavenumber, x0, y0, z0, w0, thetaGI, phiGI, alpha,      &
         beta, gamma, alphap, Mrank, Nrank, Nmax, cc)        
  end if
  sum = zero
  do m = 0, Mrank
    if (m == 0) then
      do k = 1, Nrank       
        sum = sum + cc(k) * conjg(c(k)) + cc(Nmax+k) * conjg(c(Nmax+k))        
      end do
    else        
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      do l = 1, 2
        do k = 1, Nrank - m + 1          
          sum = sum + cc(N0+k) * conjg(c(N0+k)) +                                   &
                      cc(Nmax+N0+k) * conjg(c(Nmax+N0+k))
        end do
        N0 = N0 + Nrank - m + 1
      end do      
    end if
  end do
  fct  = - Pi * real(sum,O)
  wavenumber2 = wavenumber * wavenumber
  Cext =   fct / wavenumber2
  Qext =   fct / snorm  
  deallocate (cc)      
end subroutine CQextGen
! **********************************************************************************
subroutine CQextPolGen (c, Mrank, Nrank, Nmax, thetaGI, phiGI, alpha, beta, gamma,  &
           epol_beta, epol_alpha, alphapGauss, x0, y0, z0, w0, TypeExcit,           &
           wavenumber, snorm, Cext, Qext)
!------------------------------------------------------------------------------------ 
! The routine computes the extinction cross section                                 !
!       C_ext = - Pi * SUM { RE{ Amn * conjg(Fmn) + conjg(Bmn) * Gmn } } / k**2     !
! and the extinction efficiency                                                     !
!       Q_ext = C_ext / (Pi * a**2)                                                 !
! where Fmn and Gmn are the expansion coefficients of the scattered field, and Amn  !
! and Bmn are the expansion coefficients of the incident wave. The incident wave is !
! an elliptically polarized vector plane wave or a linearly polarized Gaussian beam.! 
!------------------------------------------------------------------------------------
  use parameters
  implicit none
  integer      :: Mrank, Nrank, Nmax
  real(O)      :: thetaGI, phiGI, alpha, beta, gamma, alphapGauss, x0, y0, z0, w0,  &
                  wavenumber, snorm, Cext, Qext 
  character(5) :: TypeExcit                  
  complex(O)   :: epol_beta, epol_alpha, c(2*Nmax)
!
  integer      :: m, k, N0, l
  real(O)      :: fct, wavenumber2
  complex(O)   :: sum
  complex(O),allocatable :: cc(:)
!
  allocate (cc(2*Nmax))
  if (TypeExcit(1:5) == 'PLANE') then    
    call PWcoefficients_Polab (thetaGI, phiGI, alpha, beta, gamma, epol_beta,       &
         epol_alpha, Mrank, Nrank, Nmax, cc) 	 
  else if (TypeExcit(1:5) == 'GAUSS') then
    call GBcoefficients_ab (wavenumber, x0, y0, z0, w0, thetaGI, phiGI, alpha,      &
         beta, gamma, alphapGauss, Mrank, Nrank, Nmax, cc)        
  end if
  sum = zero
  do m = 0, Mrank
    if (m == 0) then
      do k = 1, Nrank       
        sum = sum + cc(k) * conjg(c(k)) + cc(Nmax+k) * conjg(c(Nmax+k))        
      end do
    else        
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      do l = 1, 2
        do k = 1, Nrank - m + 1          
          sum = sum + cc(N0+k) * conjg(c(N0+k)) +                                   &
                      cc(Nmax+N0+k) * conjg(c(Nmax+N0+k))
        end do
        N0 = N0 + Nrank - m + 1
      end do      
    end if
  end do
  fct  = - Pi * real(sum,O)
  wavenumber2 = wavenumber * wavenumber
  Cext =   fct / wavenumber2
  Qext =   fct / snorm  
  deallocate (cc)      
end subroutine CQextPolGen
! **********************************************************************************
subroutine CQextPolPlane (c, Mrank, Nrank, Nmax, thetaGI, phiGI, alpha, beta,       & 
           gamma, epol_beta, epol_alpha, wavenumber, snorm, Cext, Qext)
!------------------------------------------------------------------------------------ 
! The routine computes the extinction cross section                                 !
!       C_ext = - Pi * SUM { RE{ Amn * conjg(Fmn) + conjg(Bmn) * Gmn } } / k**2     !
! and the extinction efficiency                                                     !
!       Q_ext = C_ext / (Pi * a**2)                                                 !
! where Fmn and Gmn are the expansion coefficients of the scattered field, and Amn  !
! and Bmn are the expansion coefficients of the incident wave. The incident wave    !
! is an elliptically polarized vector plane wave.                                   ! 
!------------------------------------------------------------------------------------
  use parameters
  implicit none
  integer      :: Mrank, Nrank, Nmax
  real(O)      :: thetaGI, phiGI, alpha, beta, gamma, wavenumber, snorm, Cext, Qext                
  complex(O)   :: epol_beta, epol_alpha, c(2*Nmax)  
!
  integer      :: m, k, N0, l
  real(O)      :: fct, wavenumber2
  complex(O)   :: sum
  complex(O),allocatable :: cc(:)
! 
  allocate (cc(2*Nmax))     
  call PWcoefficients_Polab (thetaGI, phiGI, alpha, beta, gamma, epol_beta,         &
       epol_alpha, Mrank, Nrank, Nmax, cc) 	   
  sum = zero
  do m = 0, Mrank
    if (m == 0) then
      do k = 1, Nrank       
        sum = sum + cc(k) * conjg(c(k)) + cc(Nmax+k) * conjg(c(Nmax+k))        
      end do
    else        
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      do l = 1, 2
        do k = 1, Nrank - m + 1          
          sum = sum + cc(N0+k) * conjg(c(N0+k)) +                                   &
                      cc(Nmax+N0+k) * conjg(c(Nmax+N0+k))
        end do
        N0 = N0 + Nrank - m + 1
      end do      
    end if
  end do
  fct  = - Pi * real(sum,O)
  wavenumber2 = wavenumber * wavenumber
  Cext =   fct / wavenumber2
  Qext =   fct / snorm  
  deallocate (cc)      
end subroutine CQextPolPlane
! **********************************************************************************
subroutine PhaseMatrix (S, Z)
!-----------------------------------------------------------------------------------      
! The routine computes the phase matrix Z, where S is the amplitude matrix.        !
!-----------------------------------------------------------------------------------      
  use parameters
  implicit none  
  real(O)     :: Z(4,4)
  complex(O)  :: S(2,2)
!
  Z(1,1) =   0.5_O * (abs(S(1,1))**2 + abs(S(1,2))**2 +                             &
             abs(S(2,1))**2 + abs(S(2,2))**2)     
  Z(1,2) =   0.5_O * (abs(S(1,1))**2 - abs(S(1,2))**2 +                             &
             abs(S(2,1))**2 - abs(S(2,2))**2)
  Z(1,3) = - real (S(1,1) * conjg(S(1,2)) + S(2,2) * conjg(S(2,1)),O)
  Z(1,4) = - aimag(S(1,1) * conjg(S(1,2)) - S(2,2) * conjg(S(2,1)))      
  Z(2,1) =   0.5_O * (abs(S(1,1))**2 + abs(S(1,2))**2 -                             &
             abs(S(2,1))**2 - abs(S(2,2))**2)
  Z(2,2) =   0.5_O * (abs(S(1,1))**2 - abs(S(1,2))**2 -                             &
             abs(S(2,1))**2 + abs(S(2,2))**2)
  Z(2,3) = - real (S(1,1) * conjg(S(1,2)) - S(2,2) * conjg(S(2,1)),O)
  Z(2,4) = - aimag(S(1,1) * conjg(S(1,2)) + S(2,2) * conjg(S(2,1)))      
  Z(3,1) = - real (S(1,1) * conjg(S(2,1)) + S(2,2) * conjg(S(1,2)),O)
  Z(3,2) = - real (S(1,1) * conjg(S(2,1)) - S(2,2) * conjg(S(1,2)),O)
  Z(3,3) =   real (S(1,1) * conjg(S(2,2)) + S(1,2) * conjg(S(2,1)),O)
  Z(3,4) =   aimag(S(1,1) * conjg(S(2,2)) + S(2,1) * conjg(S(1,2)))      
  Z(4,1) = - aimag(S(2,1) * conjg(S(1,1)) + S(2,2) * conjg(S(1,2)))
  Z(4,2) = - aimag(S(2,1) * conjg(S(1,1)) - S(2,2) * conjg(S(1,2)))
  Z(4,3) =   aimag(S(2,2) * conjg(S(1,1)) - S(1,2) * conjg(S(2,1)))
  Z(4,4) =   real (S(2,2) * conjg(S(1,1)) - S(1,2) * conjg(S(2,1)),O)      
end subroutine PhaseMatrix
! **********************************************************************************
subroutine ExtinctionMatrix (wavenumber, S, K)
!-----------------------------------------------------------------------------------
! The routine computes the extinction matrix K, where S is the amplitude matrix.   !     
!-----------------------------------------------------------------------------------
  use parameters
  implicit none                      
  real(O)     :: wavenumber
  complex(O)  :: S(2,2)
  real(O)     :: K(4,4), fct      
!
  fct = 2._O * Pi / wavenumber
  K(1,1) =   aimag(S(1,1) + S(2,2))  * fct
  K(1,2) =   aimag(S(1,1) - S(2,2))  * fct
  K(1,3) = - aimag(S(1,2) + S(2,1))  * fct
  K(1,4) =   real(S(2,1) - S(1,2),O) * fct      
  K(2,1) =   K(1,2)
  K(2,2) =   K(1,1)
  K(2,3) =   aimag(S(2,1) - S(1,2))   * fct
  K(2,4) =   real (S(1,2) + S(2,1),O) * fct      
  K(3,1) =   K(1,3)
  K(3,2) = - K(2,3)
  K(3,3) =   K(1,1)
  K(3,4) =   real(S(2,2) - S(1,1),O) * fct      
  K(4,1) =   K(1,4)
  K(4,2) = - K(2,4)
  K(4,3) = - K(3,4)
  K(4,4) =   K(1,1)      
end subroutine ExtinctionMatrix                   
! **********************************************************************************
subroutine extend_vector_positive (c, cc, m, Mstart, Nrank, Nmax, Nmaxmax)
  use parameters
  implicit none
  integer    :: m, Nrank, Nmax, Nmaxmax, Mstart, N0, k
  complex(O) :: c(2*Nmax), cc(2*Nmaxmax)      
!
  if (m == Mstart) then
    do k = 1, 2*Nmaxmax
      cc(k) = zero
    end do
  end if               
  if (m == 0) then
    do k = 1, Nrank
      cc(k) = c(k)
      cc(k+Nmaxmax) = c(k+Nmax)
    end do
  else      
    N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
    do k = 1, Nrank - m + 1
      cc(k+N0) = c(k)
      cc(k+N0+Nmaxmax) = c(k+Nmax)
    end do
  end if      
end subroutine extend_vector_positive
! **********************************************************************************
subroutine extend_vector_negative (c, cc, m, Nrank, Nmax, Nmaxmax)
  use parameters
  implicit none
  integer    :: m, Nrank, Nmax, Nmaxmax, N0, k
  complex(O) :: c(2*Nmax), cc(2*Nmaxmax)      
!
  N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
  N0 = N0 + Nrank - m + 1  
  do k = 1, Nrank - m + 1
    cc(k+N0) = c(k)
    cc(k+N0+Nmaxmax) = c(k+Nmax)
  end do      
end subroutine extend_vector_negative
! **********************************************************************************
subroutine delta_DSCS (Ntheta, h, v, oldh, oldv, eps, nconv)
  use parameters
  use derived_parameters
  implicit none
  integer  :: Ntheta, i, nconv, pconv, sconv
  real(O)  :: h(Ntheta), v(Ntheta), oldh(Ntheta), oldv(Ntheta), eps, delta, deltaM
  logical  :: smallh, smallv 
!
  pconv  = 0
  sconv  = 0
  smallh = .false.
  smallv = .false.
  do i = 1, Ntheta
    delta  = abs(h(i) - oldh(i))
    deltaM = eps * h(i)    
    if (h(i) >= MachEps) then
      if (delta < deltaM) pconv = pconv + 1
    else
      pconv  = pconv + 1     
      smallh = .true.
    end if
    delta  = abs(v(i) - oldv(i))
    deltaM = eps * v(i)    
    if (v(i) >= MachEps) then
      if (delta < deltaM) sconv = sconv + 1
    else
      sconv  = sconv + 1     
      smallv = .true.
    end if
    oldh(i) = h(i)
    oldv(i) = v(i)
  end do
  nconv = min(pconv, sconv)  
  if (smallh) then
    print "(/,2x,'Warning: at least one value of the parallel component ')"     
    print "(  2x,'of the DSCS is smaller than the machine precision;')"  
    if (pconv == Ntheta) nconv = 0
  end if
  if (smallv) then
    print "(/,2x,'Warning: at least one value of the perpendicular component ')"        
    print "(  2x,'of the DSCS is smaller than the machine precision;')"  
    if (sconv == Ntheta) nconv = 0
  end if                      
end subroutine delta_DSCS 
! **********************************************************************************
subroutine NodesAlphaBetaGamma (alphamin, alphamax, Nalpha, betamin, betamax, Nbeta,&
           gammamin, gammamax, Ngamma, UseSimpson, xalpha, walpha, xbeta, wbeta,    &
           xgamma, wgamma)
  use parameters
  implicit none
  integer  :: Nalpha, Nbeta, Ngamma
  real(O)  :: alphamin, alphamax, betamin, betamax, gammamin, gammamax,             &
              xalpha(Nalpha), walpha(Nalpha), xbeta(Nbeta), wbeta(Nbeta),           &
                          xgamma(Ngamma), wgamma(Ngamma)
  logical  :: UseSimpson
!
  integer  :: ialpha, ibeta, igamma
  real(O)  :: norm, aalpha, balpha, abeta, bbeta, beta, agamma, bgamma
!  
  if (Nalpha == 1) then
    xalpha(1) = 0.5_O * (alphamax + alphamin)
    walpha(1) = 1._O
  else     
    aalpha = alphamin
    balpha = alphamax
    call Simpson (aalpha, balpha, Nalpha, xalpha, walpha) 
    norm = alphamax - alphamin
    norm = 1._O / norm     
    do ialpha = 1, Nalpha
      walpha(ialpha) = walpha(ialpha) * norm
    end do       
  end if  
  if (Nbeta == 1) then
    beta = 0.5_O * (betamax + betamin)
    if (.not. UseSimpson) then
      xbeta(1) = beta
    else
      xbeta(1) = cos(beta)
    end if
    wbeta(1) = 1._O
  else  
    if (.not. UseSimpson) then
      abeta = betamin
      bbeta = betamax   
      call Gauss_Legendre (abeta, bbeta, Nbeta, wbeta, xbeta)   
    else
      abeta = cos(betamax)
      bbeta = cos(betamin)
      call Simpson (abeta, bbeta, Nbeta, xbeta, wbeta)     
    end if
    norm = cos(betamin) - cos(betamax)
    norm = 1._O / norm     
    do ibeta = 1, Nbeta    
      if (.not. UseSimpson) then 
        beta = xbeta(ibeta)
        wbeta(ibeta) = wbeta(ibeta) * norm * sin(beta)
      else
        wbeta(ibeta) = wbeta(ibeta) * norm
      end if
    end do           
  end if    
  if (Ngamma == 1) then
    xgamma(1) = 0.5_O * (gammamax + gammamin)
    wgamma(1) = 1._O
  else      
    agamma = gammamin
    bgamma = gammamax           
    call Simpson (agamma, bgamma, Ngamma, xgamma, wgamma)
    norm = gammamax - gammamin
    norm = 1._O / norm    
    do igamma = 1, Ngamma 
      wgamma(igamma) = wgamma(igamma) * norm
    end do
  end if          
end subroutine NodesAlphaBetaGamma
! **********************************************************************************
subroutine IndexElements (Nelem, MatrixElem, IndI, IndJ, NameElem)  
  implicit none
  integer      :: Nelem, MatrixElem(16), IndI(Nelem), IndJ(Nelem), i, j, ij, k
  character(2) :: AllNameElem(16), NameElem(Nelem)
  logical      :: more 
!
  AllNameElem  = (/'11', '12', '13', '14', '21', '22', '23', '24',                  &
                   '31', '32', '33', '34', '41', '42', '43', '44'/)
  do i = 1, 4
    do j = 1, 4
      ij = 10 * i + j
      k  = 0
      more = .true.
      do while (more .and. k < Nelem)
        k = k + 1
        if (ij == MatrixElem(k)) then           
          IndI(k) = i
          IndJ(k) = j
          more    = .false.
        end if
      end do      
    end do
  end do                      
  do k = 1, Nelem
    i  = IndI(k)
    j  = IndJ(k)
    ij = j + 4 * (i - 1)
    NameElem(k) = AllNameElem(ij)  
  end do           
end subroutine IndexElements 
! **********************************************************************************
! *          SCATTERING CHARACTERISTICS: PARTICLE ON OR NEAR A SUBSTRATE           *
! **********************************************************************************
subroutine DiffScatCrossSectPARTSUB (TypeScat, k, ind_ref_s, z0, snorm, Nrank,      &
           Mrank, FileTmat)     
!------------------------------------------------------------------------------------
! Let                                                                               !
!            ES = [exp(jkR) / kR] * ESinf = [exp(jkR) / R] * FS                     !
! be the scattered field in the direction (theta,phi), where FS is the far-field    !
! pattern and ESinf = k * FS. If we consider the decomposition                      !
!            ESinf = ESinf_theta * e_theta + ESinf_phi * e_phi,                     !
! where e_theta and e_phi are the spherical unit vectors, the routine computes the  !
! differential scattering cross sections for parallel and perpendicular             !
! polarizations:                                                                    !
!            h = |ESinf_theta|**2 / k**2 = |FS_theta|**2                            !
! and                                                                               !
!            v = |ESinf_phi|**2 / k**2   = |FS_phi|**2,                             !
! respectively, or the normalized differential scattering cross sections for        !
! parallel and perpendicular polarizations:                                         ! 
!     hn = |ESinf_theta|**2 / [Pi * (k*anorm)**2] = |FS_theta|**2 / [Pi * anorm**2] !
! and                                                                               !
!     vn = |ESinf_phi|**2 / [Pi * (k*anorm)**2]   = |FS_phi|**2 / [Pi * anorm**2],  !
! respectively. The subroutine also computes the electromagnetic fields ESinf_theta !
! and ES_inf_phi.                                                                   !
!                                                                                   !
! Input parameters:                                                                 !
! - TypeScat (integer) - for TypeScat = 1 the particle is illuminated by a pane     !
!   wave traveling in the ambient medium, while for TypeScat = 2 the particle is    !
!   illuminated by a plane wave traveling in the substrate. In the first case the   !
!   direction of the wave vector is (theta = beta ,phi = Pi), while in the second   !
!   case the direction is (theta = Pi - beta,phi = 0), where beta is the incident   !
!   angle.                                                                          ! 
! - k (real) - wave number of the ambient medium.                                   !
! - ind_ref_s (complex) - refractive index of the substrate.                        !
! - z0 (real) - axial position of the substrate with respect to the particle        !
!   coordinate system.                                                              !
! - snorm (real) - normalization constant.                                          !
! - FileTmat (character(80) - name of the data file containing the T matrix.        !
! - Mrank, Nrank (integer variables) - Mrank is the number of azimuthal modes       !
!   (maximum degree), while Nrank is the maximum expansion order.                   !
!------------------------------------------------------------------------------------                       
  use parameters
  implicit none 
  integer       :: TypeScat, Nrank, Mrank
  real(O)       :: snorm, k, z0
  complex(O)    :: ind_ref_s
  character(80) :: FileTmat  
!  
  integer       :: NthetaGS, Nphi, Ntheta(NphiMax), Mstart, Nmax, Nmaxmax,          &
                   m, i, itheta, iphi, ntl, mtl, NphiAL, NthetaAL
  real(O)       :: beta, alphap, phiGS, phi(NphiMax), thetamin(Nphimax),            &
                   thetamax(NphiMax)
  character(80) :: FileDSCS, FileEMF
  logical       :: ComputeDSCS, ComputeFields, normalized, WriteInputInfo
  real(O),allocatable    :: h(:), v(:)
  complex(O),allocatable :: t(:,:), c(:), c1(:), cc(:), em(:,:), emf(:,:,:)
! -----------------------------------------------------------------------------------
!                                 Read the input file                               !
! -----------------------------------------------------------------------------------               
  call readinputPARTSUB1 ( beta, alphap, ComputeDSCS, ComputeFields, NthetaGS,      &
       phiGS, Nphi, phi, Ntheta, thetamin, thetamax, normalized, FileDSCS, FileEMF, &
       WriteInputInfo )  
  if (.not. ComputeDSCS .and. .not. ComputeFields) return
! -----------------------------------------------------------------------------------
!                                       Main                                        !
! ----------------------------------------------------------------------------------- 
  print "(/,2x,'Scattering Characteristics of a Particle on or Near a Plane Surface')"  
  print "(  2x,'-------------------------------------------------------------------')"  
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)  
  allocate (c(2*Nrank), c1(2*Nrank), cc(2*Nmaxmax))        
  open (unit = iTmat, file = FileTmat, status = "old",  position = "rewind")   
  call read_HeadFileTmat (ntl, mtl)
  call check_dimensionMat (ntl, mtl, Nrank)             
  allocate (t(2*ntl, 2*mtl))
  Mstart = 0
  do m = Mstart, Mrank    
    if (m == 0) then
      Nmax = Nrank
    else
      Nmax = Nrank - m + 1
    end if    
    call read_FileTmat (ntl, mtl, t)
    if (TypeScat == 1) then                      
      call PWcoefficientsPARTSUB (beta, alphap, m, Nrank, Nmax, c)
      call PWcoefficientsPARTSUBrefl (beta, alphap, z0, k, ind_ref_s, m, Nrank,     &
           Nmax, c1)
      call sum_vectors (c, c1, 2*Nmax)     
    else
      call PWcoefficientsPARTSUBtrans (beta, alphap, z0, k, ind_ref_s, m, Nrank,    &
           Nmax, c)
    end if
    call product_matrix_vector (2*Nmax, 2*Nmax, t, 2*ntl, 2*mtl, c, c1)
    call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)         
    if (m /= 0) then
      call matrix_m_negativ (Nmax, Nmax, t, ntl, mtl)  
      if (TypeScat == 1) then    
        call PWcoefficientsPARTSUB (beta, alphap, - m, Nrank, Nmax, c)
        call PWcoefficientsPARTSUBrefl (beta, alphap, z0, k, ind_ref_s, - m, Nrank, &
             Nmax, c1)
        call sum_vectors (c, c1, 2*Nmax)      
      else
        call PWcoefficientsPARTSUBtrans (beta, alphap, z0, k, ind_ref_s, - m, Nrank,&
             Nmax, c)
      end if
      call product_matrix_vector (2*Nmax, 2*Nmax, t, 2*ntl, 2*mtl, c, c1)
      call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)                                          
    end if
  end do
  close (unit = iTmat)
  deallocate (t, c, c1) 
  if (ComputeDSCS) then
    open (unit = iDSCS, file = FileDSCS, status = "replace")
    allocate (h(NthetaGS), v(NthetaGS), em(3,NthetaGS)) 
    do i = 1, 3
      do itheta = 1, NthetaGS
        em(i,itheta) = zero
      end do
    end do         
    call DSCSPARTSUB (1, cc, k, z0, ind_ref_s, Mrank, Nrank, Nmaxmax, NthetaGS,     &
         phiGS, snorm, em, normalized, h, v)
    call DSCSPARTSUB (2, cc, k, z0, ind_ref_s, Mrank, Nrank, Nmaxmax, NthetaGS,     &
         phiGS, snorm, em, normalized, h, v)                                
    if (WriteInputInfo) call inputDSCS_EMF (.true., TypeScat, k, ind_ref_s, z0,     &
                             Mrank, Nrank, phiGS, beta, alphap, snorm, normalized,  &
                             FileTmat)
    call write_DSCSPARTSUB1 (NthetaGS, normalized, h, v)                    
    close (unit = iDSCS) 
    deallocate (h, v, em)
  end if
  if (ComputeFields) then
    open (unit = iSCAT, file = FileEMF,  status = "replace")
    NphiAL   = Nphi
    NthetaAL = 0
    do iphi = 1, Nphi
      if (Ntheta(iphi) > NthetaAL) NthetaAL = Ntheta(iphi)
    end do
    allocate (emf(3,NphiAL,NthetaAL))
    do i = 1, 3
      do iphi = 1, NphiAL
        do itheta = 1, NthetaAL
          emf(i,iphi,itheta) = zero
        end do
      end do
    end do
    call EMFPARTSUB (1, cc, k, z0, ind_ref_s, Mrank, Nrank, Nmaxmax, Nphi, phi,     &
         Ntheta, thetamin, thetamax, emf, NphiAL, NthetaAL)
    call EMFPARTSUB (2, cc, k, z0, ind_ref_s, Mrank, Nrank, Nmaxmax, Nphi, phi,     &
         Ntheta, thetamin, thetamax, emf, NphiAL, NthetaAL)
    if (WriteInputInfo) call inputDSCS_EMF (.false., TypeScat, k, ind_ref_s, z0,    &
                             Mrank, Nrank, phiGS, beta, alphap, snorm, normalized,  &
                             FileTmat)
    call write_EMF (Nphi, phi, Ntheta, thetamin, thetamax, emf, NphiAL, NthetaAL)       
    close (unit = iSCAT) 
    deallocate (emf)
  end if      
  deallocate (cc)       
end subroutine DiffScatCrossSectPARTSUB 
! **********************************************************************************
subroutine readinputPARTSUB1 ( beta, alphap, ComputeDSCS, ComputeFields, NthetaGS,  &
           phiGS, Nphi, phi, Ntheta, thetamin, thetamax, normalized, FileDSCS,      &
           FileEMF, WriteInputInfo ) 
  use parameters
  implicit none     
  integer       :: NthetaGS, Nphi, Ntheta(NphiMax), iphi, ios
  real(O)       :: beta, alphap, phiGS, phi(NphiMax), thetamin(NphiMax),            &
                   thetamax(NphiMax), deg
  character(80) :: FileDSCS, FileEMF, string
  logical       :: ComputeDSCS, ComputeFields, normalized, WriteInputInfo, XFindPar  
! -----------------------------------------------------------------------------------
!                        Read the input file FileInputPARTSUB                       !
! -----------------------------------------------------------------------------------             
  open (unit = iInputPARTSUB, file = FileInputPARTSUB, status = "old",              &
        position = "rewind")                
  beta   = 0._O
  alphap = 0._O
  string = 'IncWave'
  if (XFindPar (iInputPARTSUB, string)) then
    read (iInputPARTSUB, *, iostat = ios) beta
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable beta;')"
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
  alphap = alphap * deg
! 
  ComputeDSCS   = .true.
  ComputeFields = .true.
  NthetaGS = 181 
  phiGS    = 0._O  
  Nphi = 1
  do iphi = 1, NphiMax
    phi(iphi)      = 0._O
    Ntheta(iphi)   = 1
    thetamin(iphi) = 0._O
    thetamax(iphi) = 0._O  
  end do 
  normalized = .true.
  FileDSCS   = '../OUTPUTFILES/DSCSpartsub.dat'
  FileEMF    = '../OUTPUTFILES/EMFields.dat'
  WriteInputInfo = .true.
  string     = 'DSCSEM'
  if (XFindPar (iInputPARTSUB, string)) then
    read (iInputPARTSUB, *, iostat = ios) ComputeDSCS 
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ComputeDSCS;')"
      stop
    end if
    read (iInputPARTSUB, *, iostat = ios) ComputeFields
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ComputeFields;')"
      stop
    end if
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
    read (iInputPARTSUB, *, iostat = ios) Nphi
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nphi;')"
      stop
    end if
    do iphi = 1, Nphi
      read (iInputPARTSUB, *, iostat = ios) phi(iphi)
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable phi;')"
        stop
      end if
      read (iInputPARTSUB, *, iostat = ios) Ntheta(iphi)
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Ntheta;')"
        stop
      end if
      read (iInputPARTSUB, *, iostat = ios) thetamin(iphi)
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable thetamin;')"
        stop
      end if      
      read (iInputPARTSUB, *, iostat = ios) thetamax(iphi)
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable thetamax;')"
        stop
      end if      
    end do    
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
    read (iInputPARTSUB, *, iostat = ios) FileEMF
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable FileEMF;')"
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
  if (ComputeDSCS) then      
    call check_azimuthal_plane (phiGS)
    phiGS = phiGS * deg                 
  end if
  if (ComputeFields) then
    call check_Nphi (Nphi)
    call check_teta_phiminmax (.false., Nphi, phi, Ntheta, thetamin, thetamax)
    do iphi = 1, Nphi    
      phi(iphi)      = phi(iphi)      * deg
      thetamin(iphi) = thetamin(iphi) * deg
      thetamax(iphi) = thetamax(iphi) * deg      
    end do
  end if
  close (unit = iInputPARTSUB)
end subroutine readinputPARTSUB1
! **********************************************************************************
subroutine DSCSPARTSUB (TypeField, c, wavenumber, z0, ind_ref, Mrank, Nrank, Nmax,  &
           NthetaGS, phiAZIMUT, snorm, em, normalized, h, v)
!------------------------------------------------------------------------------------     
! The routine computes the differential scattering cross sections of an axisymmetric!
! particle on a plane substrate in the azimuthal plane phiAZIMUTH.                  !
! -----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: TypeField, Mrank, Nrank, Nmax, NthetaGS
  real(O)    :: phiAZIMUT, snorm, h(NthetaGS), v(NthetaGS), wavenumber, z0      
  complex(O) :: ind_ref, em(3,NthetaGS), c(2*Nmax)
  logical    :: normalized
!
  integer    :: itheta, m, k, N0, l
  real(O)    :: thetaGS, phiGS, fact
  complex(O) :: sum(3)
  complex(O),allocatable :: Minf(:,:), Ninf(:,:)      
!
  if (normalized) then
    fact = snorm
  else 
    fact = wavenumber * wavenumber
  end if
  fact = 1._O / fact
  allocate (Minf(3,Nmax), Ninf(3,Nmax)) 
  do itheta = 1, NthetaGS  
    thetaGS = Pi / 2._O + real(itheta - 1,O) * Pi / real(NthetaGS - 1,O)
    phiGS   = phiAZIMUT
    if (thetaGS > Pi) then
      thetaGS = 2._O * Pi - thetaGS
      phiGS   = phiAZIMUT + Pi
    end if          
    if (TypeField == 1) then               
      call MN_infinit_complete (thetaGS, phiGS, Mrank, Nrank, Nmax, .true.,         &
           Minf, Ninf)
    else if (TypeField == 2) then            
      call MN_infinit_reflection_complete (wavenumber, z0, ind_ref, thetaGS, phiGS, &
           Mrank, Nrank, Nmax, Minf, Ninf)             
    end if                         
    sum(1) = zero
    sum(2) = zero
    sum(3) = zero
    do m = 0, Mrank
      if (m == 0) then
        do k = 1, Nrank
          sum(1) = sum(1) + (Minf(1,k) * c(k) + Ninf(1,k) * c(Nmax+k))
          sum(2) = sum(2) + (Minf(2,k) * c(k) + Ninf(2,k) * c(Nmax+k))
          sum(3) = sum(3) + (Minf(3,k) * c(k) + Ninf(3,k) * c(Nmax+k))
        end do
      else          
        N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
        do l = 1, 2
          do k = 1, Nrank - m + 1
            sum(1) = sum(1) + (Minf(1,N0+k) * c(N0+k) + Ninf(1,N0+k) * c(Nmax+N0+k))
            sum(2) = sum(2) + (Minf(2,N0+k) * c(N0+k) + Ninf(2,N0+k) * c(Nmax+N0+k))
            sum(3) = sum(3) + (Minf(3,N0+k) * c(N0+k) + Ninf(3,N0+k) * c(Nmax+N0+k))
          end do           
          N0 = N0 + Nrank - m + 1
        end do
      end if
    end do                          
    em(1,itheta) = em(1,itheta) + sum(1)
    em(2,itheta) = em(2,itheta) + sum(2)
    em(3,itheta) = em(3,itheta) + sum(3)                        
    h(itheta) = abs(em(2,itheta))**2 * fact
    v(itheta) = abs(em(3,itheta))**2 * fact      
  end do
  deallocate (Minf, Ninf)      
end subroutine DSCSPARTSUB
! **********************************************************************************
subroutine EMFPARTSUB (TypeField, c, wavenumber, z0, ind_ref, Mrank, Nrank, Nmax,  &
           Nphi, phi, Ntheta, thetamin, thetamax, emf, NphiAL, NthetaAL)
!------------------------------------------------------------------------------------     
! The routine computes the electromagnetic scattered field of an axisymmetric       !
! particle on a plane substrate.                                                    !
! -----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: TypeField, Mrank, Nrank, Nmax, Nphi, Ntheta(NphiMax), NphiAL,       &
                NthetaAL
  real(O)    :: phi(NphiMax), thetamin(NphiMax), thetamax(NphiMax), wavenumber, z0      
  complex(O) :: ind_ref, emf(3,NphiAL,NthetaAL), c(2*Nmax)
!
  integer    :: iphi, itheta, m, k, N0, l
  real(O)    :: thetaGS, phiGS, dtheta
  complex(O) :: sum(3)
  complex(O),allocatable :: Minf(:,:), Ninf(:,:)      
!
  allocate (Minf(3,Nmax), Ninf(3,Nmax)) 
  do iphi = 1, Nphi
    phiGS = phi(iphi)
    if (Ntheta(iphi) /= 1) then
      dtheta = (thetamax(iphi) - thetamin(iphi)) / (Ntheta(iphi) - 1)
    else
      dtheta = 0._O
    end if          
    do itheta = 1, Ntheta(iphi)
      thetaGS = thetamin(iphi) + (itheta - 1) * dtheta              
      if (TypeField == 1) then               
        call MN_infinit_complete (thetaGS, phiGS, Mrank, Nrank, Nmax, .true.,       &
             Minf, Ninf)
      else if (TypeField == 2) then            
        call MN_infinit_reflection_complete (wavenumber, z0, ind_ref, thetaGS,      &
             phiGS, Mrank, Nrank, Nmax, Minf, Ninf)                
      end if                       
      sum(1) = zero
      sum(2) = zero
      sum(3) = zero
      do m = 0, Mrank
        if (m == 0) then
          do k = 1, Nrank
            sum(1) = sum(1) + (Minf(1,k) * c(k) + Ninf(1,k) * c(Nmax+k))
            sum(2) = sum(2) + (Minf(2,k) * c(k) + Ninf(2,k) * c(Nmax+k))
            sum(3) = sum(3) + (Minf(3,k) * c(k) + Ninf(3,k) * c(Nmax+k))
          end do
        else          
          N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
          do l = 1, 2
            do k = 1, Nrank - m + 1
              sum(1) = sum(1) + (Minf(1,N0+k) * c(N0+k) + Ninf(1,N0+k) *            &
                                   c(Nmax+N0+k))
              sum(2) = sum(2) + (Minf(2,N0+k) * c(N0+k) + Ninf(2,N0+k) *            &
                                   c(Nmax+N0+k))
              sum(3) = sum(3) + (Minf(3,N0+k) * c(N0+k) + Ninf(3,N0+k) *            &
                                   c(Nmax+N0+k))
            end do           
            N0 = N0 + Nrank - m + 1
          end do
        end if
      end do                        
      emf(1,iphi, itheta) = emf(1,iphi, itheta) + sum(1)
      emf(2,iphi, itheta) = emf(2,iphi, itheta) + sum(2)
      emf(3,iphi, itheta) = emf(3,iphi, itheta) + sum(3)                                 
    end do
  end do
  deallocate (Minf, Ninf)      
end subroutine EMFPARTSUB
! **********************************************************************************
subroutine delta_DSCSPARTSUB (m, Ntheta, h, v, oldh, oldv, eps, nconv)
  use parameters
  use derived_parameters
  implicit none
  integer  :: m, Ntheta, i, nconv, pconv, sconv
  real(O)  :: h(Ntheta), v(Ntheta), oldh(Ntheta), oldv(Ntheta), eps, delta, deltaM
  logical  :: smallh, smallv 
!
  pconv  = 0
  sconv  = 0
  smallh = .false.
  smallv = .false.  
  do i = 1, Ntheta - 2
    delta  = abs(h(i+1) - oldh(i+1))
    deltaM = eps * h(i+1)    
    if (h(i+1) >= MachEps) then
      if (delta < deltaM) pconv = pconv + 1
    else
      pconv  = pconv + 1     
      smallh = .true.            
    end if
    delta  = abs(v(i+1) - oldv(i+1))
    deltaM = eps * v(i+1)    
    if (v(i+1) >= MachEps) then
      if (delta < deltaM) sconv = sconv + 1
    else
      sconv  = sconv + 1     
      smallv = .true.
    end if
  end do
  do i = 1, Ntheta 
    oldh(i) = h(i)
    oldv(i) = v(i)
  end do
  nconv = min(pconv,sconv)  
  if (smallh .and. m /= 0) then
    print "(/,2x,'Warning: at least one value of the parallel component ')"     
    print "(  2x,'of the DSCS is smaller than the machine precision;')"  
    if (pconv == Ntheta - 2) nconv = 0    
  end if
  if (smallv .and. m /= 0) then
    print "(/,2x,'Warning: at least one value of the perpendicular component ')"        
    print "(  2x,'of the DSCS is smaller than the machine precision;')"  
    if (sconv == Ntheta - 2) nconv = 0    
  end if                      
end subroutine delta_DSCSPARTSUB         





     
