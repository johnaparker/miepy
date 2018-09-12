subroutine SCT
!------------------------------------------------------------------------------------
!   The routine computes the scattering characteristics of a particle using the     !
!   previously calculated T matrix. The T matrix is read from the file FileTmat.    ! 
!   The input parameters are specified in the file "InputSCT.dat". The significance !
!   of the input parameters is given in the description file "Description.txt".     ! 
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
  call readinputSCT ( .true., wavelength, k, FileTmat, Nrank, Mrank,                &
       axsym, sphere, chiral, RandomOrientation, MirorSym, DoNumAvrg,               &
       UseSimpson, ReducedOrder, deltaOrder, alphamin, alphamax, Nalpha,            &
       betamin, betamax, Nbeta, gammamin, gammamax, Ngamma, anorm, TypeExcit,       &
       thetaGI, phiGI, epol_beta, epol_alpha, alphapGauss, x0, y0, z0, w0,          &
       ComputeDSCS, phiGS, NthetaGS, ExtThetaDom, normalized, FileDSCS,             &
       ComputeScatPar, NthetaRND, thetaminRND, thetamaxRND, Nphi, phi, Ntheta,      &
       thetamin, thetamax, FileScat, Nelem, MatrixElem, ComputeAsymPar,             &
       NthetaAsym, NphiAsym, PrnProgress, WriteInputInfo, snorm )   
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
end 
