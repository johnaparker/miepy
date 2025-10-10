module parameters
  implicit none
  integer,parameter       :: O = kind(1.d0)
!
  real(O),parameter       :: Pi    = 3.1415926535897932384626433832795028841971_O
  complex(O),parameter    :: im    = (0._O,1._O),                                        &
                             one   = (1._O,0._O),                                        &
                             zero  = (0._O,0._O)
!
  integer,parameter       :: NsurfPD =     10,                                           &
                             NrankPD =    200,                                           &
                             NfacePD = 100000,                                           &
                             NparPD  =     10,                                           &
                             NphiMax =    361
!
  integer,parameter       :: iTmat        =  8,                                          &
                             iTmatInfo    =  9,                                          &
                             iDSCS        = 10,                                          &
                             iSCAT        = 11,                                          &
                             iSS          = 12,                                          & 			                               
                             iFEM         = 13
!
  integer,parameter       :: iamat = 14,                                                 &
                             ibmat = 15,                                                 &
                             icmat = 16,                                                 &
                             idmat = 17,                                                 &
                             iemat = 18,                                                 &
                             ifmat = 19
!
  integer,parameter       :: iOutput           = 20,                                     &
                             iInput            = 21,                                     &
                             iInputAXSYM       = 22,                                     &
                             iInputNONAXSYM    = 23,                                     &
                             iInputNONAXSYMFEM = 24,                                     &
                             iInputCOMP        = 25,                                     &
                             iInputLAY         = 26,                                     &
!
                             iInputINHOM       = 27,                                     &
                             iInputINHOM2SPH   = 28,                                     &
                             iInputINHOMSPH    = 29,                                     &
                             iInputINHOMSPHREC = 30,                                     &
!
                             iInputMULT        = 31,                                     &
                             iInputMULT2SPH    = 32,                                     &
                             iInputMULTSPH     = 33,                                     &
                             iInputMULTSPHREC  = 34,                                     &
!
                             iInputSPHERE      = 35,                                     &
                             iInputPARTSUB     = 36,                                     &
                             iInputANIS        = 37,                                     &
                             iInputEFMED       = 38,                                     &
!
                             iInputSCT         = 39,                                     &
                             iInputSCTAVRGSPH  = 40
!
  character(80),parameter ::                                                             &
                FileOutput           = "../OUTPUTFILES/Output.dat",                      &
                FileInput            = "../INPUTFILES/Input.dat",                        &
                FileInputAXSYM       = "../INPUTFILES/InputAXSYM.dat",                   &
                FileInputNONAXSYM    = "../INPUTFILES/InputNONAXSYM.dat",                &
                FileInputNONAXSYMFEM = "../INPUTFILES/InputNONAXSYMFEM.dat",             &
                FileInputCOMP        = "../INPUTFILES/InputCOMP.dat",                    &
                FileInputLAY         = "../INPUTFILES/InputLAY.dat",                     &
!
                FileInputINHOM       = "../INPUTFILES/InputINHOM.dat",                   &
                FileInputINHOM2SPH   = "../INPUTFILES/InputINHOM2SPH.dat",               & 
                FileInputINHOMSPH    = "../INPUTFILES/InputINHOMSPH.dat" ,               &
                FileInputINHOMSPHREC = "../INPUTFILES/InputINHOMSPHREC.dat",             &
!               			
                FileInputMULT        = "../INPUTFILES/InputMULT.dat",                    &
                FileInputMULT2SPH    = "../INPUTFILES/InputMULT2SPH.dat" ,               &
                FileInputMULTSPH     = "../INPUTFILES/InputMULTSPH.dat",                 &
                FileInputMULTSPHREC  = "../INPUTFILES/InputMULTSPHREC.dat",              &
!               									
                FileInputSPHERE      = "../INPUTFILES/InputSPHERE.dat",                  &
                FileInputPARTSUB     = "../INPUTFILES/InputPARTSUB.dat",                 &
                FileInputANIS        = "../INPUTFILES/InputANIS.dat",                    &
                FileInputEFMED       = "../INPUTFILES/InputEFMED.dat",                   &
!
                FileInputSCT         = "../INPUTFILES/InputSCT.dat",                     &
                FileInputSCTAVRGSPH  = "../INPUTFILES/InputSCTAVRGSPH.dat"
!
  character(15),parameter :: PathOUTPUT = "../OUTPUTFILES/",                             &
                             PathTEMP   = "../TEMPFILES/",                               &
                             PathGEOM   = "../GEOMFILES/"                       
end module parameters
! **************************************************************************************
module derived_parameters
  use parameters
  implicit none
!  
  integer,save :: NBaseDig
  integer,save :: NIterBes
  integer,save :: NIterPol
  real(O),save :: LargestPosNumber
  real(O),save :: SmallestPosNumber
  real(O),save :: MachEps           
  real(O),save :: ZeroCoord
  real(O),save :: TolJ0Val 
  real(O),save :: TolRootPol
  real(O),save :: InitBesVal
  real(O),save :: FactNBes
  real(O),save :: LargestBesVal
  real(O),save :: MaxArgBes
  real(O),save :: UpperBoundSeq 
  real(O),save :: LowerBoundSeq   
  real(O),save :: ZeroSinXX  
  real(O),save :: ZeroLUVal 
  real(O),save :: LargestSplineVal  
end module derived_parameters
                      
