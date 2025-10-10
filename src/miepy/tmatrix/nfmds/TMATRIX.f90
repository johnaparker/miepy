program TMATRIX
  implicit none
  integer  :: TypeCode 
! -----------------------------------------------------------------------------------
!                           Select the type of T-matrix code                        ! 
! -----------------------------------------------------------------------------------   
   print "(/,2x,'T-Matrix Code for Light Scattering Calculation')"
   print "(  2x,'----------------------------------------------')"       
   print "(/,2x,'- enter an integer specifying the scattering problem:                ')"
   print "(  2x,'   1 - dielectric, perfectly conducting or chiral, axisymmetric particle;   ')"   
   print "(  2x,'   2 - dielectric, perfectly conducting or chiral, nonaxisymmetric particle;')"            
   print "(  2x,'   3 - axisymmetric, composite particle;                             ')"
   print "(  2x,'   4 - axisymmetric, layered particle;                               ')"
   print "(  2x,'   5 - inhomogeneous, axisymmetric particle with an arbitrary inclusion;    ')"  
   print "(  2x,'   6 - inhomogeneous sphere with a spherical inclusion;              ')"   
   print "(  2x,'   7 - inhomogeneous sphere with an arbitrarily shaped inclusion;    ')"
   print "(  2x,'   8 - inhomogeneous sphere with multiple spherical inclusions;      ')"
   print "(  2x,'   9 - cluster of arbitrarily shaped particles;                      ')"
   print "(  2x,'  10 - two dielectric spheres;                                       ')"
   print "(  2x,'  11 - cluster of dielectric spheres;                                ')"      
   print "(  2x,'  12 - cluster of dielectric spheres (recursive algorithm);          ')"
   print "(  2x,'  13 - concentrically, layered sphere;                               ')"
   print "(  2x,'  14 - dielectric or perfectly conducting, axisymmetric particle     ')"
   print "(  2x,'       on or near a plane surface;                                   ')"  
   print "(  2x,'  15 - scattering characteristics using the previously calculated T matrix; ')" 
   print "(  2x,'  16 - scattering characteristics of an ensemble of spherical particles;    ')" 
   print "(  2x,'  17 - effective wave number of a medium with spherical particles;          ')"                    
   call read_integer (TypeCode)   
   if (TypeCode < 1 .or. TypeCode > 17) then
     print "(/,2x,'Input error: the integer specifying the type of the T-matrix    ')"
     print "(  2x,'code is out of bounds;                                          ')"
     stop
   end if      
! -----------------------------------------------------------------------------------
!                                        Main                                       !
! -----------------------------------------------------------------------------------
  select case (TypeCode)
  case (1)
    call TAXSYM
  case (2)
    call TNONAXSYM
  case (3)    
    call TCOMP
  case (4)    
    call TLAY
  case (5)    
    call TINHOM
  case (6)    
    call TINHOM2SPH
  case (7)    
    call TINHOMSPH
  case (8)    
    call TINHOMSPHREC
  case (9)    
    call TMULT
  case (10)    
    call TMULT2SPH 
  case (11)    
    call TMULTSPH
  case (12)    
    call TMULTSPHREC 
  case (13)    
    call TSPHERE      
  case (14)    
    call TPARTSUB     
  case (15)    
    call SCT
  case (16)    
    call SCTAVRGSPH
  case (17)    
    call EFMED           
  end select
end 
