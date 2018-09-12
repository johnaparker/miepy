! **********************************************************************************
! *                       RANDOM DISTRIBUTION OF PARTICLES                         * 
! *    ------------------------------------------------------------------------    *
! *    Partial list of subroutines:                                                *
! *      SphConfig,  SeqAddMet,  RANSLATEC,  RANLAP,  RANZIG,  readinputRND        *
! **********************************************************************************
subroutine SphConfig (Ntry, Npart, Ncs, Npartcs, r, Rcs, xp, yp, zp)
!-----------------------------------------------------------------------------------
! The routine generates a random distribution of particles using the sequential    !
! addition method.                                                                 !
!                                                                                  !
! Input parameters:                                                                !
! - Ntry (integer) - maximum number of calls of the sequential addition method     !
!   routine for generating the desired random distribution of particles.           !
! - Npart (integer) - number of particles.                                         !
! - Ncs (integer) - number of concentric (auxiliary) spheres.                      !
! - r (real) - particle radius.                                                    !
! - Rcs (real array) - radii of auxiliary spheres.                                 !
!                                                                                  !
! Output parameters:                                                               !
! - Npartcs (integer array) - number of particles inside each auxiliary            !
!   sphere.                                                                        !
! - xp, yp, zp (real arrays) - Cartesian coordinates of the particles.             !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer       :: Ntry, Npart, Ncs, Npartcs(Ncs)
  real(O)       :: r, Rcs(Ncs), xp(Ncs,Npart), yp(Ncs,Npart), zp(Ncs,Npart)
!  
  integer       :: i, ics
  real(O)       :: Rcirc, d 
  character(20) :: TypeRND
  real(O),allocatable :: x(:), y(:), z(:)
!                 
  Rcirc = Rcs(Ncs) - r 
  allocate (x(Npart), y(Npart), z(Npart))
  call readinputRND (TypeRND)
  call SeqAddMet (TypeRND, Ntry, Npart, Rcirc, r, x, y, z)
  do ics = 1, Ncs
    Npartcs(ics) = 0
  end do  
  do i = 1, Npart
    d = sqrt(x(i)**2 + y(i)**2 + z(i)**2)       
    if (d < Rcs(1)) then
      Npartcs(1) = Npartcs(1) + 1      
      xp(1,Npartcs(1)) = x(i) 
      yp(1,Npartcs(1)) = y(i)
      zp(1,Npartcs(1)) = z(i)
    end if
    do ics = 2, Ncs       
      if (d >= Rcs(ics-1) .and. d < Rcs(ics)) then
        Npartcs(ics) = Npartcs(ics) + 1      
        xp(ics,Npartcs(ics)) = x(i) 
        yp(ics,Npartcs(ics)) = y(i)
        zp(ics,Npartcs(ics)) = z(i)
      end if
    end do              
  end do    
  deallocate (x, y, z)  
end subroutine SphConfig
! ************************************************************************************
subroutine readinputRND ( TypeRND ) 
  use parameters
  implicit none
  integer       :: ios   
  character(20) :: TypeRND    
  character(80) :: string
  logical       :: XFindPar
!
  open (unit = iInput, file = FileInput, status = "old", position = "rewind")  
  TypeRND = 'SLAT'
  string  = 'RandomNumbers'
  if (XFindPar (iInput, string)) then
    read (iInput, *, iostat = ios) TypeRND
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable TypeRND;')"
      stop
    end if
  else
    print "(/,2x,'Group name TypeRND not found;')"
    stop  
  end if  
  call check_Random (TypeRND)   
  close (unit = iInput)   
end subroutine readinputRND   
! **********************************************************************************
subroutine SeqAddMet (TypeRND, Ntry, Npart, Rcirc, r, x, y, z)
  use parameters  
  implicit none
  integer       :: Ntry, Npart
  real(O)       :: r, Rcirc, x(Npart), y(Npart), z(Npart)
  character(20) :: TypeRND
!  
  integer       :: i, itry, j, idum
  real(O)       :: RANSLATEC, RANLAP, RANZIG, u, v, w, diam, d 
  logical       :: tryagain, overlap, out
  integer       :: iseed(4)
!  
  if (TypeRND(1:4) == 'LPCK') then      
    iseed(1) = 11
    iseed(2) = 22
    iseed(3) = 33
    iseed(4) = 44
  else if (TypeRND(1:4) == 'ZIGG') then
    idum = - 15    
  end if       
  diam  = 2._O * r
  do i = 1, Npart
    tryagain = .true.
    itry     = 0
    do while (tryagain)
      itry = itry + 1                                               
      if (TypeRND(1:4) == 'SLAT') then                                                            
        u = RANSLATEC (0._O)                   
        v = RANSLATEC (0._O)      
        w = RANSLATEC (0._O)           
      else if (TypeRND(1:4) == 'LPCK') then
        u = RANLAP (iseed)                  
        v = RANLAP (iseed)
        w = RANLAP (iseed)
      else if (TypeRND(1:4) == 'ZIGG') then
        u = RANZIG (idum)                   
        v = RANZIG (idum)
        w = RANZIG (idum)		
      end if	            	 	      
      x(i) = Rcirc * (2._O * u - 1._O)
      y(i) = Rcirc * (2._O * v - 1._O)
      z(i) = Rcirc * (2._O * w - 1._O)
      overlap = .false.        
      if (i > 1) then
        do j = 1, i - 1
          d = sqrt((x(i) - x(j))**2 + (y(i) - y(j))**2 + (z(i) - z(j))**2)
          if (d < diam) overlap = .true.             
        end do
      end if
      out = .false.
      d   = sqrt(x(i)**2 + y(i)**2 + z(i)**2)
      if (d > Rcirc) out = .true.
      if ((.not.out .and. .not.overlap) .or. (itry == Ntry)) tryagain = .false.
    end do  
    if ((out .or. overlap) .and. (itry == Ntry)) then
      print "(/,2x,'Error in subroutine SeqAddMet:')"
      print "(  2x,'the integer Ntry is too low for the desired configuration;')"     
      stop      
    end if
  end do
end subroutine SeqAddMet
! ***********************************************************************************
! *                          RANDOM NUMBER GENERATORS                               *
! *********************************************************************************** 
function RANSLATEC (r) result (rand)
!-----------------------------------------------------------------------------------
! Linear Congruential Generator (modified routine from Slatec library)             !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  real(O) :: r, rand
!  
  integer           :: iy0, iy1  
  integer,parameter :: ia1 = 1536, ia0 = 1029, ia1ma0 = 507, ic = 1731, ib = 2048
  integer,save      :: ix1 = 0,    ix0 = 0     
!        
  if (r == 0._O) then
    iy0 = ia0 * ix0
    iy1 = ia1 * ix1 + ia1ma0 * (ix0 - ix1) + iy0
    iy0 = iy0 + ic
    ix0 = mod (iy0, ib)
    iy1 = iy1 + (iy0 - ix0) / ib
    ix1 = mod (iy1, ib)        
  else if (r > 0._O) then
    ix1 = mod(r,1._O) * 4194304._O + 0.5_O
    ix0 = mod(ix1, ib)
    ix1 = (ix1 - ix0) / ib		
  end if      
  rand = ix1 * ib + ix0
  rand = rand / 4194304._O
end function RANSLATEC
! **********************************************************************************
function RANLAP (iseed) result (rand)
!-----------------------------------------------------------------------------------
! Modified routine from Lapack library                                             !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer  :: iseed(4)
  real(O)  :: rand
!
  integer,parameter :: m1 = 494, m2 = 322, m3 = 2508, m4 = 2549, ipw2 = 4096
  real(O),parameter :: r = 1._O / ipw2
  integer           :: it1, it2, it3, it4  
!        
  it4 = iseed(4) * m4
  it3 = it4 / ipw2
  it4 = it4 - ipw2 * it3
  it3 = it3 + iseed(3) * m4 + iseed(4) * m3
  it2 = it3 / ipw2
  it3 = it3 - ipw2 * it2
  it2 = it2 + iseed(2) * m4 + iseed(3) * m3 + iseed(4) * m2
  it1 = it2 / ipw2
  it2 = it2 - ipw2 * it1
  it1 = it1 + iseed(1) * m4 + iseed(2) * m3 + iseed(3) * m2 + iseed(4) * m1
  it1 = mod(it1,ipw2)
!
  iseed(1) = it1
  iseed(2) = it2
  iseed(3) = it3
  iseed(4) = it4
!
  rand = r * (real(it1,O) + r * (real(it2,O) + r * (real(it3,O) + r * (real(it4,O)))))
end function RANLAP      
! **********************************************************************************
function RANZIG (jsr) result (rand)
!-----------------------------------------------------------------------------------
! Ziggurat method (routine of Marsaglia and Tsang                                  !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer   :: jsr, jz, ival
  real(O)   :: rand

  jz   = jsr
  jsr  = IEOR( jsr, ISHFT( jsr,  13 ) )
  jsr  = IEOR( jsr, ISHFT( jsr, -17 ) )
  jsr  = IEOR( jsr, ISHFT( jsr,   5 ) )
  ival = jz + jsr
!
  rand = 0.5_O + 0.2328306e-9_O * ival  
end function RANZIG  
  
  
  
  
  
  
  
  
  
  
  
    



