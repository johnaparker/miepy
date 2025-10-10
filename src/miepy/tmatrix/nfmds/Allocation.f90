module allocation
  use parameters
  implicit none 
!  
  integer,allocatable,save       :: Nsurf(:), Nparam(:), Nrankp(:), Nrankp1(:)
  real(O),allocatable,save       :: surf(:,:), zpart(:), zRe(:,:), zIm(:,:),        &
                                    zRe1(:,:), zIm1(:,:), lnorm(:), EpsZReIm(:)			       			       			          
  complex(O),allocatable,save    :: ind_ref(:)
  logical,allocatable,save       :: ComplexPlane(:)
!
  integer,allocatable,save       :: Mrankcs(:), Nrankcs(:)
  real(O),allocatable,save       :: Rcs(:)
!
  integer,allocatable,save       :: Mrankp(:)
  real(O),allocatable,save       :: xp(:), yp(:), zp(:),                            &
                                    alphap(:), betap(:), gammap(:)
  character(80),allocatable,save :: FileTmatp(:)
  logical,allocatable,save       :: axsymp(:), chiralp(:)
!
  real(O),allocatable,save       :: rp(:)
  complex(O),allocatable,save    :: ind_refp(:)  
!
  integer,allocatable,save       :: IndI(:), IndJ(:)
  character(2),allocatable,save  :: NameElem(:)
end module allocation  
