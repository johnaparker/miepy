! **********************************************************************************
! *                ROUTINES FOR SOLVING LINEAR ALGEBRAIC EQUATIONS                 *
! *     -----------------------------------------------------------------------    *      
! *     Partial list of subroutines:                                               *
! *        inverse_matrix ,    LU_SYSTEM,    LU_SYSTEM_DIRECT,     LUDCMP,         *
! *        LUBKSB,             LU,           equilibrate,          LUDCMP1,        *
! *        LUBKSB1,            LU1,          BICGSTAB,             Precond,        *
! *        MatSolve,           readinputMatrSolv                                   *
! ********************************************************************************** 
subroutine inverse_matrix (a, nap, map, b, nbp, mbp, n)
!------------------------------------------------------------------------------------
! The routine computes the inverse of a matrix, i.e., B(n,n) = [A(n,n)]^(-1).       ! 
!                                                                                   !
! Input parameters:                                                                 !
! - a (complex array) - square matrix.                                              !
! - nap, map (integer variables) - specify the physical dimensions of the           ! 
!   matrix a.                                                                       !
! - nbp, mbp (integer variables) - specify the physical dimensions of the           !
!   inverse matrix.                                                                 !
! - n (integer) - real dimension of the square matrix.                              !
!                                                                                   !
! Output parameter:                                                                 !  
! - b (complex array) - inverse of the matrix a.                                    !
!                                                                                   !
! The model control parameters specified in the file "../INPUTFILES/Input.dat" are  !
! given below.                                                                      !
!                                                                                   !
! - TypeMatrSolv (character array) - specifies the method for solving linear        !
!   algebraic equations. The permissive values are:                                 !
!   - 'LU1'  - LU decomposition method from Numerical Recipes,                      !
!   - 'LU2'  - LU decomposition method from Lapack library,                         !
!   - 'LU3'  - LU decomposition method from Lapack library with solution            !
!      improvement,                                                                 !
!   - 'BICG' - Bi-Conjugate-Gradients method.                                       !
!   The recommended value of the parameter TypeMatrSolv is 'LU2'.                   !
!                                                                                   !
! - itmax (integer) - maximum number of iteration for the 'LU3' method or BICG      !
!   method.                                                                         !
!                                                                                   !
! - TypePrecond (character array) - specifies the type of preconditioning for the   !
!   BICG method. The permissive values are:                                         !
!   - 'NEUMANN' - incomplete Neumann series,                                        !
!   - 'SILU'    - incomplete LU factorization,                                      !
!   - ' '       - diagonal elements of the matrix.                                  !
!                                                                                   !
! - NPrecOrder (integer) - truncation order of the Neumann series. This parameter   !
!   is used if TypePrecond = 'NEUMANN'.                                              !
!                                                                                   !
! - epsilon (real) - solution tolerance for the BICG method.                        !
!------------------------------------------------------------------------------------     
  use parameters
  implicit none
  integer    :: n, nap, map, nbp, mbp
  complex(O) :: a(nap,map), b(nbp,mbp)
!                     
  integer    :: i, j, itmax, NPrecOrder  
  real(O)    :: d, epsilon
  integer,allocatable    :: indx(:)
  real(O),allocatable    :: row(:), column(:)
  complex(O),allocatable :: c(:)
  character(20)          :: TypeMatrSolv, TypePrecond
!        
  call readinputMatrSolv ( TypeMatrSolv, itmax, TypePrecond, NPrecOrder, epsilon )    
  if (TypeMatrSolv(1:3) == 'LU1') then
    call identity_matrix (n, b, nbp, mbp)
    allocate (indx(n), c(n), row(n), column(n))      
    call equilibrate (a, nap, map, n, n, row, column)
    call LUDCMP (a, nap, map, n, indx, d)
    do j = 1, n
      do i = 1, n
        c(i) = b(i,j) * row(i)
      end do 
      call LUBKSB (a, nap, map, n, indx, c)
      do i = 1, n
        b(i,j) = c(i) * column(i)
      end do 
    end do 
    deallocate (indx, c, row, column)      
  else if (TypeMatrSolv(1:3) == 'LU2') then 
    call identity_matrix (n, b, nbp, mbp)
    allocate (indx(n), row(n), column(n))
    call equilibrate (a, nap, map, n, n, row, column)
    call LUDCMP1 (a, nap, map, n, indx)
    do j = 1, n
      do i = 1, n
        b(i,j) = row(i) * b(i,j)
      end do
    end do
    call LUBKSB1 (a, nap, map, n, n, indx, b, nbp, mbp)
    do j = 1, n
      do i = 1, n
        b(i,j) = column(i) * b(i,j)
      end do
    end do
    deallocate (indx, row, column)                                                  
  else if (TypeMatrSolv(1:3) == 'LU3') then
    call identity_matrix (n, b, nbp, mbp)
    allocate (row(n), column(n))
    call equilibrate (a, nap, map, n, n, row, column)
    do j = 1, n
      do i = 1, n
        b(i,j) = row(i) * b(i,j)
      end do
    end do      
    call LU1 (a, nap, map, b, nbp, mbp, n, n, itmax)
    do j = 1, n
      do i = 1, n
        b(i,j) = column(i) * b(i,j)
      end do
    end do      
    deallocate (row, column)
  else if (TypeMatrSolv(1:4) == 'BICG') then   
    call identity_matrix (n, b, nbp, mbp)           
    call BICGSTAB (TypePrecond, NPrecOrder, itmax, epsilon, n, n, a, nap, map,      &
         b, nbp, mbp)
  end if  
end subroutine inverse_matrix 
! ***********************************************************************************
subroutine LU_SYSTEM (a, nap, map, b, nbp, mbp, n)
!------------------------------------------------------------------------------------
! Given the square matrices a and b, the routine computes                           !
!                      x(n,n) = b(n,n) * a(n,n)^(-1)                                !
! using the following scheme:                                                       !
! - transpose the matrices a and b, i.e. compute a(n,n)^T and b(n,n)^T,             !
! - compute y(n,n) = [a(n,n)^T]^(-1) * b(n,n)^T, that is, solve the matrix system   !
!   [a(n,n)^T] * y(n,n) = b(n,n)^T, with y being stored in b,                       !
! - transpose the matrix y, i.e. compute x(n,n)= y(n,n)^T                           ! 
!                                                                                   !
! Input parameters:                                                                 !
! - a (complex array) - square matrix.                                              !
! - nap, map (integer variables) - specify the physical dimensions of the           ! 
!   matrix a.                                                                       !
! - nbp, mbp (integer variables) - specify the physical dimensions of the           !
!   matrix b.                                                                       !
! - n (integer) - real dimension of the square matrices.                            !
!                                                                                   !
! Input/Output parameter:                                                           !  
! - b (complex array) - input matrix b and output matrix x.                         !
!                                                                                   !
! The model control parameters TypeMatrSolv, itmax, TypePrecond, NPrecOrder and     !
! epsilon are specified in the file "../INPUTFILES/Input.dat".                      !
!------------------------------------------------------------------------------------ 
  use parameters
  implicit none
  integer    :: n, nap, map, nbp, mbp
  complex(O) :: a(nap,map), b(nbp,mbp)
!
  integer    :: i, j, itmax, NPrecOrder  
  real(O)    :: d, epsilon
  integer,allocatable    :: indx(:)
  real(O),allocatable    :: row(:), column(:)
  complex(O),allocatable :: c(:)
  character(20)          :: TypeMatrSolv, TypePrecond
!    
  call readinputMatrSolv ( TypeMatrSolv, itmax, TypePrecond, NPrecOrder, epsilon )      
  if (TypeMatrSolv(1:3) == 'LU1') then
    allocate (indx(n), c(n), row(n), column(n))
    call transpose (a, nap, map, n)
    call transpose (b, nbp, mbp, n)
    call equilibrate (a, nap, map, n, n, row, column)                                                   
    call LUDCMP (a, nap, map, n, indx, d)
    do j = 1, n
      do i = 1, n
        c(i) = b(i,j) * row(i)
      end do 
      call LUBKSB (a, nap, map, n, indx, c)
      do i = 1, n
        b(i,j) = c(i) * column(i)
      end do 
    end do 
    call transpose (b, nbp, mbp, n)     
    deallocate (indx, c, row, column)          
  else if (TypeMatrSolv(1:3) == 'LU2') then
    allocate (indx(n), row(n), column(n))
    call transpose (a, nap, map, n)
    call transpose (b, nbp, mbp, n)
    call equilibrate (a, nap, map, n, n, row, column)
    call LUDCMP1 (a, nap, map, n, indx)
    do j = 1, n
      do i = 1, n
        b(i,j) = row(i) * b(i,j)
      end do
    end do
    call LUBKSB1 (a, nap, map, n, n, indx, b, nbp, mbp)
    do j = 1, n
      do i = 1, n
        b(i,j) = column(i) * b(i,j)
      end do
    end do   
    call transpose (b, nbp, mbp, n) 
    deallocate (indx, row, column) 
  else if (TypeMatrSolv(1:3) == 'LU3') then
    allocate (row(n), column(n))
    call transpose (a, nap, map, n)
    call transpose (b, nbp, mbp, n)
    call equilibrate (a, nap, map, n, n, row, column)
    do j = 1, n
      do i = 1, n
        b(i,j) = row(i) * b(i,j)
      end do
    end do      
    call LU1 (a, nap, map, b, nbp, mbp, n, n, itmax)
    do j = 1, n
      do i = 1, n
        b(i,j) = column(i) * b(i,j)
      end do
    end do   
    call transpose (b, nbp, mbp, n) 
    deallocate (row, column)
  else if (TypeMatrSolv(1:4) == 'BICG') then    
    call transpose (a, nap, map, n)
    call transpose (b, nbp, mbp, n)    
    call BICGSTAB (TypePrecond, NPrecOrder, itmax, epsilon, n, n, a, nap, map,      &
         b, nbp, mbp)                
    call transpose (b, nbp, mbp, n)                                     
  end if  
end subroutine LU_SYSTEM   
! ***********************************************************************************
subroutine LU_SYSTEM_DIRECT (a, nap, map, b, nbp, mbp, n, m)
!------------------------------------------------------------------------------------
! The routine computes                                                              !
!                    x(n,m) = [a(n,n)]^(-1) * b(n,m),                               !
! that is, solves the matrix system                                                 !
!                    a(n,n) * x(n,m) = b(n,m),                                      !
! with x being stored in b.                                                         !
!                                                                                   !
! Input parameters:                                                                 !
! - a (complex array) - square matrix.                                              !
! - nap, map (integer variables) - specify the physical dimensions of the           ! 
!   matrix a.                                                                       !
! - nbp, mbp (integer variables) - specify the physical dimensions of the           !
!   matrix b.                                                                       !
! - n (integer) - real dimension of the matrix a. n also specifies the number of    !
!                 lines of the matrix b.                                            !
! - m (integer) - specifies the number of columns of the matrix b.                  ! 
!                                                                                   !
! Input/Output parameter:                                                           !  
! - b (complex array) - input matrix b and output matrix x.                         !
!                                                                                   !
! The model control parameters TypeMatrSolv, itmax, TypePrecond, NPrecOrder and     !
! epsilon are specified in the file "../INPUTFILES/Input.dat".                      !
!------------------------------------------------------------------------------------      
  use parameters
  implicit none
  integer    :: n, m, nap, map, nbp, mbp
  complex(O) :: a(nap,map), b(nbp,mbp)
!             
  integer    :: i, j, itmax, NPrecOrder  
  real(O)    :: d, epsilon
  integer,allocatable    :: indx(:)
  real(O),allocatable    :: row(:), column(:)
  complex(O),allocatable :: c(:)
  character(20)          :: TypeMatrSolv, TypePrecond
!    
  call readinputMatrSolv ( TypeMatrSolv, itmax, TypePrecond, NPrecOrder, epsilon )
  if (TypeMatrSolv(1:3) == 'LU1') then
    allocate (indx(n), c(n), row(n), column(n))
    call equilibrate (a, nap, map, n, n, row, column)
    call LUDCMP (a, nap, map, n, indx, d)       
    do j = 1, m
      do i = 1, n
        c(i) = b(i,j) * row(i)
      end do
      call LUBKSB (a, nap, map, n, indx, c)
      do i = 1, n
        b(i,j) = c(i) * column(i)
      end do
    end do
    deallocate (indx, c, row, column)                  
  else if (TypeMatrSolv(1:3) == 'LU2') then
    allocate (indx(n), row(n), column(n))
    call equilibrate (a, nap, map, n, n, row, column)
    call LUDCMP1 (a, nap, map, n, indx)
    do j = 1, m
      do i = 1, n
        b(i,j) = row(i) * b(i,j)
      end do
    end do
    call LUBKSB1 (a, nap, map, n, m, indx, b, nbp, mbp)
    do j = 1, m
      do i = 1, n
        b(i,j) = column(i) * b(i,j)
      end do
    end do
    deallocate (indx, row, column)                        
  else if (TypeMatrSolv(1:3) == 'LU3') then
    allocate (row(n), column(n))
    call equilibrate (a, nap, map, n, n, row, column)
    do j = 1, m
      do i = 1, n
        b(i,j) = row(i) * b(i,j)
      end do
    end do      
    call LU1 (a, nap, map, b, nbp, mbp, n, m, itmax)
    do j = 1, m
      do i = 1, n
        b(i,j) = column(i) * b(i,j)
      end do
    end do      
    deallocate (row, column)
  else if (TypeMatrSolv(1:4) == 'BICG') then                
  call BICGSTAB (TypePrecond, NPrecOrder, itmax, epsilon, n, m, a, nap, map,      &
       b, nbp, mbp)                                                                   
  end if      
end subroutine LU_SYSTEM_DIRECT
! ***********************************************************************************
subroutine readinputMatrSolv ( TypeMatrSolv, itmax, TypePrecond, NPrecOrder, epsilon ) 
  use parameters
  implicit none                     
  integer       :: ios, itmax, NPrecOrder  
  real(O)       :: epsilon
  character(20) :: TypeMatrSolv, TypePrecond
  character(80) :: string
  logical       :: XFindPar
!
  open (unit = iInput, file = FileInput, status = "old", position = "rewind")  
  TypeMatrSolv = 'LU1' 
  itmax        = 10
  TypePrecond  = 'NEUMANN'
  NPrecOrder   = 6
  epsilon      = 1.e-6_O
  string       = 'MatrixSolver'
  if (XFindPar (iInput, string)) then
    read (iInput, *, iostat = ios) TypeMatrSolv
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable TypeMatrSolv;')"
      stop
    end if
    read (iInput, *, iostat = ios) itmax
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable itmax;')"
      stop
    end if
    read (iInput, *, iostat = ios) TypePrecond
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable TypePrecond;')"
      stop
    end if
    read (iInput, *, iostat = ios) NPrecOrder
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable NPrecOrder;')"
      stop
    end if 
    read (iInput, *, iostat = ios) epsilon
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsilon;')"
      stop
    end if
  else
    print "(/,2x,'Group name MatrixSolver not found;')"
    stop  
  end if  
  call check_MatrixSolver (TypeMatrSolv)
  close (unit = iInput)  
end subroutine readinputMatrSolv  
! ************************************************************************************
! *               LU DECOMPOSITION - F90-VERSION FROM NUMERICAL RECIPES              *
! ************************************************************************************
subroutine LUDCMP (a, np, mp, n, indx, d)
  use parameters
  use derived_parameters
  implicit none
  integer    :: n, np, mp, indx(n)
  real(O)    :: d
  complex(O) :: a(np,mp)
!      
  integer    :: imax, i, j, k 
  real(O)    :: aamax, tampon 
  complex(O) :: sum, dum
  complex(O),allocatable :: vv(:)
!               
  allocate (vv(n))
  d = 1._O
  do i = 1, n
    aamax = 0._O
    do j = 1, n
      if (abs(a(i,j)) > aamax) aamax = abs(a(i,j))
    end do 
    if (aamax < SmallestPosNumber) then       
      print "(/,2x,'Error in subroutine LUDCMP in module MatrixSolv: singularity,')" 
      print "(  2x,'the matrix inversion can not be performed with the LU ')" 
      print "(  2x,'algorithm from Numerical Recipes;')"      
      stop
    end if 
    tampon = 1._O / aamax
    vv(i)  = cmplx(tampon,0.0,O)
  end do 
  do j = 1, n
    do i = 1, j - 1 
      sum = a(i,j)    
      do k = 1, i - 1
        sum = sum - a(i,k) * a(k,j)
      end do      
      a(i,j) = sum
    end do 
    aamax = 0._O
    do i = j, n 
      sum = a(i,j)     
      do k = 1, j - 1
        sum = sum - a(i,k) * a(k,j)
      end do      
      a(i,j) = sum
      dum    = vv(i) * abs(sum)
      if (abs(dum) >= aamax) then
        imax  = i
        aamax = abs(dum)
      end if
    end do 
    if (j /= imax) then
      do k = 1, n
        dum       = a(imax,k)
        a(imax,k) = a(j,k)
        a(j,k)    = dum
      end do 
      d = - d
      vv(imax) = vv(j)
    end if
    indx(j) = imax
    if (abs(a(j,j)) < SmallestPosNumber) a(j,j) = cmplx(ZeroLUVal,ZeroLUVal,O)
    if (j /= n) then
      dum = 1._O / a(j,j)
      do  i = j + 1, n
        a(i,j) = a(i,j) * dum
      end do 
    end if
  end do 
  deallocate (vv)      
end subroutine LUDCMP
! **********************************************************************************
subroutine LUBKSB (a, np, mp, n, indx, b)
  use parameters
  implicit none
  integer    :: n, np, mp, indx(n)
  complex(O) :: a(np,mp), b(n)
!      
  integer    :: i, j, ii, ll 
  complex(O) :: sum
!
  ii = 0
  do i = 1, n
    ll    = indx(i)
    sum   = b(ll)
    b(ll) = b(i)
    if (ii /= 0) then     
      do j = ii, i - 1
        sum = sum - a(i,j) * b(j)
      end do           
    else if (abs(sum) /= 0._O) then
      ii = i
    end if
    b(i) = sum
  end  do 
  do i = n, 1, -1
    sum = b(i)
    if (i < n) then
      do j = i + 1, n
        sum = sum - a(i,j) * b(j)
      end do          
    end if
    b(i) = sum / a(i,i)
  end do         
end subroutine LUBKSB
! **********************************************************************************
subroutine LU (a, np, mp, b, n)
  use parameters
  implicit none
  integer    :: n, np, mp
  complex(O) :: a(np,mp), b(n)
!
  real(O)    :: d
  integer,allocatable :: indx(:)
!
  allocate (indx(n))
  call LUDCMP (a, np, mp, n, indx, d)
  call LUBKSB (a, np, mp, n, indx, b)
  deallocate (indx)      
end subroutine LU
! ***********************************************************************************
subroutine equilibrate (a, np, mp, n, m, r, c)
! ----------------------------------------------------------------------------------
! This subroutine computes row and column scalings intended to equilibrate         !
! a matrix a and reduce its condition number. The equilibrated matrix anew is      !
!                   anew = r * a * c,                                              !
! where r and c are real diagonal matrices. Therefore                              !
!                   a = r^(-1) * anew * c^(-1).                                    !
! ----------------------------------------------------------------------------------        
  use parameters
  use derived_parameters
  implicit none
  integer    :: np, mp, n, m
  real(O)    :: r(m), c(n)
  complex(O) :: a(np,mp)
!
  integer    :: i, j
  real(O)    :: smlnum, bignum, norm, rcmin, rcmax, amax, rowcnd, colcnd, cj, tresh   
!
  smlnum = sqrt(SmallestPosNumber)
  bignum = 1._O / smlnum
  tresh  = 1.e-2_O
  do i = 1, m
    r(i) = 0._O
  end do
  do j = 1, n 
    do i = 1, m
      norm = abs(real(a(i,j),O)) + abs(aimag(a(i,j)))
      r(i) = max(r(i), norm)
    end do
  end do    
  rcmin = bignum 
  rcmax = 0._O 
  do i = 1, m
    rcmax = max(rcmax, r(i))
    rcmin = min(rcmin, r(i))
  end do
  amax  = rcmax
  if (rcmin == 0._O) then
    do i = 1, m
      if (r(i) == 0._O) then
        print "(/,2x,'Error in subroutine equilibrate in module MatrixSolv:')" 
        print "(  2x,'the row ',i5,' of the matrix is exactly zero;')", i             
        stop
      end if
    end do   
  else
    do i = 1, m
      r(i) = 1._O / min(max(r(i), smlnum), bignum)
    end do
    rowcnd = max(rcmin, smlnum) / min(rcmax, bignum)
  end if
!
  do j = 1, n
    c(j) = 0._O
  end do
  do j = 1, n
    do i = 1, m
      norm = abs(real(a(i,j),O)) + abs(aimag(a(i,j)))
      c(j) = max(c(j), norm * r(i))
    end do
  end do
  rcmin = bignum 
  rcmax = 0._O 
  do j = 1, n
    rcmax = max(rcmax, c(j))
    rcmin = min(rcmin, c(j))
  end do
  if (rcmin == 0._O) then
    do j = 1, n
      if (c(j) == 0._O) then
        print "(/,2x,'Error in subroutine equilibrate in module MatrixSolv:')" 
        print "(  2x,'the column ',i5,' of the matrix is exactly zero;')", j                  
        stop        
      end if
    end do   
  else
    do j = 1, n
      c(j) = 1._O / min(max(c(j), smlnum), bignum)
    end do
    colcnd = max(rcmin, smlnum) / min(rcmax, bignum)
  end if    
!
  if (rowcnd >= tresh .and. amax >= smlnum .and. amax <= bignum) then
   if (colcnd >= tresh) then
     do i = 1, m
       r(i) = 1._O
     end do
     do j = 1, n
       c(j) = 1._O
     end do      
   else    
     do j = 1, n
       cj = c(j)
       do i = 1, m
         a(i,j) = a(i,j) * cj
       end do
     end do
     do i = 1, m
       r(i) = 1._O
     end do         
   end if   
 else if (colcnd >= tresh) then  
   do j = 1, n
     do i = 1, m
       a(i,j) = r(i) * a(i,j)
     end do
   end do
   do j = 1, n
     c(j) = 1._O
   end do   
 else   
   do j = 1, n
     cj = c(j)
     do i = 1, m
       a(i,j)= r(i) * a(i,j) * cj
     end do
   end do   
 end if 
end subroutine equilibrate    
! ***********************************************************************************
! *            LU DECOMPOSITION - SIMPLIFIED VERSIONS FROM LAPACK LIBRARY           *
! ***********************************************************************************
subroutine LUDCMP1 (a, nap, map, n, ipiv)
  use parameters
  use derived_parameters
  implicit none
  integer    :: n, nap, map, ipiv(n)
  complex(O) :: a(nap,map)
!
  integer    :: j, nt, imax, jp, i, k, l
  real(O)    :: amax
  complex(O) :: temp    
!       
  do j = 1, n
    nt   = n - j + 1
    imax = 1
    if (nt > 1) then
      amax = abs(a(j,j))
      do i = 2, nt                      
        if (abs(a(j+i-1,j)) > amax) then
          imax = i
          amax = abs(a(j+i-1,j))
        end if  
      end do
    end if
    nt = n - j              
    jp = j - 1 + imax
    ipiv(j) = jp
    if (abs(a(jp,j)) > SmallestPosNumber) then
      if (jp /= j) then
        do i = 1, n
          temp    = a(j,i)
          a(j,i)  = a(jp,i)
          a(jp,i) = temp
        end do
      end if
      if (j < n) then
        temp = one / a(j,j)
        do i = 1,nt
          a(j+i,j) = temp * a(j+i,j)
        end do
      end if
    else          
      print "(/,2x,'Error in subroutine LUDCMP1 in module MatrixSolv: singularity,')" 
      print "(  2x,'the matrix inversion can not be performed with the LU ')" 
      print "(  2x,'algorithm from Lapack library;')"         
      stop                  
    end if                       
    if (j < n) then
      do k = 1, nt
        if (a(j,j+k) /= zero) then
          temp = - one * a(j,j+k)
          do l = 1, nt
            a(j+l,j+k) = a(j+l,j+k) + a(j+l,j) * temp
          end do
        end if
      end do 
    end if
  end do
end subroutine LUDCMP1
! **********************************************************************************
subroutine LUBKSB1 (a, nap, map, n, m, ipiv, b, nbp, mbp)
  use parameters
  implicit none
  integer    :: n, m, nap, map, nbp, mbp, ipiv(n)
  complex(O) :: a(nap,map), b(nbp,mbp)
!
  integer    :: i, ip, k, j
  complex(O) :: temp    
!            
  do i = 1, n
    ip = ipiv(i)
    if (ip /= i) then
      do k = 1, m
        temp    = b(i,k)
        b(i,k)  = b(ip,k)
        b(ip,k) = temp
      end do
    end if
  end do                
  do j = 1, m
    do k = 1, n
      if (b(k,j) /= zero) then
        do i = k + 1, n
          b(i,j) = b(i,j) - b(k,j) * a(i,k)
        end do
      end if
    end do
  end do        
  do j = 1, m
    do k = n, 1, - 1
      if (b(k,j) /= zero ) then
        b(k,j) = b(k,j) / a(k,k)
        do i = 1, k - 1
          b(i,j) = b(i,j) - b(k,j) * a(i,k)
        end do
      end if
    end do
  end do      
end subroutine LUBKSB1 
! **********************************************************************************
subroutine LU1 (a, nap, map, b, nbp, mbp, n, m, itmax)
  use parameters
  use derived_parameters
  implicit none
  integer    :: n, m, nap, map, nbp, mbp, itmax
  complex(O) :: a(nap,map), b(nbp,mbp)
!
  integer    :: i, j, k, nz, count
  real(O)    :: eps, safemin, safe1, safe2, lstres, s, absb, absx, absa, absw
  complex(O) :: sum
  logical    :: more
  integer,allocatable    :: ipiv(:)
  real(O),allocatable    :: rwork(:), berr(:) 
  complex(O),allocatable :: af(:,:), x(:,:), work(:,:)  
!
  allocate (ipiv(nap), af(nap,map), x(nbp,mbp), work(nbp,1), rwork(nap), berr(mbp))
  do i = 1, n
    do j = 1, n
      af(i,j) = a(i,j)
    end do
    do j = 1, m
      x(i,j) = b(i,j)
    end do
  end do
  call LUDCMP1 (af, nap, map, n, ipiv)
  call LUBKSB1 (af, nap, map, n, m, ipiv, x, nbp, mbp)
!
  nz  = n + 1
  eps = MachEps
  safemin = SmallestPosNumber
  safe1   = nz * safemin
  safe2   = safe1 / eps 
  do j = 1, m
    count  = 1
    lstres = 3._O
    more   = .true.
    do while (more)
      do i = 1, n
        work(i,1) = b(i,j)
      end do
      do i = 1, n
        sum = zero
        do k = 1, n
          sum = sum + a(i,k) * x(k,j)
        end do
        work(i,1) = work(i,1) - sum
      end do
      do i = 1, n
        absb = abs(real(b(i,j))) + abs(aimag(b(i,j)))
        rwork(i) = absb
      end do
      do k = 1, n
        absx = abs(real(x(k,j))) + abs(aimag(x(k,j)))
        do i = 1, n
          absa = abs(real(a(i,k))) + abs(aimag(a(i,k)))
          rwork(i) = rwork(i) + absa * absx
        end do
      end do
      s = 0._O
      do i = 1, n
        absw = abs(real(work(i,1))) + abs(aimag(work(i,1)))          
        if (rwork(i) > safe2 ) then
          s = max(s,absw / rwork(i))
        else                                      
          s = max(s,(absw + safe1) / (rwork(i) + safe1))
        end if
      end do
      berr(j) = s
      more = .false.
      if (berr(j) > eps .and. 2._O * berr(j) <= lstres .and. count <= itmax) then
        more = .true.
        call LUBKSB1 (af, nap, map, n, 1, ipiv, work, nbp, 1)                                   
        do i = 1, n
          x(i,j) = x(i,j) + work(i,1)
        end do
        lstres = berr(j)
        count  = count + 1        
      end if
    end do                                              
  end do              
  do i = 1, n
    do j = 1, m
      b(i,j) = x(i,j)
    end do
  end do
  deallocate (ipiv, af, x, work, rwork, berr)      
end subroutine LU1
! ***********************************************************************************
! *                                BICG Method                                      *
! ***********************************************************************************
subroutine BICGSTAB (TypePrecond, NPrecOrder, Niter, eps, n, m, a, nap, map,         &
           b, nbp, mbp)
  use parameters
  implicit none
  integer       :: NPrecOrder, Niter, n, m, nap, map, nbp, mbp
  real (O)      :: eps 
  complex(O)    :: a(nap,map), b(nbp,mbp)
  character(20) :: TypePrecond
  
  integer       :: i, j, k
  real(O)       :: error, normvct, normy, invnormy, normt, normt2
  complex(O)    :: alpha, beta, ro, roold, omega, cdot, vr0, minus_one, beta_omega,  &
                   roold_omega
  logical       :: more
  complex(O),allocatable :: r(:), r0(:), p(:), v(:), s(:), t(:), x(:), y(:),         &
                            q1(:), q2(:), pold(:), xold(:), dx(:), pm(:,:), pm0(:)
!
  minus_one = - one  
  allocate (r(n), r0(n), p(n), v(n), s(n), t(n), q1(n), q2(n), pold(n), xold(n), dx(n))
  allocate (x(n), y(n)) 
  if (TypePrecond(1:7) == 'NEUMANN' .or. TypePrecond(1:4) == 'SILU') then 
    allocate (pm(n,n)) 
    call Precond (TypePrecond, NPrecOrder, a, nap, map, pm, n, n, n) 
  else
    allocate (pm0(n))
    do i = 1, n
      pm0(i) = one / a(i,i)
    end do
  end if
  do j = 1, m
    do i = 1, n
      y(i) = b(i,j)
    end do
    normy = normvct (n, y)
    if (normy == 0._O) then
      do i = 1, n
        x(i) = zero
      end do
      go to 1
    end if
    invnormy = 1._O / normy
    call copy_vector (y, xold, n) 
    call product_matrix_vector3 (n, n, a, nap, map, xold, y, minus_one, one, r0) 
    if (normvct (n, r0) == 0._O) then
      call copy_vector (xold, x, n)
      go to 1
    end if
    call copy_vector (r0, r, n)
    more = .true.
    k    = 0
    do while (more) 
      k  = k + 1
      ro = cdot (n, r, r0)      
      if (abs(ro) == 0._O) then   
        call copy_vector (r0, r, n)   
        ro = cdot (n, r, r0)      
      end if  
      if (k == 1) then
        call copy_vector (r, p, n)  
      else
        roold_omega = roold * omega
        beta = ro * alpha / roold_omega      
        beta_omega  = - beta * omega        
        do i = 1, n
          p(i) = r(i) + beta * pold(i) + beta_omega * v(i)
        end do
      end if
      if (TypePrecond(1:7) == 'NEUMANN') then
        call product_matrix_vector  (n, n, pm, n, n, p, q1)  
      else if (TypePrecond(1:4) == 'SILU') then  
        call MatSolve (n, a, nap, map, pm, n, n, p, q1)
      else          
        do i = 1, n
          q1(i) = p(i) * pm0(i)
        end do
      end if 
      call product_matrix_vector (n, n, a, nap, map, q1, v)
      vr0 = cdot (n, v, r0)
      if (abs(vr0) /= 0._O) then
        alpha = ro / vr0
      else
        print "(/,2x,'Error in subroutine BICGSTAB: singularity,')" 
        print "(  2x,'BICG algorithm can not be used for matrix inversion;')"   
        stop 
      end if                 
      do i = 1, n
        s(i) = r(i) - alpha * v(i)
      end do
      error = normvct (n,s) * invnormy 
      if (error < eps) then             
        do i = 1, n
          x(i) = xold(i) + alpha * q1(i)
        end do
        more = .false.
      else
        if (TypePrecond(1:7) == 'NEUMANN') then
          call product_matrix_vector  (n, n, pm, n, n, s, q2)  
        else if (TypePrecond(1:4) == 'SILU') then  
          call MatSolve (n, a, nap, map, pm, n, n, s, q2)
        else           
          do i = 1, n
            q2(i) = s(i) * pm0(i)
          end do
        end if
        call product_matrix_vector (n, n, a, nap, map, q2, t)
        normt = normvct (n,t)
        if (normt /= 0._O) then
          normt2 = normt * normt
          omega  = cdot (n, s, t) / normt2
        else
          print "(/,2x,'Error in subroutine BICGSTAB: singularity,')"   
          print "(  2x,'BICG algorithm can not be used for matrix inversion;')" 
          stop 
        end if               
        do i = 1, n
          x(i) = xold(i) + alpha * q1(i) + omega * q2(i)
        end do                      
        do i = 1, n
          r(i) = s(i) - omega * t(i)
        end do         
        do i = 1, n
          dx(i) = x(i) - xold(i)
        end do   
        error = normvct (n, dx) * invnormy      
        if (error < eps .or. k == Niter) then
          more = .false.          
        else
          call copy_vector (x, xold, n)  
          call copy_vector (p, pold, n)  
          roold = ro
          if (abs(roold) == 0._O .or. abs(omega) == 0._O) then
            print "(/,2x,'Error in subroutine BICGSTAB: singularity,')"
            print "(  2x,'BICG algorithm can not be used for matrix inversion;')"    
            stop  
          end if
        end if
      end if  
    end do    
    if (k == Niter .and. error > eps) then
      print "(/,2x, a)",                                                             &
     'Error in subroutine BICGSTAB:  convergence is not achieved because'
      print "(  2x, a)",                                                             &
     'the tolerance TolConv or the iteration number NmaxIter are too low;'
      stop
    end if
1   continue
    do i = 1, n
      b(i,j) = x(i)
    end do             
  end do
  deallocate (r, r0, p, v, s, t, q1, q2, pold, xold, dx, x, y)
  if (TypePrecond(1:7) == 'NEUMANN' .or. TypePrecond(1:4) == 'SILU') then 
    deallocate (pm) 
  else
    deallocate (pm0)
  end if
end subroutine BICGSTAB    
!************************************************************************************
subroutine Precond (TypePrecond, Norder, a, nap, map, p, npp, mpp, n)
  use parameters
  implicit none
  integer       :: Norder, n, nap, map, npp, mpp
  complex(O)    :: a(nap,map), p(npp,mpp)
  character(20) :: TypePrecond
!  
  integer       :: iorder, i, j, k, imax   
  real(O)       :: maxl, absa
  complex(O)    :: sum, Yik, IP  
  complex(O),allocatable :: c(:), d(:)
!  
  call identity_matrix (n, p, npp, mpp)
  if (TypePrecond(1:7) == 'NEUMANN') then    
    allocate (c(n), d(n))  
    do j = 1, n
      maxl = 0._O
      do i = 1, n
        absa = abs(a(i,j))
        if (absa > maxl) then
          maxl = absa
          imax = i
        end if
      end do      
      if (abs(a(imax,j)) /= 0._O) then
        d(j) = one / a(imax,j) / real(n * n,O)    
      else           
        print "(/,2x,'Error in subroutine Precond: singularity,')" 
        print "(  2x,'Neumann preconditioning can not be used;')"   
        stop
      end if    
    end do                       
    do iorder = 1, Norder
      do j = 1, n
        do i = 1, n
          sum = zero
          do k = 1, n
            if (k == i) then
              Yik = one - a(i,k) * d(k)
            else
              Yik = - a(i,k) * d(k) 
            end if           
            sum = sum + Yik * p(k,j)
          end do
          c(i) = sum
        end do
        do i = 1, n
          p(i,j) = c(i)
          if (i == j) p(i,j) = p(i,j) + one
        end do        
      end do
    end do  
    do i = 1, n
      do j = 1, n
        p(i,j) = p(i,j) * d(i)
      end do
    end do    
    deallocate (c, d)   
  else if (TypePrecond(1:4) == 'SILU') then 
    do i = 1, n
      sum = zero
      do j = 1, i - 1
        sum = sum + a(i,j) * P(j,j) * a(j,i)
      end do
      IP = a(i,i) - sum
      if (abs(IP) /= 0._O) then
        P(i,i) = one / IP
      else
        print "(/,2x,'Error in subroutine Precond: singularity,')"
        print "(  2x,'SILU preconditioning can not be used;')"    
        stop
      end if    
    end do                   
  end if    
end subroutine Precond     
!************************************************************************************
subroutine MatSolve (n, a, nap, map, p, npp, mpp, y, x)
  use parameters
  implicit none
  integer     :: n, nap, map, npp, mpp  
  complex(O)  :: a(nap,map), p(npp,mpp), y(n), x(n)
!  
  integer     :: i, j
  complex(O)  :: sum 
  complex(O),allocatable :: z(:)
!  
  allocate (z(n))
  z(1) = y(1) * p(1,1)
  do i = 2, n
    sum = zero
    do j = 1, i - 1
      sum = sum + a(i,j) * z(j)
    end do
    z(i) = ( y(i) - sum ) * p(i,i)
  end do  
  x(n) = z(n)
  do i = n - 1, 1, -1
    sum = zero
    do j = i + 1, n
      sum = sum + a(i,j) * x(j)
    end do
    x(i) = z(i) - sum * p(i,i)
  end do
  deallocate (z)
end subroutine MatSolve




 
