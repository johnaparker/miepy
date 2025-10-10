! **********************************************************************************
! *                         ELEMENTARY MATRIX OPERATIONS                           *
! *    -----------------------------------------------------------------------     *
! *    Partial list of subroutines:                                                *
! *      determinant,                          transpose,                          *
! *      product_vector_vector,                product_vector_vector1,             *
! *      product_vector_vector2,               product_vector_vector_AVplus,       *
! *      product_vector_vector_AVminus,        product_matrix_vector,              *
! *      product_matrix_vector1,               product_matrix_vector2,             *
! *      product_matrix_vector3,               product_matrix_vector4,             *
! *      product_matrices,                     product_matrices1,                  *
! *      product_matrices2,                    product_matrixSUB,                  *
! *      sum_matrices,                         sum_diagonal_elements,              *
! *      sum_vectors,                          copy_matrix,                        *
! *      copy_vector,                          minus_matrix,                       *
! *      matrix_m_negativ,                     identity_matrix,                    *
! *      zero_matrix,                          cdot,                               *
! *      normvct                                                                   *
! **********************************************************************************
subroutine determinant (a, np, mp, n, det)
  use parameters
  implicit none
  integer    :: n, np, mp
  complex(O) :: a(np,mp), det
!
  integer    :: i
  real(O)    :: d
  integer,allocatable :: indx(:)
!
  allocate (indx(n))
  call LUDCMP (a, np, mp, n, indx, d)
  det = one
  do i = 1, n
    det = det * a(i,i)
  end do
  det = d * det
  deallocate (indx)
end subroutine determinant
! **********************************************************************************
subroutine transpose (a, np, mp, n)
  use parameters
  implicit none
  integer    :: n, np, mp
  complex(O) :: a(np,mp)
!      
  integer    :: i, j
  complex(O) :: tampon
!
  do i = 1, n - 1
    do j = i + 1, n
      tampon = a(i,j)
      a(i,j) = a(j,i)
      a(j,i) = tampon
    end do
  end do     
end subroutine transpose
! **********************************************************************************
subroutine product_vector_vector (N, a, b, c)
  use parameters
  implicit none
  integer    :: N, i
  complex(O) :: a(N), b(N), c(N)
!
  do i = 1, N
    c(i) = a(i) * b(i)
  end do      
end subroutine product_vector_vector 
! ***********************************************************************************
subroutine product_vector_vector1 (n, a, b, alpha, beta)
!------------------------------------------------------------------------------------
! The routine computes a = alpha * a + beta * b.                                    !
!------------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: n, i
  real(O)    :: absa, absb
  complex(O) :: a(*), b(*), alpha, beta 
!
  absa = abs(alpha)
  absb = abs(beta)  
  if (absa == 0._O) then
    if (absb == 0._O) then
      do i = 1, n
        a(i) = zero
      end do  
    else
      do i = 1, n
        a(i) = beta * b(i)
      end do  
    end if
  else
    if (absb == 0._O) then
      do i = 1, n
        a(i) = alpha * a(i)
      end do  
    else
      do i = 1, n
        a(i) = alpha * a(i) + beta * b(i)
      end do  
    end if                    
  end if         
end subroutine product_vector_vector1 
! ***********************************************************************************
subroutine product_vector_vector2 (n, a, b, c, d, alpha, beta, gamma)
!------------------------------------------------------------------------------------
! The routine computes d = alpha * a + beta * b + gama * c.                         !
!------------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: n, i
  real(O)    :: absa, absb, absg
  complex(O) :: a(*), b(*), c(*), d(*), alpha, beta, gamma 
!
  absa = abs(alpha)
  absb = abs(beta)
  absg = abs(gamma)  
  if (absa == 0._O) then
    if (absb == 0._O) then
      if (absg == 0._O) then
        do i = 1, n
          d(i) = zero
        end do  
      else
        do i = 1, n
          d(i) = gamma * c(i)
        end do  
      end if                              
    else
      if (absg == 0._O) then
        do i = 1, n
          d(i) = beta * b(i)
        end do  
      else
        do i = 1, n
          d(i) = beta * b(i) + gamma * c(i)
        end do  
      end if                            
    end if
  else
    if (absb == 0._O) then
      if (absg == 0._O) then
        do i = 1, n
          d(i) = alpha * a(i)
        end do  
      else
        do i = 1, n
          d(i) = alpha * a(i) + gamma * c(i)
        end do  
      end if                              
    else
      if (absg == 0._O) then
        do i = 1, n
          d(i) = alpha * a(i) + beta * b(i)
        end do  
      else
        do i = 1, n
          d(i) = alpha * a(i) + beta * b(i) + gamma * c(i)
        end do  
      end if                            
    end if                            
  end if          
end subroutine product_vector_vector2   
!***********************************************************************************
subroutine product_vector_vector_AVplus (Nc, Nmax, TV, nTV, c, c1)
  use parameters
  implicit none
  integer    :: Nc, Nmax, NTV, i, j
  complex(O) :: TV(nTV), c(2*Nmax), c1(2*Nmax), sum
!
  do i = 1, 2*Nmax
    sum = zero
    do j = 1, 2*Nmax
      sum = sum + TV(Nc+i+2*(j-1)*Nmax) * c(j)
    end do
    c1(i) = sum
  end do    
end subroutine product_vector_vector_AVplus
!***********************************************************************************
subroutine product_vector_vector_AVminus (Nc, Nmax, TV, nTV, c, c1)
  use parameters
  implicit none
  integer    :: Nc, Nmax, NTV, i, j
  complex(O) :: TV(nTV), c(2*Nmax), c1(2*Nmax), sum
!
  do i = 1, Nmax
    sum = zero
    do j = 1, Nmax
      sum = sum + TV(Nc+i+2*(j-1)*Nmax) * c(j)
    end do
    c1(i) = sum
    sum   = zero
    do j = Nmax + 1, 2*Nmax
      sum = sum + TV(Nc+i+2*(j-1)*Nmax) * c(j)
    end do
    c1(i) = c1(i) - sum         
  end do
  do i = Nmax + 1, 2*Nmax
    sum = zero
    do j = 1, Nmax
      sum = sum + TV(Nc+i+2*(j-1)*Nmax) * c(j)
    end do
    c1(i) = - sum
    sum   =   zero
    do j = Nmax + 1, 2 * Nmax
      sum = sum + TV(Nc+i+2*(j-1)*Nmax) * c(j)
    end do
    c1(i) = c1(i) + sum         
  end do      
end subroutine product_vector_vector_AVminus                                 
! **********************************************************************************
subroutine product_matrix_vector (N, M, a, np, mp, b, c)
!-----------------------------------------------------------------------------------     
! The routine computes c(N) = a(N,M) * b(M).                                       !
!-----------------------------------------------------------------------------------      
  use parameters
  implicit none
  integer    :: N, M, np, mp, j, i
  complex(O) :: a(np,mp), b(M), c(N)
!  
  do i = 1, N
    c(i) = zero
  end do
  do j = 1, M
    if (b(j) /= zero) then
      do i = 1, N
        c(i) = c(i) + a(i,j) * b(j)
      end do
    end if
  end do               
end subroutine product_matrix_vector
! **********************************************************************************
subroutine product_matrix_vector1 (N, M, c, a, nap, map, b, nbp, mbp)
!-----------------------------------------------------------------------------------
! The routine computes b(i,j) = c(i) * a(i,j) for i = 1,..,N; j = 1,...,M.         !  
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: N, M, nap, map, nbp, mbp, j, i
  complex(O) :: a(nap,map), b(nbp,mbp), c(N)
!
  do i = 1, N
    do j = 1, M
      b(i,j) = c(i) * a(i,j)
    end do
  end do      
end subroutine product_matrix_vector1
! **********************************************************************************
subroutine product_matrix_vector2 (N, M, c, a, np, mp)
!-----------------------------------------------------------------------------------      
! The routine computes a(i,j) = c(i) * a(i,j) for i = 1,..,N; j = 1,...,M.         !
!-----------------------------------------------------------------------------------      
  use parameters
  implicit none
  integer    :: N, M, np, mp, j, i
  complex(O) :: a(np,mp), c(N)
!
  do i = 1, N
    do j = 1, M
      a(i,j) = c(i) * a(i,j)
    end do
  end do      
end subroutine product_matrix_vector2
! ***********************************************************************************
subroutine product_matrix_vector3 (n, m, a, nap, map, u, v, alpha, beta, x)
!-----------------------------------------------------------------------------------
! The routine computes x = alpha * A * u + beta * v.                               !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: n, m, nap, map
  complex(O) :: a(nap,map), u(*), v(*), x(*), alpha, beta
!  
  integer    :: i, j
  complex(O) :: sum
!  
  do i = 1, n
    sum = zero
    do j = 1, m
      if (u(j) /= zero) sum = sum + a(i,j) * u(j)
    end do            
    x(i) = alpha * sum + beta * v(i)
  end do  
end subroutine  product_matrix_vector3  
! ***********************************************************************************
subroutine product_matrix_vector4 (n, m, a, nap, map, u, v, alpha, beta)
!-----------------------------------------------------------------------------------
! The routine computes v = alpha * A * u + beta * v.                               !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: n, m, nap, map
  complex(O) :: a(nap,map), u(*), v(*), alpha, beta
!  
  integer    :: i, j
  complex(O) :: sum
!  
  do i = 1, n
    sum = zero
    do j = 1, m
      if (u(j) /= zero) sum = sum + a(i,j) * u(j)
    end do
    v(i) = alpha * sum + beta * v(i)
  end do
end subroutine  product_matrix_vector4 
!***********************************************************************************
subroutine product_matrices (N, M, L, a, nap, map, b, nbp, mbp)
!-----------------------------------------------------------------------------------      
! The routine computes c(N,L) = a(N,M) * b(M,L), where c(N,L) is stored in         !
! a(nap,map). The natural conditions: nap >= N, map >= M  and the additional       !
! condition map >= L must be satisfied.                                            !
!----------------------------------------------------------------------------------- 
  use parameters
  implicit none       
  integer    :: N, M, L, nap, map, nbp, mbp, j, i, k
  complex(O) :: a(nap,map), b(nbp,mbp)
  complex(O), allocatable :: c(:)
!
  if (map < L) then    
    print "(/,2x,'Error in subroutine product_matrices:')"
    print "(  2x,'the real dimension exceeds the physical dimension of the matrix;')"
    stop
  end if
  allocate (c(L))
  do i = 1, N
    do j = 1, L
      c(j) = zero
    end do
    do k = 1, M
      if (a(i,k) /= zero) then
        do j = 1, L
          c(j) = c(j) + a(i,k) * b(k,j)
        end do
      end if
    end do
    do j = 1, L
      a(i,j) = c(j)
    end do      
  end do  
  deallocate (c)     
end subroutine product_matrices
!***********************************************************************************
subroutine product_matrices1 (N, M, L, a, nap, map, b, nbp, mbp, c, ncp, mcp)
!-----------------------------------------------------------------------------------      
! The routine computes c(N,L) = a(N,M) * b(M,L).                                   !
!-----------------------------------------------------------------------------------      
  use parameters
  implicit none 
  integer    :: N, M, L, nap, map, nbp, mbp, ncp, mcp, j, i, k
  complex(O) :: a(nap,map), b(nbp,mbp), c(ncp,mcp)
!
  do i = 1, N
    do j = 1, L
          c(i,j) = zero
    end do
    do k = 1, M
      if (a(i,k) /= zero ) then
        do j = 1, L
          c(i,j) = c(i,j) + a(i,k) * b(k,j)
        end do
      end if
    end do              
  end do    
end subroutine product_matrices1   
! **********************************************************************************
subroutine product_matrices2 (index, N, M, L, a, nap, map, b, nbp, mbp, c, ncp, mcp)
!-----------------------------------------------------------------------------------
! The routine computes c(N,L) = c(N,L) + a(N,M) * b(M,L).                          !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none 
  integer    :: N, M, L, j, i, k, index, nap, map, nbp, mbp, ncp, mcp
  complex(O) :: a(nap,map), b(nbp,mbp), c(ncp,mcp)
!
  if (index == 1) then
    do i = 1, N
      do j = 1, L
        c(i,j) = zero
      end do
    end do
  end if        
  do i = 1, N
    do k = 1, M
      if (a(i,k) /= zero) then
        do j = 1, L
          c(i,j) = c(i,j) + a(i,k) * b(k,j)
        end do
      end if
    end do          
  end do  
end subroutine product_matrices2    
! **********************************************************************************
subroutine product_matrixSUB (a, na, ma, b, nb, mb, n)
!-----------------------------------------------------------------------------------
! The routine computes a(n,n) = I(n,n) - a(n,n) * b(n,n).                          !
!----------------------------------------------------------------------------------- 
  use parameters
  implicit none
  integer    :: n, na, ma, nb, mb
  complex(O) :: a(na,ma), b(nb,mb)
!
  integer    :: i, j, k
  complex(O) :: sum       
  complex(O),allocatable :: c(:)      
!
  allocate (c(n))
  do i = 1, n
    do j = 1, n
      sum = zero
      do k = 1, n
        sum = sum + a(i,k) * b(k,j)
      end do 
      c(j) = sum
    end do
    do j = 1, n
      a(i,j) = c(j)
    end do 
  end do
  do i = 1, n
    do j = 1, n
      if (j == i) then
        a(i,j) = one - a(i,j)
      else
        a(i,j) = - a(i,j)
      end if
    end do
  end do        
  deallocate (c)      
end subroutine product_matrixSUB         
! **********************************************************************************
subroutine sum_matrices (N, M, a, nap, map, b, nbp, mbp)
!-----------------------------------------------------------------------------------
! The routine computes a(N,M) = a(N,M) + b(N,M).                                   !
!-----------------------------------------------------------------------------------      
  use parameters
  implicit none
  integer    :: N, M, nap, map, nbp, mbp, j, i
  complex(O) :: a(nap,map), b(nbp,mbp)
!
  do i = 1, N
    do j = 1, M
      a(i,j) = a(i,j) + b(i,j)           
    end do        
  end do      
end subroutine sum_matrices
! **********************************************************************************
subroutine sum_diagonal_elements (N, c, a, nap, map)
!-----------------------------------------------------------------------------------       
! The routine computes a(i,i) = c(i) + a(i,i) for i = 1,..,N.                      !
!-----------------------------------------------------------------------------------      
  use parameters
  implicit none
  integer    :: N, nap, map, i
  complex(O) :: c(N), a(nap,map)
!
  do i = 1, N
    a(i,i) = c(i) + a(i,i)           
  end do    
end subroutine sum_diagonal_elements
! **********************************************************************************
subroutine sum_vectors (a, b, n)
  use parameters
  implicit none
  integer    :: n, i
  complex(O) :: a(n), b(n)
!
  do i = 1, n
    a(i) = a(i) + b(i)
  end do
end subroutine sum_vectors   
! **********************************************************************************
subroutine copy_matrix (N, M, a, nap, map, b, nbp, mbp)
  use parameters
  implicit none
  integer    :: N, M, nap, map, nbp, mbp, i, j
  complex(O) :: a(nap,map), b(nbp,mbp)
!
  do i = 1, N
    do j = 1, M
      b(i,j) = a(i,j)
    end do         
  end do      
end subroutine copy_matrix  
! **********************************************************************************
subroutine copy_vector (c, c1, N)
  use parameters
  implicit none
  integer    :: N, i
  complex(O) :: c(N), c1(N)      
!
  do i = 1, N
    c1(i) = c(i)
  end do      
end subroutine copy_vector       
! **********************************************************************************
subroutine minus_matrix (N, M, a, nap, map)
!-----------------------------------------------------------------------------------      
! The routine computes a(N,M) = - a(N,M).                                          !
!-----------------------------------------------------------------------------------     
  use parameters
  implicit none
  integer    :: N, M, nap, map, j, i
  complex(O) :: a(nap,map)
!
  do i = 1, N
    do j = 1, M
      a(i,j) = - a(i,j)
    end do        
  end do      
end subroutine minus_matrix
! **********************************************************************************
subroutine matrix_m_negativ (N, M, a, nap, map)
!-----------------------------------------------------------------------------------
! The routine transforms the matrix a(2*N,2*M) as follows:                         !
!                                                                                  !
!             | a11(N,M)   a12(N,M) |      | a11(N,M)  -a12(N,M) |                 !
!             |                     |  --> |                     |.                !
!             | a21(N,M)   a22(N,M) |      |-a21(N,M)   a22(N,M) |                 !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: N, M, nap, map, i, j
  complex(O) :: a(2*nap,2*map)
!
  do i = 1, N
    do j = 1, M
      a(i,j+M) = - a(i,j+M)
      a(i+N,j) = - a(i+N,j)
    end do         
  end do      
end subroutine matrix_m_negativ
!***********************************************************************************
subroutine identity_matrix (N, a, nap, map)
  use parameters
  implicit none 
  integer    :: N, nap, map, i, j
  complex(O) :: a(nap,map)
!
  do i = 1, N
    do j = 1, N
      a(i,j) = zero
      if (j == i) a(i,j) = one
    end do        
  end do      
end subroutine identity_matrix
! **********************************************************************************
subroutine zero_matrix (N, M, a, nap, map)
  use parameters
  implicit none 
  integer    :: N, M, nap, map, i, j
  complex(O) :: a(nap,map)   
!
  do i = 1, N
    do j = 1, M
      a(i,j) = zero                           
    end do        
  end do      
end subroutine zero_matrix
!************************************************************************************
function cdot (n, a, b) result (prod)
  use parameters
  implicit none
  integer      :: n, i
  complex(O)   :: a(*), b(*), prod, sum
    
  sum = zero
  do i = 1, n    
    sum = sum + a(i) * conjg(b(i))
  end do
  prod = sum
end function cdot
! ***********************************************************************************   
function normvct (n, cx) result (norm)
  use parameters
  implicit none       
  integer     :: n, i
  real(O)     :: norm, abscx, sum
  complex(O)  :: cx(*)
!                
  sum = 0._O
  do i = 1, n         
    abscx = abs(cx(i))      
    sum   = sum + abscx * abscx
  end do       
  norm = sqrt(sum)       
end function normvct          
