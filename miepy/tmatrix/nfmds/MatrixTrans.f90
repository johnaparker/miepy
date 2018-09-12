! **********************************************************************************
! *                           MATRIX TRANSFORMATIONS                               *
! *    ------------------------------------------------------------------------    *                           
! *    Partial list of subroutines:                                                *
! *      form_Tvector,            extend_matrix1,          extend_matrix2,         *
! *      extend_matrix3,          extend_matrix4,          extract_matrix1,        *
! *      extract_matrix2,         extract_matrix3,         dimension_column,       *
! *      dimension_row,           matrix_Nrank_m,          matrix_Nrank_m_left,    *
! *      vector_Nrank_m,          matrix_Nrank_1_right,    matrix_Mrank_1_right,   *
! *      matrix_Nrank_1_left,     matrix_Mrank_1_left,     vector_Nrank_1,         *
! *      vector_Mrank_1,          internal_matrix_Nrank_m, incident_matrix_Nrank_m,*
! *      difusion_matrix_Nrank_m, internal_matrix_LAY_Nrank_m                      *
! ********************************************************************************** 
subroutine form_Tvector (Nrank, Nmax, m, tvg, tv)
  use parameters
  implicit none
  integer    :: Nrank, Nmax, m, k, n
  complex(O) :: tvg(2*Nrank), tv(2*Nmax)
!
  if (m == 0) then
    do k = 1, Nmax
      tv(k)      = tvg(k)
      tv(k+Nmax) = tvg(k+Nrank)      
    end do
  else
    do k = 1, Nmax
      n = m + k - 1                                                                                              
      tv(k)      = tvg(n)
      tv(k+Nmax) = tvg(n+Nrank)         
    end do
  endif   
end subroutine form_Tvector
!*********************************************************************************** 
subroutine extend_matrix1 (ipart, jpart, Npart, Nmaxp, Nl, Nc, a, nap, map, aa,     &
           naap, maap)
!-----------------------------------------------------------------------------------      
! Extend the matrix a(2*Nmaxp(ipart),2*Nmaxp(jpart)) into the global matrix        !
! aa(2*Nmaxmax,2*Nmaxmax).                                                         !
!-----------------------------------------------------------------------------------      
  use parameters
  implicit none 
  integer    :: ipart, jpart, Npart, Nmaxp(Npart), Nl, Nc, nap, map, naap, maap, i, j
  complex(O) :: a(2*nap,2*map), aa(2*naap,2*maap)  
!
  do i = 1, 2*Nmaxp(ipart)
    do j = 1, 2*Nmaxp(jpart)
      aa(i+Nl,j+Nc) = - a(i,j)
    end do        
  end do            
end subroutine extend_matrix1
!***********************************************************************************
subroutine extend_matrix2 (ipart, Npart, Nmaxp, Nmax, Nmaxmax, Nl, a, nap, map, aa, &
           naap, maap)
!-----------------------------------------------------------------------------------     
! Extend the matrix a(2*Nmaxp(ipart),2*Nmax) into the global matrix                !
! aa(2*Nmaxmax,2*Nmax).                                                            !
!-----------------------------------------------------------------------------------    
  use parameters
  implicit none 
  integer    :: ipart, Npart, Nmaxp(Npart), Nmax, Nmaxmax, Nl, nap, map, naap,      &
                maap, i, j  
  complex(O) :: a(2*nap,2*map), aa(2*naap,2*maap)      
!
  if (ipart == 1) then
    do i = 1, 2*Nmaxmax
      do j = 1, 2*Nmax
        aa(i,j) = zero
      end do        
    end do
  end if
  do i = 1, 2*Nmaxp(ipart)
    do j = 1, 2*Nmax
      aa(i+Nl,j) = a(i,j)
    end do        
  end do
end subroutine extend_matrix2
!***********************************************************************************
subroutine extend_matrix3 (Nmax, Nmax1, Nl, Nc, a, nap, map, aa, naap, maap)
!-----------------------------------------------------------------------------------
! Extend the matrix a(2*Nmax,2*Nmax1) into the global matrix                       !
! aa(2*Nmaxmax,2*Nmaxmax).                                                         !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none 
  integer    :: Nmax, Nmax1, Nl, Nc, nap, map, naap, maap, i, j
  complex(O) :: a(2*nap,2*map), aa(2*naap,2*maap)      
!
  do i = 1, 2*Nmax
    do j = 1, 2*Nmax1
      aa(i+Nl,j+Nc) = - a(i,j)
    end do        
  end do            
end subroutine extend_matrix3
!***********************************************************************************
subroutine extend_matrix4 (ipart, Nmax, Nmax1, Nmaxmax, Nl, a, nap, map, aa,        &
           naap, maap)
!-----------------------------------------------------------------------------------
! Extend the matrix a(2*Nmax,2*Nmax1) into the global matrix aa(2*Nmaxmax,2*Nmax1).!      
!-----------------------------------------------------------------------------------
  use parameters
  implicit none 
  integer    :: ipart, Nmax, Nmax1, Nmaxmax, Nl, nap, map, naap, maap, i, j
  complex(O) :: a(2*nap,2*map), aa(2*naap,2*maap)    
!
  if (ipart == 1) then
    do i = 1, 2*Nmaxmax
      do j = 1, 2*Nmax1
        aa(i,j) = zero
      end do        
    end do
  end if
  do i = 1, 2*Nmax
    do j = 1, 2*Nmax1
      aa(i+Nl,j) = a(i,j)
    end do        
  end do  
end subroutine extend_matrix4
!***********************************************************************************
subroutine extract_matrix1 (ipart, Npart, Nmaxp, Nmax, Nl, a, nap, map, aa,         &
           naap, maap)
!-----------------------------------------------------------------------------------
! Extract the matrix a(2*Nmaxp(ipart),2*Nmax) from aa(2*Nmaxmax,2*Nmax).           !
!-----------------------------------------------------------------------------------   
  use parameters
  implicit none 
  integer    :: ipart, Npart, Nmaxp(Npart), Nmax, nap, map, naap, maap, Nl, i, j
  complex(O) :: a(2*nap,2*map), aa(2*naap,2*maap)
!
  do i = 1, 2*Nmaxp(ipart)
    do j = 1, 2*Nmax
      a(i,j) = aa(i+Nl,j)
    end do        
  end do     
end subroutine extract_matrix1
!***********************************************************************************
subroutine extract_matrix2 (Nmax, Nmax1, Nl, a, nap, map, aa, naap, maap)
!-----------------------------------------------------------------------------------
! Extract the matrix a(2*Nmax,2*Nmax1) from aa(2*Nmaxmax,2*Nmax1).                 !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none 
  integer    :: Nmax, Nmax1, Nl, nap, map, naap, maap, i, j
  complex(O) :: a(2*nap,2*map), aa(2*naap,2*maap)
!
  do i = 1, 2*Nmax
    do j = 1, 2*Nmax1
      a(i,j) = aa(i+Nl,j)
    end do        
  end do     
end subroutine extract_matrix2
! **********************************************************************************
subroutine extract_matrix3 (index, Nmaxpmax, aa, naap, maap, a, nap, map)
!-----------------------------------------------------------------------------------
! Extract the matrix aa(2*Nmax,2*Nmax) from a(2*Nmaxpmax,2*Nmaxpmax).              !
!-----------------------------------------------------------------------------------
  use parameters
  implicit none
  integer    :: index, Nmaxpmax, naap, maap, nap, map, i, j
  complex(O) :: aa(2*naap,2*maap), a(2*nap,2*map)
!
  if (index == 1) then
    do i = 1, 2*Nmaxpmax
      do j = 1, 2*Nmaxpmax
        a(i,j) = aa(i,j)
      end do
    end do      
  else
    do i = 1, 2*Nmaxpmax
      do j = 1, 2*Nmaxpmax
        a(i,j) = aa(i+2*Nmaxpmax,j)
      end do
    end do      
  end if
end subroutine extract_matrix3 
!***********************************************************************************
function dimension_column (m, Npart, Nrankp) result (NmaxC)
  implicit none
  integer :: m, Npart, Nrankp(Npart), NmaxC, Nstart, p, Nmax
!
  Nstart = 0
  do p = 1, Npart
    if (Nrankp(p) >= abs(m)) then
      if (abs(m) == 0) then
        Nmax = Nrankp(p)
      else
        Nmax = Nrankp(p) - abs(m) + 1
      end if
      Nstart = Nstart + Nmax
    end if
  end do
  NmaxC = Nstart 
end function dimension_column
! **********************************************************************************
function dimension_row (m, Nrank) result (Nmax)
  implicit none
  integer :: m, Nrank, Nmax
!
  if (abs(m) == 0) then
    Nmax = Nrank
  else
    Nmax = Nrank - abs(m) + 1
  end if      
end function dimension_row      
!***********************************************************************************
subroutine matrix_Nrank_m (Nmax, a, nap, map)
  use parameters
  implicit none
  integer    :: Nmax, nap, map, i
  complex(O) :: a(2*nap,2*map)
!
  do i = 1, Nmax
    a(i,Nmax)        = zero
    a(Nmax,i)        = zero
    a(i,2*Nmax)      = zero
    a(Nmax,i+Nmax)   = zero
    a(i+Nmax,Nmax)   = zero
    a(2*Nmax,i)      = zero 
    a(i+Nmax,2*Nmax) = zero
    a(2*Nmax,i+Nmax) = zero
  end do             
end subroutine matrix_Nrank_m 
!***********************************************************************************
subroutine matrix_Nrank_m_left (Nmax, a, nap, map)
  use parameters
  implicit none
  integer    :: Nmax, nap, map, i
  complex(O) :: a(2*nap,2*map)
!
  do i = 1, Nmax
    a(i,Nmax)        = zero
    a(Nmax,i)        = zero
    a(i,2*Nmax)      = zero
    a(Nmax,i+Nmax)   = zero
    a(i+Nmax,Nmax)   = zero
    a(2*Nmax,i)      = zero
    a(i+Nmax,2*Nmax) = zero
    a(2*Nmax,i+Nmax) = zero
  end do
  a(Nmax,Nmax) = one
  a(Nmax,2*Nmax)   =   one
  a(2*Nmax,Nmax)   =   one
  a(2*Nmax,2*Nmax) = - one      
end subroutine matrix_Nrank_m_left
!***********************************************************************************
subroutine vector_Nrank_m (Nmax, c)
  use parameters
  implicit none
  integer    :: Nmax
  complex(O) :: c(2*Nmax)
!
  c(Nmax)  = zero
  c(2*Nmax)= zero      
end subroutine vector_Nrank_m
!***********************************************************************************
subroutine matrix_Nrank_1_right (Mrank, Nrank, Nmax, a, nap, map)
  use parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax, nap, map, m, i, l, N0
  complex(O) :: a(2*nap,2*map)
!
  do m = 0, Mrank
    if (m == 0) then
      do i = 1, Nmax
        a(i,Nrank)           = zero
        a(Nrank,i)           = zero
        a(i,Nrank+Nmax)      = zero
        a(Nrank,i+Nmax)      = zero
        a(i+Nmax,Nrank)      = zero
        a(Nrank+Nmax,i)      = zero
        a(i+Nmax,Nrank+Nmax) = zero
        a(Nrank+Nmax,i+Nmax) = zero
      end do
    else          
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      do l = 1, 2
        do i = 1, Nmax
          a(i,N0+Nrank-m+1)           = zero
          a(N0+Nrank-m+1,i)           = zero
          a(i,N0+Nrank-m+1+Nmax)      = zero
          a(N0+Nrank-m+1,i+Nmax)      = zero
          a(i+Nmax,N0+Nrank-m+1)      = zero
          a(N0+Nrank-m+1+Nmax,i)      = zero
          a(i+Nmax,N0+Nrank-m+1+Nmax) = zero
          a(N0+Nrank-m+1+Nmax,i+Nmax) = zero
        end do
        N0 = N0 + Nrank - m + 1
      end do              
    end if
  end do   
end subroutine matrix_Nrank_1_right
!***********************************************************************************
subroutine matrix_Mrank_1_right (Mrank, Nrank, Nmax, a, nap, map)
  use parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax, nap, map, m, i, j, l, N0
  complex(O) :: a(2*nap,2*map)
!
  m  = Mrank
  N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
  do l = 1, 2           
    do i = 1, Nmax
      do j = 1, Nrank - m + 1
        a(i,j+N0)           = zero
        a(j+N0,i)           = zero
        a(i,j+N0+Nmax)      = zero
        a(j+N0,i+Nmax)      = zero
        a(i+Nmax,j+N0)      = zero
        a(j+N0+Nmax,i)      = zero
        a(i+Nmax,j+N0+Nmax) = zero
        a(j+N0+Nmax,i+Nmax) = zero
      end do
    end do
    N0 = N0 + Nrank - m + 1
  end do                   
end subroutine matrix_Mrank_1_right
!***********************************************************************************
subroutine matrix_Nrank_1_left (Mrank, Nrank, Nmax, a, nap, map)
  use parameters
  implicit none 
  integer    :: Mrank, Nrank, Nmax, nap, map, m, i, l, N0
  complex(O) :: a(2*nap,2*map)   
!
  do m = 0, Mrank
    if (m == 0) then
      do i = 1, Nmax
        a(i,Nrank)           = zero
        a(Nrank,i)           = zero
        a(i,Nrank+Nmax)      = zero
        a(Nrank,i+Nmax)      = zero
        a(i+Nmax,Nrank)      = zero
        a(Nrank+Nmax,i)      = zero
        a(i+Nmax,Nrank+Nmax) = zero
        a(Nrank+Nmax,i+Nmax) = zero
      end do
    else      
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      do l = 1, 2
        do i = 1, Nmax
          a(i,N0+Nrank-m+1)           = zero
          a(N0+Nrank-m+1,i)           = zero
          a(i,N0+Nrank-m+1+Nmax)      = zero
          a(N0+Nrank-m+1,i+Nmax)      = zero
          a(i+Nmax,N0+Nrank-m+1)      = zero
          a(N0+Nrank-m+1+Nmax,i)      = zero
          a(i+Nmax,N0+Nrank-m+1+Nmax) = zero
          a(N0+Nrank-m+1+Nmax,i+Nmax) = zero
        end do
        N0 = N0 + Nrank - m + 1
      end do              
    end if
  end do
  do m = 0, Mrank
    if (m == 0) then
      a(Nrank,Nrank)           =   one
      a(Nrank,Nrank+Nmax)      =   one
      a(Nrank+Nmax,Nrank)      =   one
      a(Nrank+Nmax,Nrank+Nmax) = - one
    else      
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      do l = 1, 2
        a(N0+Nrank-m+1,N0+Nrank-m+1)           =   one
        a(N0+Nrank-m+1,N0+Nrank-m+1+Nmax)      =   one
        a(N0+Nrank-m+1+Nmax,N0+Nrank-m+1)      =   one
        a(N0+Nrank-m+1+Nmax,N0+Nrank-m+1+Nmax) = - one
        N0 = N0 + Nrank - m + 1
      end do              
    end if
  end do
end subroutine matrix_Nrank_1_left
!***********************************************************************************
subroutine matrix_Mrank_1_left (Mrank, Nrank, Nmax, a, nap, map)
  use parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax, nap, map, m, i, j, l, N0
  complex(O) :: a(2*nap,2*map)   
!
  m  = Mrank
  N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
  do l = 1, 2           
    do i = 1, Nmax
      do j = 1, Nrank - m + 1
        a(i,j+N0)           = zero
        a(j+N0,i)           = zero
        a(i,j+N0+Nmax)      = zero
        a(j+N0,i+Nmax)      = zero
        a(i+Nmax,j+N0)      = zero
        a(j+N0+Nmax,i)      = zero
        a(i+Nmax,j+N0+Nmax) = zero
        a(j+N0+Nmax,i+Nmax) = zero
      end do
    end do
    N0 = N0 + Nrank - m + 1
  end do              
  N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
  do l = 1, 2
    do j = 1, Nrank - m + 1
      a(j+N0,j+N0)           =   one
      a(j+N0,j+N0+Nmax)      =   one
      a(j+N0+Nmax,j+N0)      =   one
      a(j+N0+Nmax,j+N0+Nmax) = - one
    end do
    N0 = N0 + Nrank - m + 1
  end do                   
end subroutine matrix_Mrank_1_left
!***********************************************************************************
subroutine vector_Nrank_1 (Mrank, Nrank, Nmax, c)
  use parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax, m, l, N0
  complex(O) :: c(2*Nmax)
!
  do m = 0, Mrank
    if (m==0) then        
      c(Nrank)      = zero
      c(Nrank+Nmax) = zero
    else      
      N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
      do l = 1, 2      
        c(N0+Nrank-m+1)      = zero
        c(N0+Nrank-m+1+Nmax) = zero
        N0 = N0 + Nrank - m + 1         
      end do              
    end if
  end do     
end subroutine vector_Nrank_1
! **********************************************************************************
subroutine vector_Mrank_1 (Mrank, Nrank, Nmax, c)
  use parameters
  implicit none
  integer    :: Mrank, Nrank, Nmax, m, j, l, N0
  complex(O) :: c(2*Nmax)
!
  m  = Mrank
  N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
  do l = 1, 2
    do j = 1, Nrank - m + 1
      c(j+N0)      = zero
      c(j+N0+Nmax) = zero
    end do
    N0 = N0 + Nrank - m + 1
  end do              
end subroutine vector_Mrank_1
!***********************************************************************************
subroutine internal_matrix_Nrank_m (m, Npart, Nrankp, NmaxC, a, nap, map)
  use parameters
  implicit none
  integer    :: m, Npart, Nrankp(Npart), NmaxC, nap, map
  complex(O) :: a(2*nap,2*map)
!       
  integer    :: p, i, Nmaxp, Nstart
!
  Nstart = 0
  do p = 1, Npart
    if (m == 0) then
      Nmaxp = Nrankp(p)
    else
      Nmaxp = Nrankp(p) - abs(m) + 1
    end if
    do i = 1, 2*NmaxC
      a(i,Nmaxp+Nstart)       = zero
      a(Nmaxp+Nstart,i)       = zero
      a(i,Nmaxp+Nstart+NmaxC) = zero
      a(Nmaxp+Nstart+NmaxC,i) = zero
    end do
    a(Nmaxp+Nstart,Nmaxp+Nstart)             =   one
    a(Nmaxp+Nstart,Nmaxp+Nstart+NmaxC)       =   one
    a(Nmaxp+Nstart+NmaxC,Nmaxp+Nstart)       =   one
    a(Nmaxp+Nstart+NmaxC,Nmaxp+Nstart+NmaxC) = - one
    Nstart = Nstart + Nmaxp
  end do
end subroutine internal_matrix_Nrank_m
!***********************************************************************************
subroutine incident_matrix_Nrank_m (m, Npart, Nrankp, NmaxC, NmaxL, a, nap, map)
  use parameters
  implicit none
  integer    :: m, Npart, Nrankp(Npart), NmaxC, NmaxL, nap, map
  complex(O) :: a(2*nap,2*map)
!       
  integer    :: p, i, j, Nmaxp, Nstart
!
  Nstart = 0
  do p = 1, Npart
    if (m == 0) then
      Nmaxp = Nrankp(p)
    else
      Nmaxp = Nrankp(p) - abs(m) + 1
    end if
    do j = 1, 2*NmaxL
      a(Nmaxp+Nstart,j)       = zero
      a(Nmaxp+Nstart+NmaxC,j) = zero
    end do	
    Nstart = Nstart + Nmaxp        
  end do
  do i = 1, 2*NmaxC
    a(i,NmaxL)       = zero
    a(i,NmaxL+NmaxL) = zero
  end do
end subroutine incident_matrix_Nrank_m
!***********************************************************************************
subroutine difusion_matrix_Nrank_m (m, Npart, Nrankp, NmaxC, NmaxL, a, nap, map)
  use parameters
  implicit none
  integer    :: m, Npart, Nrankp(Npart), NmaxC, NmaxL, nap, map
  complex(O) :: a(2*nap,2*map)
!       
  integer    :: p, i, j, Nmaxp, Nstart
!
  Nstart = 0
  do p = 1, Npart
    if (m == 0) then
      Nmaxp = Nrankp(p)
    else
      Nmaxp = Nrankp(p) - abs(m) + 1
    end if
    do i = 1, 2*NmaxL
      a(i,Nmaxp+Nstart)       = zero
      a(i,Nmaxp+Nstart+NmaxC) = zero
    end do		         
    Nstart = Nstart + Nmaxp
  end do
  do j = 1, 2*NmaxC
    a(NmaxL,j)       = zero
    a(NmaxL+NmaxL,j) = zero
  end do 
end subroutine difusion_matrix_Nrank_m  
!***********************************************************************************
subroutine internal_matrix_LAY_Nrank_m (Npart, Nmaxp, aa, nap, map)
  use parameters
  implicit none
  integer    :: Npart, Nmaxp(Npart), nap, map      
  complex(O) :: aa(2*nap,2*map)
!
  integer    :: i, j, ipart, Nrankl, Nrankc, Nl, Nc
!
  do ipart = 1, Npart     
    if (ipart == 1) then
      Nrankl = Nmaxp(ipart)
      Nrankc = Nmaxp(ipart)           
      do j = 1, 2*Nrankc
        aa(  Nrankl,j) = zero
        aa(2*Nrankl,j) = zero
      end do
      do i = 1, 2*Nrankl
        aa(i,  Nrankc) = zero
        aa(i,2*Nrankc) = zero
        aa(i,3*Nrankc) = zero
        aa(i,4*Nrankc) = zero
      end do
      aa(  Nrankl,  Nrankc) =   one
      aa(  Nrankl,2*Nrankc) =   one
      aa(2*Nrankl,  Nrankc) =   one
      aa(2*Nrankl,2*Nrankc) = - one        
      Nl = 2 * Nmaxp(ipart)
      Nc = 0        
    else    
      Nrankl = Nmaxp(ipart-1)
      Nrankc = Nmaxp(ipart-1) 
      do j = 1, 2*Nrankc
        aa(  Nrankl+Nl,j+Nc) = zero
        aa(2*Nrankl+Nl,j+Nc) = zero
      end do
      do i = 1, 2*Nrankl
        aa(i+Nl,  Nrankc+Nc) = zero
        aa(i+Nl,2*Nrankc+Nc) = zero
        aa(i+Nl,3*Nrankc+Nc) = zero
        aa(i+Nl,4*Nrankc+Nc) = zero
      end do
      aa(  Nrankl+Nl,3*Nrankc+Nc) =   one
      aa(  Nrankl+Nl,4*Nrankc+Nc) =   one
      aa(2*Nrankl+Nl,3*Nrankc+Nc) =   one
      aa(2*Nrankl+Nl,4*Nrankc+Nc) = - one  
      Nc = Nc + 4 * Nmaxp(ipart-1)
      Nrankl = Nmaxp(ipart-1)
      Nrankc = Nmaxp(ipart) 
      do j = 1, 2*Nrankc
        aa(  Nrankl+Nl,j+Nc) = zero
        aa(2*Nrankl+Nl,j+Nc) = zero
      end do
      do i = 1, 2*Nrankl
        aa(i+Nl,  Nrankc+Nc) = zero
        aa(i+Nl,2*Nrankc+Nc) = zero
        if (ipart < Npart) then
          aa(i+Nl,3*Nrankc+Nc) = zero
          aa(i+Nl,4*Nrankc+Nc) = zero
        end if
      end do                     
      Nl = Nl + 2 * Nmaxp(ipart-1)
      Nc = Nc - 4 * Nmaxp(ipart-1)    
      Nrankl = Nmaxp(ipart)
      Nrankc = Nmaxp(ipart-1) 
      do j = 1, 2*Nrankc
        aa(  Nrankl+Nl,j+Nc) = zero
        aa(2*Nrankl+Nl,j+Nc) = zero
      end do
      do i = 1, 2*Nrankl
        aa(i+Nl,  Nrankc+Nc) = zero
        aa(i+Nl,2*Nrankc+Nc) = zero
        aa(i+Nl,3*Nrankc+Nc) = zero
        aa(i+Nl,4*Nrankc+Nc) = zero
      end do 
      Nc = Nc + 4 * Nmaxp(ipart-1)    
      Nrankl = Nmaxp(ipart)
      Nrankc = Nmaxp(ipart) 
      do j = 1, 2*Nrankc
        aa(  Nrankl+Nl,j+Nc) = zero
        aa(2*Nrankl+Nl,j+Nc) = zero
      end do
      do i = 1, 2*Nrankl
        aa(i+Nl,  Nrankc+Nc) = zero
        aa(i+Nl,2*Nrankc+Nc) = zero
        if (ipart < Npart) then
          aa(i+Nl,3*Nrankc+Nc) = zero
          aa(i+Nl,4*Nrankc+Nc) = zero
        end if
      end do  
      aa(  Nrankl+Nl,  Nrankc+Nc) =   one
      aa(  Nrankl+Nl,2*Nrankc+Nc) =   one
      aa(2*Nrankl+Nl,  Nrankc+Nc) =   one
      aa(2*Nrankl+Nl,2*Nrankc+Nc) = - one             
      Nl = Nl + 2 * Nmaxp(ipart)              
    end if
  end do   
end subroutine internal_matrix_LAY_Nrank_m   
