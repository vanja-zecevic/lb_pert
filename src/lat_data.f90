! Copyright (C) 2012 Vanja Zecevic
! Contact vanja.zecevic@sydney.uni.edu.au

! This file is part of lb_pert

! lb_pert is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! lb_pert is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE lat_data
IMPLICIT NONE
CONTAINS
!------------------------------------------------------------------------------
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!------------------------------------------------------------------------------
SUBROUTINE get_sizes_1 (d, q, ng)

INTEGER :: d, q, ng

d  = 2
q  = 9
ng = 7

END SUBROUTINE
!------------------------------------------------------------------------------
SUBROUTINE get_M_1 (q, M)

INTEGER          :: q
REAL(8) :: M(q, q)

M = TRANSPOSE(reshape(DBLE((/ &
  1,  1,  1,  1,  1,  1,  1,  1,  1, &
  0,  1, -1,  0,  0,  1, -1, -1,  1, &
  0,  0,  0,  1, -1,  1, -1,  1, -1, &
 -4, -1, -1, -1, -1,  2,  2,  2,  2, &
  4, -2, -2, -2, -2,  1,  1,  1,  1, &
  0, -2,  2,  0,  0,  1, -1, -1,  1, &
  0,  0,  0, -2,  2,  1, -1,  1, -1, &
  0,  1,  1, -1, -1,  0,  0,  0,  0, &
  0,  0,  0,  0,  0,  1,  1, -1, -1/)), shape(M)))

END SUBROUTINE
!------------------------------------------------------------------------------
SUBROUTINE get_C_1 (C, lambda, g, Vel, ng, d, q)

INTEGER          :: d, q, ng, i
REAL(8) :: C(q, q), lambda(q), g(ng), Vel(d)

C(:,:) = DBLE(0)
do i = d + 2, q 
    C(i, i) = -lambda(i)
enddo
C(4, 1) = DBLE(1)/DBLE(4) * lambda(4) * g(1)
C(4, 2) = DBLE(1)/DBLE(3) * lambda(4) * g(5) * Vel(1)
C(4, 3) = DBLE(1)/DBLE(3) * lambda(4) * g(5) * Vel(2)
C(5, 1) = DBLE(1)/DBLE(4) * lambda(5) * g(2)
C(5, 2) = DBLE(1)/DBLE(3) * lambda(5) * g(7) * Vel(1)
C(5, 3) = DBLE(1)/DBLE(3) * lambda(5) * g(7) * Vel(2)
C(6, 2) = DBLE(1)/DBLE(2) * lambda(6) * g(3)
C(7, 3) = DBLE(1)/DBLE(2) * lambda(7) * g(3)
C(8, 2) = DBLE(3)         * lambda(8) * g(4) * Vel(1)
C(8, 3) = DBLE(-3)        * lambda(8) * g(4) * Vel(2)
C(9, 2) = DBLE(3)/DBLE(2) * lambda(9) * g(6) * Vel(2)
C(9, 3) = DBLE(3)/DBLE(2) * lambda(9) * g(6) * Vel(1)

END SUBROUTINE
!------------------------------------------------------------------------------
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!------------------------------------------------------------------------------
SUBROUTINE get_sizes_2 (d, q, ng)

INTEGER :: d, q, ng

d  = 2
q  = 9
ng = 4

END SUBROUTINE
!------------------------------------------------------------------------------
SUBROUTINE get_M_2 (d, q, M, dx)

INTEGER          :: d, q
REAL(8) :: M(q, q), dx(d)

M = TRANSPOSE(reshape(DBLE((/ &
  1,  1,  1,  1,  1,  1,  1,  1,  1, &
  0,  1, -1,  0,  0,  1, -1, -1,  1, &
  0,  0,  0,  1, -1,  1, -1,  1, -1, &
 -2,  1,  1, -2, -2,  1,  1,  1,  1, &
 -2, -2, -2,  1,  1,  1,  1,  1,  1, &
  0,  0,  0,  0,  0,  1,  1, -1, -1, &
  0, -2,  2,  0,  0,  1, -1, -1,  1, &
  0,  0,  0, -2,  2,  1, -1,  1, -1, &
  4, -2, -2, -2, -2,  1,  1,  1,  1/)), shape(M)))

M(2,:) = M(2,:) * dx(1)
M(3,:) = M(3,:) * dx(2)
M(4,:) = M(4,:) * dx(1)*dx(1)
M(5,:) = M(5,:) * dx(2)*dx(2)
M(6,:) = M(6,:) * dx(1)*dx(2)
M(7,:) = M(7,:) * dx(1)*dx(2)*dx(2)
M(8,:) = M(8,:) * dx(2)*dx(1)*dx(1)
M(9,:) = M(9,:) * dx(1)*dx(1)*dx(2)*dx(2)

END SUBROUTINE
!------------------------------------------------------------------------------
SUBROUTINE get_C_2 (C, lambda, g, Vel, ng, d, q, dx)

INTEGER          :: d, q, ng, i
REAL(8) :: C(q, q), lambda(q), g(ng), Vel(d), dx(d)

C(:,:) = DBLE(0)
do i = d + 2, q 
    C(i, i) = -lambda(i)
enddo
C(4, 1) = lambda(4)*( DBLE(3)*g(1) - DBLE(2)*dx(1)*dx(1) )
C(4, 2) = lambda(4)*DBLE(6)*Vel(1)
C(5, 1) = lambda(5)*( DBLE(3)*g(1) - DBLE(2)*dx(2)*dx(2) )
C(5, 3) = lambda(5)*DBLE(6)*Vel(2)
C(6, 2) = lambda(6)*Vel(2)
C(6, 3) = lambda(6)*Vel(1)
C(7, 2) = lambda(7)*( DBLE(3)*g(2) - DBLE(2)*dx(2)*dx(2) )
C(8, 3) = lambda(8)*( DBLE(3)*g(2) - DBLE(2)*dx(1)*dx(1) )
C(9, 1) = lambda(9)*( DBLE(9)*g(3) - DBLE(6)*g(1)*( dx(1)*dx(1) &
  + dx(2)*dx(2) ) + 4.*dx(1)*dx(1)*dx(2)*dx(2) )
C(9, 2) = lambda(9)*( DBLE(18)*g(4) - DBLE(12)*dx(2)*dx(2) ) * Vel(1)
C(9, 3) = lambda(9)*( DBLE(18)*g(4) - DBLE(12)*dx(1)*dx(1) ) * Vel(2)

END SUBROUTINE
!------------------------------------------------------------------------------
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!------------------------------------------------------------------------------
SUBROUTINE get_sizes_3 (d, q, ng)

INTEGER :: d, q, ng

d  = 3
q  = 27
ng = 8

END SUBROUTINE
!------------------------------------------------------------------------------
SUBROUTINE get_M_3 (d, q, M, dx)

INTEGER          :: d, q
REAL(8) :: M(q, q), dx(d)

M = TRANSPOSE(reshape(DBLE((/ &
  1,  1,  1,  1,  1,  1,  1,  1,  1, &
  1,  1,  1,  1,  1,  1,  1,  1,  1, &
  1,  1,  1,  1,  1,  1,  1,  1,  1, &

  0,  1, -1,  0,  0,  0,  0,  1, -1, &
 -1,  1,  1, -1, -1,  1,  0,  0,  0, &
  0,  1, -1,  1, -1,  1, -1, -1,  1, &

  0,  0,  0,  1, -1,  0,  0,  1, -1, &
  1, -1,  0,  0,  0,  0,  1, -1, -1, &
  1,  1, -1,  1, -1, -1,  1,  1, -1, &

  0,  0,  0,  0,  0,  1, -1,  0,  0, &
  0,  0,  1, -1,  1, -1,  1, -1,  1, &
 -1,  1, -1, -1,  1,  1, -1,  1, -1, &

 -2,  1,  1, -2, -2, -2, -2,  1,  1, &
  1,  1,  1,  1,  1,  1, -2, -2, -2, &
 -2,  1,  1,  1,  1,  1,  1,  1,  1, &

 -2, -2, -2,  1,  1, -2, -2,  1,  1, &
  1,  1, -2, -2, -2, -2,  1,  1,  1, &
  1,  1,  1,  1,  1,  1,  1,  1,  1, &

 -2, -2, -2, -2, -2,  1,  1, -2, -2, &
 -2, -2,  1,  1,  1,  1,  1,  1,  1, &
  1,  1,  1,  1,  1,  1,  1,  1,  1, &

  0,  0,  0,  0,  0,  0,  0,  1,  1, &
 -1, -1,  0,  0,  0,  0,  0,  0,  0, &
  0,  1,  1,  1,  1, -1, -1, -1, -1, &

  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  0,  1,  1, -1, -1,  0,  0,  0, &
  0,  1,  1, -1, -1,  1,  1, -1, -1, &

  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  0,  0,  0,  0,  0,  1,  1, -1, &
 -1,  1,  1, -1, -1, -1, -1,  1,  1, &

  0, -2,  2,  0,  0,  0,  0,  1, -1, &
 -1,  1, -2,  2,  2, -2,  0,  0,  0, &
  0,  1, -1,  1, -1,  1, -1, -1,  1, &

  0, -2,  2,  0,  0,  0,  0, -2,  2, &
  2, -2,  1, -1, -1,  1,  0,  0,  0, &
  0,  1, -1,  1, -1,  1, -1, -1,  1, &

  0,  0,  0, -2,  2,  0,  0,  1, -1, &
  1, -1,  0,  0,  0,  0, -2,  2,  2, &
 -2,  1, -1,  1, -1, -1,  1,  1, -1, &

  0,  0,  0, -2,  2,  0,  0, -2,  2, &
 -2,  2,  0,  0,  0,  0,  1, -1, -1, &
  1,  1, -1,  1, -1, -1,  1,  1, -1, &

  0,  0,  0,  0,  0, -2,  2,  0,  0, &
  0,  0,  1, -1,  1, -1, -2,  2, -2, &
  2,  1, -1, -1,  1,  1, -1,  1, -1, &

  0,  0,  0,  0,  0, -2,  2,  0,  0, &
  0,  0, -2,  2, -2,  2,  1, -1,  1, &
 -1,  1, -1, -1,  1,  1, -1,  1, -1, &

  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  1, -1, -1,  1, -1,  1, -1,  1, &

  4, -2, -2, -2, -2,  4,  4,  1,  1, &
  1,  1, -2, -2, -2, -2, -2, -2, -2, &
 -2,  1,  1,  1,  1,  1,  1,  1,  1, &

  4, -2, -2,  4,  4, -2, -2, -2, -2, &
 -2, -2,  1,  1,  1,  1, -2, -2, -2, &
 -2,  1,  1,  1,  1,  1,  1,  1,  1, &

  4,  4,  4, -2, -2, -2, -2, -2, -2, &
 -2, -2, -2, -2, -2, -2,  1,  1,  1, &
  1,  1,  1,  1,  1,  1,  1,  1,  1, &

  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  0,  0,  0,  0,  0, -6, -6,  6, &
  6,  3,  3, -3, -3, -3, -3,  3,  3, &

  0,  0,  0,  0,  0,  0,  0,  0,  0, &
  0,  0, -6, -6,  6,  6,  0,  0,  0, &
  0,  3,  3, -3, -3,  3,  3, -3, -3, &

  0,  0,  0,  0,  0,  0,  0, -6, -6, &
  6,  6,  0,  0,  0,  0,  0,  0,  0, &
  0,  3,  3,  3,  3, -3, -3, -3, -3, &

  0,  0,  0,  0,  0,  4, -4,  0,  0, &
  0,  0, -2,  2, -2,  2, -2,  2, -2, &
  2,  1, -1, -1,  1,  1, -1,  1, -1, &

  0,  0,  0,  4, -4,  0,  0, -2,  2, &
 -2,  2,  0,  0,  0,  0, -2,  2,  2, &
 -2,  1, -1,  1, -1, -1,  1,  1, -1, &

  0,  4, -4,  0,  0,  0,  0, -2,  2, &
  2, -2, -2,  2,  2, -2,  0,  0,  0, &
  0,  1, -1,  1, -1,  1, -1, -1,  1, &

 -8,  4,  4,  4,  4,  4,  4, -2, -2, &
 -2, -2, -2, -2, -2, -2, -2, -2, -2, &
 -2,  1,  1,  1,  1,  1,  1,  1,  1  &

  /)), shape(M)))

M(2,:)  = M(2,:)  * dx(1)
M(3,:)  = M(3,:)  * dx(2)
M(4,:)  = M(4,:)  * dx(3)
M(5,:)  = M(5,:)  * dx(1)*dx(1)
M(6,:)  = M(6,:)  * dx(2)*dx(2)
M(7,:)  = M(7,:)  * dx(3)*dx(3)
M(8,:)  = M(8,:)  * dx(1)*dx(2)
M(9,:)  = M(9,:)  * dx(1)*dx(3)
M(10,:) = M(10,:) * dx(2)*dx(3)
M(11,:) = M(11,:) * dx(1)*dx(2)*dx(2)
M(12,:) = M(12,:) * dx(1)*dx(3)*dx(3)
M(13,:) = M(13,:) * dx(2)*dx(1)*dx(1)
M(14,:) = M(14,:) * dx(2)*dx(3)*dx(3)
M(15,:) = M(15,:) * dx(3)*dx(1)*dx(1)
M(16,:) = M(16,:) * dx(3)*dx(2)*dx(2)
M(17,:) = M(17,:) * dx(1)*dx(2)*dx(3)
M(18,:) = M(18,:) * dx(1)*dx(1)*dx(2)*dx(2)
M(19,:) = M(19,:) * dx(1)*dx(1)*dx(3)*dx(3)
M(20,:) = M(20,:) * dx(2)*dx(2)*dx(3)*dx(3)
M(21,:) = M(21,:) * dx(1)*dx(1)*dx(2)*dx(3)
M(22,:) = M(22,:) * dx(2)*dx(2)*dx(1)*dx(3)
M(23,:) = M(23,:) * dx(3)*dx(3)*dx(1)*dx(2)
M(24,:) = M(24,:) * dx(1)*dx(1)*dx(2)*dx(2)*dx(3)
M(25,:) = M(25,:) * dx(1)*dx(1)*dx(3)*dx(3)*dx(2)
M(26,:) = M(26,:) * dx(2)*dx(2)*dx(3)*dx(3)*dx(1)
M(27,:) = M(27,:) * dx(1)*dx(1)*dx(2)*dx(2)*dx(3)*dx(3)

END SUBROUTINE
!------------------------------------------------------------------------------
SUBROUTINE get_C_3 (C, lambda, g, Vel, ng, d, q, dx)

INTEGER          :: d, q, ng, i
REAL(8) :: C(q, q), lambda(q), g(ng), Vel(d), dx(d)

C(:,:) = DBLE(0)
do i = d + 2, q 
    C(i, i) = -lambda(i)
enddo
C(5,  1) = lambda(5)*( DBLE(3)*g(1) - DBLE(2)*dx(1)*dx(1) )
C(5,  2) = lambda(5)*DBLE(6)*Vel(1)
C(6,  1) = lambda(6)*( DBLE(3)*g(1) - DBLE(2)*dx(2)*dx(2) )
C(6,  3) = lambda(6)*DBLE(6)*Vel(2)
C(7,  1) = lambda(7)*( DBLE(3)*g(1) - DBLE(2)*dx(3)*dx(3) )
C(7,  4) = lambda(7)*DBLE(6)*Vel(3)
C(8,  2) = lambda(8)*Vel(2)
C(8,  3) = lambda(8)*Vel(1)
C(9,  2) = lambda(9)*Vel(3)
C(9,  4) = lambda(9)*Vel(1)
C(10, 3) = lambda(10)*Vel(3)
C(10, 4) = lambda(10)*Vel(2)
C(11, 2) = lambda(11)*( DBLE(3)*g(2) - DBLE(2)*dx(2)*dx(2) )
C(12, 2) = lambda(12)*( DBLE(3)*g(2) - DBLE(2)*dx(3)*dx(3) )
C(13, 3) = lambda(13)*( DBLE(3)*g(2) - DBLE(2)*dx(1)*dx(1) )
C(14, 3) = lambda(14)*( DBLE(3)*g(2) - DBLE(2)*dx(3)*dx(3) )
C(15, 4) = lambda(15)*( DBLE(3)*g(2) - DBLE(2)*dx(1)*dx(1) )
C(16, 4) = lambda(16)*( DBLE(3)*g(2) - DBLE(2)*dx(2)*dx(2) )

C(18, 1) = lambda(18)*( DBLE(9)*g(3) - DBLE(6)*g(1)*( dx(1)*dx(1) &
  + dx(2)*dx(2) ) + 4.*dx(1)*dx(1)*dx(2)*dx(2) )
C(18, 2) = lambda(18)*( DBLE(18)*g(4) - DBLE(12)*dx(2)*dx(2) )*Vel(1)
C(18, 3) = lambda(18)*( DBLE(18)*g(4) - DBLE(12)*dx(1)*dx(1) )*Vel(2)

C(19, 1) = lambda(19)*( DBLE(9)*g(3) - DBLE(6)*g(1)*( dx(1)*dx(1) &
  + dx(3)*dx(3) ) + 4.*dx(1)*dx(1)*dx(3)*dx(3) )
C(19, 2) = lambda(19)*( DBLE(18)*g(4) - DBLE(12)*dx(3)*dx(3) )*Vel(1)
C(19, 4) = lambda(19)*( DBLE(18)*g(4) - DBLE(12)*dx(1)*dx(1) )*Vel(3)

C(20, 1) = lambda(20)*( DBLE(9)*g(3) - DBLE(6)*g(1)*( dx(2)*dx(2) &
  + dx(3)*dx(3) ) + 4.*dx(2)*dx(2)*dx(3)*dx(3) )
C(20, 3) = lambda(20)*( DBLE(18)*g(4) - DBLE(12)*dx(3)*dx(3) )*Vel(2)
C(20, 4) = lambda(20)*( DBLE(18)*g(4) - DBLE(12)*dx(2)*dx(2) )*Vel(3)

C(21, 3) = lambda(21)*( DBLE(9)*g(5) - DBLE(6)*dx(1)*dx(1) )*Vel(3)
C(21, 4) = lambda(21)*( DBLE(9)*g(5) - DBLE(6)*dx(1)*dx(1) )*Vel(2)

C(22, 2) = lambda(22)*( DBLE(9)*g(5) - DBLE(6)*dx(2)*dx(2) )*Vel(3)
C(22, 4) = lambda(22)*( DBLE(9)*g(5) - DBLE(6)*dx(2)*dx(2) )*Vel(1)

C(23, 2) = lambda(23)*( DBLE(9)*g(5) - DBLE(6)*dx(3)*dx(3) )*Vel(2)
C(23, 3) = lambda(23)*( DBLE(9)*g(5) - DBLE(6)*dx(3)*dx(3) )*Vel(1)

C(24, 4) = lambda(24)*( DBLE(9)*g(6) - DBLE(6)*g(2)*( dx(1)*dx(1) &
  + dx(2)*dx(2) ) + 4.*dx(1)*dx(1)*dx(2)*dx(2) )
C(25, 3) = lambda(25)*( DBLE(9)*g(6) - DBLE(6)*g(2)*( dx(1)*dx(1) &
  + dx(3)*dx(3) ) + 4.*dx(1)*dx(1)*dx(3)*dx(3) )
C(26, 2) = lambda(26)*( DBLE(9)*g(6) - DBLE(6)*g(2)*( dx(2)*dx(2) &
  + dx(3)*dx(3) ) + 4.*dx(2)*dx(2)*dx(3)*dx(3) )

C(27, 1) = lambda(27)*( DBLE(27)*g(7) - DBLE(18)*g(3)*( dx(1)*dx(1) &
  + dx(2)*dx(2) + dx(3)*dx(3) ) &
  + DBLE(12)*g(1)*( dx(1)*dx(1)*dx(2)*dx(2) + dx(1)*dx(1)*dx(3)*dx(3) &
  + dx(2)*dx(2)*dx(3)*dx(3) ) &
  - DBLE(8)*dx(1)*dx(1)*dx(2)*dx(2)*dx(3)*dx(3) )
C(27, 2) = lambda(27)*( DBLE(54)*g(8) - DBLE(36)*g(4)*( dx(2)*dx(2) &
  + dx(3)*dx(3) ) &
  + DBLE(24)*dx(2)*dx(2)*dx(3)*dx(3) )*Vel(1)
C(27, 3) = lambda(27)*( DBLE(54)*g(8) - DBLE(36)*g(4)*( dx(1)*dx(1) &
  + dx(3)*dx(3) ) &
  + DBLE(24)*dx(1)*dx(1)*dx(3)*dx(3) )*Vel(2)
C(27, 4) = lambda(27)*( DBLE(54)*g(8) - DBLE(36)*g(4)*( dx(1)*dx(1) &
  + dx(2)*dx(2) ) &
  + DBLE(24)*dx(1)*dx(1)*dx(2)*dx(2) )*Vel(3)

END SUBROUTINE
END MODULE lat_data

