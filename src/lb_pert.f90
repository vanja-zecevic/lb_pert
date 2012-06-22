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

PROGRAM main
use tests
IMPLICIT NONE
INTEGER :: lattice, test_type, nthreads, d, q, ng, ierr
REAL(8) :: dt, dx, dy, dz
REAL(8), ALLOCATABLE :: dx_tmp(:)
NAMELIST /lb_pert_conf/ lattice, test_type, dt, nthreads, dx, dy, dz

lattice=0
test_type=0
OPEN(2, FILE='lb_pert.conf', STATUS='OLD')
READ(2, NML=lb_pert_conf) 
WRITE(*, NML=lb_pert_conf)

CALL get_sizes (d, q, ng, lattice)

ALLOCATE ( dx_tmp(d), STAT = ierr)
IF (ierr /= 0) STOP "Not enough memory"
IF (d==2) THEN
    dx_tmp = DBLE((/ dx, dy /))
ELSE IF (d==3) THEN
    dx_tmp = DBLE((/ dx, dy, dz /))
ENDIF

SELECT CASE (test_type)
    CASE (1)
    CALL ev_vs_k(lattice, d, q, ng, dx_tmp)
    CASE (2)
    CALL vmax_vs_lambda(lattice, dt, nthreads, d, q, ng, dx_tmp)
    CASE (3)
    CALL vmax_vs_lambda_aligned(lattice, dt, nthreads, d, q, ng, dx_tmp)
    CASE (4)
    CALL single_ev(lattice, d, q, ng, dx_tmp)
    CASE (5)
    CALL find_best_lambda(lattice, dt, nthreads, d, q, ng, dx_tmp)
    CASE DEFAULT
    STOP
END SELECT

STOP
END PROGRAM main
!------------------------------------------------------------------------------

