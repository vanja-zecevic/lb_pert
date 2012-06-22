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

MODULE core
use tools
use omp_lib
IMPLICIT NONE
REAL(8), PARAMETER :: PI = 3.141592653589793238d+0

CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE get_maxz_maxk_3D (k, nk, S_inv, M, q, d, L, C, M_inv, dt, z, &
  DUMMY, EWORK, LWORK, RWORK, INFO, maxz, maxk, maxz_pl, maxk_pl, nthreads)

INTEGER :: INFO, ikx, iky, ikz, tid, nk, d, q, LWORK, nthreads, nkz
COMPLEX(8) :: DUMMY(1,1), L(q, q), S_inv(q, q), z(q), EWORK(LWORK)
REAL(8) :: k(d), RWORK(2*q), M(q, q), M_inv(q, q), C(q, q), maxk(d), maxz, &
  dt, maxz_pl(nthreads), maxk_pl(nthreads, d)
EXTERNAL :: ZGEEV

IF (d.EQ.2) THEN
    nkz = 0
ELSE
    nkz = nk
ENDIF

!$OMP PARALLEL PRIVATE(ikx, iky, ikz, tid) &
!$OMP FIRSTPRIVATE(k, S_inv, L, z, EWORK, &
!$OMP   RWORK, maxz, maxk, INFO)
tid = OMP_GET_THREAD_NUM()
k(:) = DBLE(0)
!$OMP DO
DO ikx=0, nk
  DO iky=0, nk 
    DO ikz=0, nkz 
        IF (ikx.EQ.0) THEN
            k(1) = DBLE(0)
        ELSE
            k(1) = PI/DBLE(ikx)
        ENDIF
        IF (iky.EQ.0) THEN
            k(2) = DBLE(0)
        ELSE
            k(2) = PI/DBLE(iky)
        ENDIF
        IF (ikz.EQ.0) THEN
            k(3) = DBLE(0)
        ELSE
            k(3) = PI/DBLE(ikz)
        ENDIF

        ! Alternative calculation of k to give discrete wavelengths.
        !k(1) = PI*DBLE(ikx)/DBLE(nk)
        !k(2) = PI*DBLE(iky)/DBLE(nk)
        !k(3) = PI*DBLE(ikz)/DBLE(nk)

        ! Setup the inverse of the advection matrix S_inv
        CALL get_S_inv (S_inv, M, k, q, d)

        ! Get L
        CALL get_L_dt (L, C, M, M_inv, S_inv, q, dt)

        !------------------------------------------------------------------
        ! Find eigenvalues using LAPACK driver ZGEEV
        CALL ZGEEV('N', 'N', q, L, q, z, DUMMY, 1, DUMMY, 1, EWORK, &
          LWORK, RWORK, INFO)
        IF ( (ikx.NE.0) .OR. (iky.NE.0) .OR. (ikz.NE.0) ) THEN
            CALL update_max(maxz, maxk, z, k, d, q)
        ENDIF
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
maxz_pl(tid+1) = maxz
maxk_pl(tid+1, :) = maxk(:)

!$OMP END PARALLEL

CALL gather_max(maxz_pl, maxk_pl, maxz, maxk, d, nthreads)

END SUBROUTINE
!------------------------------------------------------------------------------
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!------------------------------------------------------------------------------
END MODULE core

