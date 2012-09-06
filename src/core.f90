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
use cfg
use omp_lib
IMPLICIT NONE


CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE get_maxz_maxk_3D (dx_vct, k, S_inv, M, q, d, L, C, M_inv, z, &
  DUMMY, EWORK, LWORK, RWORK, INFO, maxz, maxk, maxz_pl, maxk_pl)

INTEGER :: INFO, ikx, iky, ikz, tid, d, q, LWORK, nkz
COMPLEX(8) :: DUMMY(1,1), L(q, q), S_inv(q, q), z(q), EWORK(LWORK)
REAL(8) :: k(d), RWORK(2*q), M(q, q), M_inv(q, q), C(q, q), maxk(d), maxz, &
  maxz_pl(nthreads), maxk_pl(nthreads, d), dx_vct(d)
EXTERNAL :: ZGEEV

IF (d.EQ.2) THEN
    nkz = 1
ELSE
    nkz = nk
ENDIF

!$OMP PARALLEL PRIVATE(ikx, iky, ikz, tid) &
!$OMP FIRSTPRIVATE(k, S_inv, L, z, EWORK, &
!$OMP   RWORK, maxz, maxk, INFO)
tid = OMP_GET_THREAD_NUM()
k(:) = DBLE(0)
!$OMP DO
DO ikx=1, nk
  DO iky=1, nk 
    DO ikz=1, nkz 

        CALL get_k(ikx, iky, ikz, d, k, dx_vct)

        ! Setup the inverse of the advection matrix S_inv
        CALL get_S_inv (S_inv, M, k, q, d)

        ! Get L
        CALL get_L_dt (L, C, M, M_inv, S_inv, q)

        !------------------------------------------------------------------
        ! Find eigenvalues using LAPACK driver ZGEEV
        CALL ZGEEV('N', 'N', q, L, q, z, DUMMY, 1, DUMMY, 1, EWORK, &
          LWORK, RWORK, INFO)

        ! ikx = iky = ikz = 1 is DC
        IF ( (ikx.NE.1) .OR. (iky.NE.1) .OR. (ikz.NE.1) ) THEN
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

