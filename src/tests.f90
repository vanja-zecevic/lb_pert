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

MODULE tests
use cfg
use tools
use core
use omp_lib
IMPLICIT NONE

CONTAINS
!------------------------------------------------------------------------------
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!------------------------------------------------------------------------------
SUBROUTINE ev_vs_k(d, q, ng, dx_vct)

INTEGER :: INFO, ik, d, q, ng, LWORK
COMPLEX(8) :: DUMMY(1,1)
REAL(8) :: dx_vct(d), kbase
COMPLEX(8), ALLOCATABLE :: L(:, :), S_inv(:, :), z(:), EWORK(:)
REAL(8), ALLOCATABLE :: lambda(:), k(:), vel(:), RWORK(:), M(:, :), &
  M_inv(:, :), C(:, :), g(:)
EXTERNAL :: ZGEEV

!------------------------------------------------------------------------------
CALL alloc_all (L, S_inv, z, EWORK, lambda, k, vel, RWORK, M, M_inv, C, g, &
  LWORK, d, q, ng)
CALL setup_params(g, ng)

lambda(:) = 1.99d+0

IF (lattice == 2) THEN
    CALL adjust_lambda(lambda, dx_vct, 4, 1, q, d, g, ng)
    CALL adjust_lambda(lambda, dx_vct, 5, 2, q, d, g, ng)
    !lambda(7:8) = (0.9d+0)*lambda(7:8)
ELSE IF (lattice == 3) THEN
    CALL adjust_lambda(lambda, dx_vct, 5, 1, q, d, g, ng)
    CALL adjust_lambda(lambda, dx_vct, 6, 2, q, d, g, ng)
    CALL adjust_lambda(lambda, dx_vct, 7, 3, q, d, g, ng)
    !lambda(11:17) = (0.9d+0)*lambda(11:17)
    !lambda(24:26) = (0.9d+0)*lambda(24:26)
ENDIF

vel(:) = DBLE(0)
vel(1) = 0.2d+0

!------------------------------------------------------------------------------
! Do test
! Find M and M_inv
CALL get_M (d, q, M, dx_vct)
CALL invert_M_orth (M, M_inv, q)

CALL get_C (C, lambda, g, vel, ng, d, q, dx_vct)

k(:) = DBLE(0)
DO ik=0, nk 
    kbase = 1.5d+0 + 0.5d+0*DBLE(ik)/DBLE(nk)
    k(1) = kbase
    k(2) = kbase
    k(3) = kbase
    !k(1) = PI*DBLE(ik)/DBLE(nk)
    !k(2) = PI*DBLE(ik)/DBLE(nk)
    !k(3) = PI*DBLE(ik)/DBLE(nk)
    ! Setup the inverse of the advection matrix S_inv
    CALL get_S_inv (S_inv, M, k, q, d)

    ! Get L
    CALL get_L (L, C, M, M_inv, S_inv, q)

    !--------------------------------------------------------------------------
    ! Find eigenvalues using LAPACK driver ZGEEV
    CALL ZGEEV('N', 'N', q, L, q, z, DUMMY, 1, DUMMY, 1, EWORK, LWORK, RWORK, &
      INFO)
    CALL print_eigenvalues (z, q, INFO)
ENDDO

END SUBROUTINE
!------------------------------------------------------------------------------
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!------------------------------------------------------------------------------
SUBROUTINE vmax_vs_lambda(d, q, ng, dx_vct)

INTEGER :: INFO, d, q, ng, LWORK, ilambda, &
  ivel
COMPLEX(8) :: DUMMY(1,1)
COMPLEX(8), ALLOCATABLE :: L(:, :), S_inv(:, :), z(:), EWORK(:)
REAL(8), ALLOCATABLE :: lambda(:), k(:), vel(:), RWORK(:), M(:, :), &
  M_inv(:, :), C(:, :), g(:), unstablek(:), maxk(:)
REAL(8) :: maxz, vlow, vhgh, dx_vct(d)
REAL(8), ALLOCATABLE :: maxz_pl(:), maxk_pl(:, :)

!------------------------------------------------------------------------------
CALL alloc_all (L, S_inv, z, EWORK, lambda, k, vel, RWORK, M, M_inv, C, g, &
  LWORK, d, q, ng)

ALLOCATE ( unstablek(d) )
ALLOCATE ( maxk(d) )
ALLOCATE ( maxk_pl(nthreads, d) )
ALLOCATE ( maxz_pl(nthreads) )
CALL omp_set_num_threads(nthreads)

!------------------------------------------------------------------------------
CALL setup_params(g, ng)

! Find M and M_inv
CALL get_M (d, q, M, dx_vct)
CALL invert_M_orth (M, M_inv, q)

!CALL print_matrix(M, q, q)

lambda(:) = DBLE(0)
DO ilambda=0, nlambda
    lambda(:) = (1.8d+0) + DBLE(ilambda)*(0.2d+0)/DBLE(nlambda)

    IF (lattice == 2) THEN
        CALL adjust_lambda(lambda, dx_vct, 4, 1, q, d, g, ng)
        CALL adjust_lambda(lambda, dx_vct, 5, 2, q, d, g, ng)
        !lambda(7:8) = (0.9d+0)*lambda(7:8)
    ELSE IF (lattice == 3) THEN
        CALL adjust_lambda(lambda, dx_vct, 5, 1, q, d, g, ng)
        CALL adjust_lambda(lambda, dx_vct, 6, 2, q, d, g, ng)
        CALL adjust_lambda(lambda, dx_vct, 7, 3, q, d, g, ng)
        !!lambda(11:27) = lambda(11:27)*(1.9/2.0)
        !lambda(11:16) = lambda(11:16)*(1.8/2.0)
        !lambda(17)    = lambda(17)   *(1.8/2.0)
        !!lambda(18:20) = lambda(18:20)*(1.95/2.0)
        !!lambda(21:23) = lambda(21:23)*(1.95/2.0)
        !lambda(24:26) = lambda(24:26)*(1.8/2.0)
        !!lambda(27)    = lambda(27)*(1.9/2.0)
    ENDIF

    vel(:) = DBLE(0)
    vel(1) = vmax
    maxz = DBLE(0)
    CALL get_C (C, lambda, g, vel, ng, d, q, dx_vct)
    CALL get_maxz_maxk_3D (dx_vct, k, S_inv, M, q, d, L, C, M_inv, z, DUMMY, &
      EWORK, LWORK, RWORK, INFO, maxz, maxk, maxz_pl, maxk_pl)
    IF (maxz < DBLE(0)) THEN
        WRITE(*,*) "vmax is stable", lambda(1)
        EXIT
    ENDIF
    vlow = DBLE(0)
    vhgh = vmax
    DO ivel=0, nvel
        ! Bisect the interval
        vel(1) = 0.5d+0*(vlow+vhgh) 
        maxz = DBLE(0)
        CALL get_C (C, lambda, g, vel, ng, d, q, dx_vct)
        CALL get_maxz_maxk_3D (dx_vct, k, S_inv, M, q, d, L, C, M_inv, z, &
          DUMMY, EWORK, LWORK, RWORK, INFO, maxz, maxk, maxz_pl, maxk_pl)
        IF (maxz > DBLE(0)) THEN
            vhgh = vel(1)
            unstablek = maxk 
        ELSE
            vlow = vel(1)
        ENDIF
    ENDDO
    WRITE(*,*) vlow, lambda(1), unstablek 
ENDDO

END SUBROUTINE
!------------------------------------------------------------------------------
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!------------------------------------------------------------------------------
SUBROUTINE vmax_vs_lambda_aligned(d, q, ng, dx_vct)

INTEGER :: INFO, ik, d, q, ng, LWORK, ilambda, ivel, &
  itheta, ntheta, tid
COMPLEX(8) :: DUMMY(1,1)
COMPLEX(8), ALLOCATABLE :: L(:, :), S_inv(:, :), z(:), EWORK(:)
REAL(8), ALLOCATABLE :: lambda(:), k(:), vel(:), RWORK(:), M(:, :), &
  M_inv(:, :), C(:, :), g(:) , maxk(:)
REAL(8) :: maxz, vmag, dx_vct(d)
REAL(8), ALLOCATABLE :: maxz_pl(:), maxk_pl(:, :)
EXTERNAL :: ZGEEV

!------------------------------------------------------------------------------
CALL alloc_all (L, S_inv, z, EWORK, lambda, k, vel, RWORK, M, M_inv, C, g, &
  LWORK, d, q, ng)

ALLOCATE ( maxk(d) )

ALLOCATE ( maxk_pl(nthreads, d) )
ALLOCATE ( maxz_pl(nthreads) )
CALL omp_set_num_threads(nthreads)

!------------------------------------------------------------------------------
CALL setup_params(g, ng)

ntheta = 100

! Find M and M_inv
CALL get_M (d, q, M, dx_vct)
CALL Invert_M_Orth (M, M_inv, q)

lambda(:) = DBLE(0)
DO ilambda=0, nlambda
    lambda(:) = (1.8d+0) + DBLE(ilambda)*(0.2d+0)/DBLE(nlambda)

    IF (lattice == 2) THEN
        CALL adjust_lambda(lambda, dx_vct, 4, 1, q, d, g, ng)
        CALL adjust_lambda(lambda, dx_vct, 5, 2, q, d, g, ng)
        !lambda(7:8) = (0.9d+0)*lambda(7:8)
    ELSE IF (lattice == 3) THEN
        CALL adjust_lambda(lambda, dx_vct, 5, 1, q, d, g, ng)
        CALL adjust_lambda(lambda, dx_vct, 6, 2, q, d, g, ng)
        CALL adjust_lambda(lambda, dx_vct, 7, 3, q, d, g, ng)
        !lambda(11:17) = (0.9d+0)*lambda(11:17)
        !lambda(24:26) = (0.9d+0)*lambda(24:26)
    ENDIF

    DO ivel=0, nvel
        vmag = DBLE(ivel)*(0.4d+0)/DBLE(nvel) 
        maxz = DBLE(0)
        !!$OMP PARALLEL PRIVATE(itheta, ik, tid, INFO, vel) &
        !!$OMP FIRSTPRIVATE(C, k, S_inv, L, z, EWORK, &
        !!$OMP RWORK, maxz, maxk)
        vel(:) = DBLE(0)
        tid = OMP_GET_THREAD_NUM()
        !!$OMP DO
        DO itheta=0, ntheta
            vel(1) = vmag*cos(PI*(0.5d+0)*DBLE(itheta)/DBLE(ntheta))
            vel(2) = vmag*sin(PI*(0.5d+0)*DBLE(itheta)/DBLE(ntheta))

            CALL get_C (C, lambda, g, vel, ng, d, q, dx_vct)
            k(:) = DBLE(0)
            DO ik=1, nk 
                k(1) = PI*DBLE(ik)/DBLE(nk)*cos(PI*(0.5d+0)*DBLE(itheta) &
                  /DBLE(ntheta))
                k(2) = PI*DBLE(ik)/DBLE(nk)*sin(PI*(0.5d+0)*DBLE(itheta) &
                  /DBLE(ntheta))
                ! Setup the inverse of the advection matrix S_inv
                CALL get_S_inv (S_inv, M, k, q, d)

                ! Get L
                CALL get_L_dt (L, C, M, M_inv, S_inv, q)

                !---------------------------------------------------------------
                ! Find eigenvalues using LAPACK driver ZGEEV
                CALL ZGEEV('N', 'N', q, L, q, z, DUMMY, 1, DUMMY, 1, EWORK, &
                  LWORK, RWORK, INFO)
                CALL update_max(maxz, maxk, z, k, d, q) 
            ENDDO
        ENDDO
        !!$OMP END DO
        maxz_pl(tid+1) = maxz
        maxk_pl(tid+1, :) = maxk(:)

        !!$OMP END PARALLEL

        ! Print data if there was an unstable eigenvalue
        CALL gather_max(maxz_pl, maxk_pl, maxz, maxk, d, nthreads)

        IF (maxz > DBLE(0)) THEN
            WRITE(*,*) vmag, lambda(1), maxk
            EXIT
        ENDIF
    ENDDO
ENDDO

END SUBROUTINE
!------------------------------------------------------------------------------
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!------------------------------------------------------------------------------
SUBROUTINE single_ev(d, q, ng, dx_vct)

INTEGER :: INFO, i, d, q, ng, LWORK, maxzind, ierr
COMPLEX(8) :: DUMMY(1,1)
COMPLEX(8), ALLOCATABLE :: L(:, :), S_inv(:, :), z(:), EWORK(:), VR(:, :)
REAL(8), ALLOCATABLE :: lambda(:), k(:), vel(:), RWORK(:), M(:, :), &
  M_inv(:, :), C(:, :), g(:)
REAL(8) :: maxz, dx_vct(d)
EXTERNAL :: ZGEEV

!------------------------------------------------------------------------------
CALL alloc_all (L, S_inv, z, EWORK, lambda, k, vel, RWORK, M, M_inv, C, g, &
  LWORK, d, q, ng)

ALLOCATE ( VR(q, q), STAT = ierr)
IF (ierr /= 0) STOP "Not enough memory"

!------------------------------------------------------------------------------
CALL setup_params(g, ng)

lambda(:) = 1.9823d+0

vel(:) = DBLE(0)
vel(1) = 0.0577d+0 

!------------------------------------------------------------------------------
! Do test
! Find M and M_inv
CALL get_M (d, q, M, dx_vct)
CALL invert_M_orth (M, M_inv, q)

CALL get_C (C, lambda, g, vel, ng, d, q, dx_vct)

k(:) = DBLE(0)
k(1) = 1.7d+0
k(2) = 1.7d+0
! Setup the inverse of the advection matrix S_inv
CALL get_S_inv (S_inv, M, k, q, d)

! Get L
CALL get_L (L, C, M, M_inv, S_inv, q)

!--------------------------------------------------------------------------
! Find eigenvalues using LAPACK driver ZGEEV
CALL ZGEEV('N', 'V', q, L, q, z, DUMMY, 1, VR, q, EWORK, LWORK, RWORK, &
  INFO)

maxzind = 0
maxz = 0.
DO i=0, q
    IF (REAL(LOG(z(i))) > maxz) THEN
        maxz = REAL(LOG(z(i)))
        maxzind = i
    ENDIF
ENDDO

WRITE (*, *) maxzind
WRITE (*, *) REAL(LOG(z(maxzind)))
WRITE (*, *) VR(:, maxzind)

END SUBROUTINE
!------------------------------------------------------------------------------
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!------------------------------------------------------------------------------
SUBROUTINE find_best_lambda(d, q, ng, dx_vct)

INTEGER :: INFO, ikx, iky, d, q, ng, LWORK, &
  tid, first_val, ilambda7, ilambda8, ilambda9
COMPLEX(8) :: DUMMY(1,1)
COMPLEX(8), ALLOCATABLE :: L(:, :), S_inv(:, :), z(:), EWORK(:)
REAL(8), ALLOCATABLE :: lambda(:), k(:), vel(:), RWORK(:), M(:, :), &
  M_inv(:, :), C(:, :), g(:), maxk(:), best_lambda(:)
REAL(8) :: maxz, lowest_maxz, dx_vct(d)
REAL(8), ALLOCATABLE :: maxz_pl(:), maxk_pl(:, :)
EXTERNAL :: ZGEEV

!------------------------------------------------------------------------------
CALL alloc_all (L, S_inv, z, EWORK, lambda, k, vel, RWORK, M, M_inv, C, g, &
  LWORK, d, q, ng)

ALLOCATE ( maxk(d) )

ALLOCATE ( maxk_pl(nthreads, d) )
ALLOCATE ( maxz_pl(nthreads) )
ALLOCATE ( best_lambda(q) )
CALL omp_set_num_threads(nthreads)

!------------------------------------------------------------------------------
CALL Setup_Params(g, ng)

! Find M and M_inv
CALL get_M (d, q, M, dx_vct)
CALL invert_M_orth (M, M_inv, q)

lambda(:) = (1.99d+0)
IF (lattice == 2) THEN
    CALL adjust_lambda(lambda, dx_vct, 4, 1, q, d, g, ng)
    CALL adjust_lambda(lambda, dx_vct, 5, 2, q, d, g, ng)
ELSE IF (lattice == 3) THEN
    CALL adjust_lambda(lambda, dx_vct, 5, 1, q, d, g, ng)
    CALL adjust_lambda(lambda, dx_vct, 6, 2, q, d, g, ng)
    CALL adjust_lambda(lambda, dx_vct, 7, 3, q, d, g, ng)
ENDIF

vel(:) = DBLE(0)
vel(1) = 0.03d+0
first_val = 0
DO ilambda7=0, nlambda
DO ilambda8=0, nlambda
DO ilambda9=0, nlambda
    lambda(7) = lambda(1)*(0.5d+0)*(1.0d+0 + DBLE(ilambda7)/DBLE(nlambda))
    lambda(8) = lambda(1)*(0.5d+0)*(1.0d+0 + DBLE(ilambda8)/DBLE(nlambda))
    lambda(9) = lambda(1)*(0.5d+0)*(1.0d+0 + DBLE(ilambda9)/DBLE(nlambda))
    maxz = DBLE(0)

    CALL get_C (C, lambda, g, vel, ng, d, q, dx_vct)

    !!$OMP PARALLEL PRIVATE(ikx, iky, tid, INFO) &
    !!$OMP FIRSTPRIVATE(k, S_inv, L, z, EWORK, &
    !!$OMP RWORK, maxz, maxk)
    tid = OMP_GET_THREAD_NUM()
    k(:) = DBLE(0)
    !!$OMP DO
    DO ikx=0, nk
      DO iky=0, nk 
        k(1) = PI*DBLE(ikx)/DBLE(nk)
        k(2) = PI*DBLE(iky)/DBLE(nk)
        ! Setup the inverse of the advection matrix S_inv
        CALL get_S_inv (S_inv, M, k, q, d)

        ! Get L
        CALL get_L_dt (L, C, M, M_inv, S_inv, q)

        !---------------------------------------------------------------
        ! Find eigenvalues using LAPACK driver ZGEEV
        CALL ZGEEV('N', 'N', q, L, q, z, DUMMY, 1, DUMMY, 1, EWORK, &
          LWORK, RWORK, INFO)
        IF ( (ikx.NE.0) .OR. (iky.NE.0) ) THEN
          CALL update_max(maxz, maxk, z, k, d, q) 
        ENDIF 
     ENDDO
    ENDDO
    !!$OMP END DO
    maxz_pl(tid+1) = maxz
    maxk_pl(tid+1, :) = maxk(:)

    !!$OMP END PARALLEL

    ! Print data if there was an unstable eigenvalue
    CALL gather_max(maxz_pl, maxk_pl, maxz, maxk, d, nthreads)

    IF ((maxz <= lowest_maxz) .OR. (first_val.EQ.0)) THEN
        lowest_maxz = maxz
        best_lambda(:) = lambda(:)
    ENDIF
    first_val = 1
ENDDO
ENDDO
WRITE(*,*) ilambda7, lowest_maxz, best_lambda(7:9)
ENDDO

WRITE(*,*) "Best lambda:"
WRITE(*,*) best_lambda(7:9)

END SUBROUTINE

END MODULE tests

