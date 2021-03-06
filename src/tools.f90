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

MODULE tools
use lat_data
use cfg
IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
SUBROUTINE setup_params (g, ng)

REAL(8) :: g(ng)
INTEGER :: ng

IF (lattice == 1) THEN
    g(1) = DBLE(-8)        ! alpha 2 
    g(2) = DBLE(4)         ! alpha 3 
    g(3) = DBLE(-2)        ! c 1 
    g(4) = DBLE(2)/DBLE(3) ! gamma 1 
    g(5) = DBLE(18)        ! gamma 2
    g(6) = DBLE(2)/DBLE(3) ! gamma 3
    g(7) = DBLE(-18)       ! gamma 4
ELSE IF (lattice == 2) THEN
    g(1) = (xchar**2)/DBLE(3) ! c_s^2
    g(2) = (xchar**2)/DBLE(3) ! gamma 1
    g(3) = (xchar**4)/DBLE(9) ! gamma 2
    g(4) = (xchar**2)/DBLE(3) ! gamma 3
ELSE IF (lattice == 3) THEN
    g(1) = (xchar**2)/DBLE(3) ! c_s^2 
    g(2) = (xchar**2)/DBLE(3)
    g(3) = (xchar**4)/DBLE(9)
    g(4) = (xchar**2)/DBLE(3)
    g(5) = (xchar**2)/DBLE(3)
    g(6) = (xchar**4)/DBLE(9)
    g(7) = (xchar**6)/DBLE(27)
    g(8) = (xchar**4)/DBLE(9)
ENDIF

END SUBROUTINE
!------------------------------------------------------------------------------
SUBROUTINE get_S_inv (S_inv, M, k, q, d)

INTEGER          :: q, d, i
COMPLEX(8)      :: S_inv(q, q)
REAL(8) :: M(q, q), k(d)

S_inv(:,:) = DCMPLX(0)
do i = 1, q 
    S_inv(i, i) = exp( DCMPLX(0, -1)*DOT_PRODUCT(k, M(2:d+1, i)))
enddo
END SUBROUTINE
!------------------------------------------------------------------------------
SUBROUTINE invert_M_orth (M, M_inv, q)

INTEGER          :: q, i
REAL(8) :: M(q, q), M_inv(q, q), WORK(q, q)

M_inv = TRANSPOSE(M)
WORK = MATMUL(M, M_inv)
do i = 1, q
    WORK(i, i) = DBLE(1)/WORK(i, i) 
enddo
M_inv = MATMUL(TRANSPOSE(M), WORK)

! Uncomment to check matrix inverse
!WORK = MATMUL(M, M_inv)
!do i = 1, q
!    WRITE (*,*) (WORK(i, j),j=1,q)
!enddo

END SUBROUTINE
!------------------------------------------------------------------------------
SUBROUTINE get_L (L, C, M, M_inv, S_inv, q)

INTEGER          :: q, i
COMPLEX(8)      :: L(q, q), S_inv(q, q), WORK1(q, q), WORK2(q, q)
REAL(8) :: C(q, q), M(q, q), M_inv(q, q)

WORK1 = MATMUL(C, M)
WORK2 = MATMUL(M_inv, WORK1)
DO i = 1, q 
    WORK2(i, i) = WORK2(i, i) + DCMPLX(1)
ENDDO

L = MATMUL(S_inv, WORK2)

END SUBROUTINE
!------------------------------------------------------------------------------
SUBROUTINE get_L_dt (L, C, M, M_inv, S_inv, q)

INTEGER          :: q, i
COMPLEX(8)      :: L(q, q), S_inv(q, q), WORK1(q, q), WORK2(q, q)
REAL(8) :: C(q, q), M(q, q), M_inv(q, q)

WORK1 = MATMUL(C, M)
WORK2 = MATMUL(M_inv, WORK1)
DO i = 1, q 
    WORK2(i, i) = WORK2(i, i) + DCMPLX(1)
ENDDO

L = dt*MATMUL(S_inv, WORK2)

do i = 1, q 
    L(i, i) = L(i, i) + DCMPLX(DBLE(1)-dt)
enddo  

END SUBROUTINE
!------------------------------------------------------------------------------
SUBROUTINE get_C (C, lambda, g, vel, ng, d, q, dx_vct)

INTEGER          :: d, q, ng
REAL(8) :: C(q, q), lambda(q), g(ng), vel(d), dx_vct(d)

SELECT CASE (lattice)
    CASE (1)
    CALL get_C_1 (C, lambda, g, vel, ng, d, q)
    CASE (2)
    CALL get_C_2 (C, lambda, g, vel, ng, d, q, dx_vct)
    CASE (3)
    CALL get_C_3 (C, lambda, g, vel, ng, d, q, dx_vct)
    CASE DEFAULT
    stop
END SELECT

END SUBROUTINE
!------------------------------------------------------------------------------
SUBROUTINE print_eigenvalues (z, q, INFO)

INTEGER     :: q, INFO, i
COMPLEX(8) :: z(q)

IF (INFO.EQ.0) THEN
    DO i=1, q
        WRITE (*,'(E12.4)',advance='no') REAL(LOG(z(i)))
    ENDDO
    WRITE (*,*)
ELSE
    WRITE (*,*)
    WRITE (*,*) 'Failure in ZGEEV.  INFO = ', INFO
END IF

END SUBROUTINE
!------------------------------------------------------------------------------
SUBROUTINE get_M (d, q, M, dx_vct)

INTEGER :: d, q
REAL(8) :: M(q, q), dx_vct(d)

SELECT CASE (lattice)
    CASE (1)
    CALL get_M_1 (q, M)
    CASE (2)
    CALL get_M_2 (d, q, M, dx_vct)
    CASE (3)
    CALL get_M_3 (d, q, M, dx_vct)
    CASE DEFAULT
    stop "Lattice not supported, exiting."
END SELECT

END SUBROUTINE
!------------------------------------------------------------------------------
SUBROUTINE get_sizes (d, q, ng)

INTEGER :: d, q, ng

SELECT CASE (lattice)
    CASE (1)
    CALL get_sizes_1 (d, q, ng)
    CASE (2)
    CALL get_sizes_2 (d, q, ng)
    CASE (3)
    CALL get_sizes_3 (d, q, ng)
    CASE DEFAULT
    stop "Lattice not supported, exiting."
END SELECT

END SUBROUTINE
!------------------------------------------------------------------------------
REAL(8) FUNCTION maxreln(z, q)

INTEGER :: q, i
COMPLEX(8) :: z(q)

maxreln = REAL(LOG(z(1)))

DO i=2, q
    IF (REAL(LOG(z(i))) > maxreln) THEN
        maxreln = REAL(LOG(z(i)))
    ENDIF
ENDDO

RETURN
END FUNCTION
!------------------------------------------------------------------------------
SUBROUTINE alloc_all (L, S_inv, z, EWORK, lambda, k, vel, RWORK, M, M_inv, C, &
  g, LWORK, d, q, ng) 

INTEGER     :: LWORK, d, q, ng, ierr
COMPLEX(8), ALLOCATABLE &
            :: L(:, :), S_inv(:, :), z(:), EWORK(:)
REAL(8), ALLOCATABLE &
            :: lambda(:), k(:), vel(:), RWORK(:), M(:, :), M_inv(:, :), &
               C(:, :), g(:)

ALLOCATE ( L(q, q), STAT = ierr)
IF (ierr /= 0) STOP "Not enough memory"

ALLOCATE ( S_inv(q, q), STAT = ierr)
IF (ierr /= 0) STOP "Not enough memory"

ALLOCATE ( z(q), STAT = ierr)
IF (ierr /= 0) STOP "Not enough memory"

LWORK = (1+64)*q
ALLOCATE ( EWORK(LWORK), STAT = ierr)
IF (ierr /= 0) STOP "Not enough memory"

ALLOCATE ( lambda(q), STAT = ierr)
IF (ierr /= 0) STOP "Not enough memory"

ALLOCATE ( k(d), STAT = ierr)
IF (ierr /= 0) STOP "Not enough memory"

ALLOCATE ( vel(d), STAT = ierr)
IF (ierr /= 0) STOP "Not enough memory"

ALLOCATE ( RWORK(2*q), STAT = ierr)
IF (ierr /= 0) STOP "Not enough memory"

ALLOCATE ( M(q, q), STAT = ierr)
IF (ierr /= 0) STOP "Not enough memory"

ALLOCATE ( M_inv(q, q), STAT = ierr)
IF (ierr /= 0) STOP "Not enough memory"

ALLOCATE ( C(q, q), STAT = ierr)
IF (ierr /= 0) STOP "Not enough memory"

ALLOCATE ( g(ng), STAT = ierr)
IF (ierr /= 0) STOP "Not enough memory"

END SUBROUTINE
!------------------------------------------------------------------------------
SUBROUTINE adjust_lambda (lambda, dx_vct, iq, id, q, d, g, ng)

INTEGER :: iq, id, q, d, ng
REAL(8) :: lambda(q), dx_vct(d), g(ng)

lambda(iq) = DBLE(-1)*(dx_vct(id)*dx_vct(id) - g(2))/ &
  ( (DBLE(-2)*g(2))/lambda(iq) &
  + (0.5d+0)*( (DBLE(3)*g(2)) - dx_vct(id)*dx_vct(id)))

END SUBROUTINE
!------------------------------------------------------------------------------
SUBROUTINE update_max (maxz, maxk, z, k, d, q)

INTEGER :: d, q
COMPLEX(8) :: z(q)
REAL(8) :: tmp, maxz, maxk(d), k(d)

tmp = maxreln(z, q)
IF (tmp > maxz) THEN
    maxz = tmp
    maxk(:) = k(:)
ENDIF

END SUBROUTINE
!------------------------------------------------------------------------------
SUBROUTINE gather_max(maxz_pl, maxk_pl, maxz, maxk, d, ngather)

INTEGER :: d, ngather, igather
REAL(8) :: maxz_pl(ngather), maxk_pl(ngather, d), maxz, maxk(d)

DO igather=1, ngather
    IF (maxz_pl(igather) > maxz) THEN
        maxz = maxz_pl(igather)
        maxk(:) = maxk_pl(igather, :)
    ENDIF
ENDDO

END SUBROUTINE
!------------------------------------------------------------------------------
SUBROUTINE get_k_discrete(ikx, iky, ikz, d, k, dx_vct)

INTEGER :: ikx, iky, ikz, d
REAL(8) :: k(d), dx_vct(d)

IF (ikx.EQ.1) THEN
    k(1) = DBLE(0)
ELSE
    k(1) = 2.0*PI/(DBLE(ikx)*dx_vct(1))
ENDIF
IF (iky.EQ.1) THEN
    k(2) = DBLE(0)
ELSE
    k(2) = 2.0*PI/(DBLE(iky)*dx_vct(2))
ENDIF
IF (d.EQ.3) THEN
    IF (ikz.EQ.1) THEN
        k(3) = DBLE(0)
    ELSE
        k(3) = 2.0*PI/(DBLE(ikz)*dx_vct(3))
    ENDIF
ENDIF

END SUBROUTINE
!------------------------------------------------------------------------------
SUBROUTINE get_k(ikx, iky, ikz, d, k, dx_vct)

INTEGER :: ikx, iky, ikz, d
REAL(8) :: k(d), dx_vct(d)

IF (ikx.EQ.1) THEN
    k(1) = DBLE(0)
ELSE
    k(1) = PI*DBLE(ikx)/(DBLE(nk)*dx_vct(1))
ENDIF
IF (iky.EQ.1) THEN
    k(2) = DBLE(0)
ELSE
    k(2) = PI*DBLE(iky)/(DBLE(nk)*dx_vct(2))
ENDIF
IF (d.EQ.3) THEN
    IF (ikz.EQ.1) THEN
        k(3) = DBLE(0)
    ELSE
        k(3) = PI*DBLE(ikz)/(DBLE(nk)*dx_vct(3))
    ENDIF
ENDIF

END SUBROUTINE
!------------------------------------------------------------------------------
SUBROUTINE print_matrix(A, ni, nj)

INTEGER :: ni, nj, i, j
REAL(8) :: A(ni, nj)

DO i=1, ni
    WRITE(*, "(100F7.2)") ( A(i,j), j=1,nj )
ENDDO

END SUBROUTINE
!------------------------------------------------------------------------------
END MODULE tools

