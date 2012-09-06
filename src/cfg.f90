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

MODULE cfg
IMPLICIT NONE

REAL(8), PARAMETER :: PI = 3.141592653589793238d+0
INTEGER :: lattice, test_type, nthreads, nlambda, nvel, nk
REAL(8) :: dt, dx, dy, dz, xchar, vmax
NAMELIST /lb_pert_conf/ lattice, test_type, dt, nthreads, dx, dy, dz, &
  xchar, nlambda, nvel, nk, vmax

END MODULE cfg

