Copyright (C) 2012 Vanja Zecevic
Contact vanja.zecevic@sydney.uni.edu.au

This file is part of lb_pert

lb_pert is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

lb_pert is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

================================================================================
Introduction:
                              --------------------------------------------------
Numerical simulations performed using the lattice Boltzmann (LB) method may
become unstable depending on factors such as the relaxation rate and the
velocity of the flow. This software (lb_pert) performs a linear stability
analysis of the LB method over a range of parameters in order to predict stable
regions. LB_pert was written in Fortran and uses openMP in order to take
advantage of multi-core processors. In particular, the software investigates the
sensitivity of the method to a mean background velocity.

================================================================================
Author:
                              --------------------------------------------------
This software was written by, 

Vanja Zecevic
School of Aerospace, Mechanical and Mechatronic Engineering
The University of Sydney
Sydney, NSW, 2006
Australia

Please send any correspondence to,

vanja.zecevic@aeromech.usyd.edu.au


================================================================================
Brief description of operation:
                              --------------------------------------------------
Assuming a periodic domain with a uniform background flow and using a linear
approximation for the collision process, it is possible to form a recurrence
relation describing the evolution of a perturbation to the particle populations.
The solution to this relation is of the form,

f(t) = L^t f(0)

for some matrix L. The solution will remain bounded as time progresses if the
real component of the logarithm of the eigenvalues of L are below zero. In
general these eigenvalues will vary depending on the wavelength and direction of
the perturbation as well as the magnitude and direction of the background flow
and various collision parameters such as relaxation. 

================================================================================
Use:
                              --------------------------------------------------
To use lb_pert you will need to compile the program

cd lb_pert
make

then run,

bin/lb_pert

The program will read from the configuration file 'lb_pert.conf'. Other settings
need to be changed in the source code, more settings will be moved to the
configuration file in the future.

================================================================================
Lattices:
                              --------------------------------------------------
1 = D2Q9 Lallemand and Luo
2 = D2Q9 Zecevic et al.
3 = D3Q27 Zecevic et al.

================================================================================
Test routines:
                              --------------------------------------------------
1) ev_vs_k:
   This routine calculates all eigenvalues for a range of k.

2) vmax_vs_lambda_aligned
   This routine tests a range of relaxation rates (lambda) for the maximum
   stable velocity. The angle between the perturbation and the x-axis is varied
   and the velocity is aligned with the perturbation.

4) vmax_vs_lambda
   This routine tests a range of relaxation rates (lambda) for the maximum
   stable velocity. The angle between the perturbation and the x-axis is varied
   and the velocity is always in the x direction.

5) single_ev
   Get eigenvalues at a specified k, Vel and lambda.

6) find_best_lambda
   Undocumented, :)

================================================================================
Reference material:
                              --------------------------------------------------
Details regarding the linear perturbation analysis can be found in the following
papers,

[1] Lallemand, P. and Luo, L. S.
Theory of the lattice Boltzmann method: Dispersion, dissipation, isotropy,
  Galilean invariance, and stability
Phys. Rev. E 61 (2000)


