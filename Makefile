# Copyright (C) 2012 Vanja Zecevic
# Contact vanja.zecevic@sydney.uni.edu.au

# This file is part of lb_pert

# lb_pert is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# lb_pert is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


MAKEFLAGS += -j
FORTFLAGS = -Wall -O3 -s -fopenmp -Jobj/ -Iobj/

OBJECTS = \
  obj/lb_pert.o \
  obj/lat_data.o \
  obj/tests.o \
  obj/tools.o \
  obj/core.o \
  obj/cfg.o

bin/lb_pert: $(OBJECTS)
	gfortran $(FORTFLAGS) -llapack -o $@ $(OBJECTS)

# Objects
obj/lb_pert.o: src/lb_pert.f90
	gfortran $(FORTFLAGS) -c -o $@ $<

obj/lat_data.o: src/lat_data.f90
	gfortran $(FORTFLAGS) -c -o $@ $< 

obj/tests.o: src/tests.f90
	gfortran $(FORTFLAGS) -c -o $@ $< 

obj/tools.o: src/tools.f90
	gfortran $(FORTFLAGS) -c -o $@ $< 

obj/core.o: src/core.f90
	gfortran $(FORTFLAGS) -c -o $@ $< 

obj/cfg.o: src/cfg.f90
	gfortran $(FORTFLAGS) -c -o $@ $< 

# Depedancies
obj/core.o: obj/tools.mod

obj/tools.o: obj/lat_data.mod

obj/tests.o: obj/tools.mod

obj/tests.o: obj/core.mod

obj/tests.o: obj/cfg.mod

obj/lb_pert.o: obj/cfg.mod

obj/lb_pert.o: obj/tests.mod

# Modules
obj/cfg.mod: obj/cfg.o

obj/tests.mod: obj/tests.o

obj/lat_data.mod: obj/lat_data.o

obj/tools.mod: obj/tools.o

obj/core.mod: obj/core.o

# Clean
clean:
	rm -vf bin/*
	rm -vf obj/*

