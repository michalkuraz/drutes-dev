! Copyright 2008 Michal Kuraz, Petr Mayer, Copyright 2016  Michal Kuraz, Petr Mayer, Johanna Bloecher, Copyright 2019  Michal Kuraz, Petr Mayer, Johanna Bloecher, Juliana Arbelaez

! This file is part of DRUtES.
! DRUtES is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! DRUtES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with DRUtES. If not, see <http://www.gnu.org/licenses/>.

!> \file evap_fnc.f90
!! \brief This module contains subroutines that read input information from config files and additional input files
!<

module evap_globals
  use typy
   
  real(kind=rkind), parameter, public :: R_gas = 
  real(kind=rkind), parameter, public :: MolW = 
  real(kind=rkind), parameter, public :: gravity = 
  real(kind=rkind), parameter, public :: gamma_0 = 
  integer(kind=ikind), parameter, public ::re_order = 1
  integer(kind=ikind), parameter, public ::heat_order = 2
end module evap_globals
