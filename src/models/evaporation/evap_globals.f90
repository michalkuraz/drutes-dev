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
  
  real(kind=rkind), public:: b1,b2,b3
  integer, public :: file_vapor
  real(kind=rkind), public:: resistance
  
   
  real(kind=rkind), parameter, public :: R_gas = 8.314
  
  real(kind=rkind), parameter, public :: MolW = 0.018015
    !> gravity acceleration [m.s^-2]
  real(kind=rkind), parameter, public :: gravity = 9.18
   !> Reference surface tension at 25 ~deg C g.s^-2
  real(kind=rkind), parameter, public :: gamma_0 = 71.89
  integer(kind=ikind), parameter, public ::re_order = 1
  integer(kind=ikind), parameter, public ::heat_order = 2
  real(kind=rkind), parameter, public :: f_c = 0.02
   !> Gain factor [-]
  integer(kind=ikind), parameter, public :: GwT = 7
  
  integer(kind=ikind), parameter, public :: C_liq = 4188 
  
  integer(kind=ikind), parameter, public :: C_vap = 1800
  
  integer(kind=ikind), parameter, public :: C_soil =  1920
  integer(kind=ikind), parameter, public :: C_air =  1200
  
end module evap_globals
