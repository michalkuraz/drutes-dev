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

!> \file: evap_globals.f90
!! \brief: This module contains the global variable for the Vapour model 
!<

module evap_globals
  use typy
  
  !> Empirical regression parameters for the Thermal conductivity [Wm^-1 K^-1]
  real(kind=rkind), public:: b1,b2,b3
  
  real(kind=rkind), public:: resistance
  
  !> Universal Gas constant [J mol^-1 K^-1]
  real(kind=rkind), parameter, public :: R_gas = 8.314
  
  !> Molecular weight of water [kg mol^-1]
  real(kind=rkind), parameter, public :: MolWat = 0.018015
  
  !> gravity acceleration [m.s^-2]
  real(kind=rkind), parameter, public :: gravity = 9.18
  
  !> Reference surface tension at 25 ~deg C g.s^-2
  real(kind=rkind), parameter, public :: gamma_0 = 71.89
  
  !> Richards modified equation ID
  integer(kind=ikind), parameter, public ::re_order = 1
  
  !> Heat  modified equation ID
  integer(kind=ikind), parameter, public ::heat_order = 2
  
  !> Mass fraction of clay [-]
  real(kind=rkind), parameter, public :: f_c = 0.02
  
  !> Gain factor [-]
  integer(kind=ikind), parameter, public :: GwT = 7
  
  !> specific heat capacity of liquid  water [J/kg K]
  integer(kind=ikind), parameter, public :: C_liq = 4188 
  
  !> specific heat capacity of  water vapor [J/kg K]
  integer(kind=ikind), parameter, public :: C_vap = 1800
  
  !> specific heat capacity of  soil[J/kg K]
  integer(kind=ikind), parameter, public :: C_soil =  1920
  
  !> specific heat capacity of  air [J/kg K]
  integer(kind=ikind), parameter, public :: C_air =  1006
  
  !> density of soil [kg.m^-3]
  real(kind=rkind), parameter, public :: rho_soil = 2650
  
  !> density of air [kg.m^-3]
  real(kind=rkind), parameter, public :: rho_air = 1.29
  
  !> Inout file for vapur model
  integer, public :: file_vapor
  

  
  
end module evap_globals
