

! Copyright 2008 Michal Kuraz, Petr Mayer, Copyright 2016  Michal Kuraz, Petr Mayer, Johanna Bloecher
! Copyright 2019  Michal Kuraz, Petr Mayer, Johanna Bloecher, Thomas Heinze

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

!> \file ltne_globals.f90
!! \brief Global variables for heat equation using local thermal non-equilibrium.
!<



module ltne_globals

  use typy
  
  type, public :: ltnepars_str
    !> grain diameter of solid matrix [m] 
    real(kind=rkind) :: grainDiameter
    !> density of the solid [kg/m^2]
    real(kind=rkind) :: densityS
    !> heat capacity of the solid [J/kg/K]
    real(kind=rkind) :: heatCapS
    !> heat capacity of the fluid [J/kg/K]
    real(kind=rkind) :: heatCapL
    !> thermal conductivity of the solid [W/m/K]
    real(kind=rkind) :: CondS
    !> convection vector (if water flux specified by user and not by the solution of RE)
    real(kind=rkind), dimension(:), allocatable :: convection
    !> scalar values: initial temperature (if constant per layer)
    real(kind=rkind) :: Tinit
    !> scalar values: additional heat source term
    real(kind=rkind) :: Sh
    !> saturated water content, used only if NOT coupled with the Richards equation
    real(kind=rkind) :: satwatcont
    !> actual water content, used only if NOT coupled with the Richards equation
    real(kind=rkind) :: actwatcont
    !> volumetrix flux, used only if NOT coupled with the Richards equation
    real(kind=rkind), dimension(3) :: q
  end type ltnepars_str
      
  !> Parameters for the empirical heat transfer model (Wakao et al. 1979)
  real(kind=rkind), public, parameter :: betah = 2.4e-5
  real(kind=rkind), public, parameter :: gammah = 285.6
  real(kind=rkind), public, parameter :: deltah = 1.0_rkind/3.0_rkind
  real(kind=rkind), public, parameter :: epsilonh = 2.7

  integer(kind=ikind), public, parameter :: RE_order=3
  integer(kind=ikind), public, parameter :: Ts_order=1
  integer(kind=ikind), public, parameter :: Tl_order=2
  
  
  !> thermal conductivity of the fluid [W/m/K]
  real(kind=rkind), public :: CondL
  
  !> dynamic viscosity of the fluid [Pa s] 
  real(kind=rkind), public :: viscL

  !> heat capacity of the fluid [J/kg/K]
  real(kind=rkind) :: heatCapL 
  
  
  !> structure of solute parameters, allocatable, dimension is the number of materials
  type(ltnepars_str), dimension(:), allocatable, public :: ltnepar

  !> configuration file unit
  integer, public :: file_ltne

  !> logical, if true, convection is obtained from solving RE
  logical, public :: with_richards

end module ltne_globals

