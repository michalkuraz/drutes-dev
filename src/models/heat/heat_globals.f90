
! Copyright 2008 Michal Kuraz, Petr Mayer, Copyright 2016  Michal Kuraz, Petr Mayer, Johanna Bloecher

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

!> \file heat_globals.f90
!! \brief Global variables for heat equation.
!<



module heat_globals
  use typy
  
  
  type, public :: heatpars_str
    !> angle for anisothropy in 2D
    real(kind=rkind) :: anisoangle
    !> thermal conductivity (2nd order tensor) for global coordinates
    real(kind=rkind), dimension(:,:), allocatable :: lambda
    !> thermal conductivity (vector for local (unrotated) coordinates, in 2D just lambda_min and lambda_max. 
    real(kind=rkind), dimension(:), allocatable :: lambda_loc
    !> convection vector (if water flux specified by user and not by the solution of RE)
    real(kind=rkind), dimension(:), allocatable :: convection
    !> scalar values: specific capacity for liquid, solid, initial temperature (if constant per layer)
    real(kind=rkind) :: C_w, C, source, Tinit
  end type heatpars_str


  !> structure of solute parameters, allocatable, dimension is the number of materials
  type(heatpars_str), dimension(:), allocatable, public :: heatpar
  
  !> configuration file unit
  integer, public :: file_heat, file_heaticond
  
  !> logical, if true, convection is obtained from solving RE
  logical, public :: with_richards
  
  logical, public :: Tinitfile = .false.
  
  real(kind=rkind), dimension(:,:), allocatable :: Tinitpoints
  
  
  
  
end module heat_globals
