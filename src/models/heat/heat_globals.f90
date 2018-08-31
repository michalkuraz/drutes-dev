
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
    real(kind=rkind) :: anisoangle
    real(kind=rkind), dimension(:,:), allocatable :: lambda
    real(kind=rkind), dimension(:), allocatable :: lambda_loc, convection
    real(kind=rkind) :: C_w, C, source, Tinit
  end type heatpars_str


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!--contaminant.conf/matrix.conf variables--!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> structure of solute parameters
  type(heatpars_str), dimension(:), allocatable, public :: heatpar
  integer, public :: file_heat
  
  logical, public :: with_richards
  
end module heat_globals
