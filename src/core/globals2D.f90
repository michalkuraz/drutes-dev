! Copyright 2008 Michal Kuraz, Petr Mayer


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


!> contains global variables related to 2D problem
module globals2D
  use typy



  !> structure to define crossection points
  type, public :: crossection
    real(kind=rkind), dimension(:,:), allocatable :: xyz
    integer(kind=ikind), dimension(:), allocatable :: element
    integer :: unit
  end type crossection


  !!!!---following variables are used only if internal mesh generator used (rectangular regular meshes are created)
  !> the lenght of the domain
  real(kind=rkind), public :: length2D
  !> the width of the domain
  real(kind=rkind), public :: width2D
  !> the mesh density
  real(kind=rkind), public :: density2D 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> number of crosssections
  integer(kind=ikind), public :: cross_count
  !> stores crossection border points
  !! the last position stores the cross line density for points
  !<
  real(kind=rkind), dimension(:,:), allocatable,  public :: cross_array_border
  !> stores points on cross section line
  type(crossection), dimension(:), allocatable, public :: cross_array
  !>  point density of the 2D velocity vector plot
  real(kind=rkind), public :: vector_density






  

end module globals2D