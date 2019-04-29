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



module kinglobs
  use typy
  use global_objs
  
  type, public :: surfacend_str
    real(kind=rkind), dimension(3) :: xyz
    logical :: boundary
  end type surfacend_str
  
  type, public :: surface_el_str
    real(kind=rkind) :: sx, sy
    integer(kind=ikind) :: cover
  end type surface_el_str

  type(surfacend_str), dimension(:), allocatable :: watershed_nd
  type(surface_el_str), dimension(:), allocatable :: watershed_el
  
  real(kind=rkind), dimension(:), allocatable :: manning
  
  
  type, public :: raindata_str
    real(kind=rkind), dimension(2) :: xy
    type(smartarray_real), dimension(2) :: series
  end type raindata_str
  
  type(raindata_str), dimension(:), allocatable :: raindata
  integer(kind=ikind), dimension(:), allocatable :: el2pt

end module kinglobs


