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

!> \file boussglob.f90
!! \brief Global variables for Boussinesq equation.
!<

module boussglob
  use typy
  use global_objs
  
  type(smartarray_real), dimension(2), public :: bouss_slopes
  type(smartarray_real), dimension(2), public :: bouss_K
  type(smartarray_real), dimension(2), public :: bouss_rain
  real(kind=rkind), public :: bouss_por
  real(kind=rkind), public :: bouss_icond

end module boussglob
