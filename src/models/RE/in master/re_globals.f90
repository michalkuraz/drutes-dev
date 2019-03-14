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

!> \file re_globals.f90
!! \brief Global variables for the Richards equation.
!<


module re_globals
  use typy


  !> rcza structure
  type, public :: rcza
    logical                                         :: use
    real(kind=rkind)                                :: val
    real(kind=rkind), dimension(:,:), allocatable   :: tab
  end type rcza


  type, public :: soilpar
    real(kind=rkind) :: alpha, n, m, Thr, Ths, Ss
    !> hydraulic conductivity tensor of second order
    real(kind=rkind), dimension(:,:), allocatable :: Ks
    !> diagonal values of the hydraulic conductivity tensor in local system of coordinates
    real(kind=rkind), dimension(:), allocatable   :: Ks_local
    !> angle of the anisothrophy axes with respect to global axes
    real(kind=rkind), dimension(:), allocatable   :: anisoangle
    real(kind=rkind) :: initcond
    character(len=5) :: icondtype
    type(rcza)       :: rcza_set
    real(kind=rkind) :: top, bottom
    real(kind=rkind) :: sinkterm
  end type soilpar


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!--material parameters and methods for matrix--!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(soilpar), dimension(:), allocatable, target, public :: vgmatrix
  !> formula of the retention curve
  integer(kind=ikind), public :: retc_method

  !> soil properties defined for each layer
  type(soilpar), dimension(:), allocatable, public :: vgset

  integer, public :: file_waterm
   
end module re_globals
