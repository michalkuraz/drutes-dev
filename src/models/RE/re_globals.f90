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
  use global_objs


  !> rcza structure
  type, public :: rcza
    logical                                         :: use
    real(kind=rkind)                                :: val
    real(kind=rkind), dimension(:,:), allocatable   :: tab
  end type rcza


  type, public :: soilpar
    !> if method = vgfnc - van Genuuchten functions are used
    !! if method = table - constituve values are provided in table
    !<
    character(len=5) :: method
    !> use only if method = table 
    !! constable(1) = h
    !! constable(2) = theta
    !! constable(3) = K
    !! constable(4) = C
    !<
    real(kind=rkind), dimension(:), allocatable :: logh, Kr, theta, C, dKdh
    real(kind=rkind) :: step4fnc
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
  
  
  type, public :: root_str
    logical :: use
    !> parameter description
    !!  h_low - system is overwetted
    !!  h_start - system starts to be overwetted
    !!  h_end - system starts to be dry
    !!  h_end - system is overdried
    !<
    real(kind=rkind) :: h_wet, h_start, h_end, h_dry
    !> maximal water consumption
    real(kind=rkind) :: Smax
  end type root_str


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!--material parameters and methods for matrix--!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  logical, public :: re_transient = .true.

  !> formula of the retention curve
  integer(kind=ikind), public :: retc_method
  
  type(root_str), dimension(:), allocatable, public :: roots

  !> soil properties defined for each layer
  type(soilpar), dimension(:), allocatable, public :: vgset

  integer, public :: file_waterm
  integer, public :: file_evap
  integer, public :: file_roots
  
  !> evaporation input 
  real(kind=rkind), public :: albedo,latitude,elevation
  real(kind=rkind), public :: Tref = 273.15_rkind
  !> for evaporation
  integer(kind=ikind), public :: day_in_month, month_in_year,  init_year
  
  character(len=256) :: evap_name
  
  !> critical values
  real(kind=rkind), public :: h_crit_low=-1e5, h_crit_high=0, Kr_crit=1e-5, theta_crit=1e-5
  
   
end module re_globals
