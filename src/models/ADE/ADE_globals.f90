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

!> \file ADE_globals.f90
!! \brief Global bvariables for ADE problems.
!<


module ade_globals
  use typy
  
  !> parameters of sorption model
  type, public :: sorption_str
    logical :: kinetic=.false.
    !> name="freu" - Freundlich isoterm, "langmu" - Langmuir isoterm
    character(len=6) :: name
    real(kind=rkind) :: adsorb
    real(kind=rkind) :: desorb
    !>the third parameter in sorption model -- either n exponent in Freundlich or csmax in Langmuir
    real(kind=rkind) :: third
    !> bulk density
    real(kind=rkind) :: bd
    !> ratio  of solid media, if single solid medium ratio=1, if more solid media the sum between sorption(:)%ratio has to be equal 1.0
    real(kind=rkind) :: ratio 
    real(kind=rkind) :: csinit
  end type sorption_str


  !> ADE solute/material parameters array
  !<
  type, public :: soluteXsoil
    real(kind=rkind) :: difmol
    real(kind=rkind), dimension(:), allocatable :: diff_loc
    real(kind=rkind) :: anisoangle
    real(kind=rkind), dimension(:,:), allocatable :: diff
    real(kind=rkind), dimension(:), allocatable :: orders, lambda 
    real(kind=rkind) :: convection
    real(kind=rkind) :: water_cont
    character(len=2) :: icondtype
    real(kind=rkind) :: cmax
    real(kind=rkind) :: cinit
  end type soluteXsoil


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!--contaminant.conf variables--!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> structure of solute parameters
  type(soluteXsoil), dimension(:), allocatable, public :: adepar
  
  !> structure of sorption parameters and solid medium parameters, the layers are defined in lines, if solid medium scattered into different media, then use solumns
  type(sorption_str), dimension(:,:), allocatable, public :: sorption

  !> type of used sorption isotherm
  !! 0 - linear
  !! 1 - Friedrich exponential
  !! 2 - Langmuir
  !<
  integer(kind=ikind), public :: isotherm
  
  logical, public :: use_richards
  
  integer(kind=ikind) :: no_solids
  
  logical, public :: kinsorb
  
  logical, public :: use_sorption
  
  integer, public :: file_contaminant
  
end module ade_globals
