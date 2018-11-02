
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

!> \file re_dual_globals.f90
!! \brief Global variables for dual permeability model
!!<


module dual_globals
    use typy
    
!> defining soil parameters
  type, public :: soilpar
    !> inverse to air entry value or inverse to capillary rise [L^{-1}]
    real(kind=rkind) :: alpha
    !> shape parameter n
    real(kind=rkind) :: n
    !> shape parameter m
    real(kind=rkind) :: m
    !> residual water content 
    real(kind=rkind) :: Thr
    !> saturated water content 
    real(kind=rkind) :: Ths
    !> specific storage
    real(kind=rkind) :: Ss
    !> hydraulic conductivity tensor of second order
    real(kind=rkind), dimension(:,:), allocatable :: Ks
    !> hydraulic conductivity 
    real(kind=rkind), dimension(:), allocatable   :: Ks_local
    !> angle for anistropy of flow
    real(kind=rkind), dimension(:), allocatable   :: anisoangle
    !> initial condition
    real(kind=rkind) :: initcond
    !> type of initial condition: pressure head, water content or input read from file (also pressure head)
    character(len=5) :: icondtype
  end type soilpar
  
  !> exchange parameters 
  type, public :: exch_K
    !> inverse to air entry value or inverse to capillary rise [L^{-1}]
    real(kind=rkind) :: alpha
    !> shape parameter n
    real(kind=rkind) :: n
    !> shape parameter m
    real(kind=rkind) :: m
    !> hydraulic conductivity tensor of second order
    real(kind=rkind), dimension(:,:), allocatable :: Ks
    !> hydraulic conductivity 
    real(kind=rkind), dimension(:), allocatable   :: Ks_local
    !> angle for anistropy of flow
    real(kind=rkind), dimension(:), allocatable   :: anisoangle
  end type exch_K

!> coupling model type
  integer(kind=ikind), public :: coup_model
!> distance of top boundary to zero z-coordinate
  real(kind=rkind),public :: disttozero
!> Use-defined infiltration weight for boundary condition
  real(kind=rkind),public ::infweight
!> in progress: switching for different input files for fractrue and matrix (not implemented yet) 
  character(len=50) :: fracfile, matfile
 
!> exchange parameters according to Gerke and van Genuchten (1993)
 type,public :: expar
 !> soil geometry parameter 
  real(kind=rkind)::beta
  !> distance from fracture to center of agglomerate
  real(kind=rkind)::a
  !>geometry parameter, ususally 0.4
  real(kind=rkind)::gam_par
  !> matrix area weight
  real(kind=rkind)::weightm
  !> fracture area weight
  real(kind=rkind)::weightf
 end type expar
  
  !> soil and layer parameters
  type(soilpar), dimension(:), allocatable, public :: vgmatrix
  type(exch_K), dimension(:), allocatable, public :: vgexchange
  type(soilpar), dimension(:), allocatable, public :: vgfracture
  type(expar), dimension(:), allocatable, public :: exchange
  
  
  !> formula of the retention curve
  integer(kind=ikind), public :: retc_method
end module dual_globals
