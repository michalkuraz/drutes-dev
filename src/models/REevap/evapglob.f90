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


module evapglob
  use typy
  
  !> process IDs
  integer(kind=ikind), parameter, public :: RE_ord = 1, heat_ord = 2

  !> Molecular weight of water [kg mol^-1]
  real(kind=rkind), parameter, public :: MolWat = 0.018015
    !> gravity acceleration [m.s^-2]
  real(kind=rkind), parameter, public :: gravity = 9.80665
  
  !> Universal Gas constant [J mol^-1 K^-1]
  real(kind=rkind), parameter, public :: R_gas = 8.314
  
    !> karman constant [-]
  real(kind=rkind), parameter, public :: karman = 0.41
  
    !> Mass fraction of clay [-]
  real(kind=rkind), parameter, public :: fraction_clay = 0.02
  
    
  !> specific heat capacity of liquid  water [J/kg K]
  real(kind=rkind), parameter, public :: C_liq = 4188 
  
  !> specific heat capacity of  water vapor [J/kg K]
  real(kind=rkind), parameter, public :: C_vap = 1800
    
  !> specific heat capacity of  air [J/kg K]
  real(kind=rkind), parameter, public :: C_air =  1006
  
!  !> density of air [kg.m^-3]
!  real(kind=rkind), parameter, public :: rho_air = 1.29
  
  !> reference level fot the temperature measurement
  real(kind=rkind), public :: zref
  
  !> rainfall data step
  !! two options (1) day and (2) hrs
  !! if daily time step we sum up together evaporation (negative)  and rainfall intensity (positive)
  !! if hourly, then we either consider evaporation or rainfall
  !<
  character(len=3) :: rainfall_step
  
  type, public :: ebalance_str
    real(kind=rkind) :: day
    real(kind=rkind) :: latitude, noon, sunsent_angle, solar_decl
  end type ebalance_str
  
  type, public :: soil_heat_str
    real(kind=rkind) :: b1, b2, b3
  end type soil_heat_str
  
  type, public :: meteo4evap_str
    real(kind=rkind) :: time
    real(kind=rkind) :: inradiation
    real(kind=rkind) :: T_air
    real(kind=rkind) :: wind_speed
    real(kind=rkind) :: cloudiness
    real(kind=rkind) :: relhum
    integer(kind=ikind) :: datapos = 1
  end type meteo4evap_str
  
  type, public :: albedo_str
    !> method = 1 : read from file
    !> method = 2 : compute from formula 
    integer(kind=ikind) :: method
    real(kind=rkind) :: theta_min, theta_max, albd_min, albd_max, A
    real(kind=rkind), dimension(:,:), allocatable :: albdat
    integer(kind=ikind) :: datapos=1
  end type albedo_str
  
    
  type(meteo4evap_str), dimension(:), allocatable :: meteo4evap
  
  type(soil_heat_str), dimension(:), allocatable :: soil_heat_coef
  
  type(ebalance_str), dimension(:), allocatable, public :: pars4ebalance
  
  type(albedo_str), public :: albedo_conf
  
  real(kind=rkind), dimension(:,:), allocatable :: evap4rain
  
  integer(kind=ikind) :: evap4rain_pos

end module evapglob
