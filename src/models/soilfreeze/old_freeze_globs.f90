module freeze_globs
  use typy
  
  !> gravity acceleration [kg.s^-2]
  real(kind=rkind), parameter, public :: grav = 9.81
  
  !>latent heat of fusion [J.kg^-1]
  real(kind=rkind), parameter, public :: Lf = 333.7e3
  
  !> specific heat capacity ice [J/kg/K]
  real(kind=rkind), parameter, public :: Ci = 2117
  
  
  !> specific heat capacity water [J/kg/K]
  real(kind=rkind), parameter, public :: Cl = 4188  
  
  !> reference temperature for Clapeyron [K]
  real(kind=rkind), parameter, public :: Tref = 273.15

  !> density of water [kg.m^-3]
  real(kind=rkind), parameter, public :: rho_wat = 1000
  
  !> density of ice [kg.m^-3]
  real(kind=rkind), parameter, public :: rho_ice = 918
  
   !> Gain factor [-]
  real(kind=rkind), parameter, public :: gwt = 7
  
  !> Thermal conductivity [ W/m/K]
  real(kind=rkind), parameter, public :: thermal_cond = 0.5
  
end module freeze_globs
