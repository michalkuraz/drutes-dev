module freeze_globs
  use typy
  
  !> gravity acceleration [m.s^-2]
  real(kind=rkind), parameter, public :: grav = 9.81
  
  !>latent heat of fusion [J.kg^-1]
  real(kind=rkind), parameter, public :: Lf = 333.7e3
  
  !> specific heat capacity ice [J/kg/K]
  real(kind=rkind), parameter, public :: Ci = 2117
  
  !> specific heat capacity water [J/kg/K]
  real(kind=rkind), parameter, public :: Cl = 4188  
  
  !> specific heat capacity air[J/kg/K]
  real(kind=rkind), parameter, public :: Ca = 1006  
  
    !> specific heat capacity soil[J/kg/K]
  real(kind=rkind), parameter, public :: Cs = 800 
  
  !> reference temperature for Clapeyron [K]
  real(kind=rkind), parameter, public :: Tref = 273.15

  !> density of water [kg.m^-3]
  real(kind=rkind), parameter, public :: rho_wat = 1000
  
  !> density of ice [kg.m^-3]
  real(kind=rkind), parameter, public :: rho_ice = 1000
  
    !> density of soil [kg.m^-3]
  real(kind=rkind), parameter, public :: rho_soil = 2650
  
  !> density of air [kg.m^-3]
  real(kind=rkind), parameter, public :: rho_air = 1.2
  
   !> Gain factor [-]
  real(kind=rkind), parameter, public :: gwt = 7
  
  !> Impedance factor [-] (see Lundin 1990, tested 2,4,6,8,10)
  real(kind=rkind), parameter, public :: Omega = 6
  
  !> Thermal conductivity [ W/m/K]
  real(kind=rkind), parameter, public :: thermal_cond = 0.5
  
  !> Reference surface tension at 25 ~deg C g.s^-2
  real(kind=rkind), parameter, public :: surf_tens_ref = 71.89
end module freeze_globs
