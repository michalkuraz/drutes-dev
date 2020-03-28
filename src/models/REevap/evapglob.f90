module evapglob
  use typy
  
  integer(kind=ikind), public :: RE_ord = 1, heat_ord = 2
  !> Molecular weight of water [kg mol^-1]
  real(kind=rkind), parameter, public :: MolWat = 0.018015
    !> gravity acceleration [m.s^-2]
  real(kind=rkind), parameter, public :: gravity = 9.81
  
  !> Universal Gas constant [J mol^-1 K^-1]
  real(kind=rkind), parameter, public :: R_gas = 8.314
  
  !> Gain factor [-]
  integer(kind=ikind), parameter, public :: GwT = 7
  
  !> Reference surface tension at 25 ~deg C g.s^-2 (units will cancel out)
  real(kind=rkind), parameter, public :: gamma_0 = 71.89
  
    !> Mass fraction of clay [-]
  real(kind=rkind), parameter, public :: fraction_clay = 0.02
  
    
  !> specific heat capacity of liquid  water [J/kg K]
  integer(kind=ikind), parameter, public :: C_liq = 4188 
  
  !> specific heat capacity of  water vapor [J/kg K]
  integer(kind=ikind), parameter, public :: C_vap = 1800
  
  !> aerodynamic resistance [s/m]
  real(kind=rkind), public :: resistance = 30
  
  !> current air temperature [dg. C]
  real(kind=rkind), public :: temp_air = 25
  
  !> incoming short wave radiation [W.m^-2]
  real(kind=rkind), public :: radiation  =  100
  
  !> relative air huminidity
  real(kind=rkind), public :: rel_air_hum  = 50.0
  
  !> reflactance of short wave radiation
  real(kind=rkind), public :: wave_albedo = 0.2
  
  !> specific heat capacity of  air [J/kg K]
  integer(kind=ikind), parameter, public :: C_air =  1006
  
  !> density of air [kg.m^-3]
  real(kind=rkind), parameter, public :: rho_air = 1.29
  
 

end module evapglob
