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
  
  
  

end module evapglob
