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
  
  
  

end module evapglob
