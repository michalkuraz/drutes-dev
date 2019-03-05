module freeze_globs
  use typy
  
  !> gravity acceleration [kg.s^-2]
  real(kind=rkind), parameter, public :: grav = 9.8196
  
  !>latent heat for water [J.kg^-1]
  real(kind=rkind), parameter, public :: Lf = 334e3
  
  !> specific heat capacity ice
  real(kind=rkind), parameter, public :: Ci = 4000
  
  
  !> specific heat capacity water
  real(kind=rkind), parameter, public :: Cl = 4200  
  
  !> reference temperature for Clapeyron
  real(kind=rkind), parameter, public :: Tref = 273.15


end module freeze_globs
