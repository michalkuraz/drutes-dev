module freeze_globs
  use typy
  
  
  type, public :: freeze_sys
    real(kind=rkind) :: alpha, n, m, Thr, Ths, snow_density, diameter 
    real(kind=rkind) :: Cs, Ci, Cl, Ca
    real(kind=rkind) :: C1, C2, C3, C4, C5, F1, F2, beta

    !> hydraulic conductivity tensor of second order
    real(kind=rkind), dimension(:,:), allocatable :: Ks
    !> diagonal values of the hydraulic conductivity tensor in local system of coordinates
    real(kind=rkind), dimension(:), allocatable   :: Ks_local
    !> angle of the anisothrophy axes with respect to global axes
    real(kind=rkind), dimension(:), allocatable   :: anisoangle
    real(kind=rkind) :: initcond, Tinit
    character(len=5) :: icondtype
    character(len=5) :: icondtypeRE
    character(len=4) :: material
    real(kind=rkind) :: top, bottom
    real(kind=rkind) :: sinkterm
  end type freeze_sys
  
  type(freeze_sys), dimension(:), allocatable, target, public :: freeze_par

  !> freeze exchange conductivity
  real(kind=rkind), public :: hc, cumfilt 
  
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
  real(kind=rkind), parameter, public :: Omega = 7
  
  !> Thermal conductivity [ W/m/K]
  real(kind=rkind), parameter, public :: thermal_cond = 0.5
  
  !> Reference surface tension at 25 ~deg C g.s^-2
  real(kind=rkind), parameter, public :: surf_tens_ref = 71.89
  
  integer, public :: file_freeze

  integer, public :: frz_pnt
  
  logical, public :: clap
  
  logical, public :: qlt_log
  
  integer(kind = ikind), parameter, public :: wat = 1
  integer(kind = ikind), parameter, public :: heat_proc = 2

end module freeze_globs
