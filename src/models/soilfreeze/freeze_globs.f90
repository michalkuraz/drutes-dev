module freeze_globs
  use typy
  
  
  type, public :: freeze_sys
    real(kind=rkind) :: alpha, n, m, Thr, Ths, snow_density, diameter 
    real(kind=rkind) :: Cs, Ci, Cl, Ca, Ls, Li, Ll, La
    real(kind=rkind) :: C1, C2, C3, C4, C5, F1, F2, beta

    !> hydraulic conductivity tensor of second order
    real(kind=rkind), dimension(:,:), allocatable :: Ks
    !> diagonal values of the hydraulic conductivity tensor in local system of coordinates
    real(kind=rkind), dimension(:), allocatable   :: Ks_local
    !> angle of the anisothrophy axes with respect to global axes
    real(kind=rkind), dimension(:), allocatable   :: anisoangle
    real(kind=rkind) :: initcond, Tinit, iceini
    real(kind=rkind) :: Tinit_s, Tinit_l
    character(len=5) :: icondtype, icondtypeTs, icondtypeIce
    character(len=5) :: icondtypeRE
    character(len=4) :: material
    real(kind=rkind) :: top, bottom
    real(kind=rkind) :: sinkterm
  end type freeze_sys
  
  type(freeze_sys), dimension(:), allocatable, target, public :: freeze_par

  real(kind=rkind), dimension(:), allocatable, public :: T_air

  !> freeze exchange conductivity
  real(kind=rkind), public :: hc, cumfilt 
  
  !> gravity acceleration [m.s^-2]
  real(kind=rkind), parameter, public :: grav = 9.81
  
  !>latent heat of fusion [J.kg^-1]
  real(kind=rkind), parameter, public :: Lf = 333.7e3
  
  !> reference temperature for Clapeyron [K] (defined in RE globals now)
!   real(kind=rkind), parameter, public :: Tref = 273.15

  !> density of water [kg.m^-3]
  !real(kind=rkind), parameter, public :: rho_wat = 980
  
  !> density of ice [kg.m^-3]
  real(kind=rkind), parameter, public :: rho_ice = 910
  
    !> density of soil [kg.m^-3]
  real(kind=rkind), parameter, public :: rho_soil = 2650
  
  !> density of air [kg.m^-3]
  real(kind=rkind), parameter, public :: rho_air = 1.29
  
   !> Gain factor [-]
  real(kind=rkind), parameter, public :: gwt = 7
    
  !> Thermal conductivity [ W/m/K]
  real(kind=rkind), parameter, public :: thermal_cond = 0.5
  
  !> Reference surface tension at 25 ~deg C g.s^-2
  real(kind=rkind), parameter, public :: surf_tens_ref = 71.89
  
  !> dynamic viscosities of liquid water [Pa s] at 0 deg C
  real(kind=rkind), parameter, public :: ul = 1.787e-3
  
  !> dynamic viscosities of air [Pa s] at 0 deg C
  real(kind=rkind), parameter, public :: ua = 1.729e-5
  
  !> dynamic viscosities of ice [Pa s]
  real(kind=rkind), parameter, public :: ui = 10e12
  integer, public :: file_freeze

  integer, public :: frz_pnt
  
  logical, public :: clap, fr_rate
  
  logical, public :: qlt_log
      
  logical, public:: air
 !> Impedance factor [-] (see Lundin 1990, tested 2,4,6,8,10)
  real(kind=rkind), public :: Omega

  real(kind=rkind), public :: beta, beta_melt
  
  real(kind=rkind), public :: fac_scale, fac_add

  integer(kind = ikind), public :: wat, heat_proc, ice, heat_solid

end module freeze_globs
