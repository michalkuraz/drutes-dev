module re_globals
  use typy


  !> rcza structure
  type, public :: rcza
    logical                                         :: use
    real(kind=rkind)                                :: val
    real(kind=rkind), dimension(:,:), allocatable   :: tab
  end type rcza


  type, public :: soilpar
    real(kind=rkind) :: alpha, n, m, Thr, Ths, Ss
    !> hydraulic conductivity tensor of second order
    real(kind=rkind), dimension(:,:), allocatable :: Ks
    !> diagonal values of the hydraulic conductivity tensor in local system of coordinates
    real(kind=rkind), dimension(:), allocatable   :: Ks_local
    !> angle of the anisothrophy axes with respect to global axes
    real(kind=rkind), dimension(:), allocatable   :: anisoangle
    real(kind=rkind) :: initcond
    character(len=5) :: icondtype
    type(rcza)       :: rcza_set
    real(kind=rkind) :: top, bottom
    real(kind=rkind) :: sinkterm
  end type soilpar


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!--material parameters and methods for matrix--!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> soil and layer parameters
  type(soilpar), dimension(:), allocatable, target, public :: vgmatrix
  !> formula of the retention curve
  integer(kind=ikind), public :: retc_method



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!--material parameters and init conditions for fractures--!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> an array with values of dual richards equation coupling term and the fast domain ratio
  !! dual_prop(i,1) = layer i, parameter: omega_f \n
  !! dual_prop(i,2) = layer i, parameter: \f[ \frac{\beta}{\alpha_{DP}^2} \gamma_w ]\f
  !<
  real(kind=rkind), dimension(:,:), allocatable, public :: dual_prop
  !> soil and layer parameters
  type(soilpar), dimension(:), allocatable, target, public :: vgfractures
  !> formula of the retention curve
  integer(kind=ikind), public :: retc_method_f
  !> depending on current domain, this pointer is pointed to current parameters as
  !! if (domain == matrix) vgset => vgmatrix
  !! if (domain == fractures) vgset => vgfractures
  !<
   type(soilpar), dimension(:), pointer, public :: vgset
   
   integer, public :: file_waterm
   
end module re_globals