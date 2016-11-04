module dual_globals
    use typy

  type, public :: soilpar
    real(kind=rkind) :: alpha, n, m, Thr, Ths, Ss
    !> hydraulic conductivity tensor of second order
    real(kind=rkind), dimension(:,:), allocatable :: Ks
    real(kind=rkind), dimension(:), allocatable   :: Ks_local
    real(kind=rkind), dimension(:), allocatable   :: anisoangle
    real(kind=rkind) :: initcond
  end type soilpar

 type,public :: expar
  real(kind=rkind)::beta,a,gam_par,weightm,weightf
 end type expar
  
  !> soil and layer parameters
  type(soilpar), dimension(:), allocatable, public :: vgmatrix
  type(soilpar), dimension(:), allocatable, public :: vgexchange
  type(soilpar), dimension(:), allocatable, public :: vgfracture
  type(expar), dimension(:), allocatable, public :: exchange
  
  !> formula of the retention curve
  integer(kind=ikind), public :: retc_method
  !> if false getval returns hydraulic head H, if true getval returns pressure head h
  logical, public :: get_pressh_vals = .false.
end module dual_globals