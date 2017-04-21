module dual_globals
    use typy

  type, public :: soilpar
    real(kind=rkind) :: alpha, n, m, Thr, Ths, Ss
    !> hydraulic conductivity tensor of second order
    real(kind=rkind), dimension(:,:), allocatable :: Ks
    real(kind=rkind), dimension(:), allocatable   :: Ks_local
    real(kind=rkind), dimension(:), allocatable   :: anisoangle
    real(kind=rkind) :: initcond
    
    character(len=5) :: icondtype
  end type soilpar
  
  type, public :: exch_K
    real(kind=rkind) :: alpha, n, m
    real(kind=rkind), dimension(:,:), allocatable :: Ks
    real(kind=rkind), dimension(:), allocatable   :: Ks_local
    real(kind=rkind), dimension(:), allocatable   :: anisoangle
  end type exch_K

 integer(kind=ikind), public :: coup_model
 real(kind=rkind),public :: disttozero
 real(kind=rkind),public ::infweight
 character(len=50) :: fracfile, matfile
 
 type,public :: expar
  real(kind=rkind)::beta,a,gam_par,weightm,weightf
 end type expar
  
  !> soil and layer parameters
  type(soilpar), dimension(:), allocatable, public :: vgmatrix
  type(exch_K), dimension(:), allocatable, public :: vgexchange
  type(soilpar), dimension(:), allocatable, public :: vgfracture
  type(expar), dimension(:), allocatable, public :: exchange
  
  
  !> formula of the retention curve
  integer(kind=ikind), public :: retc_method
end module dual_globals