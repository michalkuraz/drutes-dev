module ncglobvars
  use typy
  use datetime

  integer :: netcdfID
  integer :: varid
  integer(kind=ikind) ::  ora_di_ini

  integer(kind=ikind), parameter :: missing = -9999
  type(datetime_t), public :: starttime, ncstart
  integer(kind=ikind) :: geograzone = 32
  real(kind=rkind), dimension(:), allocatable :: nodealt
  real(kind=rkind) :: Qmin = 250.0
  real(kind=rkind) :: LSdisp = 20000.0

	
	
	
	type :: ncfluxdata_type
	  logical :: initialized = .false.
	  logical :: slice_loaded = .false.

	  character(len=64) :: lon_name  = "lon"
	  character(len=64) :: lat_name  = "lat"
	  character(len=64) :: time_name = "time"
	  character(len=64) :: var_name  = "Qrouted"

	  integer :: lon_varid = -1
	  integer :: lat_varid = -1
	  integer :: time_varid = -1
	  integer :: q_varid = -1

	  integer(kind=ikind) :: nlon = 0_ikind
	  integer(kind=ikind) :: nlat = 0_ikind
	  integer(kind=ikind) :: ntime = 0_ikind
	  integer(kind=ikind) :: current_time_index = -1_ikind

	  real(kind=rkind) :: fill_value = -9999.0_rkind
	  
	  logical :: has_fill = .false.

	  logical :: has_bounds = .false.

	  real(kind=rkind), dimension(:), allocatable :: lon
	  real(kind=rkind), dimension(:), allocatable :: lat
	  integer(kind=ikind), dimension(:), allocatable :: time

	  real(kind=rkind), dimension(:,:), allocatable :: qslice

	  real(kind=rkind), dimension(:,:), allocatable :: lat_bnds
	  real(kind=rkind), dimension(:,:), allocatable :: lon_bnds
	  
	  logical :: bounds_index_ready = .false.

	  real(kind=rkind) :: lat_b0 = 0.0_rkind
	  real(kind=rkind) :: lon_b0 = 0.0_rkind
	  real(kind=rkind) :: dlat_b = 0.0_rkind
	  real(kind=rkind) :: dlon_b = 0.0_rkind
	end type ncfluxdata_type

	type(ncfluxdata_type) :: ncfluxdata
	
	
	
	real(kind=rkind), dimension(:,:), allocatable :: elslopes
	
	integer(kind=ikind) :: addedbc



end module ncglobvars
