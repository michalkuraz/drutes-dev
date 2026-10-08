module ncglobvars
  use typy
  use datetime
  use global_objs

  integer :: netcdfID
  integer :: varid
  integer(kind=ikind) ::  ora_di_ini

  integer(kind=ikind), parameter :: missing = -9999
  type(datetime_t), public :: starttime, ncstart
  integer(kind=ikind) :: geograzone = 32
  real(kind=rkind), dimension(:), allocatable :: nodealt
  real(kind=rkind) :: Qmin 
  real(kind=rkind) :: LSdisp ! longitudinal dispersivity [m]; legacy name retained
  real(kind=rkind) :: LSdisp_transverse ! transverse dispersivity [m]
  logical :: LSsupg = .false.
  real(kind=rkind) :: LSsupg_factor = 1.0_rkind
  logical :: LSshock = .false.
  real(kind=rkind) :: LSshock_factor = 1.0_rkind
  logical :: LSbank_noflow = .false.
  logical :: LSconservative=.false., LSbalance=.false., LSstep_active=.false.
  logical :: LSclock_override=.false.
  real(kind=rkind) :: LSstate_time=0, LSprevious_time=0, LStrial_time=0, LSoverride_time=0
  real(kind=rkind) :: LSstep_start=0, LSstep_dt=0
  real(kind=rkind), allocatable :: LSdepth_old(:,:),LSdepth_new(:,:)
  logical, allocatable :: open_edges(:,:) ! original exterior edges, excluding Dirichlet ports
  logical, allocatable :: bank_edges(:,:) ! local edges (1,2), (2,3), (3,1)
  real(kind=rkind) :: cinit_ls
  integer(kind=ikind) :: channel_count = 1_ikind

	
	
	
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
    logical, dimension(:), allocatable :: activeel
    real(kind=rkind), dimension(:,:), allocatable :: fluxvct
    real(kind=rkind), dimension(:), allocatable :: cellarea
  end type ncfluxdata_type

  type(ncfluxdata_type) :: ncfluxdata
  
  real(kind=rkind), dimension(:,:), allocatable :: nccellxy
	
	
	
  real(kind=rkind), dimension(:,:), allocatable :: elslopes
	
  integer(kind=ikind) :: addedbc
  
  type(node) :: ncnodes, channel_nd
  type(element) :: ncelements, channel_el
  
  integer(kind=ikind), dimension(:), allocatable :: el2ncgrid
  real(kind=rkind) :: vref=1.7, Qref=1672.0



contains
  function adenc_coefficient_time() result(t)
    use globals, only: time
    real(kind=rkind) :: t
    t=time
    if (LSconservative) then
      t=LSstate_time
      if (LSstep_active) t=LStrial_time
    end if
    if (LSclock_override) t=LSoverride_time
  end function
end module ncglobvars
