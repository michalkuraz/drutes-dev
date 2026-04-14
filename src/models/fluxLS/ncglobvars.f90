module ncglobvars
  use typy
  use datetime

  integer :: netcdfID
  integer :: varid
  integer :: dimid_time, dimid_lat, dimid_lon
  integer :: time_len, lat_len, lon_len
!  real(kind=rkind), dimension(:), allocatable :: timenc, lat, lon
  integer(kind=ikind), parameter :: missing = -9999
  type(datetime_t), public :: starttime, ncstart
  integer(kind=ikind) :: geograzone = 32
  real(kind=rkind), dimension(:), allocatable :: nodealt


end module ncglobvars
