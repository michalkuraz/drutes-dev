module ncglobvars
  use typy
  use datetime

  integer :: netcdfID
  integer :: ncDEM 
  integer :: varid
  integer :: dimid_time, dimid_lat, dimid_lon
  integer :: time_len, lat_len, lon_len
  real(kind=rkind), dimension(:), allocatable :: timenc, lat, lon
  integer(kind=ikind), parameter :: missing = -9999
  type(datetime_t), public :: starttime, ncstart


end module ncglobvars
