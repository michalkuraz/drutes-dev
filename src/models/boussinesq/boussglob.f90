module boussglob
  use typy
  use global_objs
  
  type(smartarray_real), dimension(2), public :: bouss_slopes
  type(smartarray_real), dimension(2), public :: bouss_K
  type(smartarray_real), dimension(2), public :: bouss_rain
  real(kind=rkind), public :: bouss_por
  real(kind=rkind), public :: bouss_icond

end module boussglob