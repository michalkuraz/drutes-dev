module kinglobs
  use typy
  
  type, public :: surface_str
    real(kind=rkind), dimension(3) :: xyz
    logical :: boundary
    integer(kind=ikind) :: cover
    real(kind=rkind) :: sx, sy
  end type surface_str

  type(surface_str), dimension(:), allocatable :: xyz

end module kinglobs


