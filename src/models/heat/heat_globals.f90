module heat_globals
  use typy
  

  type, public :: heatpars_str
    real(kind=rkind) :: anisoangle
    real(kind=rkind), dimension(:,:), allocatable :: lambda
    real(kind=rkind), dimension(:), allocatable :: lambda_loc, convection
    real(kind=rkind) :: C_w, C, source, Tinit
  end type heatpars_str


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!--contaminant.conf/matrix.conf variables--!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> structure of solute parameters
  type(heatpars_str), dimension(:), allocatable, public :: heatpar
  integer, public :: file_heat
  
  logical, public :: with_richards
  
end module heat_globals
