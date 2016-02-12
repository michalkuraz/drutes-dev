module ade_globals
  use typy
  
  !> parameters of sorption model
  type, public :: sorption_str
    logical :: kinetic
    !> name="freu" - Freundlich isoterm, "langmu" - Langmuir isoterm
    character(len=4) :: name
    real(kind=rkind) :: sorb
    real(kind=rkind) :: adsord
    !>the third parameter in sorption model -- either n exponent in Freundlich or csmax in Langmuir
    real(kind=rkind) :: third
  end type sorption_str


  !> ADE solute/material parameters array
  !<
  type, public :: soluteXsoil
    real(kind=rkind) :: difmol
    real(kind=rkind), dimension(:), allocatable :: diff_loc
    real(kind=rkind) :: anisoangle
    real(kind=rkind), dimension(:,:), allocatable :: diff
    real(kind=rkind), dimension(:), allocatable :: lambda 
    real(kind=rkind) :: bd
    type(sorption_str) :: sorption
    logical :: with_Richards
    real(kind=rkind) :: convection
  end type soluteXsoil


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!--contaminant.conf/matrix.conf variables--!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> structure of solute parameters
  type(soluteXsoil), dimension(:), allocatable, public :: adepar


  !> type of used sorption isotherm
  !! 0 - linear
  !! 1 - Friedrich exponential
  !! 2 - Langmuir
  !<
  integer(kind=ikind), public :: isotherm
  
  
  integer, public :: file_contaminant
end module ade_globals