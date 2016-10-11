module ade_globals
  use typy
  
  !> parameters of sorption model
  type, public :: sorption_str
    logical :: kinetic
    !> name="freu" - Freundlich isoterm, "langmu" - Langmuir isoterm
    character(len=6) :: name
    real(kind=rkind) :: adsorb
    real(kind=rkind) :: desorb
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
    real(kind=rkind), dimension(:), allocatable :: orders, lambda 
    real(kind=rkind) :: bd
    type(sorption_str) :: sorption
    real(kind=rkind) :: convection
    real(kind=rkind) :: water_cont
    character(len=2) :: icondtype
    real(kind=rkind) :: cmax
    real(kind=rkind) :: cinit
    real(kind=rkind) :: csinit
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
  
  logical, public :: with_richards
  
  
  integer, public :: file_contaminant
end module ade_globals